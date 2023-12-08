use std::arch::x86_64::*;
use std::thread;
use rand::Rng;

use crate::field::Field;
use crate::field::common::{M128, M512};
use crate::util::*;
use crate::hash::*;
use sha3::{Shake128, digest::{Update, ExtendableOutput, XofReader}};


#[derive(Clone)]
struct TapeHandle(*mut __m512i);
unsafe impl Send for TapeHandle {}
unsafe impl Sync for TapeHandle {}
impl Copy for TapeHandle {}

#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct Seed(pub [u8; 16]);

impl Seed {
    pub fn random() -> Self {
        let mut raw = [0u8; 16];
        for i in 0..16 {
            raw[i] = rand::thread_rng().gen();
        }
        Self(raw)
    }

    pub fn default() -> Self {
        Self([0; 16])
    }

    pub fn expand_seed(&self) -> [Seed; 2] {
        let mut hasher = Shake128::default();
        hasher.update(&self.0);
        // println!("input seed: {:?}", self.0);
        let mut reader = hasher.finalize_xof();
        // we do not care about endianness
        // since we never use seed for directly computing.
        let mut res = [Seed::default(); 2];
        let dst: &mut [u8; 32] = unsafe { std::mem::transmute(&mut res) };
        reader.read(dst);
        // println!("expand seed: {:?}", dst);
        res
    }

    pub fn expand_u64(&self) -> [u64; 4] {
        let mut hasher = Shake128::default();
        hasher.update(&self.0);
        let mut reader = hasher.finalize_xof();
        // we do not care about endianness
        // since we never use seed for directly computing.
        let mut res = [0u64; 4];
        let dst: &mut [u8; 32] = unsafe { std::mem::transmute(&mut res) };
        reader.read(dst);
        res
    }
}

impl Size for Seed {
    const SIZE: usize = 16;
}

#[derive(Clone)]
pub struct SeedTree {
    pub nodes: Vec<Seed>,
    pub leaves: Vec<u64>,
    pub hashes: Vec<SHA256>,
    pub log_n_party: usize,
    pub log_n_iter: usize,
    pub tapes: Vec<__m512i>,
    round_keys: Vec<Vec<Vec<__m512i>>>,
}

impl SeedTree {
    pub fn new(mseed: &Seed, n_party: usize, n_iter: usize) -> Self {
        let log_n_party = ceil_log(n_party);
        let log_n_iter = ceil_log(n_iter);
        let mut res = Self {
            nodes: vec![Seed::default(); (4 << (log_n_party + log_n_iter)) - 1],
            leaves: vec![0; 8 << (log_n_party + log_n_iter)],
            hashes: vec![SHA256::from(0); 1 << (log_n_party + log_n_iter)],
            log_n_party,
            log_n_iter,
            tapes: vec![],
            round_keys: vec![],
        };

        // TODO: when n_iter is not a power of 2
        res.nodes.push(*mseed);
        // for i in 0..3 {
        //     let next_id = 2 * i + 1;
        //     let tmp = res.nodes[i].expand_seed();
        //     res.nodes[next_id] = tmp[0];
        //     res.nodes[next_id+1] = tmp[1];
        //     println!("left: {:?}", res.nodes[next_id]);
        //     println!("right: {:?}", res.nodes[next_id+1]);
        // }
        for i in 0..((2 << (log_n_party + log_n_iter) - 1) - 1) {
            let next_id = 2 * i + 1;
            let tmp = res.nodes[i].expand_seed();
            res.nodes[next_id] = tmp[0];
            res.nodes[next_id+1] = tmp[1];
            // &mut res.nodes[next_id..(next_id + 2)].try_into().unwrap());
        }

        let ofs = (1 << (log_n_party + log_n_iter)) - 1;
        for i in 0..(1 << (log_n_party + log_n_iter)) {
            let leaf_ofs = i * 8;
            let tmp = res.nodes[ofs + i].expand_u64();
            res.leaves[leaf_ofs] = tmp[0];
            res.leaves[leaf_ofs + 1] = tmp[1];
            res.leaves[leaf_ofs + 2] = tmp[2];
            res.leaves[leaf_ofs + 3] = tmp[3];
        }

        res
    }

    pub fn get_seed(&self, iid: usize, sid: usize) -> &Seed {
        let id = ((sid << self.log_n_iter) + iid) * 8;
        // we do not care about endianness
        // since we never use seed for directly computing.
        unsafe { std::mem::transmute(&self.leaves[id]) }
    }

    fn thread_aes(tapes: TapeHandle, ofs: usize, size: usize, tid: usize, log_n_iter: usize, log_n_party: usize, round_keys: &Vec<Vec<Vec<__m512i>>>) {
        unsafe {
            let tapes = tapes;
            let ofs = (ofs + size * tid) as u128;
            let mut state0 = M512 { u128_array: [ofs, ofs, ofs, ofs] }.m512i;
            let mut tapes = tapes.0.add((tid * size) << (log_n_iter + log_n_party - 2));
            for _ in 0..size {
                for i in 0..(1 << log_n_iter) {
                    for pp in 0..(1 << (log_n_party - 2)) {
                        let mut state = state0;
                        for r in 0..10 {
                            let round_key = round_keys[i][pp][r];
                            if r != 9 {
                                state = _mm512_aesenc_epi128(state, round_key);
                            } else {
                                state = _mm512_aesenclast_epi128(state, round_key);
                            }
                        }
                        *tapes = state;
                        tapes = tapes.add(1);
                        // tapes[i][k][pp] = state;
                    }
                }
                state0 = _mm512_add_epi64(
                    state0, M512 { u128_array: [1, 1, 1, 1] }.m512i);
            }
        }
    }

    pub fn gen_tapes(&mut self, ofs: usize, size: usize) {
        unsafe {
            let log_n_party = self.log_n_party;
            let log_n_iter = self.log_n_iter;
            let round_keys = &self.round_keys;
            let tapes = TapeHandle(std::mem::transmute(self.tapes.as_mut_ptr()));
            let n_thread = 8;
            let thread_size = size / n_thread;
            // let n_iter = (1 << self.log_n_iter) / n_thread;
            thread::scope(|s| {
                for tid in 0..n_thread {
                    s.spawn(move || {
                        Self::thread_aes(tapes, ofs, thread_size, tid, log_n_iter, log_n_party, round_keys)
                    });
                }
            });
        }
    }

    // pub fn gen_tapes(&mut self, ofs: usize, size: usize) {
    //     let mut ofs = 0;
    //     for i in 0..(1 << self.log_n_iter) {
    //         for pp in 0..(1 << (self.log_n_party - 2)) {
    //             for k in 0..size {
    //                 unsafe {
    //                     let c = (k + ofs) as u64;
    //                     let mut state = M512 { u64_array: [c, 0, c, 0, c, 0, c, 0] }.m512i;
    //                     for r in 0..10 {
    //                         let round_key = self.round_keys[i][pp][r];
    //                         if r != 9 {
    //                             state = _mm512_aesenc_epi128(state, round_key);
    //                         } else {
    //                             state = _mm512_aesenclast_epi128(state, round_key);
    //                         }
    //                     }
    //                     // self.tapes[i][pp].push(state);
    //                     // self.tapes[i][pp][k] = state;
    //                     self.tapes[ofs] = state;
    //                     ofs += 1;
    //                 }
    //             }
    //         }
    //     }
    // }

    pub fn gen_round_keys(&mut self, size: usize) {
        self.round_keys = vec![
            vec![Vec::with_capacity(10); 1 << (self.log_n_party - 2)]; 1 << self.log_n_iter];
        for i in 0..(1 << self.log_n_iter) {
            for pp in 0..(1 << (self.log_n_party - 2)) {
                unsafe {
                    let mut key = [M128 { u8_array: [0; 16] }.m128i; 4];
                    for p in 0..4 {
                        key[p] = M128 { u8_array: self.get_seed(i, pp * 4 + p).0 }.m128i;
                    }
                    self.round_keys[i][pp].push(M512 { m128i: key }.m512i);
                    for r in 0..9 {
                        for p in 0..4 {
                            let tmp = match r {
                                0 => _mm_aeskeygenassist_si128(key[p], 0x01),
                                1 => _mm_aeskeygenassist_si128(key[p], 0x02),
                                2 => _mm_aeskeygenassist_si128(key[p], 0x04),
                                3 => _mm_aeskeygenassist_si128(key[p], 0x08),
                                4 => _mm_aeskeygenassist_si128(key[p], 0x10),
                                5 => _mm_aeskeygenassist_si128(key[p], 0x20),
                                6 => _mm_aeskeygenassist_si128(key[p], 0x40),
                                7 => _mm_aeskeygenassist_si128(key[p], 0x80),
                                8 => _mm_aeskeygenassist_si128(key[p], 0x1b),
                                9 => _mm_aeskeygenassist_si128(key[p], 0x36),
                                _ => unreachable!(),
                            };
                            let tmp = _mm_shuffle_epi32(tmp, 0xff);
                            key[p] = _mm_xor_si128(key[p], _mm_slli_si128(key[p], 4));
                            key[p] = _mm_xor_si128(key[p], _mm_slli_si128(key[p], 4));
                            key[p] = _mm_xor_si128(key[p], _mm_slli_si128(key[p], 4));
                            key[p] = _mm_xor_si128(key[p], tmp);
                        }
                        self.round_keys[i][pp].push(M512 { m128i: key }.m512i);
                    }
                }
            }
        }
        unsafe {
            // self.tapes = vec![vec![vec![M512 { u64_array: [0; 8] }.m512i; size]; 1 << (self.log_n_party - 2)]; 1 << self.log_n_iter];
            // self.tapes = vec![M512 { u64_array: [0; 8] }.m512i; size << (self.log_n_party - 2 + self.log_n_iter)];
            let len = size << (self.log_n_party - 2 + self.log_n_iter);
            self.tapes = Vec::with_capacity(len);
            self.tapes.set_len(len);
        }
    }

    pub fn commit(&mut self) -> SHA256 {
        unsafe {
            let mut size = 1 << (self.log_n_party + self.log_n_iter);
            let mut tmp = vec![SHA256::from(0); size];
            let mut src_ptr: *mut u8 = std::mem::transmute(self.leaves.as_ptr());
            let dst_ptr: *mut u8 = std::mem::transmute(self.hashes.as_mut_ptr());
            {
                SHA256::raw_commit(dst_ptr, src_ptr as *const u8, size as u64);
                size >>= 1;
                src_ptr = dst_ptr;
            }

            let dst_ptr: *mut u8 = std::mem::transmute(tmp.as_mut_ptr());
            for _ in 0..(self.log_n_party + self.log_n_iter) {
                tmp.truncate(size);
                SHA256::raw_commit(dst_ptr, src_ptr as *const u8, size as u64);
                size >>= 1;
                src_ptr = dst_ptr;
            }
            tmp[0]
        }
    }

    pub fn open(&self, sto: &mut SeedTreeOpen) {
        let mut ids = Vec::with_capacity(sto.missings.len() + 1);
        for i in 0..sto.missings.len() {
            let id = (sto.missings[i] << self.log_n_iter) + i;
            let pos = ids.binary_search(&(id as isize)).unwrap_or_else(|e| e);
            ids.insert(pos, id as isize);
        }
        for &id in ids.iter() {
            sto.hashes[id as usize] = self.hashes[id as usize];
        }
        ids.push(-1);

        let mut path_len = 0;
        for d in (1..=(self.log_n_iter + self.log_n_party)).rev() {
            let mut i = 0;
            let mut j = 0;
            while ids[i] != -1 {
                if (ids[i+1] == ids[i] + 1) && (ids[i] % 2 == 0) {
                    ids[j] = ids[i] / 2;
                    i += 2;
                    j += 1;
                } else {
                    let new_id =
                        if ids[i] % 2 == 1 {
                            ((1 << d) - 1) + ids[i] as usize - 1
                        } else {
                            ((1 << d) - 1) + ids[i] as usize + 1
                        };
                    sto.nodes[new_id] = self.nodes[new_id];
                    path_len += 1;
                    ids[j] = ids[i] / 2;
                    i += 1;
                    j += 1;
                }
            }
            ids[j] = -1;
        }

        sto.proof_size = Seed::SIZE * path_len + SHA256::SIZE * sto.missings.len();
    }
}

#[derive(Clone)]
pub struct SeedTreeOpen {
    pub missings: Vec<usize>,
    pub nodes: Vec<Seed>,
    pub proof_size: usize,
    pub hashes: Vec<SHA256>,
    pub leaves: Vec<u64>,
    pub log_n_party: usize,
    pub log_n_iter: usize,
}

impl SeedTreeOpen {
    pub fn new(tree: &SeedTree) -> Self {
        Self {
            missings: vec![],
            nodes: vec![Seed::default(); (2 << (tree.log_n_party + tree.log_n_iter)) - 1],
            proof_size: 0,
            hashes: vec![SHA256::from(0); 1 << (tree.log_n_party + tree.log_n_iter)],
            leaves: vec![0; 8 << (tree.log_n_party + tree.log_n_iter)],
            log_n_party: tree.log_n_party,
            log_n_iter: tree.log_n_iter,
        }
    }

    pub fn get_seed(&self, iid: usize, sid: usize) -> &Seed {
        let id = ((sid << self.log_n_iter) + iid) * 8;
        // we do not care about endianness
        // since we never use seed for directly computing.
        unsafe { std::mem::transmute(&self.leaves[id]) }
    }

    pub fn verify(&mut self, root: &SHA256) -> bool {
        let mut ids = Vec::with_capacity(self.missings.len() + 1);
        for i in 0..self.missings.len() {
            let id = (self.missings[i] << self.log_n_iter) + i;
            let pos = ids.binary_search(&(id as isize)).unwrap_or_else(|e| e);
            ids.insert(pos, id as isize);
        }
        let ids_backup = ids.clone();
        ids.push(-1);

        for d in (1..=(self.log_n_iter + self.log_n_party)).rev() {
            let mut i = 0;
            let mut j = 0;
            while ids[i] != -1 {
                if (ids[i+1] == ids[i] + 1) && (ids[i] % 2 == 0) {
                    ids[j] = ids[i] / 2;
                    i += 2;
                    j += 1;
                } else {
                    let new_id =
                        if ids[i] % 2 == 1 {
                            ((1 << d) - 1) + ids[i] as usize - 1
                        } else {
                            ((1 << d) - 1) + ids[i] as usize + 1
                        };
                    // generate inner seeds:
                    for dd in (d+1)..(self.log_n_iter + self.log_n_party + 1) {
                        let start = ((new_id + 1) << (dd - d)) - 1;
                        let end = ((new_id + 2) << (dd - d)) - 1;
                        for k in (start..end).step_by(2) {
                            let tmp = self.nodes[(k - 1) / 2].expand_seed();
                            self.nodes[k] = tmp[0];
                            self.nodes[k + 1] = tmp[1];
                        }
                    }
                    ids[j] = ids[i] / 2;
                    i += 1;
                    j += 1;
                }
            }
            ids[j] = -1;
        }

        let ofs = (1 << (self.log_n_iter + self.log_n_party)) - 1;
        let mut ids_idx = 0;
        for i in 0..(1 << (self.log_n_iter + self.log_n_party)) {
            let leaves_ofs = i * 8;
            if ids_idx >= ids_backup.len() || i != ids_backup[ids_idx] as usize {
                let tmp = self.nodes[ofs + i].expand_u64();
                self.leaves[leaves_ofs + 0] = tmp[0];
                self.leaves[leaves_ofs + 1] = tmp[1];
                self.leaves[leaves_ofs + 2] = tmp[2];
                self.leaves[leaves_ofs + 3] = tmp[3];
                // &mut self.leaves[leaves_ofs..(leaves_ofs + 4)].try_into().unwrap());
            } else {
                ids_idx += 1;
            }
        }

        let h = unsafe {
            let mut size = 1 << (self.log_n_party + self.log_n_iter);
            let mut tmp = vec![SHA256::from(0); size];
            let mut src_ptr: *mut u8 = std::mem::transmute(self.leaves.as_ptr());
            let dst_ptr: *mut u8 = std::mem::transmute(tmp.as_mut_ptr());
            {
                SHA256::raw_commit(dst_ptr, src_ptr as *const u8, size as u64);
                size >>= 1;
                src_ptr = dst_ptr;
            }
            for &id in ids_backup.iter() {
                tmp[id as usize] = self.hashes[id as usize];
            }

            let dst_ptr: *mut u8 = std::mem::transmute(tmp.as_mut_ptr());
            for _ in 0..(self.log_n_party + self.log_n_iter) {
                tmp.truncate(size);
                SHA256::raw_commit(dst_ptr, src_ptr as *const u8, size as u64);
                size >>= 1;
                src_ptr = dst_ptr;
            }
            tmp[0]
        };

        h == *root
    }
}

#[cfg(test)]
mod tests {
    use rand::Rng;
    use super::*;

    #[test]
    fn test_seedtree() {
        let n_party = 16;
        let n_iter = 33;
        let mseed = Seed::random();
        let mut st = SeedTree::new(&mseed, n_party, n_iter);
        let h = st.commit();

        let mut sto = SeedTreeOpen::new(&st);
        sto.missings = Vec::with_capacity(n_iter);
        for _ in 0..n_iter {
            let r: usize = rand::thread_rng().gen();
            sto.missings.push(r % n_party);
        }

        // test commit, open and verify
        st.open(&mut sto);
        assert!(sto.verify(&h));

        // test get seed ref
        let mut id = 0;
        for p in 0..n_party {
            for i in 0..n_iter {
                let s0 = *sto.get_seed(i, p);
                let raw = [sto.leaves[i], sto.leaves[i + 1]];
                assert_eq!(s0, Seed(unsafe { std::mem::transmute(raw) }));
                id += 8;
            }
        }
    }
}


#[derive(Clone, Debug)]
pub struct Tape {
    pub cache: Vec<u8>,
    pub pointer: usize,
    // byte_len: usize,
}

impl Tape {
    pub fn new(seed: &Seed, byte_len: usize) -> Self {
        let mut hasher = Shake128::default();
        hasher.update(&seed.0);
        let mut reader = hasher.finalize_xof();
        let mut cache = vec![0u8; byte_len];
        reader.read(&mut cache);
        Self {
            cache: cache,
            pointer: 0,
            // byte_len: byte_len,
        }
    }

    // pub fn get_byte(&mut self) -> u8 {
    //     self.pointer += 1;
    //     self.cache[self.pointer - 1]
    // }

    // pub fn get_bytes(&mut self, n: usize) -> &[u8] {
    //     self.pointer += n;
    //     &self.cache[(self.pointer - n)..self.pointer]
    // }

    pub fn get<T: Field>(&mut self) -> T {
        let size = T::SIZE;
        let n_word = 1.max(size / 8);
        let mut tmp: Vec<u64> = Vec::with_capacity(n_word);
        for _ in 0..n_word {
            let mut word = 0u64;
            for _ in 0..size.min(8) {
                word <<= 8;
                word += self.cache[self.pointer] as u64;
                self.pointer += 1;
            }
            tmp.push(word);
        }
        T::from(tmp)
    }
}
