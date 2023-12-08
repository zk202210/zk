use rug::Rational;

use std::marker::PhantomData;
use crate::field::Field;
use crate::hash::*;

#[repr(C, align(8))]
pub struct MerkleTree<T: Field, H: Hash<T>> {
    pub size: usize,
    pub leaf_degree: usize,
    pub slice_degree: usize,
    pub leaf_word_size: usize,
    pub nodes: Vec<H>,
    _leaves: Vec<u64>,
    field_type: PhantomData<T>,
}

impl<T: Field, H: Hash<T>> MerkleTree<T, H> {
    pub fn new(leaf_degree: usize, slice_degree: usize) -> Self {
        let size = 1 << leaf_degree;
        // TODO XXX for smaller field
        // let field_size = 1.max(T::SIZE / 8);
        let field_size = 1.max(T::SIZE / 8);
        let mut m = (field_size << slice_degree) / 8;
        if (field_size << slice_degree) % 8 > 0 {
            m += 1;
        }
        m *= 8;
        let nodes = vec![H::from(0); 2 * size - 1];
        Self {
            size,
            leaf_degree,
            slice_degree,
            leaf_word_size: m,
            nodes,
            _leaves: vec![0; m * size],
            field_type: PhantomData,
        }
    }

    pub fn clean(&mut self) {
        self.nodes = vec![];
        self._leaves = vec![];
    }

    pub fn set_leaf(&mut self, idx: usize, leaf: &Vec<T>) {
        // TODO XXX for smaller field
        let mut ofs = idx * self.leaf_word_size;
        for i in 0..leaf.len() {
            leaf[i].write_word(&mut self._leaves, &mut ofs);
        }
    }

    pub fn set_leaf_from_slice(&mut self, idx: usize, leaf: &[T], stride: usize, len: usize) {
        // TODO XXX for smaller field
        let mut ofs = idx * self.leaf_word_size;
        for i in (0..(len * stride)).step_by(stride) {
            leaf[i].write_word(&mut self._leaves, &mut ofs);
        }
    }

    pub fn commit(&mut self) {
        // TODO: slice_degree != 0
        // for i in 0..self.size {
        //     let tmp: &[u8] = self.leaves[i].as_ref();
        //     for j in 0..tmp.len() {
        //         self._leaves[i * 64 + j] = tmp[j];
        //     }
        // }
        unsafe {
            let mut src_ptr: *const u8 = std::mem::transmute(self._leaves.as_ptr());
            let dst_ptr: *mut u8 = std::mem::transmute(self.nodes.as_mut_ptr());

            let mut size = self.size as u64;

            {
                let leaf_size = (self.leaf_word_size as u64) * 8;
                let offset = (1 << self.leaf_degree) - 1;
                let layer_ptr = dst_ptr.add(offset * 32);
                H::raw_commit_multiblk(layer_ptr, src_ptr, size, leaf_size);

                src_ptr = layer_ptr as *const u8;
                size /= 2;

            }
            for d in (0..self.leaf_degree).rev() {
                let offset = (1 << d) - 1;
                let layer_ptr = dst_ptr.add(offset * 32);
                H::raw_commit(layer_ptr, src_ptr, size);

                src_ptr = layer_ptr as *const u8;
                size /= 2;
            }
        }

        // for i in 0..(self.size >> self.slice_degree) {
        //     self.nodes[(self.size >> self.slice_degree) - 1 + i].commit(
        //         &self.leaves[
        //         (i << self.slice_degree)..
        //         ((i + 1) << self.slice_degree)]);
        // }
        // for d in (0..(self.leaf_degree - self.slice_degree)).rev() {
        //     for i in ((1 << d) - 1)..((1 << (d + 1)) - 1) {
        //         let mut tmp = H::from(0);
        //         tmp.commit_hash(&self.nodes[(2 * i + 1)..(2 * i + 3)]);
        //         self.nodes[i] = tmp;
        //     }
        // }
    }

    pub fn open(&self, mto: &mut MerkleTreeOpen<T, H>) {
        // ids will be sorted and duplication removed
        let mut ids = Vec::with_capacity(mto.size + 1);
        for i in 0..mto.size {
            let ofs = mto.ids[i] * self.leaf_word_size;
            mto._values[i] = self._leaves[ofs..(ofs + self.leaf_word_size)].to_vec();

            let pos = ids.binary_search(&(mto.ids[i] as isize));
            match pos {
                Ok(_) => {},
                Err(pos) => {
                    ids.insert(pos, mto.ids[i] as isize);
                },
            }
        }
        let size = ids.len();
        ids.push(-1);
        let mut path_id = 0;
        // for d in (1..=(self.leaf_degree - self.slice_degree)).rev() {
        for d in (1..=self.leaf_degree).rev() {
            let mut i = 0;
            let mut j = 0;
            while ids[i] != -1 {
                if (ids[i+1] == ids[i] + 1) && (ids[i] % 2 == 0) {
                    ids[j] = ids[i] / 2;
                    i += 2;
                    j += 1;
                } else {
                    if ids[i] % 2 == 1 {
                        mto.path.push(self.nodes[((1 << d) - 1) + ids[i] as usize - 1]);
                    } else {
                        mto.path.push(self.nodes[((1 << d) - 1) + ids[i] as usize + 1]);
                    }
                    ids[j] = ids[i] / 2;
                    path_id += 1;
                    i += 1;
                    j += 1;
                }
            }
            ids[j] = -1;
        }
        mto.path_size = path_id;

        // assert!(mto.size < (1 << self.leaf_degree));
        // do not count duplicates
        mto.proof_size = (size << self.slice_degree) * T::SIZE + mto.path_size * H::SIZE;
    }
}

#[derive(Clone)]
pub struct MerkleTreeOpen<T, H> {
    pub leaf_degree: usize,
    pub leaf_word_size: usize,
    pub size: usize,
    pub slice_degree: usize,
    pub path_size: usize,
    pub ids: Vec<usize>,
    pub _values: Vec<Vec<u64>>,
    pub path: Vec<H>,
    pub proof_size: usize,
    field_type: PhantomData<T>,
}

impl<T: Field, H: Hash<T>> MerkleTreeOpen<T, H> {
    pub fn new(tree: &MerkleTree<T, H>, open_size: usize) -> Self {
        Self {
            leaf_degree: tree.leaf_degree,
            leaf_word_size: tree.leaf_word_size,
            size: open_size,
            slice_degree: tree.slice_degree,
            path_size: 0,
            ids: Vec::new(),
            _values: vec![vec![]; open_size << tree.slice_degree],
            path: vec![],
            proof_size: 0,
            field_type: PhantomData,
        }
    }

    pub fn get_value(&self, idx: usize) -> Vec<T> {
        let mut res = vec![];
        let mut ofs = 0;
        // TODO XXX: small field
        let step = 1.max(T::SIZE / 8);
        while ofs + step <= self.leaf_word_size {
            res.push(T::from(
                    self._values[idx]
                    [ofs..(ofs + step)].to_vec()));
            ofs += step;
        }
        res
    }

    pub fn get_avg_proof_size(&self, log_l: usize) -> usize {

        // avg_hash[i][j] is avg number of revealed inner nodes for a 2ⁱ⁺¹ size tree opening j points
        let mut avg_hash = vec![vec![]; self.leaf_degree];
        let mut avg_leaf = vec![vec![]; self.leaf_degree];
        let mut c = vec![vec![]; self.leaf_degree];
        let log_leaf_size = log_l - self.leaf_degree;
        let size = self.size;

        let mut c_1 = vec![Rational::from(0); size.min(1 << log_leaf_size) + 1];
        for i in 0..c_1.len() {
            c_1[i] = count_combinations(1 << log_leaf_size, i);
        }
        avg_hash[0] = vec![Rational::from(0); size.min(2 << log_leaf_size) + 1];
        avg_hash[0][0] = Rational::from(1);
        avg_hash[0][1] = Rational::from(1);
        avg_leaf[0] = vec![Rational::from(0); size.min(2 << log_leaf_size) + 1];
        avg_leaf[0][0] = Rational::from(0);
        avg_leaf[0][1] = Rational::from(1);
        c[0] = vec![Rational::from(0); size.min(2 << log_leaf_size) + 1];
        c[0][1] = Rational::from(2 << log_leaf_size);
        for j in 2..avg_hash[0].len() {
            let mut sum_c = Rational::from(0);
            c[0][j] = count_combinations(2 << log_leaf_size, j);
            for k in 0..=j.min(1 << log_leaf_size) {
                if j - k > size.min(1 << log_leaf_size) {
                    continue;
                }
                sum_c += c_1[k].clone() * c_1[j - k].clone();
            }
            let mut tmp_hash = Rational::from(0);
            let mut tmp_leaf = Rational::from(0);
            for k in 0..=j.min(1 << log_leaf_size) {
                if j - k > (1 << log_leaf_size) {
                    continue;
                }
                if k == 0 || j - k == 0 {
                    tmp_hash += c_1[k].clone() * c_1[j - k].clone() / sum_c.clone();
                    tmp_leaf += c_1[k].clone() * c_1[j - k].clone() / sum_c.clone();
                } else {
                    tmp_leaf += Rational::from(2) * c_1[k].clone() * c_1[j - k].clone() / sum_c.clone();
                }
            }
            avg_hash[0][j] = tmp_hash;
            avg_leaf[0][j] = tmp_leaf;
        }

        for i in 1..self.leaf_degree {
            avg_hash[i] = vec![Rational::from(0); size.min(2 << (i + log_leaf_size)) + 1];
            avg_leaf[i] = vec![Rational::from(0); size.min(2 << (i + log_leaf_size)) + 1];
            c[i] = vec![Rational::from(0); size.min(2 << (i + log_leaf_size)) + 1];
            avg_hash[i][0] = Rational::from(1);
            avg_hash[i][1] = Rational::from(i + 1);
            avg_leaf[i][0] = Rational::from(0);
            avg_leaf[i][1] = Rational::from(1);
            c[i][0] = Rational::from(1);
            c[i][1] = Rational::from(2 << (i + log_leaf_size));

            for j in 2..avg_hash[i].len() {
                let mut sum_c = Rational::from(0);
                c[i][j] = count_combinations(2 << (i + log_leaf_size), j);
                for k in 0..=j.min(size.min(1 << (i + log_leaf_size))) {
                    if j - k > size.min(1 << (i + log_leaf_size)) {
                        continue;
                    }
                    sum_c += c[i - 1][k].clone() * c[i - 1][j - k].clone();
                }

                let mut tmp_hash = Rational::from(0);
                let mut tmp_leaf = Rational::from(0);
                for k in 0..=j.min(size.min(1 << (i + log_leaf_size))) {
                    if j - k > size.min(1 << (i + log_leaf_size)) {
                        continue;
                    }
                    tmp_hash += (avg_hash[i - 1][k].clone() + avg_hash[i - 1][j - k].clone()) * c[i - 1][k].clone() * c[i - 1][j - k].clone() / sum_c.clone();
                    tmp_leaf += (avg_leaf[i - 1][k].clone() + avg_leaf[i - 1][j - k].clone()) * c[i - 1][k].clone() * c[i - 1][j - k].clone() / sum_c.clone();
                }
                avg_hash[i][j] = tmp_hash;
                avg_leaf[i][j] = tmp_leaf;
            }
        }
        // println!("hash {} avg {:.2} vs real {}", 1 << self.leaf_degree, avg_hash[self.leaf_degree - 1][size].to_f64(), self.path.len());
        // println!("leaf {} avg {:.2} vs real {}", 1 << self.leaf_degree, avg_leaf[self.leaf_degree - 1][size].to_f64(), self.size);

        (((T::SIZE << self.slice_degree) as f64) * avg_leaf[self.leaf_degree - 1][size].to_f64() +
         (H::SIZE as f64) * avg_hash[self.leaf_degree - 1][size].to_f64()) as usize
    }

    pub fn verify(&self, root: &H) -> bool {
        let mut hash = Vec::with_capacity(self.size);
        let mut tmp = Vec::with_capacity(self.size + 1);
        let leaf_size = (self.leaf_word_size as u64) * 8;
        for i in 0..self.size {
            let pos = tmp.binary_search(&(self.ids[i] as isize));
            match pos {
                Ok(_) => {},
                Err(pos) => {
                    tmp.insert(pos, self.ids[i] as isize);

                    let mut h = H::from(0);
                    h.commit_block(&self._values[i], leaf_size);
                    hash.insert(pos, h);
                },
            }
        }
        tmp.push(-1);

        let mut path_id = 0usize;
        // for _ in (1..=(self.leaf_degree - self.slice_degree)).rev() {
        for _ in (1..=self.leaf_degree).rev() {
            let mut i = 0;
            let mut j = 0;
            while tmp[i] != -1 {
                if (tmp[i+1] == tmp[i] + 1) && (tmp[i] % 2 == 0) {
                    tmp[j] = tmp[i] / 2;
                    let mut h = H::from(0);
                    h.commit2(&hash[i], &hash[i+1]);
                    hash[j] = h;
                    i += 2;
                    j += 1;
                } else {
                    if tmp[i] % 2 == 1 {
                        let mut h = H::from(0);
                        h.commit2(&self.path[path_id], &hash[i]);
                        hash[j] = h;
                    } else {
                        let mut h = H::from(0);
                        h.commit2(&hash[i], &self.path[path_id]);
                        hash[j] = h;
                    }
                    tmp[j] = tmp[i] / 2;
                    path_id += 1;
                    i += 1;
                    j += 1;
                }
            }
            tmp[j] = -1;
        }
        hash[0] == *root
    }
}

fn count_combinations(n: usize, r: usize) -> Rational {
    if r > n {
        Rational::from(0)
    } else {
        (1..=r.min(n - r))
            .map(|x| Rational::from(x))
            .fold(Rational::from(1), |acc, val| acc.clone() * (Rational::from(n) - val.clone() + Rational::from(1)) / val.clone())
    }
}

// fn count_combinations(n: usize, r: usize) -> f64 {
//     if r > n {
//         0.0
//     } else {
//         (1..=r.min(n - r))
//             .map(|x| x as f64)
//             .fold(1.0, |acc, val| acc * ((n as f64) - val + 1.0) / val)
//     }
// }

pub type MT<T> = MerkleTree<T, SHA256>;
pub type MTO<T> = MerkleTreeOpen<T, SHA256>;
