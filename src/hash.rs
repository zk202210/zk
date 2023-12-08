use sha2;
use sha2::Digest;

use super::field::*;

#[cfg(debug_assertions)]
use super::statistic::*;

#[link(name = "hashtree")]
extern {
    pub fn sha256_1_avx(output: *mut u8, input: *const u8, count: u64);
    pub fn sha256_shani(output: *mut u8, input: *const u8, count: u64);
    pub fn sha256_multiblk_shani(output: *mut u8, input: *const u8, count: u64, leaf_size: u64);
    pub fn sha256_hash_layer(output: *mut u8, input: *const u8, count: u64);
    pub fn sha256_hash_leaves(output: *mut u8, input: *const u8, count: u64, leaf_size: u64);
}


pub trait Commit<T> {
    fn commit(&mut self, preimages: &[T]);
}

pub trait CommitBlock {
    fn commit_block(&mut self, preimages: &Vec<u64>, leaf_size: u64);
}

pub trait RawCommit {
    fn raw_commit(output: *mut u8, input: *const u8, count: u64);
}

pub trait RawCommitMultiBlk {
    fn raw_commit_multiblk(output: *mut u8, input: *const u8, count: u64, leaf_size: u64);
}

pub trait CommitHash<H> {
    fn commit_hash(&mut self, preimages: &[H]);
}

pub trait Commit2<H> {
    fn commit2(&mut self, preimage0: &H, preimage1: &H);
}

pub trait Size {
    const SIZE: usize;
}

// pub trait Block:
//     Size + PartialEq + Clone + Copy + Size +
//     std::convert::From<GF2p16> +
//     std::convert::From<usize> {}

pub trait Hash<T: Field> :
    PartialEq + Clone + Copy + Size + std::fmt::Debug +
    std::convert::From<usize> +
    std::convert::AsRef<[u8]> +
    Size +
    Commit<T> +
    CommitBlock +
    RawCommit +
    RawCommitMultiBlk +
    CommitHash<Self> +
    Commit2<Self> {}


#[repr(C, align(8))]
#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Block512(pub [u64; 8]);

impl Size for Block512 {
    const SIZE: usize = 64;
}

impl std::convert::From<usize> for Block512 {
    fn from(x: usize) -> Self {
        Self([x as u64; 8])
    }
}

impl std::convert::From<GF2p16> for Block512 {
    fn from(x: GF2p16) -> Self {
        Self([x.0 as u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<GF2p32> for Block512 {
    fn from(x: GF2p32) -> Self {
        Self([x.0 as u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<GF2p64> for Block512 {
    fn from(x: GF2p64) -> Self {
        Self([x.0 as u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<GF2p128> for Block512 {
    fn from(x: GF2p128) -> Self {
        Self([x.0[0], x.0[1], 0u64, 0u64, 0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<GF2p192> for Block512 {
    fn from(x: GF2p192) -> Self {
        Self([x.0[0], x.0[1], x.0[2], 0u64, 0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<GF2p256> for Block512 {
    fn from(x: GF2p256) -> Self {
        Self([x.0[0], x.0[1], x.0[2], x.0[3], 0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<M61p2> for Block512 {
    fn from(x: M61p2) -> Self {
        Self([x.0[0], x.0[1], 0u64, 0u64, 0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<M127p2> for Block512 {
    fn from(x: M127p2) -> Self {
        Self([
             x.0[0] as u64,
             (x.0[0] >> 64) as u64,
             x.0[1] as u64,
             (x.0[1] >> 64) as u64,
             0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<P192> for Block512 {
    fn from(x: P192) -> Self {
        Self([x.0[0], x.0[1], x.0[2], 0u64, 0u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<P256> for Block512 {
    fn from(x: P256) -> Self {
        Self([x.0[0], x.0[1], x.0[2], x.0[3], 0u64, 0u64, 0u64, 0u64])
    }
}

// impl Block for Block512 {}


#[repr(C, align(8))]
#[derive(PartialEq, Clone, Copy, Debug)]
pub struct SHA256([u8; 32]);

impl<T: Field> Commit<T> for SHA256 {
    fn commit(&mut self, preimages: &[T]) {
        #[cfg(debug_assertions)] hash_call();
        let mut hasher = sha2::Sha256::new();
        for i in 0..preimages.len() {
            hasher.update(preimages[i]);
        }
        let hash = hasher.finalize();
        self.0 = <[u8; 32]>::from(hash);
    }
}

// use hashtree lib
impl CommitBlock for SHA256 {
    fn commit_block(&mut self, preimage: &Vec<u64>, leaf_size: u64) {
        // assert_eq! byte_len(preimages) < 64
        #[cfg(debug_assertions)] hash_call();
        unsafe {
            let src_ptr: *const u8 = std::mem::transmute(preimage.as_ptr());
            let dst_ptr: *mut u8 = std::mem::transmute(self.0.as_mut_ptr());
            sha256_multiblk_shani(dst_ptr, src_ptr, 1u64, leaf_size);
            // sha256_shani(dst_ptr, src_ptr, 1u64);
        }
    }
}

impl RawCommitMultiBlk for SHA256 {
    fn raw_commit_multiblk(output: *mut u8, input: *const u8, count: u64, leaf_size: u64) {
        // unsafe { sha256_hash_leaves(output, input, count, 64); }
        // assert!(leaf_size % 64 == 0);
        // assert!(count % 16 == 0);
        unsafe {
            sha256_hash_leaves(output, input, count, leaf_size);
            // sha256_multiblk_shani(output, input, count, leaf_size);
        }
    }
}

impl RawCommit for SHA256 {
    fn raw_commit(output: *mut u8, input: *const u8, count: u64) {
        // unsafe { sha256_hash_leaves(output, input, count, 64); }
        unsafe {
            sha256_hash_layer(output, input, count);
        }
    }
}

impl Size for SHA256 {
    const SIZE: usize = 32;
}

impl std::convert::From<usize> for SHA256 {
    fn from(x: usize) -> Self {
        Self([x as u8; 32])
    }
}

impl std::convert::AsRef<[u8]> for SHA256 {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl<H: AsRef<[u8]> + Copy> CommitHash<H> for SHA256 {
    fn commit_hash(&mut self, preimages: &[H]) {
        #[cfg(debug_assertions)] hash_call();
        let mut hasher = sha2::Sha256::new();
        for i in 0..preimages.len() {
            hasher.update(preimages[i]);
        }
        let hash = hasher.finalize();
        self.0 = <[u8; 32]>::from(hash);
    }
}

// use hashtree lib
impl<H: AsRef<[u8]> + Copy> Commit2<H> for SHA256 {
    fn commit2(&mut self, preimage0: &H, preimage1: &H) {
        #[cfg(debug_assertions)] hash_call();
        let mut src = preimage0.as_ref().to_vec();
        src.extend_from_slice(preimage1.as_ref());

        unsafe {
            let src_ptr: *const u8 = std::mem::transmute(src.as_ptr());
            let dst_ptr: *mut u8 = std::mem::transmute(self.0.as_mut_ptr());
            sha256_shani(dst_ptr, src_ptr, 1u64);
        }

        // let mut hasher = sha2::Sha256::new();
        // hasher.update(preimage0);
        // hasher.update(preimage1);
        // let hash = hasher.finalize();
        // self.0 = <[u8; 32]>::from(hash);
    }
}

impl<T: Field> Hash<T> for SHA256 {}
