use super::*;
use super::common::*;

#[derive(Copy, Clone)]
pub struct VGF2p32(pub __m256i, pub __m256i);

impl VGF2p32 {
    const MUL_LUT: &[u8] = include_bytes!(
        concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p32.mut.lut"));

    const MASKL: __m256i = unsafe {
        M256 { u8_array: [0x0fu8; 32] }.m256i
    };

    const MASKH: __m256i = unsafe {
        M256 { u8_array: [0xf0u8; 32] }.m256i
    };

    fn load_table(multiplier: u32, idx: usize, hl: usize) -> __m256i {
        unsafe {
            if hl == 0 {
                let multiplier = multiplier & 0xffff;
                _mm256_loadu_si256(std::mem::transmute(
                        &Self::MUL_LUT[(
                            (multiplier as usize) * 16 + idx) * 32]))
            } else {
                let multiplier = multiplier >> 16;
                _mm256_loadu_si256(std::mem::transmute(
                        &Self::MUL_LUT[(
                            (65536 + multiplier as usize) * 16 + idx) * 32]))

            }
        }
    }
}

impl Vect for VGF2p32 {
    type Input = GF2p32;
    type Output = Self;

    const DEGREE: usize = 16;

    fn new(raw: &[u8; Self::SIZE]) -> Self {
        unsafe {
            let mask = _mm256_set1_epi32(0xff);

            let tmp0 = _mm256_loadu_si256(std::mem::transmute(raw as *const [u8; Self::SIZE]));
            let tmp1 = _mm256_loadu_si256(std::mem::transmute(&raw[8..16] as *const _ as *const u8));
            let ll0 = _mm256_and_si256(tmp0, mask);
            let ll1 = _mm256_and_si256(tmp1, mask);

            let tmp0 = _mm256_srli_epi32(tmp0, 8);
            let tmp1 = _mm256_srli_epi32(tmp1, 8);
            let lh0 = _mm256_and_si256(tmp0, mask);
            let lh1 = _mm256_and_si256(tmp1, mask);

            let tmp0 = _mm256_srli_epi32(tmp0, 8);
            let tmp1 = _mm256_srli_epi32(tmp1, 8);
            let hl0 = _mm256_and_si256(tmp0, mask);
            let hl1 = _mm256_and_si256(tmp1, mask);

            let hh0 = _mm256_srli_epi32(tmp0, 8);
            let hh1 = _mm256_srli_epi32(tmp1, 8);

            let ll = _mm256_packus_epi32(ll0, ll1);
            let ll = _mm256_permute4x64_epi64(ll, 0b11011000);
            let lh = _mm256_packus_epi32(lh0, lh1);
            let lh = _mm256_permute4x64_epi64(lh, 0b11011000);
            let lo = _mm256_packus_epi16(ll, lh);
            let lo = _mm256_permute4x64_epi64(lo, 0b11011000);

            let hl = _mm256_packus_epi32(hl0, hl1);
            let hl = _mm256_permute4x64_epi64(hl, 0b11011000);
            let hh = _mm256_packus_epi32(hh0, hh1);
            let hh = _mm256_permute4x64_epi64(hh, 0b11011000);
            let hi = _mm256_packus_epi16(hl, hh);
            let hi = _mm256_permute4x64_epi64(hi, 0b11011000);

            Self(lo, hi)
        }
    }

    fn unpack(&self, dst: *mut [u8; Self::SIZE]) {
        unsafe {
            let ll = _mm256_permute4x64_epi64(self.0, 0b00010000);
            let lh = _mm256_permute4x64_epi64(self.0, 0b00110010);
            let lo = _mm256_unpacklo_epi8(ll, lh);

            let hl = _mm256_permute4x64_epi64(self.1, 0b00010000);
            let hh = _mm256_permute4x64_epi64(self.1, 0b00110010);
            let hi = _mm256_unpacklo_epi8(hl, hh);

            let lo0 = _mm256_permute4x64_epi64(lo, 0b00010000);
            let lo1 = _mm256_permute4x64_epi64(lo, 0b00110010);
            let hi0 = _mm256_permute4x64_epi64(hi, 0b00010000);
            let hi1 = _mm256_permute4x64_epi64(hi, 0b00110010);

            let res0 = _mm256_unpacklo_epi16(lo0, hi0);
            let res1 = _mm256_unpacklo_epi16(lo1, hi1);

            _mm256_storeu_si256(std::mem::transmute(dst as *mut [u8; Self::SIZE]), res0);
            let dst1 = &mut (*dst)[Self::SIZE / 2] as *mut _ as *mut u8;
            _mm256_storeu_si256(std::mem::transmute(dst1), res1);
            // _mm256_store_si256(std::mem::transmute(&mut dst[(Self::SIZE / 2)..Self::SIZE] as *mut _ as *mut u8), res1);
        }
    }
}

impl Add for VGF2p32 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
        unsafe {
            Self(
                _mm256_xor_si256(self.0, other.0),
                _mm256_xor_si256(self.1, other.1)
        )}
    }
}


impl AddAssign for VGF2p32 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
        unsafe {
            *self = Self(
                _mm256_xor_si256(self.0, other.0),
                _mm256_xor_si256(self.1, other.1)
        )}
    }
}

impl Sub for VGF2p32 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
        unsafe {
            Self(
                _mm256_xor_si256(self.0, other.0),
                _mm256_xor_si256(self.1, other.1)
        )}
    }
}

impl SubAssign for VGF2p32 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
        unsafe {
            *self = Self(
                _mm256_xor_si256(self.0, other.0),
                _mm256_xor_si256(self.1, other.1)
        )}
    }
}

impl Mul<GF2p32> for VGF2p32 {
    type Output = Self;

    fn mul(self, other: GF2p32) -> Self {
        self * (other.0 as usize)
    }
}

impl Mul<usize> for VGF2p32 {
    type Output = Self;

    fn mul(self, other: usize) -> Self {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
        unsafe {
            let res0 = {
                let k = 0;
                let table0 = Self::load_table(other as u32, 0, k);
                let table1 = Self::load_table(other as u32, 1, k);
                let table2 = Self::load_table(other as u32, 2, k);
                let table3 = Self::load_table(other as u32, 3, k);
                let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo));

                let table0 = Self::load_table(other as u32, 8, k);
                let table1 = Self::load_table(other as u32, 9, k);
                let table2 = Self::load_table(other as u32, 10, k);
                let table3 = Self::load_table(other as u32, 11, k);
                let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)));

                let k = 1;
                let table0 = Self::load_table(other as u32, 0, k);
                let table1 = Self::load_table(other as u32, 1, k);
                let table2 = Self::load_table(other as u32, 2, k);
                let table3 = Self::load_table(other as u32, 3, k);
                let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)));

                let table0 = Self::load_table(other as u32, 8, k);
                let table1 = Self::load_table(other as u32, 9, k);
                let table2 = Self::load_table(other as u32, 10, k);
                let table3 = Self::load_table(other as u32, 11, k);
                let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)))
            };

            let res1 = {
                let k = 0;
                let table0 = Self::load_table(other as u32, 4, k);
                let table1 = Self::load_table(other as u32, 5, k);
                let table2 = Self::load_table(other as u32, 6, k);
                let table3 = Self::load_table(other as u32, 7, k);
                let lo = _mm256_and_si256(self.1, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.1, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.1, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo));

                let table0 = Self::load_table(other as u32, 12, k);
                let table1 = Self::load_table(other as u32, 13, k);
                let table2 = Self::load_table(other as u32, 14, k);
                let table3 = Self::load_table(other as u32, 15, k);
                let lo = _mm256_and_si256(self.1, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.1, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.1, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)));

                let k = 1;
                let table0 = Self::load_table(other as u32, 4, k);
                let table1 = Self::load_table(other as u32, 5, k);
                let table2 = Self::load_table(other as u32, 6, k);
                let table3 = Self::load_table(other as u32, 7, k);
                let lo = _mm256_and_si256(self.1, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.1, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.1, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)));

                let table0 = Self::load_table(other as u32, 12, k);
                let table1 = Self::load_table(other as u32, 13, k);
                let table2 = Self::load_table(other as u32, 14, k);
                let table3 = Self::load_table(other as u32, 15, k);
                let lo = _mm256_and_si256(self.1, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.1, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.1, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)))
            };
            Self(res0, res1)
        }
    }
}

impl MulAssign<GF2p32> for VGF2p32 {
    fn mul_assign(&mut self, other: GF2p32) {
        *self *= other.0 as usize;
    }
}

impl MulAssign<usize> for VGF2p32 {
    fn mul_assign(&mut self, other: usize) {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
        unsafe {
            let res0 = {
                let k = 0;
                let table0 = Self::load_table(other as u32, 0, k);
                let table1 = Self::load_table(other as u32, 1, k);
                let table2 = Self::load_table(other as u32, 2, k);
                let table3 = Self::load_table(other as u32, 3, k);
                let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo));

                let table0 = Self::load_table(other as u32, 8, k);
                let table1 = Self::load_table(other as u32, 9, k);
                let table2 = Self::load_table(other as u32, 10, k);
                let table3 = Self::load_table(other as u32, 11, k);
                let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)));

                let k = 1;
                let table0 = Self::load_table(other as u32, 0, k);
                let table1 = Self::load_table(other as u32, 1, k);
                let table2 = Self::load_table(other as u32, 2, k);
                let table3 = Self::load_table(other as u32, 3, k);
                let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)));

                let table0 = Self::load_table(other as u32, 8, k);
                let table1 = Self::load_table(other as u32, 9, k);
                let table2 = Self::load_table(other as u32, 10, k);
                let table3 = Self::load_table(other as u32, 11, k);
                let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)))
            };

            let res1 = {
                let k = 0;
                let table0 = Self::load_table(other as u32, 4, k);
                let table1 = Self::load_table(other as u32, 5, k);
                let table2 = Self::load_table(other as u32, 6, k);
                let table3 = Self::load_table(other as u32, 7, k);
                let lo = _mm256_and_si256(self.1, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.1, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.1, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo));

                let table0 = Self::load_table(other as u32, 12, k);
                let table1 = Self::load_table(other as u32, 13, k);
                let table2 = Self::load_table(other as u32, 14, k);
                let table3 = Self::load_table(other as u32, 15, k);
                let lo = _mm256_and_si256(self.1, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.1, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.1, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)));

                let k = 1;
                let table0 = Self::load_table(other as u32, 4, k);
                let table1 = Self::load_table(other as u32, 5, k);
                let table2 = Self::load_table(other as u32, 6, k);
                let table3 = Self::load_table(other as u32, 7, k);
                let lo = _mm256_and_si256(self.1, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.1, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.1, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                let res = _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)));

                let table0 = Self::load_table(other as u32, 12, k);
                let table1 = Self::load_table(other as u32, 13, k);
                let table2 = Self::load_table(other as u32, 14, k);
                let table3 = Self::load_table(other as u32, 15, k);
                let lo = _mm256_and_si256(self.1, Self::MASKL);   // a0, a2
                let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(self.1, Self::MASKH), 4);    // a1, a3
                let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi
                let tmp = _mm256_xor_si256(hi, lo);
                let tmp0 = _mm256_permute4x64_epi64(
                    self.1, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1
                let lo = _mm256_and_si256(tmp0, Self::MASKL);     // a2, a0
                let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
                let hi = _mm256_srli_epi64(
                    _mm256_and_si256(tmp0, Self::MASKH), 4);      // a3, a1
                let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi
                _mm256_xor_si256(res, _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)))
            };
            *self = Self(res0, res1);
        }
    }
}

impl VField<GF2p32> for VGF2p32 {}
