use super::*;
use super::common::*;


// GF2p16 of pack 16
// ref: http://web.eecs.utk.edu/~jplank/plank/papers/FAST-2013-GF.pdf
//
// for a: u16 = a0 + (a1 << 4) + (a2 << 8) + (a3 << 12)
//      where a0, a1, a2, a3: u4
// MUL_LUT[y][0x00 + a0] = lo(y * (a0 << 0x0))
// MUL_LUT[y][0x10 + a2] = hi(y * (a2 << 0x8))
// MUL_LUT[y][0x20 + a2] = lo(y * (a2 << 0x8))
// MUL_LUT[y][0x30 + a0] = hi(y * (a0 << 0x0))
// MUL_LUT[y][0x40 + a1] = lo(y * (a1 << 0x4))
// MUL_LUT[y][0x50 + a3] = hi(y * (a3 << 0xc))
// MUL_LUT[y][0x60 + a3] = lo(y * (a3 << 0xc))
// MUL_LUT[y][0x70 + a1] = hi(y * (a1 << 0x4))

#[derive(Copy, Clone)]
pub struct VGF2p16(pub __m256i);

impl VGF2p16 {
    const MUL_LUT: &[u8] = include_bytes!(
        concat!(env!("CARGO_MANIFEST_DIR"), "/resources/packgf2p16.mut.lut"));

    const MASKL: __m256i = unsafe {
        M256 { u8_array: [0x0fu8; 32] }.m256i
    };

    const MASKH: __m256i = unsafe {
        M256 { u8_array: [0xf0u8; 32] }.m256i
    };

    fn load_table(multiplier: u16, idx: usize) -> __m256i {
        unsafe {
            _mm256_loadu_si256(std::mem::transmute(
                    &Self::MUL_LUT[(multiplier as usize * 4 + idx) * 32]))
        }
    }
}


impl Vect for VGF2p16 {
    type Input = GF2p16;
    type Output = Self;

    const DEGREE: usize = 16;

    fn new(raw: &[u8; Self::SIZE]) -> Self {
        Self(unsafe {
            let tmp = _mm256_loadu_si256(std::mem::transmute(raw as *const [u8; Self::SIZE]));
            let lo = _mm256_and_si256(tmp, _mm256_set1_epi16(0xff));
            let hi = _mm256_srli_epi16(tmp, 8);

            // 0lo 0 1lo 0 2lo 0 3lo 0 4lo 0 5lo 0 6lo 0 7lo 0
            // 8lo 0 9lo 0 alo 0 blo 0 clo 0 dlo 0 elo 0 flo 0
            // pack
            // 0hi 0 1hi 0 2hi 0 3hi 0 4hi 0 5hi 0 6hi 0 7hi 0
            // 8hi 0 9hi 0 ahi 0 bhi 0 chi 0 dhi 0 ehi 0 fhi 0
            // ->
            // 0lo 1lo 2lo 3lo 4lo 5lo 6lo 7lo
            // 0hi 1hi 2hi 3hi 4hi 5hi 6hi 7hi
            // 8lo 9lo alo blo clo dlo elo flo
            // 8hi 9hi ahi bhi chi dhi ehi fhi
            let tmp = _mm256_packus_epi16(lo, hi);

            _mm256_permute4x64_epi64(tmp, 0b11011000)
        })
    }

    fn unpack(&self, dst: *mut [u8; Self::SIZE]) {
        unsafe {
            // 0 1 2 3 4 5 6 7 lo, 8 9 a b c d e f lo,
            // 0 1 2 3 4 5 6 7 hi, 8 9 a b c d e f hi,
            // -->
            // 0 1 2 3 4 5 6 7 lo, x x x x x x x x
            // 8 9 a b c d e f lo, x x x x x x x x
            // &&
            // 0 1 2 3 4 5 6 7 hi, x x x x x x x x
            // 8 9 a b c d e f hi, x x x x x x x x
            let lo = _mm256_permute4x64_epi64(self.0, 0b00010000);
            let hi = _mm256_permute4x64_epi64(self.0, 0b00110010);

            // 0 1 2 3 4 5 6 7 lo, x x x x x x x x
            // 8 9 a b c d e f lo, x x x x x x x x
            // pack
            // 0 1 2 3 4 5 6 7 hi, x x x x x x x x
            // 8 9 a b c d e f hi, x x x x x x x x
            // ->
            // 0lo 0hi 1lo 1hi ...
            let res = _mm256_unpacklo_epi8(lo, hi);
            _mm256_storeu_si256(std::mem::transmute(dst as *mut u8), res);
        }
    }
}

impl Add for VGF2p16 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
        Self(unsafe {
            _mm256_xor_si256(self.0, other.0)
        })
    }
}


impl AddAssign for VGF2p16 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
        *self = Self(unsafe {
            _mm256_xor_si256(self.0, other.0)
        })
    }
}

impl Sub for VGF2p16 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
        Self(unsafe {
            _mm256_xor_si256(self.0, other.0)
        })
    }
}


impl SubAssign for VGF2p16 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
        *self = Self(unsafe {
            _mm256_xor_si256(self.0, other.0)
        })
    }
}


impl Mul<GF2p16> for VGF2p16 {
    type Output = Self;

    fn mul(self, other: GF2p16) -> Self {
        self * (other.0 as usize)
    }
}

impl Mul<usize> for VGF2p16 {
    type Output = Self;

    fn mul(self, other: usize) -> Self {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
        unsafe {
            let table0 = Self::load_table(other as u16, 0);
            let table1 = Self::load_table(other as u16, 2);
            let table2 = Self::load_table(other as u16, 1);
            let table3 = Self::load_table(other as u16, 3);

            // let res = _mm256_permute4x64_epi64(self.0, 0b01_00_11_10);  // a2a3; a0a1

            // let a = _mm512_inserti64x4(_mm512_castsi256_si512(self.0), res, 1);
            // let b = _mm512_inserti64x4(_mm512_castsi256_si512(self.0), res, 1);

            // let mut a = _mm512_and_si512(a, Self::MASKL); // a0, a2, a2, a0
            // let b = _mm512_srli_epi64(_mm512_and_si512(b, Self::MASKH), 4); // a1, a3, a3, a1

            // // let a = _mm512_shuffle_epi8(table0, a); // a0lo, a2hi, a2lo, a0hi
            // std::arch::asm!(
            //     "vpshufb zmm1, zmm1, zmm2",
            //     inout("zmm1") a,
            //     in("zmm2") table0,
            // );
            // let b = _mm512_shuffle_epi8(table1, b); // a1lo, a3hi, a3lo, a1hi

            // let res = _mm512_xor_si512(a, b);
            // let a = _mm512_castsi512_si256(res);
            // let b = _mm512_extracti64x4_epi64(res, 1);

            // Self(_mm256_xor_si256(a, b))

            let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
            let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
            let hi = _mm256_srli_epi64(
                _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
            let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi

            let tmp = _mm256_xor_si256(hi, lo);
            let res = _mm256_permute4x64_epi64(
                self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1

            let lo = _mm256_and_si256(res, Self::MASKL);      // a2, a0
            let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
            let hi = _mm256_srli_epi64(
                _mm256_and_si256(res, Self::MASKH), 4);       // a3, a1
            let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi

            Self(_mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo)))
        }
    }
}

impl MulAssign<GF2p16> for VGF2p16 {
    fn mul_assign(&mut self, other: GF2p16) {
        *self *= other.0 as usize;
    }
}

impl MulAssign<usize> for VGF2p16 {
    fn mul_assign(&mut self, other: usize) {
        #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
        unsafe {
            // let table0 = Self::load_table(other as u16, 0);
            // let table1 = Self::load_table(other as u16, 2);
            // let table2 = Self::load_table(other as u16, 1);
            // let table3 = Self::load_table(other as u16, 3);

            // let res = _mm256_permute4x64_epi64(self.0, 0b01_00_11_10);  // a2a3; a0a1

            // let a = _mm512_inserti64x4(_mm512_castsi256_si512(self.0), res, 1);
            // let b = _mm512_inserti64x4(_mm512_castsi256_si512(self.0), res, 1);

            // let a = _mm512_and_si512(a, Self::MASKL); // a0, a2, a2, a0
            // let b = _mm512_srli_epi64(_mm512_and_si512(b, Self::MASKH), 4); // a1, a3, a3, a1

            // let a = _mm512_shuffle_epi8(table0, a); // a0lo, a2hi, a2lo, a0hi
            // let b = _mm512_shuffle_epi8(table1, b); // a1lo, a3hi, a3lo, a1hi

            // let res = _mm512_xor_si512(a, b);
            // let a = _mm512_castsi512_si256(res);
            // let b = _mm512_extracti64x4_epi64(res, 1);

            // self.0 = _mm256_xor_si256(a, b);
            let table0 = Self::load_table(other as u16, 0);
            let table2 = Self::load_table(other as u16, 1);
            let table1 = Self::load_table(other as u16, 2);
            let table3 = Self::load_table(other as u16, 3);

            let lo = _mm256_and_si256(self.0, Self::MASKL);   // a0, a2
            let lo = _mm256_shuffle_epi8(table0, lo);   // a0lo, a2hi
            let hi = _mm256_srli_epi64(
                _mm256_and_si256(self.0, Self::MASKH), 4);    // a1, a3
            let hi = _mm256_shuffle_epi8(table1, hi);   // a1lo, a3hi

            let tmp = _mm256_xor_si256(hi, lo);
            self.0 = _mm256_permute4x64_epi64(
                self.0, 2 | (3 << 2) | (1 << 6));       // a2,a3; a0,a1

            let lo = _mm256_and_si256(self.0, Self::MASKL);   // a2, a0
            let lo = _mm256_shuffle_epi8(table2, lo);   // a2lo, a0hi
            let hi = _mm256_srli_epi64(
                _mm256_and_si256(self.0, Self::MASKH), 4);    // a3, a1
            let hi = _mm256_shuffle_epi8(table3, hi);   // a3lo, a1hi

            self.0 = _mm256_xor_si256(tmp, _mm256_xor_si256(hi, lo));
        }
    }
}


impl std::fmt::Display for VGF2p16 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        unsafe {
        write!(f, "{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}",
            _mm256_extract_epi8(self.0, 0),
            _mm256_extract_epi8(self.0, 1),
            _mm256_extract_epi8(self.0, 2),
            _mm256_extract_epi8(self.0, 3),
            _mm256_extract_epi8(self.0, 4),
            _mm256_extract_epi8(self.0, 5),
            _mm256_extract_epi8(self.0, 6),
            _mm256_extract_epi8(self.0, 7),
            _mm256_extract_epi8(self.0, 8),
            _mm256_extract_epi8(self.0, 9),
            _mm256_extract_epi8(self.0, 10),
            _mm256_extract_epi8(self.0, 11),
            _mm256_extract_epi8(self.0, 12),
            _mm256_extract_epi8(self.0, 13),
            _mm256_extract_epi8(self.0, 14),
            _mm256_extract_epi8(self.0, 15),
            _mm256_extract_epi8(self.0, 16),
            _mm256_extract_epi8(self.0, 17),
            _mm256_extract_epi8(self.0, 18),
            _mm256_extract_epi8(self.0, 19),
            _mm256_extract_epi8(self.0, 20),
            _mm256_extract_epi8(self.0, 21),
            _mm256_extract_epi8(self.0, 22),
            _mm256_extract_epi8(self.0, 23),
            _mm256_extract_epi8(self.0, 24),
            _mm256_extract_epi8(self.0, 25),
            _mm256_extract_epi8(self.0, 26),
            _mm256_extract_epi8(self.0, 27),
            _mm256_extract_epi8(self.0, 28),
            _mm256_extract_epi8(self.0, 29),
            _mm256_extract_epi8(self.0, 30),
            _mm256_extract_epi8(self.0, 31))
        }
    }
}

impl VField<GF2p16> for VGF2p16 {}
