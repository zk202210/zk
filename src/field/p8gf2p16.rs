use super::*;
use super::common::*;


#[derive(Copy, Clone, Debug)]
pub struct P8GF2p16(pub __m256i);

impl P8GF2p16 {
    const MASK: __m256i = unsafe {
        M256 { u64_array: [0xffffffffu64; 4] }.m256i
    };

    const MASK2: __m256i = unsafe {
        M256 { u32_array: [0xffffu32; 8] }.m256i
    };

    const PERMUTE: __m256i = unsafe {
        M256 { u16_array: [
            0x00, 0x01, 0x02, 0x03,
            0x08, 0x09, 0x0a, 0x0b,
            0x10, 0x11, 0x12, 0x13,
            0x18, 0x19, 0x1a, 0x1b,
        ] }.m256i
    };

    const ZERO: __m256i = unsafe {
        M256 { u8_array: [0; 32] }.m256i
    };

    const FULLMODULUS: __m256i = unsafe {
        M256 { u64_array: [
            0b10000000000101011u64, 0u64,
            0b10000000000101011u64, 0u64] }.m256i
    };

    const MODULUS: __m256i = unsafe {
        M256 { u64_array: [
            0b101011u64, 0u64,
            0b101011u64, 0u64] }.m256i
    };
}

impl Pack<8> for P8GF2p16 {
    type Output = Self;
    type Input = GF2p16;

    fn pack(raw: &[GF2p16; 8]) -> Self {
        unsafe {
            let ptr: *const __m128i = std::mem::transmute(raw);
            let tmp = _mm256_castsi128_si256(_mm_loadu_si128(ptr));
            Self(_mm256_unpacklo_epi16(
                    _mm256_permute4x64_epi64(tmp, 0b11_01_10_00),
                    Self::ZERO))
        }
    }

    fn unpack(&self, dst: &mut [GF2p16; 8]) {
        unsafe {
            let tmp = _mm256_packus_epi32(self.0, self.0);
            let tmp = _mm256_castsi256_si128(_mm256_permutexvar_epi16(Self::PERMUTE, tmp));
            _mm_storeu_si128(std::mem::transmute(dst as *mut _ as *mut u8), tmp);
        }
    }
}

impl Add for P8GF2p16 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..8 { add_call(); }
        unsafe {
            Self(_mm256_xor_si256(self.0, other.0))
        }
    }
}

impl AddAssign for P8GF2p16 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..8 { add_call(); }
        unsafe {
            self.0 = _mm256_xor_si256(self.0, other.0);
        }
    }
}

impl Sub for P8GF2p16 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..8 { add_call(); }
        unsafe {
            Self(_mm256_xor_si256(self.0, other.0))
        }
    }
}

impl SubAssign for P8GF2p16 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..8 { add_call(); }
        unsafe {
            self.0 = _mm256_xor_si256(self.0, other.0);
        }
    }
}

impl Mul for P8GF2p16 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..8 { mul_call(); }
        unsafe {
            let lo = _mm256_and_si256(_mm256_clmulepi64_epi128(self.0, other.0, 0x00), Self::MASK);
            let hi = _mm256_and_si256(_mm256_clmulepi64_epi128(self.0, other.0, 0x11), Self::MASK);

            let lo = _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(lo), 0b_11_01_10_00));
            let hi = _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(hi), 0b_11_01_10_00));

            let h = _mm256_srli_epi32(lo, 16);
            let r = _mm256_clmulepi64_epi128(h, Self::FULLMODULUS, 0);
            let r = _mm256_srli_epi32(r, 16);
            let lo = _mm256_xor_si256(lo, _mm256_clmulepi64_epi128(r, Self::MODULUS, 0));

            let h = _mm256_srli_epi32(hi, 16);
            let r = _mm256_clmulepi64_epi128(h, Self::FULLMODULUS, 0);
            let r = _mm256_srli_epi32(r, 16);
            let hi = _mm256_xor_si256(hi, _mm256_clmulepi64_epi128(r, Self::MODULUS, 0));

            Self(_mm256_and_si256(_mm256_unpacklo_epi64(lo, hi), Self::MASK2))
        }
    }
}

impl MulAssign for P8GF2p16 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..8 { mul_call(); }
        unsafe {
            let lo = _mm256_and_si256(_mm256_clmulepi64_epi128(self.0, other.0, 0x00), Self::MASK);
            let hi = _mm256_and_si256(_mm256_clmulepi64_epi128(self.0, other.0, 0x11), Self::MASK);

            let lo = _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(lo), 0b_11_01_10_00));
            let hi = _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(hi), 0b_11_01_10_00));

            let h = _mm256_srli_epi32(lo, 16);
            let r = _mm256_clmulepi64_epi128(h, Self::FULLMODULUS, 0);
            let r = _mm256_srli_epi32(r, 16);
            let lo = _mm256_xor_si256(lo, _mm256_clmulepi64_epi128(r, Self::MODULUS, 0));

            let h = _mm256_srli_epi32(hi, 16);
            let r = _mm256_clmulepi64_epi128(h, Self::FULLMODULUS, 0);
            let r = _mm256_srli_epi32(r, 16);
            let hi = _mm256_xor_si256(hi, _mm256_clmulepi64_epi128(r, Self::MODULUS, 0));

            self.0 = _mm256_and_si256(_mm256_unpacklo_epi64(lo, hi), Self::MASK2);
        }
    }
}

impl PField<GF2p16, 8> for P8GF2p16 {}
