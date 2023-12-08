use super::*;
use super::common::*;


#[derive(Copy, Clone, Debug)]
pub struct P4GF2p16(pub __m128i);

impl P4GF2p16 {
    const MASK: __m128i = unsafe {
        M128 { u64_array: [0xffffffffu64, 0xffffffffu64] }.m128i
    };

    const MASK2: __m128i = unsafe {
        M128 { u32_array:
            [0xffffu32, 0xffffu32,
             0xffffu32, 0xffffu32] }.m128i
    };

    const ZERO: __m128i = unsafe {
        M128 { u8_array: [0; 16] }.m128i
    };

    const FULLMODULUS: __m128i = unsafe {
        M128 { u64_array: [0b10000000000101011u64, 0u64] }.m128i
    };

    const MODULUS: __m128i = unsafe {
        M128 { u64_array: [0b101011u64, 0u64] }.m128i
    };
}

impl Pack<4> for P4GF2p16 {
    type Output = Self;
    type Input = GF2p16;

    fn pack(raw: &[GF2p16; 4]) -> Self {
        unsafe {
            let ptr: *const i64 = std::mem::transmute(raw);
            let tmp = _mm_set_epi64x(0, *ptr);
            Self(_mm_unpacklo_epi16(tmp, Self::ZERO))
        }
    }

    fn unpack(&self, dst: &mut [GF2p16; 4]) {
        unsafe {
            let tmp = _mm_packus_epi32(self.0, self.0);
            let ptr: *mut u64 = std::mem::transmute(dst);
            *ptr = _mm_cvtsi128_si64(tmp) as u64;
        }
    }
}

impl Add for P4GF2p16 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..4 { add_call(); }
        unsafe {
            Self(_mm_xor_si128(self.0, other.0))
        }
    }
}

impl AddAssign for P4GF2p16 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..4 { add_call(); }
        unsafe {
            self.0 = _mm_xor_si128(self.0, other.0);
        }
    }
}

impl Sub for P4GF2p16 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..4 { add_call(); }
        unsafe {
            Self(_mm_xor_si128(self.0, other.0))
        }
    }
}

impl SubAssign for P4GF2p16 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..4 { add_call(); }
        unsafe {
            self.0 = _mm_xor_si128(self.0, other.0);
        }
    }
}

impl Mul for P4GF2p16 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..4 { mul_call(); }
        unsafe {
            let lo = _mm_and_si128(_mm_clmulepi64_si128(self.0, other.0, 0x00), Self::MASK);
            let hi = _mm_and_si128(_mm_clmulepi64_si128(self.0, other.0, 0x11), Self::MASK);

            let lo = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(lo), 0b_11_01_10_00));
            let hi = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(hi), 0b_11_01_10_00));

            let h = _mm_srli_epi32(lo, 16);
            let r = _mm_clmulepi64_si128(h, Self::FULLMODULUS, 0);
            let r = _mm_srli_epi32(r, 16);
            let lo = _mm_xor_si128(lo, _mm_clmulepi64_si128(r, Self::MODULUS, 0));

            let h = _mm_srli_epi32(hi, 16);
            let r = _mm_clmulepi64_si128(h, Self::FULLMODULUS, 0);
            let r = _mm_srli_epi32(r, 16);
            let hi = _mm_xor_si128(hi, _mm_clmulepi64_si128(r, Self::MODULUS, 0));

            Self(_mm_and_si128(_mm_unpacklo_epi64(lo, hi), Self::MASK2))
        }
    }
}

impl MulAssign for P4GF2p16 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..4 { mul_call(); }
        unsafe {
            let lo = _mm_and_si128(_mm_clmulepi64_si128(self.0, other.0, 0x00), Self::MASK);
            let hi = _mm_and_si128(_mm_clmulepi64_si128(self.0, other.0, 0x11), Self::MASK);

            let lo = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(lo), 0b_11_01_10_00));
            let hi = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(hi), 0b_11_01_10_00));

            let h = _mm_srli_epi32(lo, 16);
            let r = _mm_clmulepi64_si128(h, Self::FULLMODULUS, 0);
            let r = _mm_srli_epi32(r, 16);
            let lo = _mm_xor_si128(lo, _mm_clmulepi64_si128(r, Self::MODULUS, 0));

            let h = _mm_srli_epi32(hi, 16);
            let r = _mm_clmulepi64_si128(h, Self::FULLMODULUS, 0);
            let r = _mm_srli_epi32(r, 16);
            let hi = _mm_xor_si128(hi, _mm_clmulepi64_si128(r, Self::MODULUS, 0));

            self.0 = _mm_and_si128(_mm_unpacklo_epi64(lo, hi), Self::MASK2);
        }
    }
}

impl PField<GF2p16, 4> for P4GF2p16 {}
