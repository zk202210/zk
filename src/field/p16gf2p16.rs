use super::*;
use super::common::*;


#[derive(Copy, Clone, Debug)]
pub struct P16GF2p16(pub __m512i);

impl P16GF2p16 {
    const MASK: __m512i = unsafe {
        M512 { u64_array: [0xffffffffu64; 8] }.m512i
    };

    const MASK2: __m512i = unsafe {
        M512 { u32_array: [0xffffu32; 16] }.m512i
    };

    const PERMUTE: __m512i = unsafe {
        M512 { u16_array: [
            0x00, 0x01, 0x02, 0x03,
            0x08, 0x09, 0x0a, 0x0b,
            0x10, 0x11, 0x12, 0x13,
            0x18, 0x19, 0x1a, 0x1b,
            0x04, 0x05, 0x06, 0x07,
            0x0c, 0x0d, 0x0e, 0x0f,
            0x14, 0x15, 0x16, 0x17,
            0x1c, 0x1d, 0x1e, 0x1f,
        ] }.m512i
    };

    const FULLMODULUS: __m512i = unsafe {
        M512 { u64_array: [
            0b10000000000101011u64, 0u64,
            0b10000000000101011u64, 0u64,
            0b10000000000101011u64, 0u64,
            0b10000000000101011u64, 0u64,
        ] }.m512i
    };

    const MODULUS: __m512i = unsafe {
        M512 { u64_array: [
            0b101011u64, 0u64,
            0b101011u64, 0u64,
            0b101011u64, 0u64,
            0b101011u64, 0u64,
        ] }.m512i
    };
}

impl Pack<16> for P16GF2p16 {
    type Output = Self;
    type Input = GF2p16;

    fn pack(raw: &[GF2p16; 16]) -> Self {
        // TODO SIMD
        unsafe {
            Self(_mm512_setr_epi32(
                    raw[0].0 as u32 as i32,
                    raw[1].0 as u32 as i32,
                    raw[2].0 as u32 as i32,
                    raw[3].0 as u32 as i32,
                    raw[4].0 as u32 as i32,
                    raw[5].0 as u32 as i32,
                    raw[6].0 as u32 as i32,
                    raw[7].0 as u32 as i32,
                    raw[8].0 as u32 as i32,
                    raw[9].0 as u32 as i32,
                    raw[10].0 as u32 as i32,
                    raw[11].0 as u32 as i32,
                    raw[12].0 as u32 as i32,
                    raw[13].0 as u32 as i32,
                    raw[14].0 as u32 as i32,
                    raw[15].0 as u32 as i32,
                    ))
        }
    }

    fn unpack(&self, dst: &mut [GF2p16; 16]) {
        unsafe {
            let tmp = _mm512_packus_epi32(self.0, self.0);
            let tmp = _mm512_castsi512_si256(_mm512_permutexvar_epi16(Self::PERMUTE, tmp));
            _mm256_storeu_si256(std::mem::transmute(dst as *mut _ as *mut u8), tmp);
        }
    }
}

impl Add for P16GF2p16 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..16 { add_call(); }
        unsafe {
            Self(_mm512_xor_si512(self.0, other.0))
        }
    }
}

impl AddAssign for P16GF2p16 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..16 { add_call(); }
        unsafe {
            self.0 = _mm512_xor_si512(self.0, other.0);
        }
    }
}

impl Sub for P16GF2p16 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..16 { add_call(); }
        unsafe {
            Self(_mm512_xor_si512(self.0, other.0))
        }
    }
}

impl SubAssign for P16GF2p16 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..16 { add_call(); }
        unsafe {
            self.0 = _mm512_xor_si512(self.0, other.0);
        }
    }
}

impl Mul for P16GF2p16 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..16 { mul_call(); }
        unsafe {
            let lo = _mm512_and_si512(_mm512_clmulepi64_epi128(self.0, other.0, 0x00), Self::MASK);
            let hi = _mm512_and_si512(_mm512_clmulepi64_epi128(self.0, other.0, 0x11), Self::MASK);

            let lo = _mm512_castps_si512(_mm512_permute_ps(_mm512_castsi512_ps(lo), 0b_11_01_10_00));
            let hi = _mm512_castps_si512(_mm512_permute_ps(_mm512_castsi512_ps(hi), 0b_11_01_10_00));

            let h = _mm512_srli_epi32(lo, 16);
            let r = _mm512_clmulepi64_epi128(h, Self::FULLMODULUS, 0);
            let r = _mm512_srli_epi32(r, 16);
            let lo = _mm512_xor_si512(lo, _mm512_clmulepi64_epi128(r, Self::MODULUS, 0));

            let h = _mm512_srli_epi32(hi, 16);
            let r = _mm512_clmulepi64_epi128(h, Self::FULLMODULUS, 0);
            let r = _mm512_srli_epi32(r, 16);
            let hi = _mm512_xor_si512(hi, _mm512_clmulepi64_epi128(r, Self::MODULUS, 0));

            Self(_mm512_and_si512(_mm512_unpacklo_epi64(lo, hi), Self::MASK2))
        }
    }
}

impl MulAssign for P16GF2p16 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..16 { mul_call(); }
        unsafe {
            let lo = _mm512_and_si512(_mm512_clmulepi64_epi128(self.0, other.0, 0x00), Self::MASK);
            let hi = _mm512_and_si512(_mm512_clmulepi64_epi128(self.0, other.0, 0x11), Self::MASK);

            let lo = _mm512_castps_si512(_mm512_permute_ps(_mm512_castsi512_ps(lo), 0b_11_01_10_00));
            let hi = _mm512_castps_si512(_mm512_permute_ps(_mm512_castsi512_ps(hi), 0b_11_01_10_00));

            let h = _mm512_srli_epi32(lo, 16);
            let r = _mm512_clmulepi64_epi128(h, Self::FULLMODULUS, 0);
            let r = _mm512_srli_epi32(r, 16);
            let lo = _mm512_xor_si512(lo, _mm512_clmulepi64_epi128(r, Self::MODULUS, 0));

            let h = _mm512_srli_epi32(hi, 16);
            let r = _mm512_clmulepi64_epi128(h, Self::FULLMODULUS, 0);
            let r = _mm512_srli_epi32(r, 16);
            let hi = _mm512_xor_si512(hi, _mm512_clmulepi64_epi128(r, Self::MODULUS, 0));

            self.0 = _mm512_and_si512(_mm512_unpacklo_epi64(lo, hi), Self::MASK2);
        }
    }
}

impl PField<GF2p16, 16> for P16GF2p16 {}
