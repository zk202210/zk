use super::*;
use super::common::*;


#[derive(Copy, Clone, Debug)]
pub struct P4M61p2(pub [__m256i; 2]);

impl P4M61p2 {
    const PACK_PRIME: __m256i = unsafe {
        M256 { u64_array: [M61p2::PRIME; 4] }.m256i
    };

    const PACK_PRIME_M1: __m256i = unsafe {
        M256 { u64_array: [M61p2::PRIME - 1; 4] }.m256i
    };
}


impl Pack<4> for P4M61p2 {
    type Output = Self;
    type Input = M61p2;

    fn pack(raw: &[M61p2; 4]) -> Self {
        // TODO SIMD pack/unpack
        unsafe {
            Self([
                _mm256_setr_epi64x(
                    raw[0].0[0] as i64, raw[1].0[0] as i64,
                    raw[2].0[0] as i64, raw[3].0[0] as i64),
                _mm256_setr_epi64x(
                    raw[0].0[1] as i64, raw[1].0[1] as i64,
                    raw[2].0[1] as i64, raw[3].0[1] as i64)
            ])
        }
    }

    fn unpack(&self, dst: &mut [M61p2; 4]) {
        unsafe {
            // let mut dst = [0u64; 8];

            // _mm256_store_si256(std::mem::transmute(&mut dst[0..2] as *mut _ as *mut u8), self.0[0]);
            // _mm256_store_si256(std::mem::transmute(&mut dst[4..8] as *mut _ as *mut u8), self.0[1]);

            // [
            //     M61p2([dst[0], dst[4]]),
            //     M61p2([dst[1], dst[5]]),
            //     M61p2([dst[2], dst[6]]),
            //     M61p2([dst[3], dst[7]])
            // ]

            let r01 = _mm256_permute4x64_epi64(self.0[0], 0b00010000);
            let r23 = _mm256_permute4x64_epi64(self.0[0], 0b00110010);
            let i01 = _mm256_permute4x64_epi64(self.0[1], 0b00010000);
            let i23 = _mm256_permute4x64_epi64(self.0[1], 0b00110010);

            let dst0 = _mm256_unpacklo_epi64(r01, i01);
            let dst1 = _mm256_unpacklo_epi64(r23, i23);

            _mm256_storeu_si256(std::mem::transmute(&mut dst[0..2] as *mut _ as *mut u8), dst0);
            _mm256_storeu_si256(std::mem::transmute(&mut dst[2..4] as *mut _ as *mut u8), dst1);
        }
    }
}

impl Add for P4M61p2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..4 { add_call(); }
        unsafe {
            let r = _mm256_add_epi64(self.0[0], other.0[0]);
            let i = _mm256_add_epi64(self.0[1], other.0[1]);
            let maskr = _mm256_cmpgt_epi64(r, Self::PACK_PRIME_M1);
            let maski = _mm256_cmpgt_epi64(i, Self::PACK_PRIME_M1);
            let r = _mm256_sub_epi64(r, _mm256_and_si256(maskr, Self::PACK_PRIME));
            let i = _mm256_sub_epi64(i, _mm256_and_si256(maski, Self::PACK_PRIME));
            Self([r, i])
        }
    }
}

impl AddAssign for P4M61p2 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..4 { add_call(); }
        unsafe {
            self.0[0] = _mm256_add_epi64(self.0[0], other.0[0]);
            self.0[1] = _mm256_add_epi64(self.0[1], other.0[1]);
            let maskr = _mm256_cmpgt_epi64(self.0[0], Self::PACK_PRIME_M1);
            let maski = _mm256_cmpgt_epi64(self.0[1], Self::PACK_PRIME_M1);
            self.0[0] = _mm256_sub_epi64(self.0[0], _mm256_and_si256(maskr, Self::PACK_PRIME));
            self.0[1] = _mm256_sub_epi64(self.0[1], _mm256_and_si256(maski, Self::PACK_PRIME));
        }
    }
}

impl Sub for P4M61p2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..4 { sub_call(); }
        unsafe {
            let tmpr = _mm256_xor_si256(other.0[0], Self::PACK_PRIME);
            let tmpi = _mm256_xor_si256(other.0[1], Self::PACK_PRIME);
            let r = _mm256_add_epi64(self.0[0], tmpr);
            let i = _mm256_add_epi64(self.0[1], tmpi);
            let maskr = _mm256_cmpgt_epi64(r, Self::PACK_PRIME_M1);
            let maski = _mm256_cmpgt_epi64(i, Self::PACK_PRIME_M1);
            let r = _mm256_sub_epi64(r, _mm256_and_si256(maskr, Self::PACK_PRIME));
            let i = _mm256_sub_epi64(i, _mm256_and_si256(maski, Self::PACK_PRIME));
            Self([r, i])
        }
    }
}

impl SubAssign for P4M61p2  {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..4 { sub_call(); }
        unsafe {
            let tmpr = _mm256_xor_si256(other.0[0], Self::PACK_PRIME);
            let tmpi = _mm256_xor_si256(other.0[1], Self::PACK_PRIME);
            self.0[0] = _mm256_add_epi64(self.0[0], tmpr);
            self.0[1] = _mm256_add_epi64(self.0[1], tmpi);
            let maskr = _mm256_cmpgt_epi64(self.0[0], Self::PACK_PRIME_M1);
            let maski = _mm256_cmpgt_epi64(self.0[1], Self::PACK_PRIME_M1);
            self.0[0] = _mm256_sub_epi64(self.0[0], _mm256_and_si256(maskr, Self::PACK_PRIME));
            self.0[1] = _mm256_sub_epi64(self.0[1], _mm256_and_si256(maski, Self::PACK_PRIME));
        }
    }
}

fn pack_mulx(x: __m256i, y: __m256i) -> __m256i {
    unsafe {
        let x_shift = _mm256_srli_epi64(x, 32);
        let y_shift = _mm256_srli_epi64(y, 32);

        let bd = _mm256_mul_epu32(x, y);
        let ac = _mm256_mul_epu32(x_shift, y_shift);
        let ad = _mm256_mul_epu32(x_shift, y);
        let bc = _mm256_mul_epu32(x, y_shift);

        let ad_bc = _mm256_add_epi64(ad, bc);
        let bd_srl32 = _mm256_srli_epi64(bd, 32);
        let ad_bc_srl32 = _mm256_srli_epi64(_mm256_add_epi64(ad_bc, bd_srl32), 32);
        let ad_bc_sll32 = _mm256_slli_epi64(ad_bc, 32);
        let hi = _mm256_add_epi64(ac, ad_bc_srl32);
        let lo = _mm256_add_epi64(bd, ad_bc_sll32);

        _mm256_add_epi64(_mm256_or_si256(
                _mm256_slli_epi64(hi, 3),
                _mm256_srli_epi64(lo, 61)),
                _mm256_and_si256(lo, P4M61p2::PACK_PRIME))
    }
}


impl Mul for P4M61p2 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] for _ in 0..4 { mul_call(); }
        unsafe {
            let tmp = pack_mulx(
                _mm256_add_epi64(self.0[0], self.0[1]),
                _mm256_add_epi64(other.0[0], other.0[1]));
            let ac = pack_mulx(self.0[0], other.0[0]);
            let bd = pack_mulx(self.0[1], other.0[1]);

            let mask = _mm256_cmpgt_epi64(bd, Self::PACK_PRIME_M1);
            let bd = _mm256_sub_epi64(bd, _mm256_and_si256(mask, Self::PACK_PRIME));
            let bd = _mm256_xor_si256(bd, Self::PACK_PRIME);

            let mask = _mm256_cmpgt_epi64(ac, Self::PACK_PRIME_M1);
            let ac = _mm256_sub_epi64(ac, _mm256_and_si256(mask, Self::PACK_PRIME));
            let nac = _mm256_xor_si256(ac, Self::PACK_PRIME);

            let i = _mm256_add_epi64(_mm256_add_epi64(tmp, nac), bd);
            let i = _mm256_add_epi64(_mm256_srli_epi64(i, 61), _mm256_and_si256(i, Self::PACK_PRIME));

            let r = _mm256_add_epi64(ac, bd);

            let mask = _mm256_cmpgt_epi64(r, Self::PACK_PRIME_M1);
            let r = _mm256_sub_epi64(r, _mm256_and_si256(mask, Self::PACK_PRIME));
            let mask = _mm256_cmpgt_epi64(i, Self::PACK_PRIME_M1);
            let i = _mm256_sub_epi64(i, _mm256_and_si256(mask, Self::PACK_PRIME));
            Self([r, i])
        }
    }
}

impl MulAssign for P4M61p2 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] for _ in 0..4 { mul_call(); }
        unsafe {
            let tmp = pack_mulx(
                _mm256_add_epi64(self.0[0], self.0[1]),
                _mm256_add_epi64(other.0[0], other.0[1]));
            let ac = pack_mulx(self.0[0], other.0[0]);
            let bd = pack_mulx(self.0[1], other.0[1]);

            let mask = _mm256_cmpgt_epi64(ac, Self::PACK_PRIME_M1);
            let ac = _mm256_sub_epi64(ac, _mm256_and_si256(mask, Self::PACK_PRIME));
            let nac = _mm256_xor_si256(ac, Self::PACK_PRIME);

            let mask = _mm256_cmpgt_epi64(bd, Self::PACK_PRIME_M1);
            let bd = _mm256_sub_epi64(bd, _mm256_and_si256(mask, Self::PACK_PRIME));
            let bd = _mm256_xor_si256(bd, Self::PACK_PRIME);

            let i = _mm256_add_epi64(_mm256_add_epi64(tmp, nac), bd);
            let i = _mm256_add_epi64(_mm256_srli_epi64(i, 61), _mm256_and_si256(i, Self::PACK_PRIME));

            let r = _mm256_add_epi64(ac, bd);

            let mask = _mm256_cmpgt_epi64(r, Self::PACK_PRIME_M1);
            self.0[0] = _mm256_sub_epi64(r, _mm256_and_si256(mask, Self::PACK_PRIME));
            let mask = _mm256_cmpgt_epi64(i, Self::PACK_PRIME_M1);
            self.0[1] = _mm256_sub_epi64(i, _mm256_and_si256(mask, Self::PACK_PRIME));
        }
    }
}

impl PField<M61p2, 4> for P4M61p2 {}
