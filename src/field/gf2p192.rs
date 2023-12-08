use super::*;
use super::common::*;
use crate::hash::Block512;


#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct GF2p192(pub [u64; 3]);

impl GF2p192 {
    // GF2p8 to GF2p16 lifting table
    // const LFT_LUT: &[u32] = include_dwords!(
    //     concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p192.lft.lut"));

    // GF2p16 inv table
    // const INV_LUT: &[u32] = include_dwords!(
    //     concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p192.inv.lut"));

    const MODULUS: __m128i = unsafe {
        // x^192 + x^7 + x^2 + x + 1
        M128 { u64_array: [0b10000111u64, 0u64] }.m128i
    };

    // fn reduce(a: [u128; 2]) -> u128 {
    //     // modolus: x^128 + x^7 + x^2 + x + 1
    //     let r: u128 = a[1] ^ (a[1] >> (128 - 7)) ^ (a[1] >> (128 - 2)) ^ (a[1] >> (128 - 1));
    //     a[0] ^ (r << 7) ^ (r << 2) ^ (r << 1) ^ (r << 0)
    // }
}

impl Inv for GF2p192 {
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();
        Self({
            let res = mod_inverse_u192(
                u256::from(self.0[0]) ^
                (u256::from(self.0[1]) << 64) ^
                (u256::from(self.0[2]) << 128),
                u256::from(0b10000111u64) ^
                (u256::from(1) << 192));
            [res.low_u64(), (res >> 64).low_u64(), (res >> 128).low_u64()]
        })
    }
}

impl Random for GF2p192 {
    type Output = GF2p192;
    // 0.inv = 0
    fn random() -> Self {
        let raw0: u64 = rand::thread_rng().gen();
        let raw1: u64 = rand::thread_rng().gen();
        let raw2: u64 = rand::thread_rng().gen();
        Self([raw0, raw1, raw2])
    }
}

impl Neg for GF2p192 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Add for GF2p192 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self([
             self.0[0] ^ other.0[0],
             self.0[1] ^ other.0[1],
             self.0[2] ^ other.0[2]])
    }
}

impl AddAssign for GF2p192 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        self.0[0] ^= other.0[0];
        self.0[1] ^= other.0[1];
        self.0[2] ^= other.0[2];
    }
}

impl Sub for GF2p192 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self([
             self.0[0] ^ other.0[0],
             self.0[1] ^ other.0[1],
             self.0[2] ^ other.0[2]])
    }
}

impl SubAssign for GF2p192 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        self.0[0] ^= other.0[0];
        self.0[1] ^= other.0[1];
        self.0[2] ^= other.0[2];
    }
}

impl Mul for GF2p192 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut res = self;
        res *= other;
        res
    }
}

// https://github.com/scipr-lab/libff/blob/develop/libff/algebra/fields/binary/gf192.cpp
impl MulAssign for GF2p192 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        unsafe {
            let ab0 = _mm_set_epi64x(self.0[0] as i64, other.0[0] as i64);
            let ab1 = _mm_set_epi64x(self.0[1] as i64, other.0[1] as i64);
            let ab2 = _mm_set_epi64x(self.0[2] as i64, other.0[2] as i64);

            let c0 = _mm_clmulepi64_si128(ab0, ab0, 0x01);
            let c4 = _mm_clmulepi64_si128(ab2, ab2, 0x01);
            let t = _mm_clmulepi64_si128(ab1, ab1, 0x01);
            let xor01 = _mm_xor_si128(ab0, ab1);
            let c1 = _mm_clmulepi64_si128(xor01, xor01, 0x01);
            let c1 = _mm_xor_si128(_mm_xor_si128(c1, c0), t);
            let xor12 = _mm_xor_si128(ab1, ab2);
            let c3 = _mm_clmulepi64_si128(xor12, xor12, 0x01);
            let c3 = _mm_xor_si128(_mm_xor_si128(c3, c4), t);
            let xor02 = _mm_xor_si128(ab0, ab2);
            let c2 = _mm_clmulepi64_si128(xor02, xor02, 0x01);
            let c2 = _mm_xor_si128(_mm_xor_si128(_mm_xor_si128(c2, c0), c4), t);
            let d0 = _mm_xor_si128(c0, _mm_slli_si128(c1, 8));
            let d2 = _mm_xor_si128(c4, _mm_srli_si128(c3, 8));
            let d1 = _mm_xor_si128(_mm_xor_si128(c2, _mm_srli_si128(c1, 8)),
                        _mm_slli_si128(c3, 8));
            let tmp = _mm_clmulepi64_si128(d2, Self::MODULUS, 0x01);
            let d1 = _mm_xor_si128(d1, tmp);

            let tmp = _mm_clmulepi64_si128(d2, Self::MODULUS, 0x00);
            let d1 = _mm_xor_si128(d1, _mm_srli_si128(tmp, 8));
            let d0 = _mm_xor_si128(d0, _mm_slli_si128(tmp, 8));

            let tmp = _mm_clmulepi64_si128(d1, Self::MODULUS, 0x01);
            let d0 = _mm_xor_si128(d0, tmp);

            _mm_storeu_si128(std::mem::transmute(&mut self.0 as *mut _ as *mut u8), d0);
            _mm_storel_epi64(std::mem::transmute(&mut self.0[2..] as *mut _ as *mut u8), d1);
        }
    }
}

impl Div for GF2p192 {
    type Output = Self;

    // x/0 = 0
    // TODO
    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl DivAssign for GF2p192 {
    // x/0 = 0
    // TODO
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}

impl Size for GF2p192 {
    const SIZE: usize = 24;

    fn from_bytes(raw: &[u8; 24]) -> Self {
        let raw: &[[u8; 8]; 3] = unsafe { std::mem::transmute(raw) };
        Self([
             u64::from_le_bytes(raw[0]),
             u64::from_le_bytes(raw[1]),
             u64::from_le_bytes(raw[2]),
        ])
    }
}

impl Characteristic for GF2p192 {
    const CHAR: usize = 2;
}

impl std::convert::From<GF2p8> for GF2p192 {
    // TODO
    fn from(x: GF2p8) -> Self {
        Self([x.0 as u64, 0u64, 0u64])
        // Self(Self::LFT_LUT[x.0 as usize])
    }
}

impl std::fmt::Display for GF2p192 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#018x}{:016x}{:016x}", self.0[2], self.0[1], self.0[0])
    }
}

impl std::convert::From<usize> for GF2p192 {
    fn from(x: usize) -> Self {
        Self([x as u64, 0u64, 0u64])
    }
}

impl std::convert::From<Block512> for GF2p192 {
    fn from(x: Block512) -> Self {
        Self([x.0[0], x.0[1], x.0[2]])
    }
}

impl WriteWord for GF2p192 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0[0];
        dst[*ptr + 1] = self.0[1];
        dst[*ptr + 2] = self.0[2];
        *ptr += 3;
    }
}

impl std::convert::From<Vec<u64>> for GF2p192 {
    fn from(x: Vec<u64>) -> Self {
        Self([x[0], x[1], x[2]])
    }
}

impl std::convert::Into<Vec<u64>> for GF2p192 {
    fn into(self) -> Vec<u64> {
        self.0.to_vec()
    }
}


impl std::convert::AsRef<[u8]> for GF2p192 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const [u64; 3], &[u8; 24]>
                (&self.0 as *const [u64; 3])
                as &[u8]
        }
    }
}

impl Pow for GF2p192 {}
impl Sqrt for GF2p192 {}
impl Xor for GF2p192 {}
impl _Field for GF2p192 {}

impl_pack_field!(GF2p192, P4GF2p192, 4);
impl_pack_field!(GF2p192, P8GF2p192, 8);
impl_pack_field!(GF2p192, P16GF2p192, 16);
impl_vect_field!(GF2p192, VGF2p192, 16);

impl FieldPack for GF2p192 {
    type P4 = P4GF2p192;
    type P8 = P8GF2p192;
    type P16 = P16GF2p192;
    type V16 = VGF2p192;
}

impl Field for GF2p192 {}
