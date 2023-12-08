use super::*;
use super::common::*;
use crate::hash::Block512;


#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct GF2p128(pub [u64; 2]);

impl GF2p128 {
    // GF2p8 to GF2p16 lifting table
    // const LFT_LUT: &[u32] = include_dwords!(
    //     concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p128.lft.lut"));

    // GF2p16 inv table
    // const INV_LUT: &[u32] = include_dwords!(
    //     concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p128.inv.lut"));

    const MODULUS: __m128i = unsafe {
        M128 { u64_array: [0b10000111u64, 0u64] }.m128i
    };

    // fn reduce(a: [u128; 2]) -> u128 {
    //     // modolus: x^128 + x^7 + x^2 + x + 1
    //     let r: u128 = a[1] ^ (a[1] >> (128 - 7)) ^ (a[1] >> (128 - 2)) ^ (a[1] >> (128 - 1));
    //     a[0] ^ (r << 7) ^ (r << 2) ^ (r << 1) ^ (r << 0)
    // }
}

impl Inv for GF2p128 {
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();
        Self({
            let res = mod_inverse_u128(self.0[0] as u128 ^ ((self.0[1] as u128) << 64), 0b10000111);
            [res as u64, (res >> 64) as u64]
        })
    }
}

impl Random for GF2p128 {
    type Output = GF2p128;
    // 0.inv = 0
    fn random() -> Self {
        let raw0: u64 = rand::thread_rng().gen();
        let raw1: u64 = rand::thread_rng().gen();
        Self([raw0, raw1])
    }
}

impl Neg for GF2p128 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Add for GF2p128 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self([self.0[0] ^ other.0[0], self.0[1] ^ other.0[1]])
    }
}

impl AddAssign for GF2p128 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        self.0[0] ^= other.0[0];
        self.0[1] ^= other.0[1];
    }
}

impl Sub for GF2p128 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self([self.0[0] ^ other.0[0], self.0[1] ^ other.0[1]])
    }
}

impl SubAssign for GF2p128 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        self.0[0] ^= other.0[0];
        self.0[1] ^= other.0[1];
    }
}

impl Mul for GF2p128 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut res = self;
        res *= other;
        res
    }
}

impl MulAssign for GF2p128 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        unsafe {

            let a = _mm_loadu_si128(std::mem::transmute(&self.0 as *const _ as *const u8));
            let b = _mm_loadu_si128(std::mem::transmute(&other.0 as *const _ as *const u8));

            let hi = _mm_clmulepi64_si128(a, b, 0x11);
            let lo = _mm_clmulepi64_si128(a, b, 0x00);
            let hl = _mm_clmulepi64_si128(a, b, 0x10);
            let lh = _mm_clmulepi64_si128(a, b, 0x01);

            let mid = _mm_xor_si128(hl, lh);
            let hi = _mm_xor_si128(hi, _mm_srli_si128(mid, 8));
            let lo = _mm_xor_si128(lo, _mm_slli_si128(mid, 8));

            // reduce
            let tmp = _mm_clmulepi64_si128(hi, Self::MODULUS, 0x01);
            let lo = _mm_xor_si128(lo, _mm_slli_si128(tmp, 8));
            let hi = _mm_xor_si128(hi, _mm_srli_si128(tmp, 8));

            let tmp = _mm_clmulepi64_si128(hi, Self::MODULUS, 0x00);
            let lo = _mm_xor_si128(lo, tmp);

            _mm_storeu_si128(std::mem::transmute(&mut self.0 as *mut _ as *mut u8), lo);
        }
    }
}

impl Div for GF2p128 {
    type Output = Self;

    // x/0 = 0
    // TODO
    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl DivAssign for GF2p128 {
    // x/0 = 0
    // TODO
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}

impl Size for GF2p128 {
    const SIZE: usize = 16;

    fn from_bytes(raw: &[u8; 16]) -> Self {
        let raw: &[[u8; 8]; 2] = unsafe { std::mem::transmute(raw) };
        Self([
             u64::from_le_bytes(raw[0]),
             u64::from_le_bytes(raw[1])
        ])
    }
}

impl Characteristic for GF2p128 {
    const CHAR: usize = 2;
}

impl std::convert::From<GF2p8> for GF2p128 {
    // TODO
    fn from(x: GF2p8) -> Self {
        Self([x.0 as u64, 0u64])
        // Self(Self::LFT_LUT[x.0 as usize])
    }
}

impl std::fmt::Display for GF2p128 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#018x}{:016x}", self.0[1], self.0[0])
    }
}

impl std::convert::From<usize> for GF2p128 {
    fn from(x: usize) -> Self {
        Self([x as u64, 0u64])
    }
}

impl std::convert::From<Block512> for GF2p128 {
    fn from(x: Block512) -> Self {
        Self([x.0[0], x.0[1]])
    }
}

impl WriteWord for GF2p128 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0[0];
        dst[*ptr + 1] = self.0[1];
        *ptr += 2;
    }
}

impl std::convert::From<Vec<u64>> for GF2p128 {
    fn from(x: Vec<u64>) -> Self {
        Self([x[0], x[1]])
    }
}

impl std::convert::Into<Vec<u64>> for GF2p128 {
    fn into(self) -> Vec<u64> {
        self.0.to_vec()
    }
}


impl std::convert::AsRef<[u8]> for GF2p128 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const [u64; 2], &[u8; 16]>
                (&self.0 as *const [u64; 2])
                as &[u8]
        }
    }
}

impl Pow for GF2p128 {}
impl Sqrt for GF2p128 {}
impl Xor for GF2p128 {}
impl _Field for GF2p128 {}

impl_pack_field!(GF2p128, P4GF2p128, 4);
impl_pack_field!(GF2p128, P8GF2p128, 8);
impl_pack_field!(GF2p128, P16GF2p128, 16);
impl_vect_field!(GF2p128, VGF2p128, 16);

impl FieldPack for GF2p128 {
    type P4 = P4GF2p128;
    type P8 = P8GF2p128;
    type P16 = P16GF2p128;
    type V16 = VGF2p128;
}

impl Field for GF2p128 {}
