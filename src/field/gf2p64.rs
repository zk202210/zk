use super::*;
use super::common::*;
use crate::hash::Block512;


#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct GF2p64(pub u64);

impl GF2p64 {
    // GF2p8 to GF2p16 lifting table
    // const LFT_LUT: &[u32] = include_dwords!(
    //     concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p64.lft.lut"));

    // GF2p16 inv table
    // const INV_LUT: &[u32] = include_dwords!(
    //     concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p64.inv.lut"));

    const MODULUS: __m128i = unsafe {
        M128 { u64_array: [0b11011u64, 0u64] }.m128i
    };

    // fn reduce(a: [u128; 2]) -> u128 {
    //     // modolus: x^128 + x^7 + x^2 + x + 1
    //     let r: u128 = a[1] ^ (a[1] >> (128 - 7)) ^ (a[1] >> (128 - 2)) ^ (a[1] >> (128 - 1));
    //     a[0] ^ (r << 7) ^ (r << 2) ^ (r << 1) ^ (r << 0)
    // }
}

impl Inv for GF2p64 {
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();
        Self(mod_inverse_u64(self.0, 0b11011))
    }
}

impl Random for GF2p64 {
    type Output = GF2p64;
    // 0.inv = 0
    fn random() -> Self {
        let raw: u64 = rand::thread_rng().gen();
        Self(raw)
    }
}

impl Neg for GF2p64 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Add for GF2p64 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self(self.0 ^ other.0)
    }
}

impl AddAssign for GF2p64 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        self.0 ^= other.0;
    }
}

impl Sub for GF2p64 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self(self.0 ^ other.0)
    }
}

impl SubAssign for GF2p64 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        self.0 ^= other.0;
    }
}

impl Mul for GF2p64 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut res = self;
        res *= other;
        res
    }
}

impl MulAssign for GF2p64 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        unsafe {
            let a = _mm_clmulepi64_si128(
                _mm_set_epi64x(0, self.0 as i64),
                _mm_set_epi64x(0, other.0 as i64), 0);

            let r = _mm_clmulepi64_si128(a, Self::MODULUS, 0x01);
            let r = _mm_xor_si128(a, r);
            let r = _mm_clmulepi64_si128(r, Self::MODULUS, 0x01);
            let r = _mm_xor_si128(a, r);

            self.0 = _mm_cvtsi128_si64(r) as u64;
        }
    }
}

impl Div for GF2p64 {
    type Output = Self;

    // x/0 = 0
    // TODO
    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl DivAssign for GF2p64 {
    // x/0 = 0
    // TODO
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}

impl Size for GF2p64 {
    const SIZE: usize = 8;

    fn from_bytes(raw: &[u8; 8]) -> Self {
        Self(u64::from_le_bytes(*raw))
    }
}

impl Characteristic for GF2p64 {
    const CHAR: usize = 2;
}

impl std::convert::From<GF2p8> for GF2p64 {
    // TODO
    fn from(x: GF2p8) -> Self {
        Self(x.0 as u64)
        // Self(Self::LFT_LUT[x.0 as usize])
    }
}

impl std::fmt::Display for GF2p64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#018x}", self.0)
    }
}

impl std::convert::From<usize> for GF2p64 {
    fn from(x: usize) -> Self {
        Self(x as u64)
    }
}

impl std::convert::From<Block512> for GF2p64 {
    fn from(x: Block512) -> Self {
        Self(x.0[0] as u64)
    }
}

impl WriteWord for GF2p64 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0;
        *ptr += 1;
    }
}

impl std::convert::From<Vec<u64>> for GF2p64 {
    fn from(x: Vec<u64>) -> Self {
        Self(x[0])
    }
}

impl std::convert::Into<Vec<u64>> for GF2p64 {
    fn into(self) -> Vec<u64> {
        vec![self.0 as u64]
    }
}

impl std::convert::AsRef<[u8]> for GF2p64 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const u64, &[u8; 8]>
                (&self.0 as *const u64)
                as &[u8]
        }
    }
}

impl Pow for GF2p64 {}
impl Sqrt for GF2p64 {}
impl Xor for GF2p64 {}
impl _Field for GF2p64 {}

impl_pack_field!(GF2p64, P4GF2p64, 4);
impl_pack_field!(GF2p64, P8GF2p64, 8);
impl_pack_field!(GF2p64, P16GF2p64, 16);
impl_vect_field!(GF2p64, VGF2p64, 16);

impl FieldPack for GF2p64 {
    type P4 = P4GF2p64;
    type P8 = P8GF2p64;
    type P16 = P16GF2p64;
    type V16 = VGF2p64;
}

impl Field for GF2p64 {}
