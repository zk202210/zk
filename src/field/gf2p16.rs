use super::*;
use super::common::*;
use crate::hash::Block512;

#[derive(PartialEq, Eq, Copy, Clone, Debug)]
#[repr(transparent)]
pub struct GF2p16(pub u16);


impl GF2p16 {
    // GF2p8 to GF2p16 lifting table
    const LFT_LUT: &[u16] = include_words!(
        concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p16.lft.lut"));

    // GF2p16 inv table
    const INV_LUT: &[u16] = include_words!(
        concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p16.inv.lut"));

    // GF2p16 mut table
    // const MUL_LUT: &[u16] = include_words!(
    //     concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p16.mut.lut"));

    // const FULLMODULUS: __m128i = unsafe {
    //     M128 { u64_array: [0b10000000000101011u64, 0u64] }.m128i
    // };

    // const MODULUS: __m128i = unsafe {
    //     M128 { u64_array: [0b101011u64, 0u64] }.m128i
    // };

    fn reduce(a: u64) -> u16 {
        // modolus: x^16 + x^5 + x^3 + x + 1
        let h: u64 = a >> 16;
        let r: u64 = h ^ (h >> 11) ^ (h >> 13) ^ (h >> 15);
        (a ^ (r << 5) ^ (r << 3) ^ (r << 1) ^ r) as u16
    }
}

impl Inv for GF2p16 {
    // 0.inv = 0
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();
        Self(Self::INV_LUT[self.0 as usize])
    }
}

// impl Pow for GF2p16 {
//     fn sqrt(&self, p: usize) -> Self {
//         *self
//     }
// }

impl Random for GF2p16 {
    type Output = GF2p16;
    // 0.inv = 0
    fn random() -> Self {
        let raw: u16 = rand::thread_rng().gen();
        Self(raw)
    }
}

impl Neg for GF2p16 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Add for GF2p16 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self(self.0 ^ other.0)
    }
}

impl AddAssign for GF2p16 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        *self = Self(self.0 ^ other.0);
    }
}

impl Sub for GF2p16 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self(self.0 ^ other.0)
    }
}

impl SubAssign for GF2p16 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        *self = Self(self.0 ^ other.0);
    }
}

impl Mul for GF2p16 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        // Self(Self::reduce(clmul(self.0 as u64, other.0 as u64)))
        let mut res = self;
        res *= other;
        res
    }
}

impl MulAssign for GF2p16 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        *self = Self(Self::reduce(clmul(self.0 as u64, other.0 as u64)))
    }
}

// impl MulAssign for GF2p16 {
//     fn mul_assign(&mut self, other: Self) {
//         unsafe {
//             let a = _mm_clmulepi64_si128(
//                 _mm_set_epi64x(0, self.0 as i64),
//                 _mm_set_epi64x(0, other.0 as i64), 0);
//             let h = _mm_srli_epi64(a, 16);
//             let r = _mm_clmulepi64_si128(h, Self::FULLMODULUS, 0);
//             let r = _mm_srli_epi64(r, 16);
//             let res = _mm_xor_si128(a, _mm_clmulepi64_si128(r, Self::MODULUS, 0));
//             self.0 = _mm_cvtsi128_si64(res) as u16;
//         }
//     }
// }

// impl MulAssign for GF2p16 {
//     fn mul_assign(&mut self, other: Self) {
//         let idx = other.0 as usize;
//         let idy = self.0 as usize;
//         self.0 =
//             Self::MUL_LUT[((idx & 0xff) << 9) + (idy & 0xff)] ^
//             Self::MUL_LUT[((idx & 0xff) << 9) + (idy >> 8) + 256] ^
//             Self::MUL_LUT[(((idx >> 8) + 256) << 9) + (idy & 0xff)] ^
//             Self::MUL_LUT[(((idx >> 8) + 256) << 9) + (idy >> 8) + 256];
//     }
// }


impl Div for GF2p16 {
    type Output = Self;

    // x/0 = 0
    fn div(self, other: Self) -> Self {
        Self(Self::reduce(clmul(
                    self.0 as u64,
                    Self::INV_LUT[other.0 as usize] as u64)))
    }
}

impl DivAssign for GF2p16 {
    // x/0 = 0
    fn div_assign(&mut self, other: Self) {
        *self = Self(Self::reduce(clmul(
                    self.0 as u64,
                    Self::INV_LUT[other.0 as usize] as u64)))
    }
}

impl Size for GF2p16 {
    const SIZE: usize = 2;

    fn from_bytes(raw: &[u8; Self::SIZE]) -> Self {
        Self(u16::from_le_bytes(*raw))
    }
}

impl Characteristic for GF2p16 {
    const CHAR: usize = 2;
}


impl std::convert::From<GF2p8> for GF2p16 {
    fn from(x: GF2p8) -> Self {
        Self(Self::LFT_LUT[x.0 as usize])
    }
}

impl std::fmt::Display for GF2p16 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#06x}", self.0)
    }
}

impl std::convert::From<usize> for GF2p16 {
    fn from(x: usize) -> Self {
        Self(x as u16)
    }
}

impl std::convert::From<Block512> for GF2p16 {
    fn from(x: Block512) -> Self {
        Self(x.0[0] as u16)
    }
}

impl WriteWord for GF2p16 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0 as u64;
        *ptr += 1;
    }
}

impl std::convert::From<Vec<u64>> for GF2p16 {
    fn from(x: Vec<u64>) -> Self {
        Self(x[0] as u16)
    }
}

impl std::convert::Into<Vec<u64>> for GF2p16 {
    fn into(self) -> Vec<u64> {
        vec![self.0 as u64]
    }
}

impl std::convert::AsRef<[u8]> for GF2p16 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const u16, &[u8; 2]>
                (&self.0 as *const u16)
                as &[u8]
        }
    }
}

impl Pow for GF2p16 {}
impl Sqrt for GF2p16 {}
impl Xor for GF2p16 {}
impl _Field for GF2p16 {}

pub use super::p4gf2p16::P4GF2p16;
pub use super::p8gf2p16::P8GF2p16;
pub use super::p16gf2p16::P16GF2p16;
pub use super::vgf2p16::VGF2p16;

impl FieldPack for GF2p16 {
    type P4 = P4GF2p16;
    type P8 = P8GF2p16;
    type P16 = P16GF2p16;
    type V16 = VGF2p16;
}

impl Field for GF2p16 {}
