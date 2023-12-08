use super::*;
use super::common::*;
use crate::hash::Block512;

#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct GF2p32(pub u32);

impl GF2p32 {
    // GF2p8 to GF2p16 lifting table
    const LFT_LUT: &[u32] = include_dwords!(
        concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p32.lft.lut"));

    fn reduce(a: u64) -> u32 {
        // modulus = x^32 + x^7 + x^3 + x^2 + 1
        let h: u64 = a >> 32;
        let r: u64 = h ^ (h >> 25) ^ (h >> 29) ^ (h >> 30);
        (a ^ (r << 7) ^ (r << 3) ^ (r << 2) ^ r) as u32
    }
}

impl Inv for GF2p32 {
    // 0.inv = 0
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();
        Self(mod_inverse(self.0 as u64, (1u64 << 32) ^ 0b10001101) as u32)
    }
}

impl Random for GF2p32 {
    type Output = GF2p32;
    // 0.inv = 0
    fn random() -> Self {
        let raw: u32 = rand::thread_rng().gen();
        Self(raw)
    }
}

impl Neg for GF2p32 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
}

impl Add for GF2p32 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self(self.0 ^ other.0)
    }
}

impl AddAssign for GF2p32 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        *self = Self(self.0 ^ other.0);
    }
}

impl Sub for GF2p32 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self(self.0 ^ other.0)
    }
}

impl SubAssign for GF2p32 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        *self = Self(self.0 ^ other.0);
    }
}

impl Mul for GF2p32 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] mul_call();
        Self(Self::reduce(clmul(self.0 as u64, other.0 as u64)))
    }
}

impl MulAssign for GF2p32 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        *self = Self(Self::reduce(clmul(self.0 as u64, other.0 as u64)))
    }
}

impl Div for GF2p32 {
    type Output = Self;

    // x/0 = 0
    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl DivAssign for GF2p32 {
    // x/0 = 0
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}

impl Size for GF2p32 {
    const SIZE: usize = 4;

    fn from_bytes(raw: &[u8; 4]) -> Self {
        Self(u32::from_le_bytes(*raw))
    }
}

impl Characteristic for GF2p32 {
    const CHAR: usize = 2;
}

impl std::convert::From<GF2p8> for GF2p32 {
    fn from(x: GF2p8) -> Self {
        Self(Self::LFT_LUT[x.0 as usize])
    }
}

impl std::fmt::Display for GF2p32 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#010x}", self.0)
    }
}


impl std::convert::From<usize> for GF2p32 {
    fn from(x: usize) -> Self {
        Self(x as u32)
    }
}

impl std::convert::From<Block512> for GF2p32 {
    fn from(x: Block512) -> Self {
        Self(x.0[0] as u32)
    }
}

impl WriteWord for GF2p32 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0 as u64;
        *ptr += 1;
    }
}

impl std::convert::From<Vec<u64>> for GF2p32 {
    fn from(x: Vec<u64>) -> Self {
        Self(x[0] as u32)
    }
}

impl std::convert::Into<Vec<u64>> for GF2p32 {
    fn into(self) -> Vec<u64> {
        vec![self.0 as u64]
    }
}


impl std::convert::AsRef<[u8]> for GF2p32 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const u32, &[u8; 4]>
                (&self.0 as *const u32)
                as &[u8]
        }
    }
}

impl Pow for GF2p32 {}
impl Sqrt for GF2p32 {}
impl Xor for GF2p32 {}
impl _Field for GF2p32 {}

impl_pack_field!(GF2p32, P4GF2p32, 4);
impl_pack_field!(GF2p32, P8GF2p32, 8);
impl_pack_field!(GF2p32, P16GF2p32, 16);
impl_vect_field!(GF2p32, VGF2p32, 16);

impl FieldPack for GF2p32 {
    type P4 = P4GF2p32;
    type P8 = P8GF2p32;
    type P16 = P16GF2p32;
    type V16 = VGF2p32;
}

impl Field for GF2p32 {}
