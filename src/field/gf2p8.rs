use super::*;
use super::common::*;

#[derive(PartialEq, Eq, Copy, Clone)]
pub struct GF2p8(pub u8);

impl GF2p8 {
    const INV_LUT: &[u8] = include_bytes!(
        concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p8.inv.lut"));
    const MUL_LUT: &[u8] = include_bytes!(
        concat!(env!("CARGO_MANIFEST_DIR"), "/resources/gf2p8.mut.lut"));
}

impl Pow for GF2p8 {}

impl Inv for GF2p8 {
    // 0.inv = 0
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();
        Self(Self::INV_LUT[self.0 as usize])
    }
}

impl Add for GF2p8 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self(self.0 ^ other.0)
    }
}

impl AddAssign for GF2p8 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        *self = Self(self.0 ^ other.0);
    }
}

impl Sub for GF2p8 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        Self(self.0 ^ other.0)
    }
}

impl SubAssign for GF2p8 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        *self = Self(self.0 ^ other.0);
    }
}

impl Mul for GF2p8 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] mul_call();
        Self(Self::MUL_LUT[(self.0 as usize) * 256 + (other.0 as usize)])
    }
}

impl MulAssign for GF2p8 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        *self = Self(Self::MUL_LUT[(self.0 as usize) * 256 + (other.0 as usize)]);
    }
}

impl Div for GF2p8 {
    type Output = Self;

    // x/0 = 0
    fn div(self, other: Self) -> Self {
        Self(Self::MUL_LUT[(self.0 as usize) * 256 + (Self::INV_LUT[other.0 as usize] as usize)])
    }
}

impl DivAssign for GF2p8 {
    // x/0 = 0
    fn div_assign(&mut self, other: Self) {
        *self = Self(Self::MUL_LUT[(self.0 as usize) * 256 + (Self::INV_LUT[other.0 as usize] as usize)])
    }
}

impl std::convert::From<usize> for GF2p8 {
    fn from(x: usize) -> Self {
        Self(x as u8)
    }
}

impl std::fmt::Debug for GF2p8 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#04x}", self.0)
    }
}
