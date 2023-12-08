use super::*;
use super::common::*;
use crate::hash::Block512;


#[derive(Copy, Clone)]
pub struct M61p2(pub [u64; 2]);

impl M61p2 {
    pub const PRIME: u64 = (1u64 << 61) - 1;
}

impl Size for M61p2 {
    const SIZE: usize = 16;

    fn from_bytes(raw: &[u8; 16]) -> Self {
        let raw: &[[u8; 8]; 2] = unsafe { std::mem::transmute(raw) };
        Self([
             u64::from_le_bytes(raw[0]),
             u64::from_le_bytes(raw[1])
        ])
    }
}

impl Characteristic for M61p2 {
    const CHAR: usize = (1 << 61) - 1;
}

impl Inv for M61p2 {
    // TODO
    // 0.inv = 0
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();
        let re = Self([self.0[0], 0]);
        let im = Self([self.0[1], 0]);
        let sumsq = Self::pow(re, 2) + Self::pow(im, 2);
        let denom = Self::pow(sumsq, (Self::PRIME - 2) as usize);
        let re = re * denom;
        let im = - im * denom;
        Self([re.0[0], im.0[0]])
    }
}

impl Pow for M61p2 {}

impl Sqrt for M61p2 {
    fn sqrt(x: Self) -> Self {
        // mod M61 has half quadratic residues and half non-residues (except 0)
        //
        // Fermat's little theorem
        // => a = a^{M61}
        // => a^2 = a^{M61 + 1} = a^{2^60} × a^{2^60}
        // => a = a^{2^60} or a = - a^{2^60}
        // => sqrt(a) = a^{2^59} or sqrt(a) does not exist
        //
        // sqrt(a) exists iff a^{(M61 - 1) / 2} = 1
        // sqrt(a) does not exist iff a^{(M61 - 1) / 2} = -1
        //
        // a + bi = (c + di)^2
        // => c = sqrt((delta + a) / 2)
        //    d = sqrt((delta - a) / 2) or - sqrt((delta - a) / 2)
        //    where delta = sqrt(a^2 + b^2)
        let inv_2 = Self::from(1 << 60);
        let re0 = Self::from(x.0[0] as usize);
        let im0 = Self::from(x.0[1] as usize);
        let tmp = re0 * re0 + im0 * im0;
        let delta = Self::pow(tmp, 1 << 59);

        if delta * delta != tmp {
            // not exist
            return Self::from(0);
        }

        let re_sq = (delta + re0) * inv_2;
        let im_sq = (delta - re0) * inv_2;

        let mut re = Self::pow(re_sq, 1 << 59);
        let mut im = Self::pow(im_sq, 1 << 59);
        if re * re != re_sq || im * im != im_sq {
            let delta = - delta;
            let re_sq = (delta + re0) * inv_2;
            let im_sq = (delta - re0) * inv_2;
            re = Self::pow(re_sq, 1 << 59);
            im = Self::pow(im_sq, 1 << 59);
        }
        if im0 != re * im * Self::from(2) {
            im = - im;
        }
        let res = Self([re.0[0], im.0[0]]);
        // assert_eq!(res * res, x);

        res
    }
}

impl Xor for M61p2 {
    #[inline]
    fn xor(x: Self, y: Self) -> Self {
        let xy = x * y;
        Self::from(1) - x - y + xy + xy
    }
}

impl Random for M61p2 {
    type Output = M61p2;

    fn random() -> Self {
        let res0: u64 =
            loop {
                let r: u64 = rand::thread_rng().gen();
                if (r & Self::PRIME) < Self::PRIME {
                    break r & Self::PRIME;
                }
            };
        let res1: u64 =
            loop {
                let r: u64 = rand::thread_rng().gen();
                if (r & Self::PRIME) < Self::PRIME {
                    break r & Self::PRIME;
                }
            };
        Self([res0, res1])
    }
}

#[inline]
fn m61p2_equal(x: u64, y: u64) -> bool {
    (x == y) || (x == 0 && y == M61p2::PRIME) || (x == M61p2::PRIME && y == 0)
}

impl PartialEq for M61p2 {
    fn eq(&self, other: &Self) -> bool {
        m61p2_equal(self.0[0], other.0[0]) && m61p2_equal(self.0[1], other.0[1])
    }
}

impl Eq for M61p2 {}

impl Neg for M61p2 {
    type Output = Self;

    fn neg(self) -> Self {
        let res0 = Self::PRIME ^ self.0[0];
        let res1 = Self::PRIME ^ self.0[1];
        Self([res0, res1])
    }
}

impl Add for M61p2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        let res0 = self.0[0] + other.0[0];
        let res0 = (res0 + (res0 >> 61)) & Self::PRIME;
        let res1 = self.0[1] + other.0[1];
        let res1 = (res1 + (res1 >> 61)) & Self::PRIME;
        Self([res0, res1])
    }
}

impl AddAssign for M61p2 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        self.0[0] += other.0[0];
        self.0[1] += other.0[1];
        self.0[0] = (self.0[0] + (self.0[0] >> 61)) & Self::PRIME;
        self.0[1] = (self.0[1] + (self.0[1] >> 61)) & Self::PRIME;
    }
}

impl Sub for M61p2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] sub_call();
        let res0 =
            if self.0[0] >= other.0[0] {
                self.0[0] - other.0[0]
            } else {
                self.0[0] + (Self::PRIME ^ other.0[0])
            };
        let res1 =
            if self.0[1] >= other.0[1] {
                self.0[1] - other.0[1]
            } else {
                self.0[1] + (Self::PRIME ^ other.0[1])
            };
        Self([res0, res1])
    }
}

impl SubAssign for M61p2 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] sub_call();
        if self.0[0] >= other.0[0] {
            self.0[0] -= other.0[0];
        } else {
            self.0[0] += Self::PRIME ^ other.0[0];
        }
        if self.0[1] >= other.0[1] {
            self.0[1] -= other.0[1];
        } else {
            self.0[1] += Self::PRIME ^ other.0[1];
        }
    }
}

fn mulx(a: u64, b: u64) -> u64 {
    let mut hi = a;
    unsafe {
        let res = _mulx_u64(a, b, &mut hi);
        ((hi << 3) | (res >> 61)) + (res & M61p2::PRIME)
    }
    // use std::arch::asm;
    // unsafe {
    //     asm!(
    //         "mov rdx, {a}",
    //         "mulx {d}, {c}, {b}",
    //         "shld {d}, {c}, 3",
    //         "and {c}, {modulus}",
    //         "add {c}, {d}",
    //         "mov {hi}, {c}",
    //         a = in(reg) a,
    //         b = in(reg) b,
    //         modulus = in(reg) M61p2::PRIME,
    //         hi = out(reg) hi,
    //         c = out(reg) _,
    //         d = out(reg) _,
    //         out("rdx") _,
    //     );
    //     hi
    // }
}

impl Mul for M61p2 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] mul_call();
        let mut res = self;
        res *= other;
        res
    }
}

impl MulAssign for M61p2 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        let tmp = mulx(self.0[0] + self.0[1], other.0[0] + other.0[1]);
        // (a+b)*(c+d)  in [0, 2p]
        self.0[0] = mulx(self.0[0], other.0[0]); // a*c   in [0, 2p]
        self.0[1] = mulx(self.0[1], other.0[1]); // b*d   in [0, 2p]
        for i in 0..2 {
            if self.0[i] >= Self::PRIME {
                self.0[i] -= Self::PRIME; // a*c, b*d   in [0, p]
            }
        }
        self.0[1] ^= Self::PRIME; // (a*c, -b*d) in [0, p]^2
        let tmp = tmp + (self.0[0] ^ Self::PRIME) + self.0[1]; // (a+b)*(c+d)-ac-bd in [0, 4p]
        let tmp = (tmp >> 61) + (tmp & Self::PRIME);
        self.0[0] += self.0[1];
        self.0[1] = tmp;
        for i in 0..2 {
            if self.0[i] >= Self::PRIME {
               self.0[i] -= Self::PRIME;
            }
        }
    }
}

impl Div for M61p2 {
    type Output = Self;

    // x/0 = 0
    // TODO
    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl DivAssign for M61p2 {
    // x/0 = 0
    // TODO
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}

impl std::convert::From<GF2p8> for M61p2 {
    // TODO
    fn from(x: GF2p8) -> Self {
        Self([x.0 as u64, 0u64])
        // Self(Self::LFT_LUT[x.0 as usize])
    }
}

impl std::fmt::Display for M61p2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let re = self.0[0];
        let im = self.0[1];
        let re: isize = if re > (Self::PRIME >> 1) { - ((Self::PRIME ^ re) as isize) } else { re as isize };
        let im: isize = if im > (Self::PRIME >> 1) { - ((Self::PRIME ^ im) as isize) } else { im as isize };
        write!(f, "({:#018x}, {:#018x})", re, im)
    }
}

impl std::fmt::Debug for M61p2 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let re = self.0[0];
        let im = self.0[1];
        let re: isize = if re > (Self::PRIME >> 1) { - ((Self::PRIME ^ re) as isize) } else { re as isize };
        let im: isize = if im > (Self::PRIME >> 1) { - ((Self::PRIME ^ im) as isize) } else { im as isize };
        write!(f, "({:#018x}, {:#018x})", re, im)
    }
}

impl std::convert::From<usize> for M61p2 {
    fn from(x: usize) -> Self {
        Self([(x as u64) & Self::PRIME, 0u64])
    }
}

impl std::convert::From<Block512> for M61p2 {
    fn from(x: Block512) -> Self {
        Self([x.0[0], x.0[1]])
    }
}

impl WriteWord for M61p2 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0[0];
        dst[*ptr + 1] = self.0[1];
        *ptr += 2;
    }
}

impl std::convert::From<Vec<u64>> for M61p2 {
    fn from(x: Vec<u64>) -> Self {
        Self([x[0] & Self::PRIME, x[1] & Self::PRIME])
    }
}

impl std::convert::Into<Vec<u64>> for M61p2 {
    fn into(self) -> Vec<u64> {
        self.0.to_vec()
    }
}


impl std::convert::AsRef<[u8]> for M61p2 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const [u64; 2], &[u8; 16]>
                (&self.0 as *const [u64; 2])
                as &[u8]
        }
    }
}

impl _Field for M61p2 {}

pub use super::p4m61p2::P4M61p2;
impl_pack_field!(M61p2, P8M61p2, 8);
impl_pack_field!(M61p2, P16M61p2, 16);
impl_vect_field!(M61p2, VM61p2, 16);

impl FieldPack for M61p2 {
    type P4 =  P4M61p2;
    type P8 =  P8M61p2;
    type P16 = P16M61p2;
    type V16 = VM61p2;
}

impl Field for M61p2 {}
