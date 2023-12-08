use super::*;
use super::common::*;
use crate::hash::Block512;


#[derive(Copy, Clone)]
pub struct M127p2(pub [u128; 2]);

impl M127p2 {
    pub const PRIME: u128 = (1 << 127) - 1;

    pub fn pow_u128(self, p: u128) -> Self {
        let mut p = p;
        let mut res = Self::from(1);
        let mut tmp = self;
        while p != 0u128 {
            if p & 1u128 == 1u128 {
                res *= tmp
            }
            tmp = tmp * tmp;
            p >>= 1;
        }
        res
    }
}

impl Size for M127p2 {
    const SIZE: usize = 32;

    fn from_bytes(raw: &[u8; 32]) -> Self {
        let raw: &[[u8; 16]; 2] = unsafe { std::mem::transmute(raw) };
        Self([
             u128::from_le_bytes(raw[0]),
             u128::from_le_bytes(raw[1])
        ])
    }
}

impl Characteristic for M127p2 {
    const CHAR: usize = 0;
}

impl Inv for M127p2 {
    // TODO
    // 0.inv = 0
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();
        let re = Self([self.0[0], 0]);
        let im = Self([self.0[1], 0]);
        let sumsq = Self::pow(re, 2) + Self::pow(im, 2);
        let denom = Self::pow_u128(sumsq, Self::PRIME - 2);
        let re = re * denom;
        let im = - im * denom;
        Self([re.0[0], im.0[0]])
    }
}

impl Pow for M127p2 {}

impl Sqrt for M127p2 {
    fn sqrt(x: Self) -> Self {
        // mod M127 has half quadratic residues and half non-residues (except 0)
        //
        // Fermat's little theorem
        // => a = a^{M127}
        // => a^2 = a^{M127 + 1} = a^{2^126} × a^{2^126}
        // => a = a^{2^126} or a = - a^{2^126}
        // => sqrt(a) = a^{2^125} or sqrt(a) does not exist
        //
        // sqrt(a) exists iff a^{(M127 - 1) / 2} = 1
        // sqrt(a) does not exist iff a^{(M127 - 1) / 2} = -1
        //
        // a + bi = (c + di)^2
        // => c = sqrt((delta + a) / 2)
        //    d = sqrt((delta - a) / 2) or - sqrt((delta - a) / 2)
        //    where delta = sqrt(a^2 + b^2)
        let inv_2 = Self([1u128 << 126, 0]);
        let re0 = Self([x.0[0], 0]);
        let im0 = Self([x.0[1], 0]);
        let tmp = re0 * re0 + im0 * im0;
        let delta = Self::pow_u128(tmp, 1u128 << 125);

        if delta * delta != tmp {
            // not exist
            return Self::from(0);
        }

        let re_sq = (delta + re0) * inv_2;
        let im_sq = (delta - re0) * inv_2;

        let mut re = Self::pow_u128(re_sq, 1u128 << 125);
        let mut im = Self::pow_u128(im_sq, 1u128 << 125);
        if re * re != re_sq || im * im != im_sq {
            let delta = - delta;
            let re_sq = (delta + re0) * inv_2;
            let im_sq = (delta - re0) * inv_2;
            re = Self::pow_u128(re_sq, 1u128 << 125);
            im = Self::pow_u128(im_sq, 1u128 << 125);
        }
        if im0 != re * im + re * im {
            im = - im;
        }
        let res = Self([re.0[0], im.0[0]]);
        assert_eq!(res * res, x);

        res
    }
}

impl Xor for M127p2 {
    #[inline]
    fn xor(x: Self, y: Self) -> Self {
        let xy = x * y;
        Self::from(1) - x - y + xy + xy
    }
}

impl Random for M127p2 {
    type Output = M127p2;

    fn random() -> Self {
        let res0: u128 =
            loop {
                let r: u128 = rand::thread_rng().gen();
                if (r & Self::PRIME) < Self::PRIME {
                    break r & Self::PRIME;
                }
            };
        let res1: u128 =
            loop {
                let r: u128 = rand::thread_rng().gen();
                if (r & Self::PRIME) < Self::PRIME {
                    break r & Self::PRIME;
                }
            };
        Self([res0, res1])
    }
}

#[inline]
fn m127p2_equal(x: u128, y: u128) -> bool {
    (x == y) || (x == 0 && y == M127p2::PRIME) || (x == M127p2::PRIME && y == 0)
}

impl PartialEq for M127p2 {
    fn eq(&self, other: &Self) -> bool {
        m127p2_equal(self.0[0], other.0[0]) && m127p2_equal(self.0[1], other.0[1])
    }
}

impl Eq for M127p2 {}

impl Neg for M127p2 {
    type Output = Self;

    fn neg(self) -> Self {
        let res0 = Self::PRIME ^ self.0[0];
        let res1 = Self::PRIME ^ self.0[1];
        let res = Self([res0, res1]);
        res
    }
}

impl Add for M127p2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        let res0 = self.0[0] + other.0[0];
        let res0 = (res0 + (res0 >> 127)) & Self::PRIME;
        let res1 = self.0[1] + other.0[1];
        let res1 = (res1 + (res1 >> 127)) & Self::PRIME;
        Self([res0, res1])
    }
}

impl AddAssign for M127p2 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();
        self.0[0] += other.0[0];
        self.0[1] += other.0[1];
        self.0[0] = (self.0[0] + (self.0[0] >> 127)) & Self::PRIME;
        self.0[1] = (self.0[1] + (self.0[1] >> 127)) & Self::PRIME;
    }
}

impl Sub for M127p2 {
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

impl SubAssign for M127p2 {
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

#[inline]
fn mulx(a: u128, b: u128) -> u128 {
    let mut y0: u64;
    let mut y1: u64;
    // let mut carry: u64 = 0;
    unsafe {
        std::arch::asm!(
            "mov rdx, {b0}",
            "mulx {y1}, {y0}, {a0}",
            "mulx {y5}, {y4}, {a1}",
            "mov rdx, {b1}",
            "mulx {y3}, {y2}, {a0}",
            "mulx {y7}, {y6}, {a1}",

            // y7, (y6+y5+y3), (y4+y2+y1), y0 => y7,y3,y1,y0
            "add {y1}, {y2}",
            "adc {y3}, {y5}",
            "adc {y7}, 0",
            "add {y1}, {y4}",
            "adc {y3}, {y6}",
            "adc {y7}, 0",

            // (hi << 1 | lo >> 127) => (carry, y7, y3)
            "shld {y7}, {y3}, 1",
            "mov rdx, 0",
            "adc rdx, 0",
            "shld {y3}, {y1}, 1",

            // (y1, y0) | PRIME => (y1, y0)
            "shl {y1}, 1",
            "shr {y1}, 1",

            // (carry, y7, y3) + (y1, y0) => (carry, y1, y0)
            "add {y0}, {y3}",
            "adc {y1}, {y7}",
            "adc rdx, 0",

            "shld rdx, {y1}, 1",
            "shl {y1}, 1",
            "shr {y1}, 1",
            "add {y0}, rdx",
            "adc {y1}, 0",
            // "mov {carry}, 0",
            // "adc {carry}, 0",

            out("rdx") _,
            a0 = in(reg) a as u64,
            a1 = in(reg) (a >> 64) as u64,
            b0 = in(reg) b as u64,
            b1 = in(reg) (b >> 64) as u64,
            y0 = out(reg) y0,
            y1 = out(reg) y1,
            y2 = out(reg) _,
            y3 = out(reg) _,
            y4 = out(reg) _,
            y5 = out(reg) _,
            y6 = out(reg) _,
            y7 = out(reg) _,
        );
    }
    (y0 as u128) | (y1 as u128) << 64

    // assert_eq!(carry, 0);
    // if carry != 0 {
    //     println!("carry {}, quotient {}", carry, u256::from(res) / u256::from(M127p2::PRIME));
    // }
    // assert_eq!(u256::from(res) % u256::from(M127p2::PRIME), u256::from(a) * u256::from(b) % u256::from(M127p2::PRIME));

    // use primitive_types::U128;
    // let res: u256 = U128::full_mul(U128::from(a), U128::from(b));

    // let hi = (res >> 128).as_u128();
    // let lo = res.low_u128();

    // let res = ((hi << 1) | (lo >> 127)) + (lo & M127p2::PRIME);
    // let res = (res >> 127) + (res & M127p2::PRIME);
    // res

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
    //         modulus = in(reg) M127p2::PRIME,
    //         hi = out(reg) hi,
    //         c = out(reg) _,
    //         d = out(reg) _,
    //         out("rdx") _,
    //     );
    //     hi
    // }
}

impl Mul for M127p2 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] mul_call();
        let mut res = self;
        res *= other;
        res
    }
}

impl MulAssign for M127p2 {
    fn mul_assign(&mut self, other: Self) {

        // // a.re * b.re - a.im * b.im
        // // a.re * b.im + a.re * b.im
        // let ar = u256::from(self.0[0]);
        // let ai = u256::from(self.0[1]);
        // let br = u256::from(other.0[0]);
        // let bi = u256::from(other.0[1]);

        // let p = u256::from(Self::PRIME);

        // let rr = (ar * br) % p;
        // let ii = (ai * bi) % p;
        // let r = if rr > ii {
        //     rr - ii
        // } else {
        //     rr + p - ii
        // };

        // let ri = (ar * bi) % p;
        // let ir = (ai * br) % p;
        // let i = (ri + ir) % p;

        // self.0[0] = r.as_u128();
        // self.0[1] = i.as_u128();

        #[cfg(debug_assertions)] mul_call();
        let a = self.0[0] + self.0[1];
        let b = other.0[0] + other.0[1];
        // let a = (a >> 127) + (a & M127p2::PRIME);
        // let b = (b >> 127) + (b & M127p2::PRIME);
        let mut tmp = mulx(a, b);
        // if tmp >= Self::PRIME {
        //     tmp -= Self::PRIME;
        // }
        // (a+b)*(c+d)  in [0, 2p]
        self.0[0] = mulx(self.0[0], other.0[0]); // a*c   in [0, 2p]
        self.0[1] = mulx(self.0[1], other.0[1]); // b*d   in [0, 2p]
        // for i in 0..2 {
        //     if self.0[i] >= Self::PRIME {
        //         self.0[i] -= Self::PRIME; // a*c, b*d   in [0, p]
        //     }
        // }
        self.0[1] ^= Self::PRIME; // (a*c, -b*d) in [0, p]^2
        let tmp2 = (self.0[0] ^ Self::PRIME) + self.0[1]; // -ac-bd in [0, 2p]
        let tmp2 = (tmp2 >> 127) + (tmp2 & Self::PRIME);
        tmp += tmp2;
        tmp = (tmp >> 127) + (tmp & Self::PRIME);
        self.0[0] += self.0[1];
        self.0[0] = (self.0[0] >> 127) + (self.0[0] & Self::PRIME);
        self.0[1] = tmp;
        // for i in 0..2 {
        //     if self.0[i] >= Self::PRIME {
        //        self.0[i] -= Self::PRIME;
        //     }
        // }
    }
}

impl Div for M127p2 {
    type Output = Self;

    // x/0 = 0
    // TODO
    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl DivAssign for M127p2 {
    // x/0 = 0
    // TODO
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}

impl std::convert::From<GF2p8> for M127p2 {
    // TODO
    fn from(x: GF2p8) -> Self {
        Self([x.0 as u128, 0u128])
        // Self(Self::LFT_LUT[x.0 as usize])
    }
}

impl std::fmt::Display for M127p2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let re = self.0[0];
        let im = self.0[1];
        write!(f, "({:#018x}{:016x}, {:#018x}{:016x})",
            (re >> 64) as usize, re as usize, (im >> 64) as usize, im as usize)
    }
}

impl std::fmt::Debug for M127p2 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let re = self.0[0];
        let im = self.0[1];
        write!(f, "({:#018x}{:016x}, {:#018x}{:016x})",
            (re >> 64) as usize, re as usize, (im >> 64) as usize, im as usize)
    }
}

impl std::convert::From<usize> for M127p2 {
    fn from(x: usize) -> Self {
        Self([(x as u128) & Self::PRIME, 0u128])
    }
}

impl std::convert::From<Block512> for M127p2 {
    fn from(x: Block512) -> Self {
        Self([
             ((x.0[0] as u128) + (x.0[1] as u128) << 64) & Self::PRIME,
             ((x.0[2] as u128) + (x.0[3] as u128) << 64) & Self::PRIME])
    }
}

impl WriteWord for M127p2 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0[0] as u64;
        dst[*ptr + 1] = (self.0[0] >> 64) as u64;
        dst[*ptr + 2] = self.0[1] as u64;
        dst[*ptr + 3] = (self.0[1] >> 64) as u64;
        *ptr += 4;
    }
}

impl std::convert::From<Vec<u64>> for M127p2 {
    fn from(x: Vec<u64>) -> Self {
        // println!("from vector {:016x} {:016x} {:016x} {:016x}", x[0], x[1], x[2], x[3]);
        let res =
        Self([
             ((x[0] as u128) + ((x[1] as u128) << 64)) & Self::PRIME,
             ((x[2] as u128) + ((x[3] as u128) << 64)) & Self::PRIME]);
        // println!("into {}", res);
        res
    }
}

impl std::convert::Into<Vec<u64>> for M127p2 {
    fn into(self) -> Vec<u64> {
        unsafe {
            let ptr = std::mem::transmute::<*const [u128; 2], *const [u64; 4]>
                (&self.0 as *const [u128; 2]);
            (*ptr).to_vec()
        }
    }
}


impl std::convert::AsRef<[u8]> for M127p2 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const [u128; 2], &[u8; 32]>
                (&self.0 as *const [u128; 2])
                as &[u8]
        }
    }
}

impl _Field for M127p2 {}

impl_pack_field!(M127p2, P4M127p2, 4);
impl_pack_field!(M127p2, P8M127p2, 8);
impl_pack_field!(M127p2, P16M127p2, 16);
impl_vect_field!(M127p2, VM127p2, 16);

impl FieldPack for M127p2 {
    type P4 =  P4M127p2;
    type P8 =  P8M127p2;
    type P16 = P16M127p2;
    type V16 = VM127p2;
}

impl Field for M127p2 {}
