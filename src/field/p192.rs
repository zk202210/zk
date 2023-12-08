use super::*;
use super::common::*;
use crate::hash::Block512;


#[derive(Copy, Clone, PartialEq, Eq)]
pub struct P192(pub [u64; 3]);

/* prime candidates:
 *      2**192 -  19 * (2**30) + 1
 *      2**192 - 453 * (2**31) + 1
 *      2**192 -  49 * (2**32) + 1
 *      2**192 - 209 * (2**33) + 1
 *      2**192 -  75 * (2**34) + 1
 *      2**192 - 297 * (2**35) + 1
 *      2**192 -  13 * (2**36) + 1
 *      2**192 - 165 * (2**37) + 1
 *      2**192 -  21 * (2**38) + 1
 *      2**192 -  81 * (2**39) + 1
 *      2**192 -  75 * (2**40) + 1
 *      2**192 -   5 * (2**41) + 1
 *      2**192 - 295 * (2**42) + 1
 *      2**192 -  41 * (2**43) + 1
 *      2**192 -  99 * (2**44) + 1
 *      2**192 - 227 * (2**45) + 1
 *      2**192 - 375 * (2**46) + 1
 *      2**192 -  23 * (2**47) + 1
 *      2**192 -  21 * (2**48) + 1
 *      2**192 -  95 * (2**49) + 1
 *      2**192 -  37 * (2**50) + 1
 *      2**192 -  15 * (2**51) + 1
 *      2**192 -  19 * (2**52) + 1
 *      2**192 - 345 * (2**53) + 1
 *      2**192 - 321 * (2**54) + 1
 *      2**192 - 527 * (2**55) + 1
 *      2**192 -  15 * (2**56) + 1
 *      2**192 -  27 * (2**57) + 1
 *      2**192 -  49 * (2**58) + 1
 *      2**192 - 233 * (2**59) + 1
 *      2**192 -  25 * (2**60) + 1
 *      2**192 -  17 * (2**61) + 1
 *      2**192 - 139 * (2**62) + 1
 *      2**192 -   5 * (2**63) + 1
 *      2**192 - 129 * (2**64) + 1
 */

impl P192 {
    pub const K: u64 = 15 * (1 << 56) - 1;
    // 2**192 - K
    pub const PRIME: u256 = u256([
        0xf100000000000001,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0,
    ]);
    // (prime - 1) / 2
    pub const SQRT_POW: u256 = u256([
        0xf880000000000000,
        0xffffffffffffffff,
        0x7fffffffffffffff,
        0,
    ]);

    pub fn to_u256(&self) -> u256 {
        u256([self.0[0], self.0[1], self.0[2], 0])
    }
}

impl Size for P192 {
    const SIZE: usize = 24;

    fn from_bytes(raw: &[u8; 24]) -> Self {
        let raw: &[[u8; 8]; 3] = unsafe { std::mem::transmute(raw) };
        Self([
             u64::from_le_bytes(raw[0]),
             u64::from_le_bytes(raw[1]),
             u64::from_le_bytes(raw[2])
        ])
    }
}

impl Characteristic for P192 {
    const CHAR: usize = 0;
}

impl Inv for P192 {
    // TODO
    // 0.inv = 0
    fn inv(&self) -> Self {
        #[cfg(debug_assertions)] inv_call();

        let mut old_t = u256::from(0);
        let mut t = u256::from(1);
        let mut old_r = Self::PRIME;
        let mut r = self.to_u256();
        let mut positive = true;
        let mut old_positive = true;

        while r != u256::from(0) {
            let (q, new_r) = u256::div_mod(old_r, r);

            old_r = r;
            r = new_r;

            let tmp = old_t;
            let tmp_positive = old_positive;
            old_t = t;
            old_positive = positive;
            if old_positive && tmp_positive {
                if tmp < q * old_t {
                    t = q * old_t - tmp;
                    positive = false;
                } else {
                    t = tmp - q * old_t;
                }
            } else if !old_positive && tmp_positive {
                t = tmp + q * old_t;
                positive = true;
            } else if old_positive && !tmp_positive {
                t = tmp + q * old_t;
                positive = false;
            } else {
                if tmp < q * old_t {
                    t = q * old_t - tmp;
                    positive = true;
                } else {
                    t = tmp - q * old_t;
                }
            }
        }
        if !old_positive {
            old_t = Self::PRIME - old_t;
        }
        Self([old_t.0[0], old_t.0[1], old_t.0[2]])
    }
}

impl Pow for P192 {}


fn pow_u256(x: P192, pow: u256) -> P192 {
    let mut p = pow;
    let mut res = P192::from(1);
    let mut tmp = x;
    while p != u256::from(0) {
        if p & u256::from(1) == u256::from(1) {
            res *= tmp
        }
        tmp = tmp * tmp;
        p >>= 1;
    }
    res
}

fn log2_ord(x: P192) -> usize {
    let mut res = 0;
    let mut x = x;
    while x != P192::from(1) {
        x = x * x;
        res += 1;
    }
    res
}

impl Sqrt for P192 {
    fn sqrt(a: Self) -> Self {
        // https://www.math.canterbury.ac.nz/~j.booher/expos/sqr_qnr.pdf
        let r = 56;
        let s = (Self::PRIME - u256::from(1)) >> r;
        let v = Self::from(2);
        let w = pow_u256(v, s);
        let a_inv = a.inv();

        let mut w_power = Vec::with_capacity(r);
        w_power.push(w);
        for i in 0..55 {
            let w = w_power[i];
            w_power.push(w * w);
        }

        let p = (s + u256::from(1)) >> 1;
        let mut x = pow_u256(a, p);
        let mut t = log2_ord(x * x * a_inv);

        while t != 0 {
            x = x * w_power[r - 2 - t];
            t = log2_ord(x * x * a_inv);
        }
        x
    }
}

impl Xor for P192 {
    #[inline]
    fn xor(x: Self, y: Self) -> Self {
        let xy = x * y;
        Self::from(1) - x - y + xy + xy
    }
}

impl Random for P192 {
    type Output = P192;

    fn random() -> Self {
        let (x0, x1, x2): (u64, u64, u64) =
            loop {
                let x0: u64 = rand::thread_rng().gen();
                let x1: u64 = rand::thread_rng().gen();
                let x2: u64 = rand::thread_rng().gen();
                if u256([x0, x1, x2, 0]) < Self::PRIME {
                    break (x0, x1, x2);
                }
            };
        Self([x0, x1, x2])
    }
}

impl Neg for P192 {
    type Output = Self;

    fn neg(self) -> Self {
        if self == Self::from(0) {
            return self
        }
        let mut res = self;
        unsafe {
            std::arch::asm!(
                "add {res0}, {k}",
                "adc {res1}, 0",
                "adc {res2}, 0",
                "not {res0}",
                "not {res1}",
                "not {res2}",
                "add {res0}, 1",
                "adc {res1}, 0",
                "adc {res2}, 0",
                res0 = inout(reg) res.0[0],
                res1 = inout(reg) res.0[1],
                res2 = inout(reg) res.0[2],
                k = in(reg) Self::K,
            );
        }
        res
    }
}

impl Add for P192 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        let mut res = self;
        res += other;
        res
    }
}

impl AddAssign for P192 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();

        let mut carry = 0u64;
        unsafe {
            std::arch::asm!(
                "add {self0}, {other0}",
                "adc {self1}, {other1}",
                "adc {self2}, {other2}",
                "adc {carry}, 0",
                self0 = inout(reg) self.0[0],
                self1 = inout(reg) self.0[1],
                self2 = inout(reg) self.0[2],
                carry = inout(reg) carry,
                other0 = in(reg) other.0[0],
                other1 = in(reg) other.0[1],
                other2 = in(reg) other.0[2],
            );
        }

        // let mut c = res.to_u256();
        if carry == 1 {
            // a+b = 2^192 + res
            // a+b mod p = res + k
            unsafe {
                std::arch::asm!(
                    "add {self0}, {k}",
                    "adc {self1}, 0",
                    "adc {self2}, 0",
                    self0 = inout(reg) self.0[0],
                    self1 = inout(reg) self.0[1],
                    self2 = inout(reg) self.0[2],
                    k = in(reg) Self::K,
                );
            }
        } else if self.0[2] == -1i64 as u64 &&
                  self.0[1] == -1i64 as u64 &&
                  self.0[0] >= Self::PRIME.0[0] {
            self.0[2] = 0;
            self.0[1] = 0;
            self.0[0] -= Self::PRIME.0[0];
        }
    }
}

impl Sub for P192 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] sub_call();
        let mut res = - other;
        res += self;
        res
    }
}

impl SubAssign for P192 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] sub_call();
        *self += - other;
    }
}

impl Mul for P192 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] mul_call();
        let a = self.to_u256();
        let b = other.to_u256();
        let c = u256::full_mul(a, b);

        let mut carry = 0u64;
        let mut res = Self([c.0[0], c.0[1], c.0[2]]);
        unsafe {
            std::arch::asm!(
                // xhigh * k + xlow = res + (carry << 192)
                // where res is =192 bits
                // and carry is <64 bits
                "mov rdx, {k}",
                "mulx {y1}, {y0}, {xh0}",
                "mulx {y3}, {y2}, {xh1}",
                "mulx {y5}, {y4}, {xh2}",
                "add {res0}, {y0}",
                "adc {res1}, {y1}",
                "adc {res2}, {y3}",
                "adc {carry}, {y5}",
                "add {res1}, {y2}",
                "adc {res2}, {y4}",
                "adc {carry}, 0",
                // res + carry * k => res + (carry << 192)
                "mulx {y1}, {y0}, {carry}",
                "mov {carry}, 0",
                "add {res0}, {y0}",
                "adc {res1}, {y1}",
                "adc {res2}, 0",
                "adc {carry}, 0",
                out("rdx") _,
                k = in(reg) Self::K,
                xh0 = in(reg) c.0[3],
                xh1 = in(reg) c.0[4],
                xh2 = in(reg) c.0[5],
                y0 = out(reg) _,
                y1 = out(reg) _,
                y2 = out(reg) _,
                y3 = out(reg) _,
                y4 = out(reg) _,
                y5 = out(reg) _,
                res0 = inout(reg) res.0[0],
                res1 = inout(reg) res.0[1],
                res2 = inout(reg) res.0[2],
                carry = inout(reg) carry,
            );
        }

        if carry == 1 {
            // (2^192 + res) mod p = res + k
            unsafe {
                std::arch::asm!(
                    "add {self0}, {k}",
                    "adc {self1}, 0",
                    "adc {self2}, 0",
                    self0 = inout(reg) res.0[0],
                    self1 = inout(reg) res.0[1],
                    self2 = inout(reg) res.0[2],
                    k = in(reg) Self::K,
                );
            }
        } else if res.0[2] == -1i64 as u64 &&
                  res.0[1] == -1i64 as u64 &&
                  res.0[0] >= Self::PRIME.0[0] {
            res.0[2] = 0;
            res.0[1] = 0;
            res.0[0] -= Self::PRIME.0[0];
        }

        res
    }
}

impl MulAssign for P192 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        *self = *self * other;
    }
}

impl Div for P192 {
    type Output = Self;

    // x/0 = 0
    // TODO
    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl DivAssign for P192 {
    // x/0 = 0
    // TODO
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}

impl std::convert::From<GF2p8> for P192 {
    // TODO
    fn from(x: GF2p8) -> Self {
        Self([x.0 as u64, 0u64, 0u64])
    }
}

impl std::fmt::Display for P192 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#018x}{:016x}{:016x}", self.0[2], self.0[1], self.0[0])
    }
}

impl std::fmt::Debug for P192 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:#018x}{:016x}{:016x}", self.0[2], self.0[1], self.0[0])
    }
}

impl std::convert::From<usize> for P192 {
    fn from(x: usize) -> Self {
        Self([x as u64, 0u64, 0u64])
    }
}

impl std::convert::From<Block512> for P192 {
    fn from(x: Block512) -> Self {
        Self([x.0[0], x.0[1], x.0[2]])
    }
}

impl WriteWord for P192 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0[0];
        dst[*ptr + 1] = self.0[1];
        dst[*ptr + 2] = self.0[2];
        *ptr += 3;
    }
}

impl std::convert::From<Vec<u64>> for P192 {
    fn from(x: Vec<u64>) -> Self {
        Self([x[0], x[1], x[2]])
    }
}

impl std::convert::Into<Vec<u64>> for P192 {
    fn into(self) -> Vec<u64> {
        self.0.to_vec()
    }
}


impl std::convert::AsRef<[u8]> for P192 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const [u64; 3], &[u8; 24]>
                (&self.0 as *const [u64; 3])
                as &[u8]
        }
    }
}

impl _Field for P192 {}

impl_pack_field!(P192, P4P192, 4);
impl_pack_field!(P192, P8P192, 8);
impl_pack_field!(P192, P16P192, 16);
impl_vect_field!(P192, VP192, 16);

impl FieldPack for P192 {
    type P4 =  P4P192;
    type P8 =  P8P192;
    type P16 = P16P192;
    type V16 = VP192;
}

impl Field for P192 {}
