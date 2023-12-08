use super::*;
use super::common::*;
use crate::hash::Block512;


#[derive(Copy, Clone, PartialEq, Eq)]
pub struct P256(pub [u64; 4]);

/* prime candidates:
 * 2**256 -  37 * (2**30) + 1
 * 2**256 -  23 * (2**31) + 1
 * 2**256 -  43 * (2**32) + 1
 * 2**256 - 299 * (2**33) + 1
 * 2**256 - 127 * (2**34) + 1
 * 2**256 - 105 * (2**35) + 1
 * 2**256 -  19 * (2**36) + 1
 * 2**256 -  23 * (2**37) + 1
 * 2**256 -  75 * (2**38) + 1
 * 2**256 - 281 * (2**39) + 1
 * 2**256 - 261 * (2**40) + 1
 * 2**256 -  87 * (2**41) + 1
 * 2**256 - 109 * (2**42) + 1
 * 2**256 -  83 * (2**43) + 1
 * 2**256 - 375 * (2**44) + 1
 * 2**256 - 209 * (2**45) + 1
 * 2**256 - 135 * (2**46) + 1
 * 2**256 - 243 * (2**47) + 1
 * 2**256 - 363 * (2**48) + 1
 * 2**256 - 237 * (2**49) + 1
 * 2**256 -  19 * (2**50) + 1
 * 2**256 -  21 * (2**51) + 1
 * 2**256 - 139 * (2**52) + 1
 * 2**256 -  39 * (2**53) + 1
 * 2**256 - 301 * (2**54) + 1
 * 2**256 - 165 * (2**55) + 1
 * 2**256 - 115 * (2**56) + 1
 * 2**256 - 693 * (2**57) + 1
 * 2**256 -  85 * (2**58) + 1
 * 2**256 - 113 * (2**59) + 1
 * 2**256 - 343 * (2**60) + 1
 * 2**256 - 239 * (2**61) + 1
 * 2**256 -  91 * (2**62) + 1
 * 2**256 -  51 * (2**63) + 1
 * 2**256 - 139 * (2**64) + 1
 */

impl P256 {
    pub const K: u64 = 19 * (1 << 50) - 1;
    // 2**256 - K
    pub const PRIME: u256 = u256([
        0xffb4000000000001,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
    ]);
    // (prime - 1) / 2
    pub const SQRT_POW: u256 = u256([
        0xffda000000000000,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x7fffffffffffffff,
    ]);

    pub fn to_u256(&self) -> u256 {
        u256([self.0[0], self.0[1], self.0[2], self.0[3]])
    }
}

impl Size for P256 {
    const SIZE: usize = 32;

    fn from_bytes(raw: &[u8; 32]) -> Self {
        let raw: &[[u8; 8]; 4] = unsafe { std::mem::transmute(raw) };
        Self([
             u64::from_le_bytes(raw[0]),
             u64::from_le_bytes(raw[1]),
             u64::from_le_bytes(raw[2]),
             u64::from_le_bytes(raw[3])
        ])
    }
}

impl Characteristic for P256 {
    const CHAR: usize = 0;
}

impl Inv for P256 {
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
        Self([old_t.0[0], old_t.0[1], old_t.0[2], old_t.0[3]])
    }
}

impl Pow for P256 {}


fn pow_u256(x: P256, pow: u256) -> P256 {
    let mut p = pow;
    let mut res = P256::from(1);
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

fn log2_ord(x: P256) -> usize {
    let mut res = 0;
    let mut x = x;
    while x != P256::from(1) {
        x = x * x;
        res += 1;
    }
    res
}

// TODO XXX
impl Sqrt for P256 {
    fn sqrt(a: Self) -> Self {
        // https://www.math.canterbury.ac.nz/~j.booher/expos/sqr_qnr.pdf
        let r = 50;
        let s = (Self::PRIME - u256::from(1)) >> r;
        let v = Self::from(3);
        let w = pow_u256(v, s);
        let a_inv = a.inv();

        let mut w_power = Vec::with_capacity(r);
        w_power.push(w);
        for i in 0..49 {
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

impl Xor for P256 {
    #[inline]
    fn xor(x: Self, y: Self) -> Self {
        let xy = x * y;
        Self::from(1) - x - y + xy + xy
    }
}

impl Random for P256 {
    type Output = P256;

    fn random() -> Self {
        let (x0, x1, x2, x3): (u64, u64, u64, u64) =
            loop {
                let x0: u64 = rand::thread_rng().gen();
                let x1: u64 = rand::thread_rng().gen();
                let x2: u64 = rand::thread_rng().gen();
                let x3: u64 = rand::thread_rng().gen();
                if u256([x0, x1, x2, x3]) < Self::PRIME {
                    break (x0, x1, x2, x3);
                }
            };
        Self([x0, x1, x2, x3])
    }
}

impl Neg for P256 {
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
                "adc {res3}, 0",
                "not {res0}",
                "not {res1}",
                "not {res2}",
                "not {res3}",
                "add {res0}, 1",
                "adc {res1}, 0",
                "adc {res2}, 0",
                "adc {res3}, 0",
                res0 = inout(reg) res.0[0],
                res1 = inout(reg) res.0[1],
                res2 = inout(reg) res.0[2],
                res3 = inout(reg) res.0[3],
                k = in(reg) Self::K,
            );
        }
        res
    }
}

impl Add for P256 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        #[cfg(debug_assertions)] add_call();
        let mut res = self;
        res += other;
        res
    }
}

impl AddAssign for P256 {
    fn add_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] add_call();

        let mut carry = 0u64;
        unsafe {
            std::arch::asm!(
                "add {self0}, {other0}",
                "adc {self1}, {other1}",
                "adc {self2}, {other2}",
                "adc {self3}, {other3}",
                "adc {carry}, 0",
                self0 = inout(reg) self.0[0],
                self1 = inout(reg) self.0[1],
                self2 = inout(reg) self.0[2],
                self3 = inout(reg) self.0[3],
                carry = inout(reg) carry,
                other0 = in(reg) other.0[0],
                other1 = in(reg) other.0[1],
                other2 = in(reg) other.0[2],
                other3 = in(reg) other.0[3],
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
                    "adc {self3}, 0",
                    self0 = inout(reg) self.0[0],
                    self1 = inout(reg) self.0[1],
                    self2 = inout(reg) self.0[2],
                    self3 = inout(reg) self.0[3],
                    k = in(reg) Self::K,
                );
            }
        } else if self.0[3] == -1i64 as u64 &&
                  self.0[2] == -1i64 as u64 &&
                  self.0[1] == -1i64 as u64 &&
                  self.0[0] >= Self::PRIME.0[0] {
            self.0[3] = 0;
            self.0[2] = 0;
            self.0[1] = 0;
            self.0[0] -= Self::PRIME.0[0];
        }
    }
}

impl Sub for P256 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        #[cfg(debug_assertions)] sub_call();
        let mut res = - other;
        res += self;
        res
    }
}

impl SubAssign for P256 {
    fn sub_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] sub_call();
        *self += - other;
    }
}

// TODO XXX
impl Mul for P256 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        #[cfg(debug_assertions)] mul_call();
        let a = self.to_u256();
        let b = other.to_u256();
        let c = u256::full_mul(a, b);

        let mut carry = 0u64;
        let mut res = Self([c.0[0], c.0[1], c.0[2], c.0[3]]);
        unsafe {
            std::arch::asm!(
                // xhigh * k + xlow = res + (carry << 256)
                // where res is =256 bits
                // and carry is <64 bits
                "mov rdx, {k}",
                "mulx {y1}, {y0}, {xh0}",
                "mulx {y5}, {y4}, {xh2}",
                "add {res0}, {y0}",
                "adc {res1}, {y1}",
                "adc {res2}, {y4}",
                "adc {res3}, {y5}",
                "adc {carry}, 0",
                "mulx {y5}, {y4}, {xh1}",
                "mulx {y1}, {y0}, {xh3}",
                "add {res1}, {y4}",
                "adc {res2}, {y5}",
                "adc {res3}, {y0}",
                "adc {carry}, {y1}",
                // res + carry * k => res + (carry << 256)
                "mulx {y1}, {y0}, {carry}",
                "mov {carry}, 0",
                "add {res0}, {y0}",
                "adc {res1}, {y1}",
                "adc {res2}, 0",
                "adc {res3}, 0",
                "adc {carry}, 0",
                out("rdx") _,
                k = in(reg) Self::K,
                xh0 = in(reg) c.0[4],
                xh1 = in(reg) c.0[5],
                xh2 = in(reg) c.0[6],
                xh3 = in(reg) c.0[7],
                y0 = out(reg) _,
                y1 = out(reg) _,
                // y2 = out(reg) _,
                // y3 = out(reg) _,
                y4 = out(reg) _,
                y5 = out(reg) _,
                // y6 = out(reg) _,
                // y7 = out(reg) _,
                res0 = inout(reg) res.0[0],
                res1 = inout(reg) res.0[1],
                res2 = inout(reg) res.0[2],
                res3 = inout(reg) res.0[3],
                carry = inout(reg) carry,
            );
        }

        if carry == 1 {
            // (2^256 + res) mod p = res + k
            unsafe {
                std::arch::asm!(
                    "add {self0}, {k}",
                    "adc {self1}, 0",
                    "adc {self2}, 0",
                    "adc {self3}, 0",
                    self0 = inout(reg) res.0[0],
                    self1 = inout(reg) res.0[1],
                    self2 = inout(reg) res.0[2],
                    self3 = inout(reg) res.0[3],
                    k = in(reg) Self::K,
                );
            }
        } else if res.0[3] == -1i64 as u64 &&
                  res.0[2] == -1i64 as u64 &&
                  res.0[1] == -1i64 as u64 &&
                  res.0[0] >= Self::PRIME.0[0] {
            res.0[3] = 0;
            res.0[2] = 0;
            res.0[1] = 0;
            res.0[0] -= Self::PRIME.0[0];
        }

        res
    }
}

impl MulAssign for P256 {
    fn mul_assign(&mut self, other: Self) {
        #[cfg(debug_assertions)] mul_call();
        *self = *self * other;
    }
}

impl Div for P256 {
    type Output = Self;

    // x/0 = 0
    // TODO
    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl DivAssign for P256 {
    // x/0 = 0
    // TODO
    fn div_assign(&mut self, other: Self) {
        *self *= other.inv();
    }
}

impl std::convert::From<GF2p8> for P256 {
    // TODO
    fn from(x: GF2p8) -> Self {
        Self([x.0 as u64, 0u64, 0u64, 0u64])
    }
}

impl std::fmt::Display for P256 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#018x}{:016x}{:016x}{:016x}", self.0[3], self.0[2], self.0[1], self.0[0])
    }
}

impl std::fmt::Debug for P256 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:#018x}{:016x}{:016x}{:016x}", self.0[3], self.0[2], self.0[1], self.0[0])
    }
}

impl std::convert::From<usize> for P256 {
    fn from(x: usize) -> Self {
        Self([x as u64, 0u64, 0u64, 0u64])
    }
}

impl std::convert::From<Block512> for P256 {
    fn from(x: Block512) -> Self {
        Self([x.0[0], x.0[1], x.0[2], x.0[3]])
    }
}

impl WriteWord for P256 {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize) {
        dst[*ptr] = self.0[0];
        dst[*ptr + 1] = self.0[1];
        dst[*ptr + 2] = self.0[2];
        dst[*ptr + 3] = self.0[3];
        *ptr += 4;
    }
}

impl std::convert::From<Vec<u64>> for P256 {
    fn from(x: Vec<u64>) -> Self {
        Self([x[0], x[1], x[2], x[3]])
    }
}

impl std::convert::Into<Vec<u64>> for P256 {
    fn into(self) -> Vec<u64> {
        self.0.to_vec()
    }
}


impl std::convert::AsRef<[u8]> for P256 {
    fn as_ref(&self) -> &[u8] {
        unsafe {
            std::mem::transmute::<*const [u64; 4], &[u8; 32]>
                (&self.0 as *const [u64; 4])
                as &[u8]
        }
    }
}

impl _Field for P256 {}

impl_pack_field!(P256, P4P256, 4);
impl_pack_field!(P256, P8P256, 8);
impl_pack_field!(P256, P16P256, 16);
impl_vect_field!(P256, VP256, 16);

impl FieldPack for P256 {
    type P4 =  P4P256;
    type P8 =  P8P256;
    type P16 = P16P256;
    type V16 = VP256;
}

impl Field for P256 {}
