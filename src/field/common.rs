use super::*;
use crate::hash::Block512;

#[derive(Copy, Clone)]
pub union M128 {
    pub m128i: __m128i,
    pub u8_array: [u8; 16],
    // u16_array: [u16; 8],
    pub u32_array: [u32; 4],
    pub u64_array: [u64; 2],
}

#[derive(Copy, Clone)]
pub union M256 {
    pub m256i: __m256i,
    pub u8_array: [u8; 32],
    pub u16_array: [u16; 16],
    pub u32_array: [u32; 8],
    pub u64_array: [u64; 4],
    pub m128i: [__m128i; 2],
}

#[derive(Copy, Clone)]
pub union M512 {
    pub m512i: __m512i,
    pub u8_array: [u8; 64],
    pub u16_array: [u16; 32],
    pub u32_array: [u32; 16],
    pub u64_array: [u64; 8],
    pub u128_array: [u128; 4],
    pub m128i: [__m128i; 4],
    pub m256i: [__m256i; 2],
}

pub trait Inv {
    fn inv(&self) -> Self;
}

pub trait Xor {
    #[inline]
    fn xor(x: Self, y: Self) -> Self
        where Self: Sized + Copy +
            MulAssign + Mul<Output = Self> +
            AddAssign + Add<Output = Self> +
            std::convert::From<usize> {
        Self::from(1) + x + y
    }
}

pub trait Pow {
    fn pow(x: Self, p: usize) -> Self
        where Self: Sized + Copy +
            MulAssign + Mul<Output = Self> +
            std::convert::From<usize> {
        let mut p = p;
        let mut res = Self::from(1);
        let mut tmp = x;
        while p != 0 {
            if p & 1 == 1 {
                res *= tmp
            }
            tmp = tmp * tmp;
            p >>= 1;
        }
        res
    }
}

pub trait WriteWord {
    fn write_word(&self, dst: &mut Vec<u64>, ptr: &mut usize);
}

pub trait Sqrt {
    fn sqrt(x: Self) -> Self
        where Self: Sized + Pow + Size + Copy +
            MulAssign + Mul<Output = Self> +
            std::convert::From<usize> {
        Self::pow(x, 1 << (Self::SIZE - 1))
    }
}

pub trait Random {
    type Output;
    fn random() -> Self::Output;
}

pub trait Size {
    const SIZE: usize;

    fn from_bytes(_: &[u8; Self::SIZE]) -> Self;
}

pub trait Characteristic {
    const CHAR: usize;
}

pub trait Pack<const DEGREE: usize> {
    type Input: _Field;
    type Output;

    // nightly build of rust only:
    // #![feature(generic_const_exprs)]
    //
    // for now, force DEGREE to be 4
    //
    // const DEGREE: usize = 4;
    // const DEGREE: usize;

    const SIZE: usize = DEGREE * Self::Input::SIZE;

    // fn pack(_: &[Self::Input; 4]) -> Self::Output;
    fn pack(src: &[Self::Input; DEGREE]) -> Self::Output;

    // fn unpack(&self) -> [Self::Input; 4];
    fn unpack(&self, dst: &mut [Self::Input; DEGREE]);
}

pub trait Vect {
    type Input: _Field;
    type Output;

    // nightly build of rust only:
    // #![feature(generic_const_exprs)]
    //
    // for now, force DEGREE to be 4
    //
    // const DEGREE: usize = 4;
    const DEGREE: usize;

    const SIZE: usize = Self::DEGREE * Self::Input::SIZE;

    // fn pack(_: &[Self::Input; 4]) -> Self::Output;
    fn new(src: &[u8; Self::SIZE]) -> Self::Output;

    // fn unpack(&self) -> [Self::Input; 4];
    fn unpack(&self, dst: *mut [u8; Self::SIZE]);
}

pub trait FieldPack where Self: _Field {
    // type P4;
    // type P8;
    // type P16;
    // type V16;
    type P4: PField<Self, 4>;
    type P8: PField<Self, 8>;
    type P16: PField<Self, 16>;
    type V16: VField<Self>;
}

pub trait _Field :
    Clone + Copy + Debug +
    PartialEq + Eq +
    'static + Send + Sync +
    Size + Characteristic +
    std::convert::From<usize> +
    std::convert::From<Block512> +
    std::convert::Into<Block512> +
    std::convert::From<Vec<u64>> +
    std::convert::Into<Vec<u64>> +
    std::convert::AsRef<[u8]> +
    WriteWord +
    AddAssign + SubAssign + MulAssign + DivAssign +
    Neg<Output = Self> +
    Add<Output = Self> +
    Sub<Output = Self> +
    Mul<Output = Self> +
    Div<Output = Self> +
    Inv + Pow + Sqrt + Xor +
    Random<Output = Self> +
    std::fmt::Display {}

pub trait Field : _Field + FieldPack {}

pub trait PField<T: _Field, const DEGREE: usize> :
    Clone + Copy + Debug +
    Send + Sync +
    Pack<DEGREE, Input = T, Output = Self> +
    AddAssign + SubAssign + MulAssign +
    Add<Output = Self> +
    Sub<Output = Self> +
    Mul<Output = Self> {}

pub trait VField<T: _Field> :
    Clone + Copy +
    // Debug +
    Vect<Input = T, Output = Self> +
    // std::convert::From<usize> +
    // std::convert::AsRef<[u8]> +
    AddAssign + SubAssign +
    Add<Output = Self> +
    Sub<Output = Self> +
    Mul<usize, Output = Self> +
    Mul<T, Output = Self> +
    MulAssign<usize> +
    MulAssign<T> {}
    // std::fmt::Display {}


pub fn clmul(a: u64, b: u64) -> u64 {
    return unsafe {
        _mm_cvtsi128_si64(
            _mm_clmulepi64_si128(
                _mm_set_epi64x(0, a as i64),
                _mm_set_epi64x(0, b as i64), 0)
            ) as u64
    }
}


pub fn clmul_u128(a: u128, b: u128) -> [u128; 2] {
    return unsafe {
        let x = _mm_set_epi64x((a >> 64) as i64, a as i64);
        let y = _mm_set_epi64x((b >> 64) as i64, b as i64);

        let hh = _mm_clmulepi64_si128(x, y, 17);
        let ll = _mm_clmulepi64_si128(x, y, 0);
        let hl = _mm_xor_si128(
            _mm_clmulepi64_si128(x, y, 16),
            _mm_clmulepi64_si128(x, y, 1));

        let hi = ((_mm_extract_epi64(hh, 1) as u64 as u128) << 64) ^ (
            (_mm_extract_epi64(hh, 0) as u64 as u128) ^
            (_mm_extract_epi64(hl, 1) as u64 as u128));
        let lo = (_mm_extract_epi64(ll, 0) as u64 as u128) ^ ((
            (_mm_extract_epi64(ll, 1) as u64 as u128) ^
            (_mm_extract_epi64(hl, 0) as u64 as u128)) << 64);
        [lo, hi]
    }
}


pub fn clmul_u128_lo(a: u128, b: u128) -> u128 {
    return unsafe {
        let x = _mm_set_epi64x((a >> 64) as i64, (a & ((1u128 << 64) - 1)) as i64);
        let y = _mm_set_epi64x((b >> 64) as i64, (b & ((1u128 << 64) - 1)) as i64);
        let ll = _mm_clmulepi64_si128(x, y, 0);
        let hl = _mm_xor_si128(
            _mm_clmulepi64_si128(x, y, 16),
            _mm_clmulepi64_si128(x, y, 1));

        let lo = (_mm_extract_epi64(ll, 0) as u64 as u128) ^ ((
            (_mm_extract_epi64(ll, 1) as u64 as u128) ^
            (_mm_extract_epi64(hl, 0) as u64 as u128)) << 64);
        lo
    }
}

// pub fn clmul_u256(a: u256, b: u256) -> [u256; 2] {
//     let a_lo = a.low_u128();
//     let a_hi = (a >> 128).as_u128();
//     let b_lo = b.low_u128();
//     let b_hi = (b >> 128).as_u128();
//
//     let ll = clmul_u128(a_lo, b_lo);
//     let lh = clmul_u128(a_lo, b_hi);
//     let hl = clmul_u128(a_hi, b_lo);
//     let hh = clmul_u128(a_hi, b_hi);
//
//     let lo: u256 = u256::from(ll[0]) ^ (u256::from(lh[0] ^ hl[0] ^ ll[1]) << 128);
//     let hi: u256 = (u256::from(hh[1]) << 128) ^ u256::from(lh[1] ^ hl[1] ^ hh[0]);
//
//     [lo, hi]
// }


pub fn clmul_u256_lo(a: u256, b: u256) -> u256 {
    let a_lo = a.low_u128();
    let a_hi = (a >> 128).as_u128();
    let b_lo = b.low_u128();
    let b_hi = (b >> 128).as_u128();

    let ll = clmul_u128(a_lo, b_lo);
    let lh = clmul_u128(a_lo, b_hi);
    let hl = clmul_u128(a_hi, b_lo);

    let lo: u256 = u256::from(ll[0]) ^ (u256::from(lh[0] ^ hl[0] ^ ll[1]) << 128);

    lo
}

pub fn div_quotient(a: u64, b: u64) -> u64 {
    let mut a = a;
    let mut q: u64 = 0;
    let mut diff: i64 = (b.leading_zeros() as i64) - (a.leading_zeros() as i64);

    while diff >= 0 && a != 0 {
        q |= (1u128 << diff) as u64;
        a ^= (b << diff) as u64;
        diff = (b.leading_zeros() as i64) - (a.leading_zeros() as i64);
    }
    q
}

pub fn mod_inverse(a: u64, m: u64) -> u64 {
    let mut old_t: u64 = 0;
    let mut t: u64 = 1;
    let mut old_r = m;
    let mut r = a;

    while r != 0 {
        let q = div_quotient(old_r, r);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul(q, old_t);
    }

    old_t
}

pub fn mod_inverse_u64(a: u64, m: u64) -> u64 {
    if a == 1 {
        return 1;
    }

    if a == 0 {
        return 0;
    }

    let mut old_t: u64 = 0;
    let mut t: u64 = 1;
    let mut old_r = m;
    let mut r = a;

    {
        let q = div_quotient_u64(old_r, r, 1);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul(q, old_t);
    }

    while r != 0 {
        let q = div_quotient_u64(old_r, r, 0);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul(q, old_t);
    }

    old_t
}

pub fn div_quotient_u64(a: u64, b: u64, flag: usize) -> u64 {
    let mut a = a;
    let mut q: u64 = 0;
    let mut diff: i64 =
        if flag == 1 {
            (b.leading_zeros() as i64) + 1
        } else {
            (b.leading_zeros() as i64) - (a.leading_zeros() as i64)
        };

    while diff >= 0 && a != 0 {
        q |= 1u64 << diff;
        a ^= b << diff;
        diff = (b.leading_zeros() as i64) - (a.leading_zeros() as i64);
    }
    q
}


pub fn div_quotient_u128(a: u128, b: u128, flag: usize) -> u128 {
    let mut a = a;
    let mut q: u128 = 0;
    let mut diff: i64 =
        if flag == 1 {
            (b.leading_zeros() as i64) + 1
        } else {
            (b.leading_zeros() as i64) - (a.leading_zeros() as i64)
        };

    while diff >= 0 && a != 0 {
        q |= 1u128 << diff;
        a ^= b << diff;
        diff = (b.leading_zeros() as i64) - (a.leading_zeros() as i64);
    }
    q
}


pub fn mod_inverse_u128(a: u128, m: u128) -> u128 {
    if a == 1 {
        return 1;
    }

    if a == 0 {
        return 0;
    }

    let mut old_t: u128 = 0;
    let mut t: u128 = 1;
    let mut old_r = m;
    let mut r = a;

    {
        let q = div_quotient_u128(old_r, r, 1);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul_u128_lo(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul_u128_lo(q, old_t);
    }

    while r != 0 {
        let q = div_quotient_u128(old_r, r, 0);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul_u128_lo(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul_u128_lo(q, old_t);
    }

    old_t
}

pub fn div_quotient_u192(a: u256, b: u256, flag: usize) -> u256 {
    let mut a = a;
    let mut q: u256 = u256::from(0u8);
    let mut diff: i64 =
        if flag == 1 {
            (b.leading_zeros() as i64) + 1
        } else {
            (b.leading_zeros() as i64) - (a.leading_zeros() as i64)
        };

    while diff >= 0 && a != u256::from(0u8) {
        q |= u256::from(1u8) << diff;
        a ^= b << diff;
        diff = (b.leading_zeros() as i64) - (a.leading_zeros() as i64);
    }
    q
}


pub fn mod_inverse_u192(a: u256, m: u256) -> u256 {
    if a == u256::from(1u8) {
        return u256::from(1u8);
    }

    if a == u256::from(0u8) {
        return u256::from(0u8);
    }

    let mut old_t: u256 = u256::from(0u8);
    let mut t: u256 = u256::from(1u8);
    let mut old_r = m;
    let mut r = a;

    {
        let q = div_quotient_u256(old_r, r, 0);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul_u256_lo(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul_u256_lo(q, old_t);
    }

    while r != u256::from(0u8) {
        let q = div_quotient_u256(old_r, r, 0);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul_u256_lo(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul_u256_lo(q, old_t);
    }

    old_t
}


pub fn div_quotient_u256(a: u256, b: u256, flag: usize) -> u256 {
    let mut a = a;
    let mut q: u256 = u256::from(0u8);
    let mut diff: i64 =
        if flag == 1 {
            (b.leading_zeros() as i64) + 1
        } else {
            (b.leading_zeros() as i64) - (a.leading_zeros() as i64)
        };

    while diff >= 0 && a != u256::from(0u8) {
        q |= u256::from(1u8) << diff;
        a ^= b << diff;
        diff = (b.leading_zeros() as i64) - (a.leading_zeros() as i64);
    }
    q
}


pub fn mod_inverse_u256(a: u256, m: u256) -> u256 {
    if a == u256::from(1u8) {
        return u256::from(1u8);
    }

    if a == u256::from(0u8) {
        return u256::from(0u8);
    }

    let mut old_t: u256 = u256::from(0u8);
    let mut t: u256 = u256::from(1u8);
    let mut old_r = m;
    let mut r = a;

    {
        let q = div_quotient_u256(old_r, r, 1);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul_u256_lo(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul_u256_lo(q, old_t);
    }

    while r != u256::from(0u8) {
        let q = div_quotient_u256(old_r, r, 0);

        let mut tmp = old_r;
        old_r = r;
        r = tmp ^ clmul_u256_lo(q, old_r);

        tmp = old_t;
        old_t = t;
        t = tmp ^ clmul_u256_lo(q, old_t);
    }

    old_t
}

// pub fn mulx_u256(a: u256, b: u256) -> (u256, u256) {
//     let a_lo = u256::from(a.low_u128());
//     let a_hi = a >> 128;
//     let b_lo = u256::from(b.low_u128());
//     let b_hi = b >> 128;
// 
//     let lo = a_lo * b_lo;
//     let mi = a_lo * b_hi + a_hi * b_lo;
//     let hi = a_hi * b_hi;
// 
//     let mi_lo = mi << 128;
//     let mut mi_hi = mi >> 128;
// 
//     let lo = lo + mi_lo;
//     if lo < mi_lo {
//         mi_hi += u256::from(1);
//     }
//     let hi = hi + mi_hi;
//     (lo, hi)
// }

#[macro_export]
macro_rules! impl_pack_field {
    ($U:tt, $T:tt, $degree:tt) => {
        #[derive(Copy, Clone, Debug)]
        pub struct $T(pub [$U; $degree]);

        impl $T {}

        impl Pack<$degree> for $T {
            type Output = Self;
            type Input = $U;

            // const DEGREE: usize = $degree;

            fn pack(raw: &[$U; $degree]) -> Self {
                Self(*raw)
            }

            fn unpack(&self, dst: &mut [$U; $degree]) {
                *dst = self.0;
            }
        }

        impl Add for $T {
            type Output = Self;

            fn add(self, other: Self) -> Self {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
                let mut res = self;
                for i in 0..$degree {
                    res.0[i] = self.0[i] + other.0[i];
                }
                res
            }
        }

        impl AddAssign for $T {
            fn add_assign(&mut self, other: Self) {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
                for i in 0..$degree {
                    self.0[i] += other.0[i];
                }
            }
        }

        impl Sub for $T {
            type Output = Self;

            fn sub(self, other: Self) -> Self {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { sub_call(); }
                let mut res = self;
                for i in 0..$degree {
                    res.0[i] = self.0[i] - other.0[i];
                }
                res
            }
        }

        impl SubAssign for $T {
            fn sub_assign(&mut self, other: Self) {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { sub_call(); }
                for i in 0..$degree {
                    self.0[i] -= other.0[i];
                }
            }
        }

        impl Mul for $T {
            type Output = Self;

            fn mul(self, other: Self) -> Self {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
                let mut res = self;
                for i in 0..$degree {
                    res.0[i] = self.0[i] * other.0[i];
                }
                res
            }
        }

        impl MulAssign for $T {
            fn mul_assign(&mut self, other: Self) {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
                for i in 0..$degree {
                    self.0[i] *= other.0[i];
                }
            }
        }

        impl PField<$U, $degree> for $T {}
    }
}

#[macro_export]
macro_rules! impl_vect_field {
    ($U:tt, $T:tt, $degree:tt) => {
        #[derive(Copy, Clone)]
        pub struct $T(pub [$U; $degree]);

        impl $T {}

        impl Vect for $T {
            type Output = Self;
            type Input = $U;

            const DEGREE: usize = $degree;

            fn new(raw: &[u8; Self::SIZE]) -> Self {
                unsafe {
                    Self(*std::mem::transmute::<&[u8; Self::SIZE], &[$U; $degree]>(raw))
                }
            }

            fn unpack(&self, dst: *mut [u8; Self::SIZE]) {
                unsafe {
                    let ptr: *mut [$U; $degree] = std::mem::transmute(dst);
                    *ptr = self.0;
                }
            }
        }

        impl Add for $T {
            type Output = Self;

            fn add(self, other: Self) -> Self {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
                let mut res = self;
                for i in 0..$degree {
                    res.0[i] = self.0[i] + other.0[i];
                }
                res
            }
        }

        impl AddAssign for $T {
            fn add_assign(&mut self, other: Self) {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { add_call(); }
                for i in 0..$degree {
                    self.0[i] += other.0[i];
                }
            }
        }

        impl Sub for $T {
            type Output = Self;

            fn sub(self, other: Self) -> Self {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { sub_call(); }
                let mut res = self;
                for i in 0..$degree {
                    res.0[i] = self.0[i] - other.0[i];
                }
                res
            }
        }

        impl SubAssign for $T {
            fn sub_assign(&mut self, other: Self) {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { sub_call(); }
                for i in 0..$degree {
                    self.0[i] -= other.0[i];
                }
            }
        }

        impl Mul<$U> for $T {
            type Output = Self;

            fn mul(self, other: $U) -> Self {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
                let mut res = self;
                for i in 0..$degree {
                    res.0[i] = self.0[i] * other;
                }
                res
            }
        }

        impl MulAssign<$U> for $T {
            fn mul_assign(&mut self, other: $U) {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
                for i in 0..$degree {
                    self.0[i] *= other;
                }
            }
        }

        impl Mul<usize> for $T {
            type Output = Self;

            fn mul(self, other: usize) -> Self {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
                let mut res = self;
                for i in 0..$degree {
                    res.0[i] = self.0[i] * $U::from(other);
                }
                res
            }
        }

        impl MulAssign<usize> for $T {
            fn mul_assign(&mut self, other: usize) {
                #[cfg(debug_assertions)] for _ in 0..Self::DEGREE { mul_call(); }
                for i in 0..$degree {
                    self.0[i] *= $U::from(other);
                }
            }
        }

        impl VField<$U> for $T {}
    }
}

pub(crate) use impl_pack_field;
pub(crate) use impl_vect_field;
