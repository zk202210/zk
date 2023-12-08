use crate::field::Field;

#[inline]
pub fn ceil_log(x: usize) -> usize {
    let mut res: usize = 63 - ((x as u64).leading_zeros() as usize);
    if res != ((x as u64).trailing_zeros() as usize) {
        res += 1;
    }
    res
}

// https://stackoverflow.com/questions/2602823/in-c-c-whats-the-simplest-way-to-reverse-the-order-of-bits-in-a-byte
pub fn swap_bits(b: usize, m: usize) -> usize {
    // assert!(b < (1 << 32));
    let b = (b & 0xFFFF0000) >> 16 | (b & 0x0000FFFF) << 16;
    let b = (b & 0xFF00FF00) >> 8 | (b & 0x00FF00FF) << 8;
    let b = (b & 0xF0F0F0F0) >> 4 | (b & 0x0F0F0F0F) << 4;
    let b = (b & 0xCCCCCCCC) >> 2 | (b & 0x33333333) << 2;
    let b = (b & 0xAAAAAAAA) >> 1 | (b & 0x55555555) << 1;
    b >> (32 - m)

    // let mut res: usize = 0;
    // let mut k = k;
    // for _ in 0..m {
    //     res <<= 1;
    //     res += k & 1;
    //     k >>= 1;
    // }
    // res
}

pub fn expand_binary<T: Field>(x: usize, bit_len: usize) -> Vec<T> {
    let mut res = Vec::with_capacity(bit_len);
    let mut x = x;
    for _ in 0..bit_len {
        res.push(T::from(x & 1));
        x >>= 1;
    }
    res
}

pub fn expand_binary_eq<T: Field>(x: usize, tau: &Vec<T>) -> T {
    let mut res = T::from(1);
    let mut x = x;
    for t in tau.iter() {
        res *= if (x & 1) == 1 { *t } else { T::from(1) - *t };
        x >>= 1;
    }
    res
}

pub fn expand_binary_eq2<T: Field>(x: usize, tau: &Vec<T>, rho: &Vec<T>) -> T {
    let mut res = T::from(1);
    let mut x = x;
    for i in 0..tau.len() {
        res *=
            if (x & 1) == 1 {
                one(tau[i]) * one(rho[i])
            } else {
                zer(tau[i]) * zer(rho[i])
            };
        x >>= 1;
    }
    res
}

#[inline]
pub fn one<T: Field>(a: T) -> T { a }

#[inline]
pub fn zer<T: Field>(a: T) -> T { T::from(1) - a }

#[inline]
pub fn eq3<T: Field>(a: T, b: T, c: T) -> T {
    zer(a) * zer(b) * zer(c) + one(a) * one(b) * one(c)
}

// #[inline]
// fn neq3<T: Field>(a: T, b: T, c: T) -> T {
//     zer(a) * one(b) * one(c) + one(a) * zer(b) * zer(c)
// }

#[inline]
pub fn eq2<T: Field>(a: T, b: T) -> T {
    T::xor(a, b)
}

#[inline]
pub fn neq2<T: Field>(a: T, b: T) -> T {
    // char-2 field only
    a + b
}


// gkr security: (1/16) ** 33 + ... ≈ 2^{-132}
// RS code rate:
//      1/32 for original inputs
//      1/16 for l(x) * q(x)
// target security: 128 bit
pub fn get_n_ldt_iter(log_input: usize, field_size: usize, rho: usize) -> usize {
    // field_size in bits

    // first we need to get ~129 bit security with εc. ref: 2020/654 p40
    let m7 = field_size - 6 + 1 - (log_input + 4) * 2 - 129;
    let m = (m7 as f64) / 7.0;
    let m = f64::powf(2.0, m) - 0.5;
    let m = m.floor();
    println!("m: {}", m);
    assert!(m > 2.0);

    let n_ldt_iter = 129.0 / ((rho as f64).sqrt().log2() - (1.0 + 1.0 / (2.0 * m)).log2());
    n_ldt_iter.ceil() as usize
}
