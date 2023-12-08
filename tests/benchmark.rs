#![allow(incomplete_features)]
#![allow(dead_code)]
#![feature(generic_const_exprs)]

use rand::Rng;
use std::time::Instant;
use std::ops::*;

use rustzk::field::*;
use rustzk::hash::*;
use rustzk::field_env;

const LOG_LEN: i32 = 25;

fn type_name<'a, T>() -> &'a str {
    let s = std::any::type_name::<T>();
    let parts: Vec<&str> = s.rsplit(":").collect();
    parts[0]
}

fn add<T: Add<U, Output=T>, U>(a: T, b: U) -> T {
    a + b
}

fn mul<T: Mul<U, Output=T>, U>(a: T, b: U) -> T {
    a * b
}

fn div<T: Div<U, Output=T>, U>(a: T, b: U) -> T {
    a / b
}

// fn add_assign<T: AddAssign<U>, U>(a: &mut T, b: U) {
//     *a += b;
// }
// 
// fn mul_assign<T: MulAssign<U>, U>(a: &mut T, b: U) {
//     *a *= b;
// }
// 
// fn div_assign<T: DivAssign<U>, U>(a: &mut T, b: U) {
//     *a /= b;
// }


/*
fn bench_vect_add<U, T>()
    where
        U: Field,
        T: VField<U>,
        [(); T::SIZE]:
{
    let mut a: Vec<T> = Vec::new();
    let mut b: Vec<T> = Vec::new();
    let mut c: Vec<T> = Vec::new();

    for _ in 0..((1 << LOG_LEN) / T::DEGREE) {
        let mut tmp_a: [u8; T::SIZE] = [0u8; T::SIZE];
        let mut tmp_b: [u8; T::SIZE] = [0u8; T::SIZE];
        let mut tmp_c: [u8; T::SIZE] = [0u8; T::SIZE];
        for i in 0..T::DEGREE {
            tmp_a[i] = rand::thread_rng().gen();
            tmp_b[i] = rand::thread_rng().gen();
        }
        a.push(T::new(&tmp_a));
        b.push(T::new(&tmp_b));
        c.push(T::new(&tmp_c));
    }

    let start = Instant::now();
    for i in 0..((1 << LOG_LEN) / T::DEGREE) {
        c[i] = a[i] + b[i];
    }
    let elapsed = start.elapsed();
    let in_ms = elapsed.as_millis();
    println!("{:<10} vect {}: {:>18.15} ns", type_name::<U>(), "add",
             in_ms as f64 * 1000000.0 / f64::powi(2.0, LOG_LEN));

    let mut d: Vec<[u8; T::SIZE]> = vec![[0u8; T::SIZE]; (1 << LOG_LEN) / T::DEGREE];
    let start = Instant::now();
    for i in 0..((1 << LOG_LEN) / T::DEGREE) {
        #[allow(unused_assignments)]
        let mut c = T::new(&d[i]);
        c = a[i] + b[i];
        c.unpack(&mut d[i]);
    }
    let elapsed = start.elapsed();
    let in_ms = elapsed.as_millis();
    println!("{:<10} vect {}: {:>18.15} ns", type_name::<U>(), "a&c",
             in_ms as f64 * 1000000.0 / f64::powi(2.0, LOG_LEN));
}

fn bench_vect_mul<U, T>()
    where
        U: Field,
        T: VField<U>,
        [(); T::SIZE]:
{
    let mut a: Vec<T> = Vec::new();
    let mut b: Vec<U> = Vec::new();
    let mut c: Vec<T> = Vec::new();

    for _ in 0..((1 << LOG_LEN) / T::DEGREE) {
        let mut tmp_a: [u8; T::SIZE] = [0u8; T::SIZE];
        let tmp_b = U::random();
        let mut tmp_c: [u8; T::SIZE] = [0u8; T::SIZE];
        for i in 0..T::DEGREE {
            tmp_a[i] = rand::thread_rng().gen();
        }
        a.push(T::new(&tmp_a));
        b.push(tmp_b);
        c.push(T::new(&tmp_c));
    }

    let start = Instant::now();
    for i in 0..((1 << LOG_LEN) / T::DEGREE) {
        c[i] = a[i] * b[i];
    }
    let elapsed = start.elapsed();
    let in_ms = elapsed.as_millis();
    println!("{:<10} vect {}: {:>18.15} ns", type_name::<U>(), "mul",
             in_ms as f64 * 1000000.0 / f64::powi(2.0, LOG_LEN));

    let mut d: Vec<[u8; T::SIZE]> = vec![[0u8; T::SIZE]; (1 << LOG_LEN) / T::DEGREE];
    let start = Instant::now();
    for i in 0..((1 << LOG_LEN) / T::DEGREE) {
        #[allow(unused_assignments)]
        let mut c = T::new(&d[i]);
        c = a[i] * b[i];
        c.unpack(&mut d[i]);
    }
    let elapsed = start.elapsed();
    let in_ms = elapsed.as_millis();
    println!("{:<10} vect {}: {:>18.15} ns", type_name::<U>(), "m&c",
             in_ms as f64 * 1000000.0 / f64::powi(2.0, LOG_LEN));
}


fn bench_pack_ops<U, T, F>(op: F, name: &str)
    where
        U: Field,
        T: PField<U>,
        F: Fn(T, T) -> T,
        [(); T::DEGREE]:
{
    let mut a: Vec<T> = Vec::new();
    let mut b: Vec<T> = Vec::new();
    let mut c: Vec<T> = Vec::new();

    for _ in 0..((1 << LOG_LEN) / T::DEGREE) {
        // NIGHTLY, replace 4 with T::DEGREE
        let mut tmp_a: [U; T::DEGREE] = [U::from(0); T::DEGREE];
        let mut tmp_b: [U; T::DEGREE] = [U::from(0); T::DEGREE];
        let mut tmp_c: [U; T::DEGREE] = [U::from(0); T::DEGREE];
        for i in 0..T::DEGREE {
            tmp_a[i] = U::random();
            tmp_b[i] = U::random();
        }
        a.push(T::pack(&tmp_a));
        b.push(T::pack(&tmp_b));
        c.push(T::pack(&tmp_c));
    }

    let mut d: Vec<[U; T::DEGREE]> = vec![[U::from(0); T::DEGREE]; (1 << LOG_LEN) / T::DEGREE];

    let start = Instant::now();
    for i in 0..((1 << LOG_LEN) / T::DEGREE) {
        // c[i] = op(a[i], b[i]);
        c[i] = op(a[i], b[i]);
    }
    let elapsed = start.elapsed();
    let in_ms = elapsed.as_millis();
    println!("{:<10} pack 1/{:<2} {:<5}: {:>18.15} ns", type_name::<T>(), T::DEGREE, name,
             in_ms as f64 * 1000000.0 / f64::powi(2.0, LOG_LEN));

    let start = Instant::now();
    for i in 0..((1 << LOG_LEN) / T::DEGREE) {
        #[allow(unused_assignments)]
        let mut c = T::pack(&d[i]);
        c = op(a[i], b[i]);
        c.unpack(&mut d[i]);
    }
    let elapsed = start.elapsed();
    let in_ms = elapsed.as_millis();
    println!("{:<10} pack 1/{:<2} {}&c: {:>18.15} ns", type_name::<T>(), T::DEGREE, name,
             in_ms as f64 * 1000000.0 / f64::powi(2.0, LOG_LEN));
}

fn bench_pack_add<U: Field, T: PField<U>>()
    where [(); T::DEGREE]:
{
    bench_pack_ops::<U, T, _>(add::<T, T>, "add");
}

fn bench_pack_mul<U: Field, T: PField<U>>()
    where [(); T::DEGREE]:
{
    bench_pack_ops::<U, T, _>(mul::<T, T>, "mul");
}
*/

fn bench_ops<T: Field, F: Fn(T, T) -> T>(op: F, name: &str) {
    let mut a: Vec<T> = Vec::new();

    let mut i = 0;
    while i < (1 << LOG_LEN) {
        let tmp = T::random();
        if tmp != T::from(0) {
            i += 1;
            a.push(tmp);
        }
    }

    let mut b: Vec<T> = Vec::new();

    let mut i = 0;
    while i < (1 << LOG_LEN) {
        let tmp = T::random();
        if tmp != T::from(0) {
            i += 1;
            b.push(tmp);
        }
    }

    let mut c: Vec<T> = vec![T::from(0); 1 << LOG_LEN];

    let start = Instant::now();
    for i in 0..(1 << LOG_LEN) {
        c[i] = op(a[i], b[i]);
    }
    let elapsed = start.elapsed();
    let millis = elapsed.as_millis();
    println!("{:<10} {}: {:>18.15} ns", type_name::<T>(), name,
             millis as f64 * 1000000.0 / f64::powi(2.0, LOG_LEN));
}

fn bench_add<T: Field>() {
    bench_ops::<T, _>(add::<T, T>, "add");
}

fn bench_mul<T: Field>() {
    bench_ops::<T, _>(mul::<T, T>, "mul");
}

fn bench_div<T: Field>() {
    bench_ops::<T, _>(div::<T, T>, "div");
}


/*
fn bench_pack_add_t<T, P4, P8, P16>()
    where
        T: Field,
        P4: PField<T>,
        P8: PField<T>,
        P16: PField<T>,
        [(); P4::DEGREE]:,
        [(); P8::DEGREE]:,
        [(); P16::DEGREE]:
{
    bench_pack_add::<T, P4>();
    bench_pack_add::<T, P8>();
    bench_pack_add::<T, P16>();
}

fn bench_pack_mul_t<T, P4, P8, P16>()
    where
        T: Field,
        P4: PField<T>,
        P8: PField<T>,
        P16: PField<T>,
        [(); P4::DEGREE]:,
        [(); P8::DEGREE]:,
        [(); P16::DEGREE]:
{
    bench_pack_mul::<T, P4>();
    bench_pack_mul::<T, P8>();
    bench_pack_mul::<T, P16>();
}
*/


fn bench_add_all(name: &str) {
    println!("{}", name);
    println!("-------------------------------------------------");
    bench_add::<GF2p16>();
    bench_add::<GF2p32>();
    bench_add::<GF2p64>();
    bench_add::<GF2p128>();
    bench_add::<GF2p192>();
    bench_add::<GF2p256>();
    bench_add::<M61p2>();
    bench_add::<M127p2>();
    bench_add::<P192>();
    println!("-------------------------------------------------");
    println!("");
}

fn bench_mul_all(name: &str) {
    println!("{}", name);
    println!("-------------------------------------------------");
    bench_mul::<GF2p16>();
    bench_mul::<GF2p32>();
    bench_mul::<GF2p64>();
    bench_mul::<GF2p128>();
    bench_mul::<GF2p192>();
    bench_mul::<GF2p256>();
    bench_mul::<M61p2>();
    bench_mul::<M127p2>();
    bench_mul::<P192>();
    println!("-------------------------------------------------");
    println!("");
}

fn bench_div_all(name: &str) {
    println!("{}", name);
    println!("-------------------------------------------------");
    bench_div::<GF2p16>();
    bench_div::<GF2p32>();
    bench_div::<GF2p64>();
    bench_div::<GF2p128>();
    bench_div::<GF2p192>();
    bench_div::<GF2p256>();
    bench_div::<M61p2>();
    bench_div::<M127p2>();
    bench_div::<P192>();
    println!("-------------------------------------------------");
    println!("");
}

/*
fn bench_pack_add_all(name: &str) {
    println!("{}", name);
    println!("-------------------------------------------------");
    field_env! { bench_pack_add_t::<GF2p16>() };
    println!("");
    field_env! { bench_pack_add_t::<GF2p32>() };
    println!("");
    field_env! { bench_pack_add_t::<GF2p64>() };
    println!("");
    field_env! { bench_pack_add_t::<GF2p128>() };
    println!("");
    field_env! { bench_pack_add_t::<GF2p256>() };
    println!("");
    field_env! { bench_pack_add_t::<M61p2>() };
    println!("-------------------------------------------------");
    println!("");
}

fn bench_pack_mul_all(name: &str) {
    println!("{}", name);
    println!("-------------------------------------------------");
    field_env! { bench_pack_mul_t::<GF2p16>() };
    println!("");
    field_env! { bench_pack_mul_t::<GF2p32>() };
    println!("");
    field_env! { bench_pack_mul_t::<GF2p64>() };
    println!("");
    field_env! { bench_pack_mul_t::<GF2p128>() };
    println!("");
    field_env! { bench_pack_mul_t::<GF2p256>() };
    println!("");
    field_env! { bench_pack_mul_t::<M61p2>() };
    println!("-------------------------------------------------");
    println!("");
}

fn bench_vect_add_all(name: &str) {
    println!("{}", name);
    println!("-------------------------------------------------");
    bench_vect_add::<GF2p16, VGF2p16>();
    println!("-------------------------------------------------");
    println!("");
}

fn bench_vect_mul_all(name: &str) {
    println!("{}", name);
    println!("-------------------------------------------------");
    bench_vect_mul::<GF2p16, VGF2p16>();
    println!("-------------------------------------------------");
    println!("");
}
*/

fn bench_sha(name: &str) {
    println!("{}", name);
    println!("-------------------------------------------------");
    let mut a: Vec<GF2p128> = Vec::new();
    for _ in 0..(1 << LOG_LEN) {
        a.push(GF2p128::random());
    }

    let mut h: Vec<SHA256> = vec![SHA256::from(0); 1 << LOG_LEN];

    let start = Instant::now();
    for i in 0..(1 << LOG_LEN) {
        h[i].commit(&a[i..i+1]);
    }
    let elapsed = start.elapsed();
    let in_ms = elapsed.as_millis();
    println!("SHA256: {:?} ns", in_ms as f64 * 1000000.0 / f64::powi(2.0, LOG_LEN));
    println!("-------------------------------------------------");
    println!("");
}


fn run<F: Fn(&str)>(f: F) {
    f(type_name::<F>());
}


#[test]
fn bench_all() {
    println!("");
    println!("");

    run(bench_add_all);
    run(bench_mul_all);
    // bench_div_all();
    // run(bench_pack_add_all);
    // run(bench_pack_mul_all);
    // run(bench_vect_add_all);
    // run(bench_vect_mul_all);

    // run(bench_sha);
}
