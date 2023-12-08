use primitive_types::U256 as u256;

use rustzk::field::*;

fn pow_u256<T: Field>(x: T, p: u256) -> T {
    let mut p = p;
    let mut res = T::from(1);
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

#[test]
fn find_root() {
    println!("");
    type T = GF2p192;

    let max_size = (u256::from(1) << 192) - u256::from(1);
    let pow255 = max_size / u256::from(255);

    let mut i = 2;
    let r = loop {
        let r = pow_u256(T::from(i), pow255);
        if r != T::from(1) &&
            T::pow(r, 85) != T::from(1) &&
            T::pow(r, 3) != T::from(1) &&
            T::pow(r, 5) != T::from(1) &&
            T::pow(r, 17) != T::from(1) {
                println!("found root @{}: {}", i, r);
                break r;
        }
        i += 1;
    };

    let mut gen = T::from(1);
    for i in 0..255 {
        gen *= r;
        if T::pow(gen, 8) + T::pow(gen, 4) + T::pow(gen, 3) + gen == T::from(1) {
            println!("found gen=r^{:<3}: {}", i + 1, gen);
        }
    }
}
