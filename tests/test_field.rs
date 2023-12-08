use rand::Rng;
use rustzk::field::*;
use rustzk::poly::*;

fn test_field_arith<T: Field>() {
    test_field_add_sub::<T>();
    test_field_mul_add::<T>();
    test_field_mul_div::<T>();
    test_field_neg::<T>();
}

fn test_field_neg<T: Field>() {
    for _ in 0..1000 {
        let a = T::random();
        assert_eq!(- a, T::from(0) - a);
    }
}

fn test_field_mul_add<T: Field>() {
    for _ in 0..1000 {
        let a0 = T::random();
        let a1 = T::random();
        let b: T = loop {
            let tmp = T::random();
            if tmp != T::from(0) {
                break tmp;
            }
        };

        let c = (a0 + a1) * b;
        let c0 = a0 * b;
        let c1 = a1 * b;

        assert_eq!(c - c0, c1);

    }
}

fn test_field_add_sub<T: Field>() {
    for _ in 0..1000 {
        let a = T::random();
        let b: T = loop {
            let tmp = T::random();
            if tmp != T::from(0) {
                break tmp;
            }
        };

        let c = a + b;

        assert_eq!(c - a, b);

    }
}

fn test_field_mul_div<T: Field>() {
    let a = - T::from(1);
    let b = - T::from(2);
    let c = - T::from(0);
    assert_eq!(a * a, T::from(1));
    assert_eq!(b * a, T::from(2));
    assert_eq!(c * c, T::from(0));
    assert_eq!(a * c, T::from(0));
    assert_eq!(b * c, T::from(0));
    for _ in 0..1000 {
        let a = T::random();
        let b: T = loop {
            let tmp = T::random();
            if tmp != T::from(0) {
                break tmp;
            }
        };

        let c = a / b;

        assert_eq!(c * b, a);

    }
}

// #[test]
// fn test_fft_bit_reversal() {
//     type T = M61p2;
// 
//     let log_len = 4;
//     let log_k = 2;
//     let log_l = log_len + log_k;
//     let coset = Coset::<T>::init(log_l, log_k, 1);
// 
//     let mut values = vec![];
//     for i in 0..(1 << log_len) {
//         values.push(T::from(i))
//     }
// 
//     let coeff = coset.multiplicative_ifft(&values, log_len);
//     let values_ex = coset.multiplicative_fft(&coeff, log_len, log_l);
//     println!("{:?}", values_ex);
// 
//     let values_ex2 = coset.interpolate(&values, log_len, log_l);
// 
//     assert_eq!(values_ex, values_ex2);
// }


fn test_afft<T: Field>() {
    let log_len = 4;
    let log_k = 1;
    let log_l = log_len + log_k;
    // let coset = Coset::<T>::init(log_l, log_k);
    let coset = Coset::<T>::init(log_l, log_k, 1);

    let mut values = vec![T::from(0usize); 1 << log_len];
    for i in 0..(1 << log_len) {
        values[i] = T::random();
    }

    let coeffs = (coset.ifft)(&coset, &values, log_len);
    let values_ex = (coset.fft)(&coset, &coeffs, log_len, log_l);

    for i in 0..(1 << log_len) {
        assert_eq!(values[i], values_ex[i]);
    }

    // let new_values = (coset.fft)(&coset, &coeffs, log_len, log_len);
    // for i in 0..(1 << log_len) {
    //     assert_eq!(values[i], new_values[i]);
    // }

    // let new_coeffs = (coset.ifft)(&coset, &values_ex, log_l);
    // for i in 0..(1 << log_len) {
    //     assert_eq!(coeffs[i], new_coeffs[i]);
    // }
    // for i in coeffs.len()..new_coeffs.len() {
    //     assert_eq!(T::from(0), new_coeffs[i]);
    // }
}


fn test_mfft<T: Field>() {
    let log_len = 4;
    let log_k = 1;
    let log_l = log_len + log_k;
    // let coset = Coset::<T>::init(log_l, log_k);
    let coset = Coset::<T>::init(log_l, log_k, 1);

    let mut values = vec![T::from(0usize); 1 << log_len];
    for i in 0..(1 << log_len) {
        values[i] = T::random();
    }

    let coeffs = (coset.ifft)(&coset, &values, log_len);
    let values_ex = (coset.fft)(&coset, &coeffs, log_len, log_l);

    for i in 0..(1 << log_len) {
        assert_eq!(values[i], values_ex[i << log_k]);
    }

    // let new_values = (coset.fft)(&coset, &coeffs, log_len, log_len);
    // for i in 0..(1 << log_len) {
    //     assert_eq!(values[i], new_values[i]);
    // }

    // let new_coeffs = (coset.ifft)(&coset, &values_ex, log_l);
    // for i in 0..(1 << log_len) {
    //     assert_eq!(coeffs[i], new_coeffs[i]);
    // }
    // for i in coeffs.len()..new_coeffs.len() {
    //     assert_eq!(T::from(0), new_coeffs[i]);
    // }
}

// fn test_poly_over_vanish<T: Field>() {
//     let log_len = 10;
//     let log_k = 5;
//     let log_l = log_len + log_k;
//     let coset = Coset::<T>::init(log_l, log_k);
//     let mut inputs = vec![T::from(0usize); 1 << log_len];
//     let mut inputs_sq = vec![T::from(0usize); 1 << log_len];
//     for i in 0..(1 << log_len) {
//         inputs[i] = T::random();
//         inputs_sq[i] = inputs[i] * inputs[i];
//     }
// 
//     let coeff = coset.ifft_slow(&inputs, log_len);
//     let coeff2 = poly_mul(&coeff, &coeff);
// 
//     let outputs = coset.fft_slow(&coeff, log_l);
//     let outputs_sq = coset.fft_slow(&coeff2, log_l);
// 
//     for i in 0..outputs.len() {
//         assert_eq!(outputs_sq[i], outputs[i] * outputs[i]);
//     }
// 
//     let (q_coeff, r_coeff) = coset.poly_novel_over_vanish(&coeff2);
// 
//     let r = coset.fft_slow(&r_coeff, log_l);
//     let q = coset.fft_slow(&q_coeff, log_l);
// 
//     for i in 0..inputs.len() {
//         assert_eq!(T::from(0), poly_eval(&coset.vanishing_coeff, T::from(i)));
//     }
// 
//     for i in 0..inputs.len() {
//         assert_eq!(outputs_sq[i], r[i] + q[i] * poly_eval(&coset.vanishing_coeff, T::from(i)));
//     }
// 
// }

#[test]
fn test_fft_all() {
    test_mfft::<M61p2>();
    test_mfft::<P192>();
    test_mfft::<P256>();
    // test_afft::<GF2p16>();
    // test_afft::<GF2p32>();
    // test_afft::<GF2p64>();
    test_afft::<GF2p128>();
    test_afft::<GF2p192>();
    test_afft::<GF2p256>();
}

// fn test_poly_over_vanish_all() {
//     // test_poly_over_vanish::<GF2p16>();
//     // test_poly_over_vanish::<GF2p32>();
//     // test_poly_over_vanish::<GF2p64>();
//     // test_poly_over_vanish::<GF2p128>();
//     // test_poly_over_vanish::<GF2p256>();
// }

// #[test]
fn test_coset() {
    type T = M61p2;
    let log_len = 20;
    let log_k = 5;
    let log_l = log_len + log_k;
    let coset = Coset::<T>::init(log_l, log_k, 1);
}


#[test]
fn test_field() {
    test_field_arith::<GF2p16>();
    test_field_arith::<GF2p32>();
    test_field_arith::<GF2p64>();
    test_field_arith::<GF2p128>();
    test_field_arith::<GF2p192>();
    test_field_arith::<GF2p256>();
    test_field_arith::<P192>();
    test_field_arith::<P256>();
    test_field_arith::<M61p2>();
    test_field_arith::<M127p2>();
    // test_fft_all();
    // test_poly_over_vanish_all();
}
