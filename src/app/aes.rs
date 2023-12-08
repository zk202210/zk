use std::sync::Arc;

use crate::field::*;
use crate::field::common::*;
use crate::circuit::*;
use crate::gkr::*;
use crate::util::*;

pub fn test_bits<T: Field>(a: T, n: usize) -> T {
    let mut res = T::from(1);
    for i in 0..(1 << n) {
        res *= a - T::from(i);
    }
    res
}

pub fn aes_enc_test<T: Field>(state: &[T; 16], inv: &[T; 16], aff: &[T; 16], mix: &[T; 16], last: bool) -> bool {
    let inv_11b = T::from(0x11b).inv();
    let inv_101 = T::from(0x101).inv();

    // omit 8-bit tests for inv, aff, mix
    let mut res: Vec<T> = vec![];
    for k in 0..16 {
        // shift row
        let shift_k = (((k / 4) - k % 4) % 4) * 4 + (k % 4);
        let tmp = inv_11b * (T::from(1) + state[k] * inv[shift_k]);
        let test = test_bits(tmp, 7);
        res.push(test * state[k]);
        res.push(test * inv[shift_k]);

        let tmp = (inv[k] * T::from(0x1f) + aff[k] + T::from(0x63)) * inv_101;
        res.push(test_bits(tmp, 4));

        if !last {
            let tmp = (aff[k] * T::from(2) + mix[k]) * inv_11b;
            res.push(test_bits(tmp, 7));
        }
    }

    for i in 0..res.len() {
        if res[i] != T::from(0) {
            return false;
        }
    }
    true
}

pub fn fieldize<T: Field>(u8array: &[u8; 16]) -> [T; 16] {
    let mut res = [T::from(0); 16];
    for i in 0..16 {
        res[i] = T::from(u8array[i] as usize);
    }
    res
}

pub fn aes_enc(state: &[u8], roundkey: &[u8], last: bool) -> [[u8; 16]; 4] {
    let mut inv = [0u8; 16];
    let mut affine = [0u8; 16];
    let mut time2 = [0u8; 16];
    let mut res = [0u8; 16];

    // inverse
    for k in 0..16 {
        let s = GF2p8(state[k]);
        let i = s.inv();
        inv[k] = i.0;
    }

    // shiftrow
    affine[0] = inv[0];  affine[1] = inv[5];  affine[2] = inv[10]; affine[3] = inv[15];
    affine[4] = inv[4];  affine[5] = inv[9];  affine[6] = inv[14]; affine[7] = inv[3];
    affine[8] = inv[8];  affine[9] = inv[13]; affine[10]= inv[2];  affine[11]= inv[7];
    affine[12]= inv[12]; affine[13]= inv[1];  affine[14]= inv[6];  affine[15]= inv[11];

    let inv = affine.clone();

    for k in 0..16 {
        let a = clmul(affine[k] as u64, 0x1fu64) ^ 0x63u64;
        let a = a ^ (a >> 8);
        affine[k] = a as u8;
        // let a = affine[k] as u16;
        // let a = a ^ (a << 1) ^ (a << 2) ^ (a << 3) ^ (a << 4) ^ 0x63u16;
        // let a = a ^ (a >> 8);
        // affine[k] = a as u8;
        let a = (a & 0xff) * 2;
        if a < 0x100 {
            time2[k] = a as u8;
        } else {
            time2[k] = (a ^ 0x11b) as u8;
        }
        // assert_eq!(GF2p8(affine[k]) * GF2p8(2), GF2p8(time2[k]));
    }

    if last {
        for k in 0..16 {
            res[k] = roundkey[k] ^ affine[k];
        }
    } else {
        // mixcolumn and add roundkey
        for k in 0..4 {
            res[k * 4 + 0] =
                roundkey[k * 4 + 0] ^ time2[k * 4 + 0] ^ time2[k * 4 + 1] ^
                affine[k * 4 + 1] ^ affine[k * 4 + 2] ^ affine[k * 4 + 3];
            res[k * 4 + 1] =
                roundkey[k * 4 + 1] ^ time2[k * 4 + 1] ^ time2[k * 4 + 2] ^
                affine[k * 4 + 2] ^ affine[k * 4 + 3] ^ affine[k * 4 + 0];
            res[k * 4 + 2] =
                roundkey[k * 4 + 2] ^ time2[k * 4 + 2] ^ time2[k * 4 + 3] ^
                affine[k * 4 + 3] ^ affine[k * 4 + 0] ^ affine[k * 4 + 1];
            res[k * 4 + 3] =
                roundkey[k * 4 + 3] ^ time2[k * 4 + 3] ^ time2[k * 4 + 0] ^
                affine[k * 4 + 0] ^ affine[k * 4 + 1] ^ affine[k * 4 + 2];
        }
    }

    [res, inv, affine, time2]
}

pub fn gen_ctraes_witness<T: Field>(key: &[u8; 16], nonce: &[u8; 16], plain: &[u8; 16], log_n_aes: usize) -> (Vec<T>, Vec<T>) {
    let mut wit = vec![T::from(0); 512 << log_n_aes];
    let mut fin = vec![T::from(0); 16 << log_n_aes];

    for p in 0..(1 << log_n_aes) {
        let mut block_nonce = nonce.clone();
        block_nonce[0] ^= p as u8;
        block_nonce[1] ^= (p >> 8) as u8;
        let (new_wit, new_fin) = gen_aes_witness(key, &block_nonce, plain);
        for i in 0..464 {
            wit[i + 512 * p] = new_wit[i];
        }
        for i in 0..16 {
            fin[i + 16 * p] = new_fin[i];
        }

        if p == 0 {
            for i in 0..32 {
                wit[512 - 48 + i] = new_wit[512 + i];
            }
            for i in 0..48 {
                wit[512 * 2 - 48 + i] = new_wit[512 + 32 + i];
                wit[512 * 3 - 48 + i] = new_wit[512 + 32 + 48 + i];
                wit[512 * 4 - 48 + i] = new_wit[512 + 32 + 96 + i];
            }
            for i in 0..32 {
                wit[512 * 5 - 48 + i] = new_wit[512 + 32 + 144 + i];
            }

            for i in 0..40 {
                wit[512 * 7 - 48 + i] = new_wit[512 + 32 + 11 * 16 + i];
            }
            for i in 0..40 {
                wit[512 * 8 - 48 + i] = new_wit[512 + 32 + 11 * 16 + 10 * 4 + i];
            }
        }
    }

    (wit, fin)
}

pub fn gen_aes_witness<T: Field>(key: &[u8; 16], nonce: &[u8; 16], plain: &[u8; 16]) -> (Vec<T>, Vec<T>) {
    // I A I A I A ... I A M M M M M M M M M ...(512)
    let mut wit = vec![T::from(0); 1024];
    // key schedule
    let rcon: [u8; 10] = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36];
    // let w = fieldize(&key).to_vec();
    let mut w: Vec<u8> = vec![];
    let mut wi: Vec<u8> = vec![];
    let mut wa: Vec<u8> = vec![];
    for i in 0..16 {
        w.push(key[i]);
    }
    for r in 0..10 {
        for k in 0..4 {
            let inv = GF2p8(w[16 * r + k + 12]).inv().0;
            wi.push(inv);
            let a = clmul(inv as u64, 0x1fu64) ^ 0x63u64;
            let a = a ^ (a >> 8);
            let aff = a as u8;
            wa.push(aff);
        }
        for k in 0..4 {
            let mut tmp = w[16 * r + k] ^ wa[r * 4 + (k + 1) % 4];
            if k == 0 {
                tmp ^= rcon[r];
            }
            w.push(tmp);
        }
        for k in 4..16 {
            let tmp = w[16 * r + k] ^ w[16 * r + k + 16 - 4];
            w.push(tmp);
        }
    }

    let mut state = nonce.clone();
    for i in 0..16 {
        state[i] ^= w[i];
    }
    for r in 0..10 {
        let last = r == 9;
        let round_wit = aes_enc(&state, &w[16 * (r + 1)..], last);
        state = round_wit[0];

        let mut ofs = r * 32;
        for i in 0..16 {
            wit[ofs + i] = T::from(round_wit[1][i] as usize);
        }
        ofs += 16;
        for i in 0..16 {
            wit[ofs + i] = T::from(round_wit[2][i] as usize);
        }
        // ofs += 16;
        ofs = 10 * 32 + r * 16;
        if !last {
            for i in 0..16 {
                wit[ofs + i] = T::from(round_wit[3][i] as usize);
            }
            // ofs += 16;
        }
    }

    for i in 0..16 {
        state[i] ^= plain[i];
    }

    let mut ofs = 512;
    for i in ofs..(ofs + 16) {
        wit[i] = T::from(nonce[i - ofs] as usize);
    }
    ofs += 16;
    for i in ofs..(ofs + 16) {
        wit[i] = T::from(plain[i - ofs] as usize);
    }
    ofs += 16;
    // round key
    for i in ofs..(ofs + 11 * 16) {
        wit[i] = T::from(w[i - ofs] as usize);
    }
    // inversed round key
    ofs += 11 * 16;
    for i in ofs..(ofs + 10 * 4) {
        wit[i] = T::from(wi[i - ofs] as usize);
    }
    // affine transformed round key
    ofs += 10 * 4;
    for i in ofs..(ofs + 10 * 4) {
        wit[i] = T::from(wa[i - ofs] as usize);
    }

    (wit, fieldize(&state).to_vec())
}

    // method 1
    // input layer
    //
    // Each var is 4 gates,
    // round keys: W0, W1, ..., W43,
    // i-th round witness, 0 ≤ i ≤ 10
    //             Iᵢ0, Iᵢ1, Iᵢ2, Iᵢ3, (inverse)
    //             Jᵢ0, Jᵢ1, Jᵢ2, Jᵢ3, (witness for inverse)
    //             Aᵢ0, Aᵢ1, Aᵢ2, Aᵢ3, (affine transform)
    //             Bᵢ0, Bᵢ1, Bᵢ2, Bᵢ3, (witness for affine transform)
    //             Mᵢ0, Mᵢ1, Mᵢ2, Mᵢ3, (witness for mix-column, i ≠ 10)
    // Inner States Sᵢ0, Sᵢ1, Sᵢ2, Sᵢ3
    //              Pᵢ0, Pᵢ1, Pᵢ2, Pᵢ3
    //
    // where Sᵢ₋₁0 × Iᵢ0 = Pᵢ0 + Jᵢ0 × 0x11b
    //       Pᵢ0 × (1 - Sᵢ₋₁0) = 0
    //       Pᵢ0 × (1 - Iᵢ0) = 0
    //       Iᵢ0 × 0x1f = Aᵢ0 + 0x63 + Bᵢ0 × 0x101 ...... (1)
    //       Mᵢ0 × 0x11b⁻¹ + Aᵢ0 × 2 × 0x11b⁻¹ = 1 or 0
    //
    //  Ai0..Ai3 :=
    //      a0, a1, ..., af
    //  Apply ShiftRow:
    //      a0, a5, aa, af; a4, a9, ae, a3; a8, ad, aa, a7; ac, a1, a6, ab
    //  Apply MixColumn
    //      s0, s1, s2, s3; s4, s5, s6, s7; ...
    //  where
    //      s0 = 2a0 + 3a5 + aa + af = m0 + m5 + a5 + aa + af
    //      s1 = a0 + 2a5 + 3aa + af = a0 + m5 + ma + aa + af
    //      s2 = a0 + a5 + 2aa + 3af = a0 + a5 + ma + mf + af
    //      s3 = 3a0 + a5 + aa + 2af = a0 + a5 + aa + mf + m0
    //      ...
    //      Sᵢ0 = (s0, s1, s2, s3) + W₄ᵢ
    //      Sᵢ1 = (s4, s5, s6, s7) + W₄ᵢ₊₁
    //      Sᵢ2 = (s8, s9, sa, sb) + W₄ᵢ₊₂
    //      Sᵢ3 = (sc, sd, se, sf) + W₄ᵢ₊₃
    // Rearrange label for A, B, M to bring ShiftRow forward in eq (1)
    //
    // 1st layer:
    //      256 × 16 × 3 gates: 8-bit test for I, A, M
    //      128 × 16 × 2 gates: 7-bit test for J, B
    //      16 rel gates: i
    //      16 scl gates: k := 0x11b × j
    //      16 lin gates: c := (0x63 + a + 0x11 × b) / 0x1f
    //      16 lin gates: z := m × 0x11b⁻¹ + a × 2 × 0x11b⁻¹
    //      16 ad4 gates: r(j) := w(j) + a(j) + m(j) + m((j + 1) % 4 + j/4 * 4)
    //      4  sum gates: sum0 := sum(a0…b3),
    //                    sum1 := sum(a4…b7),
    //                    sum2 := sum(a8…bb),
    //                    sum3 := sum(ac…bf)
    // 2nd layer:
    //      128 × 16 × 3 gates: 8-bit test for I, A, M
    //       64 × 16 × 2 gates: 7-bit test for J, B
    //      16 rel gates: i
    //      16 rel gates: k = 0x11b × j
    //      16 add gates: c+i = 0
    //          // ShiftRow here
    //          // c0 + i0, c1 + i5, c2 + ia, c3 + if
    //          // c4 + i4, c5 + i9, c6 + ie, c7 + i3
    //          // c8 + i8, c9 + id, ca + i2, cb + i7
    //          // cc + ic, cd + i1, ce + i6, cf + ib
    //      16 tst gates: z(z+1) = 0
    //      16 add gates: sᵢ := r + sum (needed by next round)
    //
    // 3rd layer:
    //      64 × 16 × 3 gates: 8-bit test for I, A, M
    //      32 × 16 × 2 gates: 7-bit test for J, B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 rel gates: k = 0x11b × j
    //      16 mul gates: si := sᵢ₋₁ × i
    //      16 rel gates: c+i = 0
    //      16 rel gates: z(z+1) = 0
    //
    // 4th layer:
    //      32 × 16 × 3 gates: 8-bit test for I, A, M
    //      16 × 16 × 2 gates: 7-bit test for J, B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 add gates: p := k + si = 0x11b × j + sᵢ₋₁ × i
    //      16 rel gates: c+i = 0
    //      16 rel gates: z(z+1) = 0
    //
    // 5th layer:
    //      16 × 16 × 3 gates: 8-bit test for I, A, M
    //       8 × 16 × 2 gates: 7-bit test for J, B
    //      16 nab gates: (1-p)sᵢ₋₁
    //      16 nab gates: (1-p)i
    //      16 rel gates: c+i = 0
    //      16 rel gates: z(z+1) = 0
    //
    //
    //
    // method 2
    // input layer
    //
    // Each var is 4 gates,
    // round keys: W0, W1, W2, W3, W7, W11, .. W43
    // i-th round witness, 0 ≤ i ≤ 10
    //             Iᵢ0, Iᵢ1, Iᵢ2, Iᵢ3, (inverse)
    //             Aᵢ0, Aᵢ1, Aᵢ2, Aᵢ3, (affine transform)
    //             Mᵢ0, Mᵢ1, Mᵢ2, Mᵢ3, (witness for mix-column, M = 2A mod 0x11b, i ≠ 10)
    // Inner States Sᵢ0, Sᵢ1, Sᵢ2, Sᵢ3
    //              Jᵢ0, Jᵢ1, Jᵢ2, Jᵢ3, (witness for inverse)
    //              Bᵢ0, Bᵢ1, Bᵢ2, Bᵢ3, (witness for affine transform)
    //
    // where Sᵢ₋₁0 × Iᵢ0 = 1 + Jᵢ0 × 0x11b or (Sᵢ₋₁0 = 0 and Iᵢ0 = 0)
    //       J ∈ [0, 127] or (Sᵢ₋₁0 = 0 and Iᵢ0 = 0)
    //       Iᵢ0 × 0x1f = Aᵢ0 + 0x63 + Bᵢ0 × 0x101
    //       B ∈ [0, 127]
    //       Mᵢ0 × 0x11b⁻¹ + Aᵢ0 × 2 × 0x11b⁻¹ = 1 or 0
    //
    //  Ai0..Ai3 :=
    //      a0, a1, ..., af
    //  Apply ShiftRow:
    //      a0, a5, aa, af; a4, a9, ae, a3; a8, ad, aa, a7; ac, a1, a6, ab
    //  Apply MixColumn
    //      s0, s1, s2, s3; s4, s5, s6, s7; ...
    //  where
    //      s0 = 2a0 + 3a5 + aa + af = m0 + m5 + a5 + aa + af
    //      s1 = a0 + 2a5 + 3aa + af = a0 + m5 + ma + aa + af
    //      s2 = a0 + a5 + 2aa + 3af = a0 + a5 + ma + mf + af
    //      s3 = 3a0 + a5 + aa + 2af = a0 + a5 + aa + mf + m0
    //      ...
    //      Sᵢ0 = (s0, s1, s2, s3) + W₄ᵢ
    //      Sᵢ1 = (s4, s5, s6, s7) + W₄ᵢ₊₁
    //      Sᵢ2 = (s8, s9, sa, sb) + W₄ᵢ₊₂
    //      Sᵢ3 = (sc, sd, se, sf) + W₄ᵢ₊₃
    //
    // 1st layer:
    //      128 × 16 × 3 gates: 8-bit test for I, A, M which are row-shifted
    //      16 rel gates: i (reverse row shift)
    //      16 lin gates: b := (0x63 + a + 0x1f × i) / 0x101
    //      16 lin gates: z := m × 0x11b⁻¹ + a × 2 × 0x11b⁻¹
    //      16 add6 gates:
    //                    sᵢ(j) := w(j) + m(j) + m(j+1) + a(j+1) + a(j+2) + a(j+3)
    //                    where j+k := ((j%4 + k) % 4 + j/4 × 4)
    //      16 add6 and mul and linear gates:
    //                    j := (sᵢ₋₁ × i_(reverse row shift) + 1) × 0x11b⁻¹
    // 2nd layer:
    //      64 × 16 × 3 gates: 8-bit test for I, A, M
    //      64 × 16 × 1 gates: 7-bit test for J
    //       8 × 16 × 1 gates: 4-bit test for B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 tst gates: z(z+1) = 0
    // 3rd layer:
    //      32 × 16 × 3 gates: 8-bit test for I, A, M
    //      32 × 16 × 1 gates: 7-bit test for J
    //       4 × 16 × 1 gates: 4-bit test for B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 rel gates: z(z+1) = 0
    // 4th layer:
    //      16 × 16 × 3 gates: 8-bit test for I, A, M
    //      16 × 16 × 1 gates: 7-bit test for J
    //       2 × 16 × 1 gates: 4-bit test for B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 rel gates: z(z+1) = 0
    // 5th layer:
    //      8 × 16 × 3 gates: 8-bit test for I, A, M
    //      8 × 16 × 1 gates: 7-bit test for J
    //      1 × 16 × 1 gates: 4-bit test for B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 rel gates: z(z+1) = 0
    // 6th layer:
    //      4 × 16 × 3 gates: 8-bit test for I, A, M
    //      4 × 16 × 1 gates: 7-bit test for J
    //      1 × 16 × 1 gates: 4-bit test for B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 rel gates: z(z+1) = 0
    // 7th layer:
    //      2 × 16 × 3 gates: 8-bit test for I, A, M
    //      2 × 16 × 1 gates: 7-bit test for J
    //      1 × 16 × 1 gates: 4-bit test for B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 rel gates: z(z+1) = 0
    // 8th layer:
    //      1 × 16 × 3 gates: 8-bit test for I, A, M
    //      1 × 16 × 1 gates: 7-bit test for J =: p
    //      1 × 16 × 1 gates: 4-bit test for B
    //      16 rel gates: sᵢ₋₁
    //      16 rel gates: i
    //      16 rel gates: z(z+1) = 0
    // 9th layer:
    //      1 × 16 × 3 gates: 8-bit test for I, A, M
    //      1 × 16 × 1 gates: 4-bit test for B
    //      16 mul gates: p × sᵢ₋₁
    //      16 mul gates: p × i
    //      16 rel gates: z(z+1) = 0

macro_rules! eq1 {
    (-$x:expr) => {
        (T::from(1) - $x)
    };
    (+$x:expr) => {
        ($x)
    };
    (-$x:expr, $($s:tt $y:expr),+) => {
        (T::from(1) - $x) * eq1!($($s $y),+)
    };
    (+$x:expr, $($s:tt $y:expr),+) => {
        ($x) * eq1!($($s $y),+)
    };
}

macro_rules! eq0 {
    (+$x:expr) => {
        (T::from(1) - $x)
    };
    (-$x:expr) => {
        ($x)
    };
    (+$x:expr, $($s:tt $y:expr),+) => {
        (T::from(1) - $x) * eq0!($($s $y),+)
    };
    (-$x:expr, $($s:tt $y:expr),+) => {
        ($x) * eq0!($($s $y),+)
    };
}

macro_rules! eq {
    (-$x:expr, $($s:tt $y:expr),+) => {
        (((T::from(1) - $x) * eq1!($($s $y),+)) + (($x) * eq0!($($s $y),+)))
    };
    (+$x:expr, $($s:tt $y:expr),+) => {
        ((($x) * eq1!($($s $y),+)) + ((T::from(1) - $x) * eq0!($($s $y),+)))
    }
}

pub fn pred_ctraes_mid_layer<T: Field>(log_n_aes: usize, utils: &Vec<T>, d: usize, tau: &Vec<T>, rho: &Vec<T>, _sgm: &Vec<T>) -> [T; 4] {
    let mut res = [T::from(0); 4];

    let mut block_eq2 = T::from(1);
    for i in 0..log_n_aes {
        block_eq2 *= eq2(tau[i + 11], rho[i + 11]);
    }

    let gate_bit_2_eq2 = eq2(tau[0], rho[0]) * eq2(tau[1], rho[1]);
    let gate_bit_3_eq2 = gate_bit_2_eq2 * eq2(tau[2], rho[2]);
    let gate_bit_4_eq2 = gate_bit_3_eq2 * eq2(tau[3], rho[3]);

    let mut round_0_7_eq2 = zer(tau[10]) * zer(rho[10]);
    round_0_7_eq2 *= eq2(tau[9], rho[9]);
    round_0_7_eq2 *= eq2(tau[8], rho[8]);
    let mut round_8_9_eq2 = one(tau[10]) * one(rho[10]);
    round_8_9_eq2 *= zer(tau[9]) * zer(rho[9]);
    round_8_9_eq2 *= zer(tau[8]) * zer(rho[8]);
    let round_0_9_eq2 =
        (round_0_7_eq2 + round_8_9_eq2) *
        eq2(tau[7], rho[7]);

    let round_15_eq2 = one(tau[10]) * one(rho[10]) *
        one(tau[9]) * one(rho[9]) *
        one(tau[8]) * one(rho[8]) *
        one(tau[7]) * one(rho[7]);

    // 3/4 gate of first half of round 0-9 is addmul util[1]
    let pred_p1d2c0 = block_eq2 * gate_bit_4_eq2 * round_0_9_eq2 *
        zer(tau[6]) * zer(rho[6]) *
        (zer(tau[5]) * zer(rho[5]) *
         eq2(tau[4], rho[4]) +
         one(tau[5]) * one(rho[5]) *
         zer(tau[4]) * zer(rho[4]));
    res[0] += pred_p1d2c0;
    res[1] += pred_p1d2c0 * utils[0];

    // let pred_relay = block_eq2 * zer_sgm * gate_bit_4_eq2 *
    //     (round_0_9_eq2 *
    //      (zer(tau[6]) * zer(rho[6]) *
    //       eq2(tau[5], rho[5]) * eq2(tau[4], rho[4]) +
    //       one(tau[6]) * one(rho[6]) *
    //       one(tau[5]) * one(rho[5]) *
    //       one(tau[4]) * one(rho[4])) +
    // 16-47 gates, second half of round 0-9 is relay
    // 0-16 of round 15 is relay
    let mut pred_relay = block_eq2 * gate_bit_4_eq2 *
        (round_0_9_eq2 * one(tau[6]) * one(rho[6]) *
         (zer(tau[5]) * zer(rho[5]) * one(tau[4]) * one(rho[4]) +
          one(tau[5]) * one(rho[5]) * zer(tau[4]) * zer(rho[4])) +
         round_15_eq2 *
         zer(tau[6]) * zer(rho[6]) *
         zer(tau[5]) * zer(rho[5]) *
         zer(tau[4]) * zer(rho[4]));

    // 0-15 gates, second half of round 0-9 is addmul
    let mut pred_addmul = block_eq2 * gate_bit_4_eq2 *
        round_0_9_eq2 * one(tau[6]) * one(rho[6]) *
        zer(tau[5]) * zer(rho[5]) *
        zer(tau[4]) * zer(rho[4]);

    // last 1/4 gate of first half is addmul/relay
    if d < 5 {
        pred_addmul += block_eq2 * gate_bit_4_eq2 *
            round_0_9_eq2 * zer(tau[6]) * zer(rho[6]) *
            one(tau[5]) * one(rho[5]) *
            one(tau[4]) * one(rho[4]);
    } else {
        pred_relay += block_eq2 * gate_bit_4_eq2 *
            round_0_9_eq2 * zer(tau[6]) * zer(rho[6]) *
            one(tau[5]) * one(rho[5]) *
            one(tau[4]) * one(rho[4]);
    }

    // last 1/4 gate of second half is addmul/relay
    if d < 2 {
        pred_addmul += block_eq2 * gate_bit_4_eq2 *
            round_0_9_eq2 * one(tau[6]) * one(rho[6]) *
            one(tau[5]) * one(rho[5]) *
            one(tau[4]) * one(rho[4]);
    } else {
        pred_relay += block_eq2 * gate_bit_4_eq2 *
            round_0_9_eq2 * one(tau[6]) * one(rho[6]) *
            one(tau[5]) * one(rho[5]) *
            one(tau[4]) * one(rho[4]);
    }

    res[0] += pred_addmul;
    res[1] += pred_relay + pred_addmul * utils[1];


    // block 0 smallround 0-9:
    // 10.0 = 0b1010 | 0
    // 10.1 = 0b1010 | 1
    // 10.2 = 0b1011 | 0
    // 10.3 = 0b1011 | 1
    //
    // 10.4 = 0b1100 | 0
    // 10.5 = 0b1100 | 1
    // 10.6 = 0b1101 | 0
    // 10.7 = 0b1101 | 1
    //
    // 10.8 = 0b1110 | 0
    // 10.9 = 0b1110 | 1
    let mut block_0_eq2 = T::from(1);
    for i in 0..log_n_aes {
        block_0_eq2 *= zer(tau[i + 11]) * zer(rho[i + 11]);
    }

    // TODO merge
    let smallround_0_3_eq2 = one(tau[10]) * one(rho[10]) * zer(tau[9]) * zer(rho[9]) * one(tau[8]) * one(rho[8]);

    let smallround_4_7_eq2 = one(tau[10]) * one(rho[10]) * one(tau[9]) * one(rho[9]) * zer(tau[8]) * zer(rho[8]);

    let smallround_8_9_eq2 = one(tau[10]) * one(rho[10]) * one(tau[9]) * one(rho[9]) * one(tau[8]) * one(rho[8]) * zer(tau[7]) * zer(rho[7]);

    let smallround_0_9_eq2 = eq2(tau[6], rho[6]) * (eq2(tau[7], rho[7]) * (smallround_0_3_eq2 + smallround_4_7_eq2) + smallround_8_9_eq2);

    // gate 0-8 addmul
    let pred_addmul = block_0_eq2 * smallround_0_9_eq2 * gate_bit_3_eq2 *
        zer(tau[5]) * zer(rho[5]) *
        zer(tau[4]) * zer(rho[4]) *
        zer(tau[3]) * zer(rho[3]);
    res[0] += pred_addmul;
    res[1] += pred_addmul * utils[0];

    // gate 20-23 relay
    // gate 24-27 relay
    // gate 32-47 relay
    let mut pred_relay = block_0_eq2 * smallround_0_9_eq2 *
        (gate_bit_4_eq2 * zer(tau[4]) * zer(rho[4]) * one(tau[5]) * one(rho[5]) +
         gate_bit_2_eq2 * zer(tau[5]) * zer(rho[5]) * one(tau[4]) * one(rho[4]) *
         (zer(tau[3]) * zer(rho[3]) * one(tau[2]) * one(rho[2]) +
          one(tau[3]) * one(rho[3]) * zer(tau[2]) * zer(rho[2])));

    // gate 16-19 addmul
    let mut pred_addmul = block_0_eq2 * smallround_0_9_eq2 * gate_bit_2_eq2 *
        zer(tau[5]) * zer(rho[5]) *
        one(tau[4]) * one(rho[4]) *
        zer(tau[3]) * zer(rho[3]) *
        zer(tau[2]) * zer(rho[2]);

    // gate 28-31 addmul/relay
    if d < 5 {
        pred_addmul += block_0_eq2 * smallround_0_9_eq2 * gate_bit_2_eq2 *
            zer(tau[5]) * zer(rho[5]) *
            one(tau[4]) * one(rho[4]) *
            one(tau[3]) * one(rho[3]) *
            one(tau[2]) * one(rho[2]);
    } else {
        pred_relay += block_0_eq2 * smallround_0_9_eq2 *
            gate_bit_2_eq2 *
            zer(tau[5]) * zer(rho[5]) *
            one(tau[4]) * one(rho[4]) *
            one(tau[3]) * one(rho[3]) *
            one(tau[2]) * one(rho[2]);
    }
    res[0] += pred_addmul;
    res[1] += pred_addmul * utils[1] + pred_relay;

    res
}

pub fn pred_ctraes_last_layer<T: Field>(log_n_aes: usize, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
    let mut res = [T::from(0); 4];

    let mut zer_sgm = T::from(1);
    for i in 0..sgm.len() {
        zer_sgm *= zer(sgm[i]);
    }

    let mut block_eq3 = T::from(1);
    for i in 0..log_n_aes {
        block_eq3 *= eq3(tau[i + 11], rho[i + 11], sgm[i + 11]);
    }

    let mut block_eq2 = T::from(1);
    for i in 0..log_n_aes {
        block_eq2 *= eq2(tau[i + 11], rho[i + 11]);
    }

    let gate_bit_2_eq2 = eq2(tau[0], rho[0]) * eq2(tau[1], rho[1]);
    let gate_bit_3_eq2 = gate_bit_2_eq2 * eq2(tau[2], rho[2]);
    let gate_bit_4_eq2 = gate_bit_3_eq2 * eq2(tau[3], rho[3]);

    let gate_bit_2_eq3 = eq3(tau[0], rho[0], sgm[0]) * eq3(tau[1], rho[1], sgm[1]);
    let gate_bit_3_eq3 = gate_bit_2_eq3 * eq3(tau[2], rho[2], sgm[2]);
    let gate_bit_4_eq3 = gate_bit_3_eq3 * eq3(tau[3], rho[3], sgm[3]);

    let mut round_0_7_eq2 = zer(tau[10]) * zer(rho[10]);
    let mut round_0_7_eq3 = round_0_7_eq2 * zer(sgm[10]);
    round_0_7_eq2 *= eq2(tau[9], rho[9]);
    round_0_7_eq2 *= eq2(tau[8], rho[8]);
    round_0_7_eq3 *= eq3(tau[9], rho[9], sgm[9]);
    round_0_7_eq3 *= eq3(tau[8], rho[8], sgm[8]);
    let mut round_8_9_eq2 = one(tau[10]) * one(rho[10]);
    round_8_9_eq2 *= zer(tau[9]) * zer(rho[9]);
    round_8_9_eq2 *= zer(tau[8]) * zer(rho[8]);
    let round_8_9_eq3 = round_8_9_eq2 * one(sgm[10]) * zer(sgm[9]) * zer(sgm[8]);
    let round_0_9_eq2 =
        (round_0_7_eq2 + round_8_9_eq2) *
        eq2(tau[7], rho[7]);
    let round_0_9_eq3 =
        (round_0_7_eq3 + round_8_9_eq3) *
        eq3(tau[7], rho[7], sgm[7]);

    let round_15_eq2 = one(tau[10]) * one(rho[10]) *
        one(tau[9]) * one(rho[9]) *
        one(tau[8]) * one(rho[8]) *
        one(tau[7]) * one(rho[7]);

    // first half of each 0-9 round is relay
    // 48-63 gates, second half of round 0-9 is relay
    // 0-16 of round 15 is relay
    let pred_relay = block_eq2 * zer_sgm * gate_bit_4_eq2 *
        (round_0_9_eq2 *
         (zer(tau[6]) * zer(rho[6]) *
          eq2(tau[5], rho[5]) * eq2(tau[4], rho[4]) +
          one(tau[6]) * one(rho[6]) *
          one(tau[5]) * one(rho[5]) *
          one(tau[4]) * one(rho[4])) +
         round_15_eq2 *
         zer(tau[6]) * zer(rho[6]) *
         zer(tau[5]) * zer(rho[5]) *
         zer(tau[4]) * zer(rho[4]));

    res[1] += pred_relay;

    // 0-15 gates, second half of each 0-9 round is mul
    // 16-31 gates, second half of each 0-9 round is mul
    let pred_mul = block_eq3 * round_0_9_eq3 * gate_bit_4_eq3 *
        one(tau[6]) * one(rho[6]) * one(sgm[6]) * // second half
        zer(tau[5]) * zer(rho[5]) * zer(rho[4]) *
        (zer(sgm[5]) * zer(tau[4]) * one(sgm[4]) +
         one(sgm[5]) * one(tau[4]) * zer(sgm[4]));
    res[0] += pred_mul;

    // block 0 smallround 0-9:
    // 10.0 = 0b1010 | 0
    // 10.1 = 0b1010 | 1
    // 10.2 = 0b1011 | 0
    // 10.3 = 0b1011 | 1
    //
    // 10.4 = 0b1100 | 0
    // 10.5 = 0b1100 | 1
    // 10.6 = 0b1101 | 0
    // 10.7 = 0b1101 | 1
    //
    // 10.8 = 0b1110 | 0
    // 10.9 = 0b1110 | 1
    let mut block_0_eq2 = T::from(1);
    for i in 0..log_n_aes {
        block_0_eq2 *= zer(tau[i + 11]) * zer(rho[i + 11]);
    }
    let mut block_0_eq3 = block_0_eq2;
    for i in 0..log_n_aes {
        block_0_eq3 *= zer(sgm[i + 11]);
    }

    // TODO merge
    let smallround_0_3_eq2 = one(tau[10]) * one(rho[10]) * zer(tau[9]) * zer(rho[9]) * one(tau[8]) * one(rho[8]);
    let smallround_0_3_eq3 = smallround_0_3_eq2 * one(sgm[10]) * zer(sgm[9]) * one(sgm[8]);

    let smallround_4_7_eq2 = one(tau[10]) * one(rho[10]) * one(tau[9]) * one(rho[9]) * zer(tau[8]) * zer(rho[8]);
    let smallround_4_7_eq3 = smallround_4_7_eq2 * one(sgm[10]) * one(sgm[9]) * zer(sgm[8]);

    let smallround_8_9_eq2 = one(tau[10]) * one(rho[10]) * one(tau[9]) * one(rho[9]) * one(tau[8]) * one(rho[8]) * zer(tau[7]) * zer(rho[7]);
    let smallround_8_9_eq3 = smallround_8_9_eq2 * one(sgm[10]) * one(sgm[9]) * one(sgm[8]) * zer(sgm[7]);

    let smallround_0_9_eq2 = eq2(tau[6], rho[6]) * (eq2(tau[7], rho[7]) * (smallround_0_3_eq2 + smallround_4_7_eq2) + smallround_8_9_eq2);
    let smallround_0_9_eq3 = eq3(tau[6], rho[6], sgm[6]) * (eq3(tau[7], rho[7], sgm[7]) * (smallround_0_3_eq3 + smallround_4_7_eq3) + smallround_8_9_eq3);

    // gate 0-7 relay
    // gate 28-31
    // gate 32-47 relay
    let pred_relay = block_0_eq2 * smallround_0_9_eq2 * zer_sgm *
        (zer(tau[4]) * zer(rho[4]) *
         (gate_bit_3_eq2 *
          zer(tau[5]) * zer(rho[5]) *
          zer(tau[3]) * zer(rho[3]) +
          gate_bit_4_eq2 *
          one(tau[5]) * one(rho[5])) +
         gate_bit_2_eq2 * zer(tau[5]) * zer(rho[5]) *
         one(tau[4]) * one(rho[4]) *
         one(tau[3]) * one(rho[3]) *
         one(tau[2]) * one(rho[2]));
    res[1] += pred_relay;

    // gate 16-19 mul
    // gate 20-23 mul
    let pred_mul = block_0_eq3 * smallround_0_9_eq3 * gate_bit_2_eq3 *
        zer(tau[5]) * zer(rho[5]) * zer(sgm[5]) *
        one(tau[4]) * one(rho[4]) * one(sgm[4]) *
        zer(tau[3]) * zer(rho[3]) * zer(rho[2]) *
        (zer(sgm[3]) * zer(tau[2]) * one(sgm[2]) +
         one(sgm[3]) * one(tau[2]) * zer(sgm[2]));
    res[0] += pred_mul;

    res
}

pub fn pred_ctraes_first_layer<T: Field>(log_n_aes: usize, utils: &Vec<T>, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
    // input:  log_n_aes + 9 bits
    //         xxxxxxxxx | x x x x x | x x x x
    //
    // output: log_n_aes + 11 bits
    // output: log_n_aes + 4 bits + 7 bits
    //         xxxxxxxxx | x x x x | x x x x x x x
    //         xxxxxxxxx | x x x x | x x x x x x x
    let mut res = [T::from(0); 4];

    let mut block_eq3 = T::from(1);
    for i in 0..log_n_aes {
        block_eq3 *= eq3(tau[i + 11], rho[i + 9], sgm[i + 9]);
    }

    let mut block_eq2 = T::from(1);
    for i in 0..log_n_aes {
        block_eq2 *= eq2(tau[i + 11], rho[i + 9]);
    }

    let mut block_eq2_sgm = T::from(1);
    for i in 0..log_n_aes {
        block_eq2_sgm *= eq2(tau[i + 11], sgm[i + 9]);
    }

    let mut block_0_hi_eq2 = zer(tau[11]) * zer(tau[12]) * zer(tau[13]);
    for i in 0..(log_n_aes - 3) {
        block_0_hi_eq2 *= zer(tau[i + 14]) * zer(rho[i + 12]);
    }

    let mut block_0_hi_eq3 = block_0_hi_eq2;
    for i in 0..(log_n_aes - 3) {
        block_0_hi_eq3 *= zer(sgm[i + 12]);
    }
    // let block_0_eq2 = block_0_hi_eq2 *
    //     zer(rho[9]) * zer(rho[10]) * zer(rho[11]);
    let block_0_eq3 = block_0_hi_eq3 *
        zer(rho[9]) * zer(rho[10]) * zer(rho[11]) *
        zer(sgm[9]) * zer(sgm[10]) * zer(sgm[11]);

    let gate_bit_2_eq2 = eq2(tau[0], rho[0]) * eq2(tau[1], rho[1]);
    let gate_bit_3_eq2 = gate_bit_2_eq2 * eq2(tau[2], rho[2]);
    let gate_bit_4_eq2 = gate_bit_3_eq2 * eq2(tau[3], rho[3]);

    let gate_bit_2_eq3 = eq3(tau[0], rho[0], sgm[0]) * eq3(tau[1], rho[1], sgm[1]);
    let gate_bit_3_eq3 = gate_bit_2_eq3 * eq3(tau[2], rho[2], sgm[2]);
    let gate_bit_4_eq3 = gate_bit_3_eq3 * eq3(tau[3], rho[3], sgm[3]);

    let mut round_0_7_eq2 = zer(tau[10]) * zer(rho[8]);
    let mut round_0_7_eq3 = round_0_7_eq2 * zer(sgm[8]);
    round_0_7_eq2 *= eq2(tau[9], rho[7]);
    round_0_7_eq2 *= eq2(tau[8], rho[6]);
    round_0_7_eq3 *= eq3(tau[9], rho[7], sgm[7]);
    round_0_7_eq3 *= eq3(tau[8], rho[6], sgm[6]);
    let mut round_8_9_eq2 = one(tau[10]) * one(rho[8]);
    round_8_9_eq2 *= zer(tau[9]) * zer(rho[7]);
    round_8_9_eq2 *= zer(tau[8]) * zer(rho[6]);
    let round_8_9_eq3 = round_8_9_eq2 * one(sgm[8]) * zer(sgm[7]) * zer(sgm[6]);
    let round_0_9_eq2 =
        (round_0_7_eq2 + round_8_9_eq2) *
        eq2(tau[7], rho[5]);
    let round_0_9_eq3 =
        (round_0_7_eq3 + round_8_9_eq3) *
        eq3(tau[7], rho[5], sgm[5]);

    // let round_9_eq2 = round_8_9_eq2 * one(tau[7]) * one(rho[5]);
    // let round_9_eq3 = round_8_9_eq3 * one(tau[7]) * one(rho[5]) * one(sgm[5]);
    // let round_0_8_eq2 = round_0_9_eq2 - round9_eq2;
    // let round_0_8_eq3 = round_0_9_eq3 - round9_eq3;

    let mut round_0_7_eq2_sgm = zer(tau[10]) * zer(sgm[8]);
    round_0_7_eq2_sgm *= eq2(tau[9], sgm[7]) * eq2(tau[8], sgm[6]);
    let mut round_8_9_eq2_sgm = one(tau[10]) * one(sgm[8]);
    round_8_9_eq2_sgm *= zer(tau[9]) * zer(sgm[7]) * zer(tau[8]) * zer(sgm[6]);
    // let round_0_9_eq2_sgm = (round_0_7_eq2_sgm + round_8_9_eq2_sgm) * eq2(tau[7], sgm[5]);

    let round_15 = one(tau[10]) * one(tau[9]) * one(tau[8]) * one(tau[7]);

    let mut zer_rho = T::from(1);
    for i in 0..rho.len() {
        zer_rho *= zer(rho[i]);
    }

    let mut zer_sgm = T::from(1);
    for i in 0..sgm.len() {
        zer_sgm *= zer(sgm[i]);
    }

    let mut zer_hi_rho = T::from(1);
    for i in 12..rho.len() {
        zer_hi_rho *= zer(rho[i]);
    }

    let mut zer_hi_tau = T::from(1);
    for i in 11..tau.len() {
        zer_hi_tau *= zer(tau[i]);
    }

    // first 32 gates of first 10 rounds of each block is vanishing gate for I A
    // x(x + 1)
    // input:     xxxxxxx  x x | x x   x x | q q q q q
    //                    10 9   8 7   6 5   4
    // output: xxxxxxxxx | x x   x x | 0 0 | q q q q q
    let pred_vanish_i_a =
        round_0_9_eq3 * gate_bit_4_eq3 * block_eq3 *
        eq3(tau[4], rho[4], sgm[4]) *
        zer(tau[5]) * zer(tau[6]);
    res[0] += pred_vanish_i_a;
    res[1] += pred_vanish_i_a;


    // 32~47 gates of first 9 rounds of each block is vanishing gate for M
    // tau                                          rho/sgm
    //
    // 0000 010 qqqq <-> round 0 <-> 20 * (16 B) <-> 10100qqqq
    // 0001 010 qqqq <-> round 1 <-> 21 * (16 B) <-> 10101qqqq
    // 0010 010 qqqq <-> round 2 <-> 22 * (16 B) <-> 10110qqqq
    // 0011 010 qqqq <-> round 3 <-> 23 * (16 B) <-> 10111qqqq
    //
    // 0100 010 qqqq <-> round 4 <-> 24 * (16 B) <-> 11000qqqq
    // 0101 010 qqqq <-> round 5 <-> 25 * (16 B) <-> 11001qqqq
    // 0110 010 qqqq <-> round 6 <-> 26 * (16 B) <-> 11010qqqq
    // 0111 010 qqqq <-> round 7 <-> 27 * (16 B) <-> 11011qqqq
    //
    // 1000 10 qqqq <-> round 8 <-> 28 * (16 B) <-> 11100qqqq

    let mut round_first_8 = eq!(+tau[9], +rho[7], +sgm[7], -rho[6], -sgm[6]);
    round_first_8 *= eq3(tau[8], rho[5], sgm[5]);
    round_first_8 *= eq3(tau[7], rho[4], sgm[4]);
    round_first_8 *= zer(tau[10]);
    let mut round_9 = one(tau[10]) * zer(tau[9]) * zer(tau[8]) * zer(tau[7]);
    round_9 *= one(rho[7]) * one(rho[6]) * zer(rho[5]) * zer(rho[4]);
    round_9 *= one(sgm[7]) * one(sgm[6]) * zer(sgm[5]) * zer(sgm[4]);
    let round_first_9 = (round_first_8 + round_9) * one(rho[8]) * one(sgm[8]);
    let pred_vanish_m = round_first_9 * gate_bit_4_eq3 * block_eq3 *
        zer(tau[6]) * one(tau[5]) * zer(tau[4]);
    res[0] += pred_vanish_m;
    res[1] += pred_vanish_m;

    // 32-39 vanishing for Nonce/PriKey
    // n_ofs =    0b111010000
    // w_ofs = 0b001111010000
    let mut pred_vanish_nk = block_0_eq3 * gate_bit_2_eq3 *
        one(tau[10]) * zer(tau[9]) * zer(tau[8]) * one(tau[7]) *
        zer(tau[6]) * one(tau[5]) * zer(tau[4]) * zer(tau[3]) * zer(tau[2]) *
        zer(rho[2]) * zer(rho[3]) * one(rho[4]) * zer(rho[5]) * one(rho[6]) * one(rho[7]) * one(rho[8]) *
        zer(sgm[2]) * zer(sgm[3]) * one(sgm[4]) * zer(sgm[5]) * one(sgm[6]) * one(sgm[7]) * one(sgm[8]);
    pred_vanish_nk += block_0_hi_eq3 * gate_bit_2_eq3 *
        one(tau[10]) * zer(tau[9]) * zer(tau[8]) * one(tau[7]) *
        zer(tau[6]) * one(tau[5]) * zer(tau[4]) * zer(tau[3]) * one(tau[2]) *
        one(rho[9]) * zer(rho[10]) * zer(rho[11]) * zer(rho[2]) * zer(rho[3]) * one(rho[4]) * zer(rho[5]) * one(rho[6]) * one(rho[7]) * one(rho[8]) *
        one(sgm[9]) * zer(sgm[10]) * zer(sgm[11]) * zer(sgm[2]) * zer(sgm[3]) * one(sgm[4]) * zer(sgm[5]) * one(sgm[6]) * one(sgm[7]) * one(sgm[8]);
    res[0] += pred_vanish_nk;
    res[1] += pred_vanish_nk;


    // 48~63 gates of first 10 rounds of each block is
    // linear gates: b := (0x63 + a + 0x1f × i) / 0x101
    // xxxx 011 qqqq <-> xxxx 0qqqq & xxxx 1qqqq
    let pred_linear_b = round_0_9_eq3 * gate_bit_4_eq3 * block_eq3 *
        zer(tau[6]) * one(tau[5]) * one(tau[4]) *
        zer(rho[4]) * one(sgm[4]);
    res[1] += pred_linear_b * utils[2];
    res[2] += pred_linear_b * utils[3];
    // XXX pred should for sgm = 0
    // res[3] += pred_linear_b * utils[4];
    res[3] += zer_sgm * zer_rho *
        (zer(tau[10]) + one(tau[10]) * zer(tau[9]) * zer(tau[8])) *
        zer(tau[6]) * one(tau[5]) * one(tau[4]) * utils[4];

    // 64-79   J gates
    // 96-111 Si gates
    let mut pred_w_i = T::from(0);
    let mut pred_ma_i = T::from(0);
    let mut pred_w = T::from(0);
    let mut pred_ma = T::from(0);
    let mut pred_ctr_0 = T::from(0);
    let mut pred_ctr_1 = T::from(0);
    let mut pred_ctr_0_s = T::from(0);
    let mut pred_ctr_1_s = T::from(0);
    for r in 0..10 {
        // let i_ofs = r * 16 * 2;
        let a_ofs = r * 16 * 2 + 16;
        let m_ofs = r * 16 + 10 * 32;
        let w_ofs = 512 * (r / 3 + 2) - 48 + (r % 3) * 16;
        // let k_ofs = r * 128 + 64;

        // let mut pred_w_i_r = T::from(0);
        for k in 0..16 {
            let r0 = r;
            let a_ofs = a_ofs - 32;
            let m_ofs = m_ofs - 16;
            let k0 = k;
            let k1 = (k0 % 4 + 1) % 4 + k0 / 4 * 4;
            let k2 = (k0 % 4 + 2) % 4 + k0 / 4 * 4;
            let k3 = (k0 % 4 + 3) % 4 + k0 / 4 * 4;
            let shift_k = (((k / 4) - k % 4) % 4) * 4 + (k % 4);

            let w = expand_binary(w_ofs + k0, 12);
            let i = expand_binary(shift_k, 4);
            let k = expand_binary(k, 4);
            let r = expand_binary(r, 4);

            // assert_eq!(k[4], T::from(0));
            // assert_eq!(k[5], T::from(0));
            // assert_eq!(k[6], T::from(1));
            // assert_eq!(i[4], T::from(0));

            let mut tmp_k = T::from(1);
            for j in 0..4 {
                tmp_k *= T::from(1) + k[j] + tau[j];
            }

            let mut tmp_w = T::from(1);
            for j in 0..12 {
                tmp_w *= T::from(1) + w[j] + rho[j];
            }

            let mut tmp_i = T::from(1);
            for j in 0..4 {
                tmp_i *= T::from(1) + i[j] + sgm[j];
            }

            for j in 0..4 {
                tmp_k *= T::from(1) + r[j] + tau[7 + j];
                tmp_i *= T::from(1) + r[j] + sgm[5 + j];
            }

            if r0 == 0 {
                let mut tmp_n = T::from(1);
                let n = expand_binary(512 - 48 + k0, 12);
                for j in 0..12 {
                    tmp_n *= T::from(1) + n[j] + rho[j];
                }
                tmp_w += tmp_n;

            } else {
                let bits = [m_ofs + k0, m_ofs + k1, a_ofs + k1, a_ofs + k2, a_ofs + k3];
                let mut tmp_ma = T::from(0);
                for b in bits.iter() {
                    let mut tmp = T::from(1);
                    let b = expand_binary(*b, 9);
                    for j in 0..9 {
                        tmp *= T::from(1) + b[j] + rho[j];
                    }
                    tmp_ma += tmp;
                }
                tmp_ma *= tmp_k;
                pred_ma += tmp_ma;
                pred_ma_i += tmp_ma * tmp_i;
            }

            tmp_w *= tmp_k;
            pred_w += tmp_w;
            pred_w_i += tmp_w * tmp_i;

            if k0 == 0 && r0 == 0 {
                pred_ctr_0_s += tmp_k;
                pred_ctr_0 += tmp_i * tmp_k;
            }
            if k0 == 1 && r0 == 0 {
                pred_ctr_1_s += tmp_k;
                pred_ctr_1 += tmp_i * tmp_k;
            }
        }
    }

    // Si 96-111, 96 = 64 + 32
    pred_w *= zer_sgm * zer_hi_rho * one(tau[6]) * one(tau[5]) * zer(tau[4]);
    pred_ma *= zer_sgm * block_eq2 * zer(tau[4]) * one(tau[5]) * one(tau[6]);
    res[1] += pred_w + pred_ma;
    // Si pred_ctr
    let mut tmp_0 = T::from(0);
    let mut tmp_1 = T::from(0);
    for i in 0..log_n_aes.min(8) {
        tmp_0 += T::from(1 << i) * one(tau[11 + i]);
    }
    for i in 8..log_n_aes {
        tmp_1 += T::from(1 << (i - 8)) * one(tau[11 + i]);
    }
    res[3] += zer_sgm * zer_rho * zer(tau[4]) * one(tau[5]) * one(tau[6]) *
        (tmp_0 * pred_ctr_0_s + tmp_1 * pred_ctr_1_s);

    // J 64-79
    pred_w_i *= zer(sgm[4]) *
        zer(tau[4]) * zer(tau[5]) * one(tau[6]) *
        block_eq2_sgm * zer_hi_rho * utils[6];
    pred_ma_i *= zer(sgm[4]) *
        zer(tau[4]) * zer(tau[5]) * one(tau[6]) *
        block_eq3 * utils[6];
    res[0] += pred_w_i + pred_ma_i;
    res[3] += zer_sgm * zer_rho *
        one(tau[6]) * zer(tau[5]) * zer(tau[4]) * utils[6] *
        (zer(tau[10]) + one(tau[10]) * zer(tau[9]) * zer(tau[8]));

    // J pred_ctr
    let mut tmp_0 = T::from(0);
    let mut tmp_1 = T::from(1);
    // TODO O(log_n_aes) alg
    for i in 0..log_n_aes.min(8) {
        let mut tt = T::from(1 << i) * one(tau[11 + i]) * one(sgm[9 + i]);
        for j in 0..log_n_aes.min(8) {
            if i != j {
                tt *= eq2(tau[11 + j], sgm[9 + j]);
            }
        }
        tmp_0 += tt;
        tmp_1 *= eq2(tau[11 + i], sgm[9 + i]);
    }
    let mut ttmp_1 = T::from(0);
    for i in 8..log_n_aes {
        tmp_0 *= eq2(tau[11 + i], sgm[9 + i]);
        let mut tt = T::from(1 << (i - 8)) * one(tau[11 + i]) * one(sgm[9 + i]);
        for j in 8..log_n_aes {
            if i != j {
                tt *= eq2(tau[11 + j], sgm[9 + j]);
            }
        }
        ttmp_1 += tt;
    }
    res[2] += (pred_ctr_0 * tmp_0 + pred_ctr_1 * tmp_1 * ttmp_1) *
        zer_rho * zer(sgm[4]) *
        one(tau[6]) * zer(tau[5]) * zer(tau[4]) * utils[6];


    let mut shiftrow_rho = T::from(0);
    // shiftrow_rho += zer(tau[0]) * zer(tau[1]) * zer(rho[0]) * zer(rho[1]) *
    //     eq2(tau[2], rho[2]) * eq2(tau[3], rho[3]);
    // shiftrow_rho += zer(tau[0]) * one(tau[1]) * zer(rho[0]) * one(rho[1]) *
    //     eq2(tau[2], rho[2]) * neq2(tau[3], rho[3]);
    shiftrow_rho += zer(tau[0]) * zer(rho[0]) * eq2(tau[2], rho[2]) *
        (zer(tau[1]) * zer(rho[1]) * eq2(tau[3], rho[3]) +
         one(tau[1]) * one(rho[1]) * neq2(tau[3], rho[3]));
    shiftrow_rho += one(tau[0]) * zer(tau[1]) * one(rho[0]) * zer(rho[1]) *
        (zer(tau[3]) * zer(tau[2]) * one(rho[3]) * one(rho[2]) +
         zer(tau[3]) * one(tau[2]) * zer(rho[3]) * zer(rho[2]) +
         one(tau[3]) * zer(tau[2]) * zer(rho[3]) * one(rho[2]) +
         one(tau[3]) * one(tau[2]) * one(rho[3]) * zer(rho[2]));
    shiftrow_rho += one(tau[0]) * one(tau[1]) * one(rho[0]) * one(rho[1]) *
        (zer(tau[3]) * zer(tau[2]) * zer(rho[3]) * one(rho[2]) +
         zer(tau[3]) * one(tau[2]) * one(rho[3]) * zer(rho[2]) +
         one(tau[3]) * zer(tau[2]) * one(rho[3]) * one(rho[2]) +
         one(tau[3]) * one(tau[2]) * zer(rho[3]) * zer(rho[2]));
    // 80~95 gates of first 10 rounds of each block is
    // relay for I (rev shift row)
    // rrrr 101 0000 <-> rrrr 0 0000
    // rrrr 101 0100 <-> rrrr 0 0100
    // rrrr 101 1000 <-> rrrr 0 1000
    // rrrr 101 1100 <-> rrrr 0 1100
    // rrrr 101 0010 <-> rrrr 0 1010
    // rrrr 101 0110 <-> rrrr 0 1110
    // rrrr 101 1010 <-> rrrr 0 0010
    // rrrr 101 1110 <-> rrrr 0 0110
    //
    // rrrr 101 0001 <-> rrrr 0 1101
    // rrrr 101 0101 <-> rrrr 0 0001
    // rrrr 101 1001 <-> rrrr 0 0101
    // rrrr 101 1101 <-> rrrr 0 1001
    //
    // rrrr 101 0011 <-> rrrr 0 0111
    // rrrr 101 0111 <-> rrrr 0 1011
    // rrrr 101 1011 <-> rrrr 0 1111
    // rrrr 101 1111 <-> rrrr 0 0011
    let pred_relay_i = block_eq2 * shiftrow_rho * zer_sgm *
        round_0_9_eq2 *
        zer(rho[4]) * one(tau[6]) * zer(tau[5]) * one(tau[4]);
    res[1] += pred_relay_i;

    // 112-128
    // linear gates: z := m × 0x11b⁻¹ + a × 2 × 0x11b⁻¹
    // a_ofs = r * 16 * 2 + 16;
    // m_ofs = r * 16 + 10 * 32;
    // z: zzzz111xxxx, a: aaaa1xxxx, m: 1mmmmxxxx
    //    0000            0000           0100
    //    0001            0001           0101
    //    0010            0010           0110
    //    0011            0011           0111
    //    0100            0100           1000
    //    0101            0101           1001
    //    0110            0110           1010
    //    0111            0111           1011
    //    1000            1000           1100
    let pred_linear_z = block_eq3 * gate_bit_4_eq3 *
        one(tau[4]) * one(tau[5]) * one(tau[6]) * one(rho[4]) * one(sgm[8]) *
        (eq3(tau[7], rho[5], sgm[4]) * eq3(tau[8], rho[6], sgm[5]) *
         (eq0!(+tau[10], +tau[9], +rho[8], +rho[7], +sgm[7], -sgm[6]) +
          eq0!(+tau[10], -tau[9], +rho[8], -rho[7], -sgm[7], +sgm[6])) +
         eq0!(-tau[10], +tau[9], +tau[8], +tau[7],
              -rho[8], +rho[7], +rho[6], +rho[5],
              -sgm[7], -sgm[6], +sgm[5], +sgm[4]));
    res[1] += pred_linear_z * utils[5];
    res[2] += pred_linear_z * utils[6];

    // key schedule
    let mut pred_wi_wa = T::from(0);
    let mut pred_wj = T::from(0);
    let mut pred_const = T::from(0);
    let mut pred_relay = T::from(0);
    let mut pred_linear = T::from(0);
    let mut pred_sum = T::from(0);
    let mut pred_sum_const = T::from(0);
    for r in 0..10 {
        let mut ofs = 10 * 8 * 4 + (r << 4);
        let k = expand_binary(ofs, 9);
        let wi = expand_binary(((512 * 7 - 48) >> 2) + r, 10);
        let mut tmp = T::from(1);
        for j in 0..9 {
            tmp *= T::from(1) + k[j] + tau[j + 2];
        }
        for j in 0..10 {
            tmp *= T::from(1) + wi[j] + rho[j + 2];
            tmp *= T::from(1) + wi[j] + sgm[j + 2];
        }
        pred_wi_wa += tmp;

        ofs += 1;
        let k = expand_binary(ofs, 9);
        let wa = expand_binary(((512 * 8 - 48) >> 2) + r, 10);
        let mut tmp = T::from(1);
        for j in 0..9 {
            tmp *= T::from(1) + k[j] + tau[j + 2];
        }
        for j in 0..10 {
            tmp *= T::from(1) + wa[j] + rho[j + 2];
            tmp *= T::from(1) + wa[j] + sgm[j + 2];
        }
        pred_wi_wa += tmp;

        ofs += 3;
        let k = expand_binary(ofs, 9);
        let w3 = expand_binary(128 * (r / 3 + 2) - 12 + (r % 3) * 4 + 3, 10);
        let mut tmp = T::from(1);
        for j in 0..9 {
            tmp *= T::from(1) + k[j] + tau[j + 2];
        }
        tmp *= utils[6];
        pred_const += tmp;
        for j in 0..10 {
            tmp *= T::from(1) + wi[j] + rho[j + 2];
            tmp *= T::from(1) + w3[j] + sgm[j + 2];
        }
        pred_wj += tmp;

        ofs += 1;
        let k = expand_binary(ofs, 9);
        let mut tmp = T::from(1);
        for j in 0..9 {
            tmp *= T::from(1) + k[j] + tau[j + 2];
        }
        for j in 0..10 {
            tmp *= T::from(1) + wi[j] + rho[j + 2];
        }
        pred_relay += tmp;

        ofs += 1;
        let k = expand_binary(ofs, 9);
        let mut tmp = T::from(1);
        for j in 0..9 {
            tmp *= T::from(1) + k[j] + tau[j + 2];
        }
        for j in 0..10 {
            tmp *= T::from(1) + w3[j] + rho[j + 2];
        }
        pred_relay += tmp;

        ofs += 1;
        let k = expand_binary(ofs, 9);
        let mut tmp = T::from(1);
        for j in 0..9 {
            tmp *= T::from(1) + k[j] + tau[j + 2];
        }
        pred_const += tmp * utils[4];
        for j in 0..10 {
            tmp *= T::from(1) + wi[j] + rho[j + 2];
            tmp *= T::from(1) + wa[j] + sgm[j + 2];
        }
        pred_linear += tmp;

        ofs += 1;
        for k in 0..4 {
            let k0 = k;
            let k = expand_binary(ofs + k0, 9);
            let wn = expand_binary(128 * ((r + 1) / 3 + 2) - 12 + ((r + 1) % 3) * 4 + k0, 10);
            let w = expand_binary(128 * (r / 3 + 2) - 12 + (r % 3) * 4 + k0, 10);

            let mut tmp = T::from(1);
            for j in 0..9 {
                tmp *= T::from(1) + k[j] + tau[j + 2];
            }
            let mut tmp1 = T::from(1);
            for j in 0..10 {
                tmp1 *= T::from(1) + wn[j] + rho[j + 2];
            }
            let mut tmp2 = T::from(1);
            for j in 0..10 {
                tmp2 *= T::from(1) + w[j] + rho[j + 2];
            }

            if k0 == 0 {
                pred_sum_const += tmp * utils[7 + r];
                pred_relay += tmp * (tmp1 + tmp2);
                let wa = expand_binary(((512 * 8 - 48) >> 2) + r, 10);
                let mut tmp3 = T::from(1);
                for j in 0..10 {
                    tmp3 *= T::from(1) + wa[j] + rho[j + 2];
                }
                pred_sum += tmp * tmp3 *
                    (eq0!(+tau[0], +tau[1], -rho[0], +rho[1]) +
                     eq0!(-tau[0], +tau[1], +rho[0], -rho[1]) +
                     eq0!(+tau[0], -tau[1], -rho[0], -rho[1]) +
                     eq0!(-tau[0], -tau[1], +rho[0], +rho[1]));
            } else {
                let wn = expand_binary(128 * ((r + 1) / 3 + 2) - 12 + ((r + 1) % 3) * 4 + k0 - 1, 10);
                let mut tmp3 = T::from(1);
                for j in 0..10 {
                    tmp3 *= T::from(1) + wn[j] + rho[j + 2];
                }
                pred_relay += tmp * (tmp1 + tmp2 + tmp3);
            }
        }
    }
    pred_wi_wa *= block_0_hi_eq3 * gate_bit_2_eq3;
    pred_linear *= block_0_hi_eq3 * gate_bit_2_eq3;
    pred_wj *= block_0_hi_eq3 * gate_bit_2_eq3;
    pred_relay *= gate_bit_2_eq2;
    pred_relay += pred_sum;
    pred_relay *= block_0_hi_eq2 * zer_sgm;
    res[0] += pred_wi_wa + pred_wj;
    res[1] += pred_wi_wa + pred_relay + pred_linear * utils[2];
    res[2] += pred_linear * utils[3];
    res[3] += zer_sgm * zer_rho * zer_hi_tau *
        (pred_const + pred_sum_const * zer(tau[0]) * zer(tau[1]));


    // fin values
    // a_ofs = 9 * 16 * 2 + 16 =      0b100110000 per block
    // w_ofs = 512 * 5 - 48 + 16 = 0b100111100000 fix block
    // p_ofs = 512 - 48 + 16 =     0b000111100000 fix block
    let mut block_0_7_rho = T::from(1);
    for i in 0..(log_n_aes - 3) {
        block_0_7_rho *= zer(rho[12 + i]);
    }
    let pred_sum = zer_sgm * gate_bit_4_eq2 * round_15 *
        zer(tau[6]) * zer(tau[5]) * zer(tau[4]) * one(rho[8]) * one(rho[5]) *
        (block_eq2 *
         one(rho[4]) * zer(rho[6]) * zer(rho[7]) +
         zer(rho[4]) * one(rho[6]) * one(rho[7]) *
         zer(rho[9]) * zer(rho[10]) * block_0_7_rho);
    res[1] += pred_sum;
    res
}


pub fn build_circuit<'a, T: Field>(log_n_aes: usize) -> Circuit<'a, T> {
    assert!(log_n_aes >= 3);
    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);
    {
        c.layers.push(Layer::<T>::default());
        // 1st 512 gates: IAM0,IAM1,IAM2,IAM3,...,IAM9,IA10
        // 2nd 512 gates: 16 nonce, 16 plaintext, 44 × 4 W, 10 × 4 A, 10 × 4 M
        c.layers[depth].bit_len = 9 + log_n_aes;
        // TODO: XXX useless gates at input layer
        for i in 0..(512 << log_n_aes) {
            c.layers[depth].gates.push(Gate::default());
            c.layers[depth].gates[i].gtype = GateType::IN;
            c.layers[depth].gates[i].inputs.push(i as u32);
            c.layers[depth].values.push(T::from(0));
        }
        depth += 1;
    }

    let mut vanishing = [
        T::from(1 << 0), T::from(1 << 1), T::from(1 << 2), T::from(1 << 3),
        T::from(1 << 4), T::from(1 << 5), T::from(1 << 6), T::from(1 << 7)];

    // first layer

    let round_log_size = 3 + 4;
    let block_log_size = round_log_size + 4;
    let log_size = block_log_size + log_n_aes;
    // length 8 * 16 * 10 + 4 * 16 * 10 + 16 < 8 * 16 * 16
    //
    // aes round: 8 × 16 gates,
    // 16 × 3 gates for I A M(N) 8-bit tests, N for nonce+plaintext
    // 16 lin gates: b := (0x63 + a + 0x1f × i) / 0x101
    // 16 sum mul linear gates:
    //               j := (sᵢ₋₁ × i_(reverse row shift) + 1) × 0x11b⁻¹
    //               where j+k := ((j%4 + k) % 4 + j/4 × 4)
    // 16 rel gates: i (reverse row shift)
    // 16 sum gates: sᵢ₋₁(j) := w(j) + m(j) + m(j+1) + a(j+1) + a(j+2) + a(j+3)
    // 16 lin gates: z := m × 0x11b⁻¹ + a × 2 × 0x11b⁻¹
    //
    // key schedule round: 4 * 16 gates
    //  4 × 2 gates for WI WA 8-bit tests
    //  4 rel gates: wi
    //  4 rel gates: w₄ₖ₊₃
    //  4 lin gates: wj := (w₄ₖ₊₃ × wi + 1) × 0x11b⁻¹
    //  4 lin gates: wb := (0x63 + wa + 0x1f × wi) / 0x101
    // 16 sum gates: w₄ₖ w₄ₖ₊₁ w₄ₖ₊₂ w₄ₖ₊₃
    c.layers.push(Layer::<T>::default());
    let layer = &mut c.layers[depth];
    layer.bit_len = log_size;
    layer.gates = vec![Gate::default(); 1 << log_size];
    layer.values = vec![T::from(0); 1 << log_size];

    // constants
    layer.utils.push(T::from(0));
    layer.utils.push(vanishing[0]);
    for i in 1..8 {
        vanishing[i] = vanishing[i] * (vanishing[i] - vanishing[0]);
    }
    let inv_101 = T::from(0x101).inv();
    layer.utils.push(inv_101 * T::from(0x1f));  // # 2
    layer.utils.push(inv_101);                  // # 3
    layer.utils.push(inv_101 * T::from(0x63));  // # 4
    let inv_11b = T::from(0x11b).inv();
    layer.utils.push(inv_11b * T::from(2));     // # 5
    layer.utils.push(inv_11b);                  // # 6
    let rcon = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36];
    for i in 0..10 {                            // # 7 + i
        layer.utils.push(T::from(rcon[i]));
    }
    for p in 0..(1 << log_n_aes) {              // # 17 + 4 * p + 0/1/2/3
        layer.utils.push(T::from(p & 0xff));
        layer.utils.push(T::from(p >> 8));
    }

    for p in 0..(1 << log_n_aes) {
        for r in 0..10 {
            let i_ofs: u32 = (p as u32) * 512 + (r as u32) * 16 * 2;
            let a_ofs: u32 = (p as u32) * 512 + (r as u32) * 16 * 2 + 16;
            let m_ofs: u32 = (p as u32) * 512 + (r as u32) * 16 + 10 * 32;

            let n_ofs: u32 = 512 - 48;
            // let p_ofs: u32 = 512 - 48 + 16;
            // let n_ofs: u32 = (p as u32) * 512 + (r as u32) * 16 + 512;                 // nonce
            // let w_ofs: u32 = (p as u32) * 512 + (r as u32) * 16 + 512 + 32;            // w0,w1,w2,w3
            // let w3_ofs: u32 = (p as u32) * 512 + (r as u32) * 16 + 512 + 32 + 12;      // w4r+3
            let w_ofs: u32 =
                512 * (r / 3 + 2) as u32 - 48 + (r % 3) as u32 * 16;
            let w_next_ofs: u32 = 512 * ((r + 1) / 3 + 2) as u32 - 48 + ((r + 1) % 3) as u32 * 16;
            let w3_ofs: u32 = w_ofs + 12;
            // let wi_ofs: u32 = (p as u32) * 512 + (r as u32) * 4 + 512 + 32 + 11 * 16;
            // let wa_ofs: u32 = (p as u32) * 512 + (r as u32) * 4 + 512 + 32 + 11 * 16 + 10 * 4;
            let wi_ofs: u32 = 512 * 7 - 48 + (r as u32) * 4;
            let wa_ofs: u32 = 512 * 8 - 48 + (r as u32) * 4;

            // 8-bit test for I, A
            let mut round_ofs: usize = (r << round_log_size) + (p << block_log_size);
            for k in 0..(16 * 2) {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::ADDMUL;
                let k = k as u32;
                g.inputs = vec![i_ofs + k, i_ofs + k];
                g.utils = vec![0, 1];
            }
            round_ofs += 16 * 2;

            // 8-bit test for M or key+nonce
            if r != 9 {
                for k in 0..16 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    g.inputs = vec![m_ofs + k, m_ofs + k];
                    g.utils = vec![0, 1];
                }
            } else if p == 0 {
                // nonce
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    g.inputs = vec![n_ofs + k, n_ofs + k];
                    g.utils = vec![0, 1];
                }
                let round_ofs = round_ofs + 4;
                // key: w0, w1, w2, w3
                let w_ofs: u32 = 512 * 2 - 48;
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    g.inputs = vec![w_ofs + k, w_ofs + k];
                    g.utils = vec![0, 1];
                }
            }
            round_ofs += 16;

            // linear gates: b := (0x63 + a + 0x1f × i) / 0x101
            for k in 0..16 {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::FLINEAR;
                let k = k as u32;
                g.inputs = vec![i_ofs + k, a_ofs + k];
                g.utils = vec![2, 3, 4];
            }
            round_ofs += 16;

            // sum-linear gates:
            // J := (sᵢ₋₁ × i_(reverse row shift) + 1) × 0x11b⁻¹
            if r != 0 {
                for k in 0..16 {
                    // let w_ofs = w_ofs - 16u32;
                    let a_ofs = a_ofs - 32u32;
                    let m_ofs = m_ofs - 16u32;
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUMLIN;
                    let k0 = k as u32;
                    let k1 = (k0 % 4 + 1) % 4 + k0 / 4 * 4;
                    let k2 = (k0 % 4 + 2) % 4 + k0 / 4 * 4;
                    let k3 = (k0 % 4 + 3) % 4 + k0 / 4 * 4;
                    let shift_k = ((((k / 4) - k % 4) % 4) * 4 + (k % 4)) as u32;
                    g.inputs = vec![
                        w_ofs + k0, m_ofs + k0, m_ofs + k1,
                        a_ofs + k1, a_ofs + k2, a_ofs + k3,
                        i_ofs + shift_k
                    ];
                    g.utils = vec![6, 6, 0];
                }
            }
            else {
                for k in 0..16 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUMLIN;
                    let k = k as u32;
                    let shift_k = ((((k / 4) - k % 4) % 4) * 4 + (k % 4)) as u32;
                    g.inputs = vec![
                        w_ofs + k,  // w
                        n_ofs + k,  // nonce
                        i_ofs + shift_k
                    ];
                    g.utils = if k == 0 {
                            vec![6, 6, 17 + 2 * (p as u32)]
                        } else if k == 1 {
                            vec![6, 6, 17 + 2 * (p as u32) + 1]
                        } else {
                            vec![6, 6, 0]
                    };
                }
            }
            round_ofs += 16;

            // relay for I (rev shift row)
            for k in 0..16 {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let shift_k = ((((k / 4) - k % 4) % 4) * 4 + (k % 4)) as u32;
                g.inputs = vec![i_ofs + shift_k];
            }
            round_ofs += 16;

            // sum gates:
            // sᵢ(j) := w(j) + m(j) + m(j+1) + a(j+1) + a(j+2) + a(j+3)
            // where j+k := ((j%4 + k) % 4 + j/4 × 4)
            if r != 0 {
                for k in 0..16 {
                    // let w_ofs = w_ofs - 16u32;
                    let a_ofs = a_ofs - 32u32;
                    let m_ofs = m_ofs - 16u32;
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUM;
                    let k0 = k as u32;
                    let k1 = (k0 % 4 + 1) % 4 + k0 / 4 * 4;
                    let k2 = (k0 % 4 + 2) % 4 + k0 / 4 * 4;
                    let k3 = (k0 % 4 + 3) % 4 + k0 / 4 * 4;
                    g.inputs = vec![
                        w_ofs + k0, m_ofs + k0, m_ofs + k1,
                        a_ofs + k1, a_ofs + k2, a_ofs + k3
                    ];
                }
            } else {
                {
                    let k = 0;
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUMCONS;
                    let k = k as u32;
                    g.inputs = vec![
                        n_ofs + k,      // nonce
                        w_ofs + k,      // w
                    ];
                    g.utils = vec![17 + 2 * (p as u32) + 0];
                }
                {
                    let k = 1;
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUMCONS;
                    let k = k as u32;
                    g.inputs = vec![
                        n_ofs + k,      // nonce
                        w_ofs + k,      // w
                    ];
                    g.utils = vec![17 + 2 * (p as u32) + 1];
                }
                for k in 2..16 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUM;
                    let k = k as u32;
                    g.inputs = vec![
                        n_ofs + k,      // nonce
                        w_ofs + k,      // w
                    ];
                }
            }
            round_ofs += 16;

            // linear gates: z := m × 0x11b⁻¹ + a × 2 × 0x11b⁻¹
            if r != 9 {
                for k in 0..16 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::LINEAR;
                    let k = k as u32;
                    g.inputs = vec![a_ofs + k, m_ofs + k];
                    g.utils = vec![5, 6];
                }
            }
            // round_ofs += 16;

            if p == 0 {
                let mut round_ofs: usize = r << (round_log_size - 1);
                round_ofs += 10 * 8 * 16;
                // key schedule
                //
                // 8-bit test for WI, WA
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    g.inputs = vec![wi_ofs + k, wi_ofs + k];
                    g.utils = vec![0, 1];
                }
                round_ofs += 4;
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    g.inputs = vec![wa_ofs + k, wa_ofs + k];
                    g.utils = vec![0, 1];
                }
                round_ofs += 12;

                // wj := (w₄ᵣ₊₃ × wi + 1) × 0x11b⁻¹
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUMLIN;
                    let k = k as u32;
                    g.inputs = vec![wi_ofs + k, w3_ofs + k];
                    g.utils = vec![6, 6, 0];
                }
                round_ofs += 4;

                // relay for WI
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::RELAY;
                    let k = k as u32;
                    // let shift_k = (k / 4) * 4 + (k % 4 + 3) % 4 as u32;
                    g.inputs = vec![wi_ofs + k];
                }
                round_ofs += 4;

                // relay for W₄ᵣ₊₃
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::RELAY;
                    let k = k as u32;
                    // let shift_k = (k / 4) * 4 + (k % 4 + 3) % 4 as u32;
                    g.inputs = vec![w3_ofs + k];
                }
                round_ofs += 4;

                // lin gates: wb := (0x63 + wa + 0x1f × wi) / 0x101
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::FLINEAR;
                    let k = k as u32;
                    g.inputs = vec![wi_ofs + k, wa_ofs + k];
                    g.utils = vec![2, 3, 4];
                }
                round_ofs += 4;

                // roundkey
                for k in 0..16 {
                    // w4r = w4r-4 + sub(rot(w4r-1))
                    // w4r+i = w4r+i-4 + w4r+i-1
                    let g = &mut layer.gates[round_ofs + k];
                    let k = k as u32;
                    if k == 0 {
                        g.gtype = GateType::SUMCONS;
                        g.inputs = vec![w_next_ofs, wa_ofs + 1, w_ofs + k];
                        g.utils = vec![(7 + r) as u32];
                    } else if k < 4 {
                        g.gtype = GateType::SUM;
                        g.inputs = vec![w_next_ofs + k, wa_ofs + (k + 1) % 4, w_ofs + k];
                    } else {
                        g.gtype = GateType::SUM;
                        g.inputs = vec![w_next_ofs + k, w_next_ofs + k - 4, w_ofs + k];
                    }
                }
            }
            // round_ofs += 16;
        }

        // fin
        let round_ofs = (15 << round_log_size) + (p << block_log_size);
        let a_ofs: u32 = 9 * 16 * 2 + 16 + (p as u32) * 512;
        let w_ofs: u32 = 512 * 5 - 48 + 16;
        // let w_ofs: u32 = 10 * 16 + 512 + 32;
        // let p_ofs: u32 = 512 + 16;
        let p_ofs: u32 = 512 - 48 + 16;
        for k in 0..16 {
            let g = &mut layer.gates[round_ofs + k];
            g.gtype = GateType::SUM;
            let k = k as u32;
            g.inputs = vec![w_ofs + k, a_ofs + k, p_ofs + k];
            // a + w + plaintext = ciphertext
        }
    }
    let utils = layer.utils.clone();
    layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
            pred_ctraes_first_layer(log_n_aes, &utils, tau, rho, sgm)));
    depth += 1;

    // mid layers
    for d in 1..8 {
        let round_log_size = 3 + 4;
        let block_log_size = round_log_size + 4;
        let log_size = block_log_size + log_n_aes;
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_size;
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];

        // constants
        layer.utils.push(vanishing[d]);
        layer.utils.push(vanishing[d - 1]);
        for i in (d + 1)..8 {
            vanishing[i] = vanishing[i] * (vanishing[i] - vanishing[d]);
        }

        for p in 0..(1 << log_n_aes) {
            for r in 0..10 {
                let mut ofs: usize = (r << round_log_size) + (p << block_log_size);
                for k in 0..(16 * 3) {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::P1D2C0;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k];
                    g.utils = vec![0];
                }
                ofs += 16 * 3;

                if d < 5 {
                    for k in 0..16 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::P1D2C0;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                        g.utils = vec![1];
                    }
                } else {
                    for k in 0..16 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::RELAY;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                    }
                }
                ofs += 16;

                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::P1D2C0;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k];
                    g.utils = vec![1];
                }
                ofs += 16;

                for k in 0..(16 * 2) {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::RELAY;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k];
                }
                ofs += 16 * 2;

                if d < 2 {
                    for k in 0..16 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::P1D2C0;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                        g.utils = vec![1];
                    }
                } else {
                    for k in 0..16 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::RELAY;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                    }
                }
                // ofs += 16;

                // key schedule
                if p == 0 {
                    let mut ofs: usize = r << (round_log_size - 1);
                    ofs += 10 * 8 * 16;

                    // 8-bit test for WI, WA
                    for k in 0..(4 * 2) {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::P1D2C0;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                        g.utils = vec![0];
                    }
                    // round_ofs += 4 * 2;
                    ofs += 16;

                    // wj := (w₄ᵣ₊₃ × wi + 1) × 0x11b⁻¹
                    for k in 0..4 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::P1D2C0;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                        g.utils = vec![1];
                    }
                    ofs += 4;

                    for k in 0..(4 * 2) {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::RELAY;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                    }
                    ofs += 8;

                    if d < 5 {
                        for k in 0..16 {
                            let g = &mut layer.gates[ofs + k];
                            g.gtype = GateType::P1D2C0;
                            let k = k as u32;
                            let ofs = ofs as u32;
                            g.inputs = vec![ofs + k];
                            g.utils = vec![1];
                        }
                    } else {
                        for k in 0..16 {
                            let g = &mut layer.gates[ofs + k];
                            g.gtype = GateType::RELAY;
                            let k = k as u32;
                            let ofs = ofs as u32;
                            g.inputs = vec![ofs + k];
                        }
                    }
                    ofs += 4;

                    for k in 0..16 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::RELAY;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                    }
                    // ofs += 16
                }
            }
            let ofs = (15 << round_log_size) + (p << block_log_size);
            for k in 0..16 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
        }
        layer.val_degree = vec![2];
        layer.num_table = vec![2];
        layer.num_phase = 1;
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_d2_no_const);
        let utils = layer.utils.clone();
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
            pred_ctraes_mid_layer(log_n_aes, &utils, d, tau, rho, sgm)));
        depth += 1;
    }

    {
        let round_log_size = 3 + 4;
        let block_log_size = round_log_size + 4;
        let log_size = block_log_size + log_n_aes;
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_size;
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];

        for p in 0..(1 << log_n_aes) {
            for r in 0..10 {
                let mut ofs: usize = (r << round_log_size) + (p << block_log_size);
                for k in 0..(16 * 4) {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::RELAY;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k];
                }
                ofs += 16 * 4;

                // J * I
                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::MUL;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k, ofs + k + 16];
                }
                ofs += 16;

                // J * S
                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::MUL;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k - 16, ofs + k + 16];
                }
                ofs += 16;
                ofs += 16;

                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::RELAY;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k];
                }
                // ofs += 16;

                // if p == 0 {
                if p == 0 {
                    let mut ofs: usize = r << (round_log_size - 1);
                    ofs += 10 * 8 * 16;
                    for k in 0..8 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::RELAY;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                    }
                    // ofs += 8;
                    ofs += 16;

                    // J * I
                    for k in 0..4 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::MUL;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k, ofs + k + 4];
                    }
                    ofs += 4;

                    // J * S
                    for k in 0..4 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::MUL;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k - 4, ofs + k + 4];
                    }
                    ofs += 4;
                    ofs += 4;
                    for k in 0..4 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::RELAY;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                    }
                    ofs += 4;

                    for k in 0..16 {
                        let g = &mut layer.gates[ofs + k];
                        g.gtype = GateType::RELAY;
                        let k = k as u32;
                        let ofs = ofs as u32;
                        g.inputs = vec![ofs + k];
                    }
                    // // ofs += 16;
                }
            }

            let ofs = (15 << round_log_size) + (p << block_log_size);
            for k in 0..16 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
        }

        c.layers[depth].predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
            pred_ctraes_last_layer(log_n_aes, tau, rho, sgm)));
        depth += 1;
    }

    c.depth = depth;
    c
}

pub fn build_single_circuit<'a, T: Field>() -> Circuit<'a, T> {
    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);
    {
        c.layers.push(Layer::<T>::default());
        // 1st 512 gates: IAM0,IAM1,IAM2,IAM3,...,IAM9,IA10
        // 2nd 512 gates: 16 nonce, 16 plaintext, 44 × 4 W, 10 × 4 A, 10 × 4 M
        c.layers[depth].bit_len = 10;
        // TODO: XXX useless gates at input layer
        for i in 0..1024 {
            c.layers[depth].gates.push(Gate::default());
            c.layers[depth].gates[i].gtype = GateType::IN;
            c.layers[depth].gates[i].inputs.push(i as u32);
            c.layers[depth].values.push(T::from(0));
        }
        depth += 1;
    }

    let mut vanishing = [
        T::from(1 << 0), T::from(1 << 1), T::from(1 << 2), T::from(1 << 3),
        T::from(1 << 4), T::from(1 << 5), T::from(1 << 6), T::from(1 << 7)];

    // first layer
    {
        let round_log_size = 3 + 4;
        let log_size = round_log_size + 4;
        // length 8 * 16 * 10 + 4 * 16 * 10 + 16 < 8 * 16 * 16
        //
        // aes round: 8 × 16 gates,
        // 16 × 3 gates for I A M(N) 8-bit tests, N for nonce+plaintext
        // 16 lin gates: b := (0x63 + a + 0x1f × i) / 0x101
        // 16 sum mul linear gates:
        //               j := (sᵢ₋₁ × i_(reverse row shift) + 1) × 0x11b⁻¹
        //               where j+k := ((j%4 + k) % 4 + j/4 × 4)
        // 16 rel gates: i (reverse row shift)
        // 16 sum gates: sᵢ₋₁(j) := w(j) + m(j) + m(j+1) + a(j+1) + a(j+2) + a(j+3)
        // 16 lin gates: z := m × 0x11b⁻¹ + a × 2 × 0x11b⁻¹
        //
        // key schedule round: 4 * 16 gates
        //  4 × 2 gates for WI WA 8-bit tests
        //  4 rel gates: wi
        //  4 rel gates: w₄ₖ₊₃
        //  4 lin gates: wj := (w₄ₖ₊₃ × wi + 1) × 0x11b⁻¹
        //  4 lin gates: wb := (0x63 + wa + 0x1f × wi) / 0x101
        // 16 sum gates: w₄ₖ w₄ₖ₊₁ w₄ₖ₊₂ w₄ₖ₊₃
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_size;
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];

        // constants
        layer.utils.push(T::from(0));
        layer.utils.push(vanishing[0]);
        for i in 1..8 {
            vanishing[i] = vanishing[i] * (vanishing[i] - vanishing[0]);
        }
        let inv_101 = T::from(0x101).inv();
        layer.utils.push(inv_101 * T::from(0x1f));  // # 2
        layer.utils.push(inv_101);                  // # 3
        layer.utils.push(inv_101 * T::from(0x63));  // # 4
        let inv_11b = T::from(0x11b).inv();
        layer.utils.push(inv_11b * T::from(2));     // # 5
        layer.utils.push(inv_11b);                  // # 6
        let rcon = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36];
        for i in 0..10 {                            // # 7 + i
            layer.utils.push(T::from(rcon[i]));
        }

        for r in 0..10 {
            let i_ofs: u32 = (r as u32) * 16 * 2;
            let a_ofs: u32 = (r as u32) * 16 * 2 + 16;
            // let m_ofs: u32 = (r as u32) * 16 * 3 + (16 * 2);
            let m_ofs: u32 = (r as u32) * 16 + 10 * 32;
            let n_ofs: u32 = (r as u32) * 16 + 512;                 // nonce
            let w_ofs: u32 = (r as u32) * 16 + 512 + 32;            // w0,w1,w2,w3
            let w3_ofs: u32 = (r as u32) * 16 + 512 + 32 + 12;      // w4r+3
            let wi_ofs: u32 = (r as u32) * 4 + 512 + 32 + 11 * 16;
            let wa_ofs: u32 = (r as u32) * 4 + 512 + 32 + 11 * 16 + 10 * 4;

            // 8-bit test for I, A
            let mut round_ofs: usize = r << round_log_size;
            for k in 0..(16 * 2) {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::ADDMUL;
                let k = k as u32;
                g.inputs = vec![i_ofs + k, i_ofs + k];
                g.utils = vec![0, 1];
            }
            round_ofs += 16 * 2;

            // 8-bit test for M or key+nonce
            if r != 9 {
                for k in 0..16 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    g.inputs = vec![m_ofs + k, m_ofs + k];
                    g.utils = vec![0, 1];
                }
            } else {
                // nonce
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    g.inputs = vec![n_ofs + k, n_ofs + k];
                    g.utils = vec![0, 1];
                }
                let round_ofs = round_ofs + 4;
                // key: w0, w2, w2, w3
                for k in 0..4 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    g.inputs = vec![w_ofs + k, w_ofs + k];
                    g.utils = vec![0, 1];
                }
            }
            round_ofs += 16;

            // linear gates: b := (0x63 + a + 0x1f × i) / 0x101
            for k in 0..16 {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::FLINEAR;
                let k = k as u32;
                g.inputs = vec![i_ofs + k, a_ofs + k];
                g.utils = vec![2, 3, 4];
            }
            round_ofs += 16;

            // sum-linear gates:
            // J := (sᵢ₋₁ × i_(reverse row shift) + 1) × 0x11b⁻¹
            if r != 0 {
                for k in 0..16 {
                    // let w_ofs = w_ofs - 16u32;
                    let a_ofs = a_ofs - 32u32;
                    let m_ofs = m_ofs - 16u32;
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUMLIN;
                    let k0 = k as u32;
                    let k1 = (k0 % 4 + 1) % 4 + k0 / 4 * 4;
                    let k2 = (k0 % 4 + 2) % 4 + k0 / 4 * 4;
                    let k3 = (k0 % 4 + 3) % 4 + k0 / 4 * 4;
                    let shift_k = ((((k / 4) - k % 4) % 4) * 4 + (k % 4)) as u32;
                    g.inputs = vec![
                        w_ofs + k0, m_ofs + k0, m_ofs + k1,
                        a_ofs + k1, a_ofs + k2, a_ofs + k3,
                        i_ofs + shift_k
                    ];
                    g.utils = vec![6, 6, 0];
                }
            } else {
                for k in 0..16 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUMLIN;
                    let k = k as u32;
                    let shift_k = ((((k / 4) - k % 4) % 4) * 4 + (k % 4)) as u32;
                    g.inputs = vec![
                        512 + k,        // nonce
                        512 + 32 + k,   // w
                        i_ofs + shift_k
                    ];
                    g.utils = vec![6, 6, 0];
                }
            }
            round_ofs += 16;

            // relay for I (rev shift row)
            for k in 0..16 {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let shift_k = ((((k / 4) - k % 4) % 4) * 4 + (k % 4)) as u32;
                g.inputs = vec![i_ofs + shift_k];
            }
            round_ofs += 16;

            // sum gates:
            // sᵢ(j) := w(j) + m(j) + m(j+1) + a(j+1) + a(j+2) + a(j+3)
            // where j+k := ((j%4 + k) % 4 + j/4 × 4)
            if r != 0 {
                for k in 0..16 {
                    // let w_ofs = w_ofs - 16u32;
                    let a_ofs = a_ofs - 32u32;
                    let m_ofs = m_ofs - 16u32;
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUM;
                    let k0 = k as u32;
                    let k1 = (k0 % 4 + 1) % 4 + k0 / 4 * 4;
                    let k2 = (k0 % 4 + 2) % 4 + k0 / 4 * 4;
                    let k3 = (k0 % 4 + 3) % 4 + k0 / 4 * 4;
                    g.inputs = vec![
                        w_ofs + k0, m_ofs + k0, m_ofs + k1,
                        a_ofs + k1, a_ofs + k2, a_ofs + k3
                    ];
                }
            } else {
                for k in 0..16 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::SUM;
                    let k = k as u32;
                    g.inputs = vec![
                        n_ofs + k,      // nonce
                        w_ofs + k,      // w
                    ];
                }
            }
            round_ofs += 16;

            // linear gates: z := m × 0x11b⁻¹ + a × 2 × 0x11b⁻¹
            if r != 9 {
                for k in 0..16 {
                    let g = &mut layer.gates[round_ofs + k];
                    g.gtype = GateType::LINEAR;
                    let k = k as u32;
                    g.inputs = vec![a_ofs + k, m_ofs + k];
                    g.utils = vec![5, 6];
                }
            }
            // round_ofs += 16;


            let mut round_ofs: usize = r << (round_log_size - 1);
            round_ofs += 10 * 8 * 16;
            // key schedule
            //
            // 8-bit test for WI, WA
            for k in 0..(4 * 2) {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::ADDMUL;
                let k = k as u32;
                g.inputs = vec![wi_ofs + k, wi_ofs + k];
                g.utils = vec![0, 1];
            }
            // round_ofs += 4 * 2;
            round_ofs += 16;

            // wj := (w₄ᵣ₊₃ × wi + 1) × 0x11b⁻¹
            for k in 0..4 {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::SUMLIN;
                let k = k as u32;
                g.inputs = vec![wi_ofs + k, w3_ofs + k];
                g.utils = vec![6, 6, 0];
            }
            round_ofs += 4;

            // relay for WI
            for k in 0..4 {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                // let shift_k = (k / 4) * 4 + (k % 4 + 3) % 4 as u32;
                g.inputs = vec![wi_ofs + k];
            }
            round_ofs += 4;

            // relay for W₄ᵣ₊₃
            for k in 0..4 {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                // let shift_k = (k / 4) * 4 + (k % 4 + 3) % 4 as u32;
                g.inputs = vec![w3_ofs + k];
            }
            round_ofs += 4;

            // lin gates: wb := (0x63 + wa + 0x1f × wi) / 0x101
            for k in 0..4 {
                let g = &mut layer.gates[round_ofs + k];
                g.gtype = GateType::FLINEAR;
                let k = k as u32;
                g.inputs = vec![wi_ofs + k, wa_ofs + k];
                g.utils = vec![2, 3, 4];
            }
            round_ofs += 4;

            // roundkey
            for k in 0..16 {
                // w4r = w4r-4 + sub(rot(w4r-1))
                // w4r+i = w4r+i-4 + w4r+i-1
                let g = &mut layer.gates[round_ofs + k];
                let k = k as u32;
                if k == 0 {
                    g.gtype = GateType::SUMCONS;
                    g.inputs = vec![w_ofs + 16, wa_ofs + 1, w_ofs + k];
                    g.utils = vec![(7 + r) as u32];
                } else if k < 4 {
                    g.gtype = GateType::SUM;
                    g.inputs = vec![w_ofs + 16 + k, wa_ofs + (k + 1) % 4, w_ofs + k];
                } else {
                    g.gtype = GateType::SUM;
                    g.inputs = vec![w_ofs + 16 + k, w_ofs + 16 + k - 4, w_ofs + k];
                }
            }
            // round_ofs += 16;
        }

        let round_ofs = (15 << round_log_size) + 0;
        let a_ofs: u32 = 9 * 16 * 2 + 16;
        let w_ofs: u32 = 10 * 16 + 512 + 32;
        let p_ofs: u32 = 512 + 16;
        for k in 0..16 {
            let g = &mut layer.gates[round_ofs + k];
            g.gtype = GateType::SUM;
            let k = k as u32;
            g.inputs = vec![w_ofs + k, a_ofs + k, p_ofs + k];
            // a + w + plaintext = ciphertext
        }
        depth += 1;
    }

    // mid layers
    for d in 1..8 {
        let round_log_size = 3 + 4;
        let log_size = round_log_size + 4;
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_size;
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];

        // constants
        layer.utils.push(T::from(0));
        layer.utils.push(vanishing[d]);
        layer.utils.push(vanishing[d - 1]);
        for i in (d + 1)..8 {
            vanishing[i] = vanishing[i] * (vanishing[i] - vanishing[d]);
        }

        for r in 0..10 {
            let mut ofs: usize = r << round_log_size;
            for k in 0..(16 * 3) {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::ADDMUL;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k, ofs + k];
                g.utils = vec![0, 1];
            }
            ofs += 16 * 3;

            if d < 5 {
                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k, ofs + k];
                    g.utils = vec![0, 2];
                }
            } else {
                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::RELAY;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k];
                }
            }
            ofs += 16;

            for k in 0..16 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::ADDMUL;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k, ofs + k];
                g.utils = vec![0, 2];
            }
            ofs += 16;

            for k in 0..(16 * 2) {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
            ofs += 16;

            if d < 2 {
                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k, ofs + k];
                    g.utils = vec![0, 2];
                }
            } else {
                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::RELAY;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k];
                }
            }
            // ofs += 16;

            // key schedule
            let mut ofs: usize = r << (round_log_size - 1);
            ofs += 10 * 8 * 16;
            //
            // 8-bit test for WI, WA
            for k in 0..(4 * 2) {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::ADDMUL;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k, ofs + k];
                g.utils = vec![0, 1];
            }
            // round_ofs += 4 * 2;
            ofs += 16;

            // wj := (w₄ᵣ₊₃ × wi + 1) × 0x11b⁻¹
            for k in 0..4 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::ADDMUL;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k, ofs + k];
                g.utils = vec![0, 2];
            }
            ofs += 4;

            for k in 0..(4 * 2) {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
            ofs += 8;

            if d < 5 {
                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::ADDMUL;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k, ofs + k];
                    g.utils = vec![0, 2];
                }
            } else {
                for k in 0..16 {
                    let g = &mut layer.gates[ofs + k];
                    g.gtype = GateType::RELAY;
                    let k = k as u32;
                    let ofs = ofs as u32;
                    g.inputs = vec![ofs + k];
                }
            }
            ofs += 4;

            for k in 0..16 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
            // ofs += 16
        }
        let ofs = (15 << round_log_size) + 0;
        for k in 0..16 {
            let g = &mut layer.gates[ofs + k];
            g.gtype = GateType::RELAY;
            let k = k as u32;
            let ofs = ofs as u32;
            g.inputs = vec![ofs + k];
        }

        depth += 1;
    }

    {
        let round_log_size = 3 + 4;
        let log_size = round_log_size + 4;
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_size;
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];

        for r in 0..10 {
            let mut ofs: usize = r << round_log_size;
            for k in 0..(16 * 4) {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
            ofs += 16 * 4;

            // J * I
            for k in 0..16 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::MUL;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k, ofs + k + 16];
            }
            ofs += 16;

            // J * S
            for k in 0..16 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::MUL;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k - 16, ofs + k + 16];
            }
            ofs += 16;
            ofs += 16;

            for k in 0..16 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
            // ofs += 16;

            let mut ofs: usize = r << (round_log_size - 1);
            ofs += 10 * 8 * 16;
            for k in 0..8 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
            // ofs += 8;
            ofs += 16;

            // J * I
            for k in 0..4 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::MUL;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k, ofs + k + 4];
            }
            ofs += 4;

            // J * S
            for k in 0..4 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::MUL;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k - 4, ofs + k + 4];
            }
            ofs += 4;
            ofs += 4;
            for k in 0..4 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
            ofs += 4;

            for k in 0..16 {
                let g = &mut layer.gates[ofs + k];
                g.gtype = GateType::RELAY;
                let k = k as u32;
                let ofs = ofs as u32;
                g.inputs = vec![ofs + k];
            }
            // ofs += 16;
        }
        let ofs = (15 << round_log_size) + 0;
        for k in 0..16 {
            let g = &mut layer.gates[ofs + k];
            g.gtype = GateType::RELAY;
            let k = k as u32;
            let ofs = ofs as u32;
            g.inputs = vec![ofs + k];
        }

        depth += 1;
    }

    c.depth = depth;
    c
}
