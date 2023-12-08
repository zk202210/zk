use std::sync::Arc;

use crate::field::*;
use crate::circuit::*;
use crate::gkr::*;
use crate::util::*;

pub const ROTATE: [[usize; 5]; 5] = [
    [ 0, 36,  3, 41, 18],
    [ 1, 44, 10, 45,  2],
    [62,  6, 43, 15, 61],
    [28, 55, 25, 21, 56],
    [27, 20, 39,  8, 14],
];

pub const RC: [usize; 24] = [
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808A,
    0x8000000080008000,
    0x000000000000808B,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008A,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000A,
    0x000000008000808B,
    0x800000000000008B,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800A,
    0x800000008000000A,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008,
];



// sha3-256 merkle tree without padding
// leaf inputs are 512-bit each
pub fn build_circuit<'a, T: Field>(log_leaf_size: usize) -> Circuit<'a, T> {
    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);

    /*
    Round[b](A,RC) {
        # θ step
        C[x] = A[x,0] xor A[x,1] xor A[x,2] xor A[x,3] xor A[x,4],   for x in 0…4
        D[x] = C[x-1] xor rot(C[x+1],1),                             for x in 0…4
        A[x,y] = A[x,y] xor D[x],                           for (x,y) in (0…4,0…4)

        # ρ and π steps
        B[y,2*x+3*y] = rot(A[x,y], r[x,y]),                 for (x,y) in (0…4,0…4)

        # χ step
        A[x,y] = B[x,y] xor ((not B[x+1,y]) and B[x+2,y]),  for (x,y) in (0…4,0…4)

        # ι step
        A[0,0] = A[0,0] xor RC

        return A
        }
    */

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_leaf_size + 1 + 9;
        layer.bit_len = log_size;
        // TODO: XXX useless gates at input layer
        for i in 0..(1 << log_size) {
            layer.gates.push(Gate::default());
            layer.gates[i].gtype = GateType::IN;
            layer.gates[i].inputs.push(i as u32);
            layer.values.push(T::from(0));
        }
        depth += 1;
    }

    for d in 0..24 {
        // linear layers
        {
            c.layers.push(Layer::<T>::default());
            let layer = &mut c.layers[depth];
            // 1600 states => 2048
            let log_size = log_leaf_size + 1 + 11 + 1;
            layer.bit_len = log_size;
            // first round
            if d == 0 {
                layer.utils = vec![T::from(1)];
                layer.num_phase = 1;
                layer.num_table = vec![2];
                layer.val_degree = vec![2];
                layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_d2_no_const);
                layer.sumcheck_init[0] = Arc::new(
                    move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                        tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size - 3)]; 2];
                        let tbl_prev_p = &tbl.tbl_prev_point[0];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 9;
                            let mut i = k << 12;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let a_idx = (15 + y - 3 * x) * 3 % 5;
                                        let a_idy = x;
                                        let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;
                                        let c0_idx = (a_idx + 4) % 5;
                                        let c0_idz = a_idz;
                                        let c1_idx = (a_idx + 1) % 5;
                                        let c1_idz = (a_idz + 63) % 64;
                                        let mut tmp = vec![];
                                        if a_idx * 5 + a_idy < 8 {
                                            tbl.tbl_coeff[0][ofs + (a_idx * 5 + a_idy) * 64 + a_idz] += tbl_prev_p[i];
                                            tmp.push((ofs + (a_idx * 5 + a_idy) * 64 + a_idz) as u32);
                                        }
                                        for c_idy in 0..=4 {
                                            if c0_idx * 5 + c_idy < 8 {
                                                tbl.tbl_coeff[0][ofs + (c0_idx * 5 + c_idy) * 64 + c0_idz] += tbl_prev_p[i];
                                                tmp.push((ofs + (c0_idx * 5 + c_idy) * 64 + c0_idz) as u32);
                                            }
                                        }
                                        for c_idy in 0..=4 {
                                            if c1_idx * 5 + c_idy < 8 {
                                                tbl.tbl_coeff[0][ofs + (c1_idx * 5 + c_idy) * 64 + c1_idz] += tbl_prev_p[i];
                                                tmp.push((ofs + (c1_idx * 5 + c_idy) * 64 + c1_idz) as u32);
                                            }
                                        }
                                        i += 1;
                                    }
                                }
                            }
                            for z in 0..512 {
                                let u = (k << 9) + z;
                                let i = (k << 12) + 2048 + z;
                                if k < (1 << log_leaf_size) {
                                    tbl.tbl_coeff[1][u] += tbl_prev_p[i];
                                    tbl.tbl_coeff[0][u] -= tbl_prev_p[i];
                                } else {
                                    tbl.tbl_coeff[0][u] += tbl_prev_p[i];
                                }
                            }
                        }
                    });
                layer.evaluate[0] = Arc::new(
                    move |_: &Layer<T>, inputs: &Vec<T>| {
                        let mut res = vec![T::from(0); 1 << log_size];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 9;
                            let mut i = k << 12;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let a_idx = (15 + y - 3 * x) * 3 % 5;
                                        let a_idy = x;
                                        let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;
                                        let c0_idx = (a_idx + 4) % 5;
                                        let c0_idz = a_idz;
                                        let c1_idx = (a_idx + 1) % 5;
                                        let c1_idz = (a_idz + 63) % 64;
                                        if a_idx * 5 + a_idy < 8 {
                                            res[i] += inputs[ofs + (a_idx * 5 + a_idy) * 64 + a_idz];
                                        }
                                        for c_idy in 0..=4 {
                                            if c0_idx * 5 + c_idy < 8 {
                                                res[i] += inputs[ofs + (c0_idx * 5 + c_idy) * 64 + c0_idz];
                                            }
                                        }
                                        for c_idy in 0..=4 {
                                            if c1_idx * 5 + c_idy < 8 {
                                                res[i] += inputs[ofs + (c1_idx * 5 + c_idy) * 64 + c1_idz];
                                            }
                                        }
                                        i += 1;
                                    }
                                }
                            }
                            for z in 0..512 {
                                let u = (k << 9) + z;
                                let i = (k << 12) + 2048 + z;
                                if k < (1 << log_leaf_size) {
                                    res[i] = inputs[u] * (inputs[u] - T::from(1));
                                } else {
                                    res[i] += inputs[u];
                                }
                            }
                        }
                        res
                    });
                // TODO XXX, reuse tbl0 & tbl1
                layer.predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                        let mut sum = T::from(0);
                        let mut tbl0 = vec![T::from(0); 1 << 12];
                        let mut tbl1 = vec![T::from(0); 1 << 9];
                        tbl0[0] = T::from(1);
                        tbl1[0] = T::from(1);
                        for k in 0..12 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl0[i + mask] = tbl0[i] * tau[k];
                                tbl0[i] = tbl0[i] - tbl0[i + mask];
                            }
                        }
                        for k in 0..9 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl1[i + mask] = tbl1[i] * rho[k];
                                tbl1[i] = tbl1[i] - tbl1[i + mask];
                            }
                        }
                        let mut i = 0;
                        for x in 0..=4 {
                            for y in 0..=4 {
                                for z in 0..64 {
                                    let a_idx = (15 + y - 3 * x) * 3 % 5;
                                    let a_idy = x;
                                    let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;
                                    let c0_idx = (a_idx + 4) % 5;
                                    let c0_idz = a_idz;
                                    let c1_idx = (a_idx + 1) % 5;
                                    let c1_idz = (a_idz + 63) % 64;
                                    if a_idx * 5 + a_idy < 8 {
                                        sum += tbl0[i] * tbl1[(a_idx * 5 + a_idy) * 64 + a_idz];
                                    }
                                    for c_idy in 0..=4 {
                                        if c0_idx * 5 + c_idy < 8 {
                                            sum += tbl0[i] * tbl1[(c0_idx * 5 + c_idy) * 64 + c0_idz];
                                        }
                                    }
                                    for c_idy in 0..=4 {
                                        if c1_idx * 5 + c_idy < 8 {
                                            sum += tbl0[i] * tbl1[(c1_idx * 5 + c_idy) * 64 + c1_idz];
                                        }
                                    }
                                    i += 1;
                                }
                            }
                        }
                        let mut relay = one(tau[11]) * zer(tau[10]) * zer(tau[9]);
                        for i in 0..9 {
                            relay *= eq2(tau[i], rho[i]);
                        }
                        let mut block_eq = T::from(1);
                        for i in 12..(tau.len() - 1) {
                            block_eq *= eq2(tau[i], rho[i - 3]);
                        }
                        relay *= block_eq;
                        [relay * zer(tau[log_leaf_size + 12]) * zer(rho[log_leaf_size + 9]),
                         (relay + sum * block_eq) * eq2(tau[log_leaf_size + 12], rho[log_leaf_size + 9]),
                         T::from(0), T::from(0)]
                    }));

            // rest rounds
            } else {
                layer.num_phase = 1;
                layer.num_table = vec![1];
                layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
                layer.sumcheck_init[0] = Arc::new(
                    move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                        tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_size]];
                        let tbl_prev_p = &tbl.tbl_prev_point[0];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 12;
                            let mut i = ofs;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let a_idx = (15 + y - 3 * x) * 3 % 5;
                                        let a_idy = x;
                                        let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;
                                        let c0_idx = (a_idx + 4) % 5;
                                        let c0_idz = a_idz;
                                        let c1_idx = (a_idx + 1) % 5;
                                        let c1_idz = (a_idz + 63) % 64;
                                        tbl.tbl_coeff[0][ofs + (a_idx * 5 + a_idy) * 64 + a_idz] += tbl_prev_p[i];
                                        for c_idy in 0..=4 {
                                            tbl.tbl_coeff[0][ofs + (c0_idx * 5 + c_idy) * 64 + c0_idz] += tbl_prev_p[i];
                                        }
                                        for c_idy in 0..=4 {
                                            tbl.tbl_coeff[0][ofs + (c1_idx * 5 + c_idy) * 64 + c1_idz] += tbl_prev_p[i];
                                        }
                                        i += 1;
                                    }
                                }
                            }
                            for i in 0..512 {
                                let u = ofs + 2048 + i;
                                tbl.tbl_coeff[0][u] += tbl_prev_p[u];
                            }
                        }
                    });
                layer.evaluate[0] = Arc::new(
                    move |_: &Layer<T>, inputs: &Vec<T>| {
                        let mut res = vec![T::from(0); 1 << log_size];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 12;
                            let mut i = ofs;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let a_idx = (15 + y - 3 * x) * 3 % 5;
                                        let a_idy = x;
                                        let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;
                                        let c0_idx = (a_idx + 4) % 5;
                                        let c0_idz = a_idz;
                                        let c1_idx = (a_idx + 1) % 5;
                                        let c1_idz = (a_idz + 63) % 64;
                                        res[i] += inputs[ofs + (a_idx * 5 + a_idy) * 64 + a_idz];
                                        for c_idy in 0..=4 {
                                            res[i] += inputs[ofs + (c0_idx * 5 + c_idy) * 64 + c0_idz];
                                        }
                                        for c_idy in 0..=4 {
                                            res[i] += inputs[ofs + (c1_idx * 5 + c_idy) * 64 + c1_idz];
                                        }
                                        i += 1;
                                    }
                                }
                            }
                            for i in 0..512 {
                                let u = ofs + 2048 + i;
                                res[u] += inputs[u];
                            }
                        }
                        res
                    });
                layer.predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                        let mut sum = T::from(0);
                        let mut tbl0 = vec![T::from(0); 1 << 12];
                        let mut tbl1 = vec![T::from(0); 1 << 12];
                        // TODO XXX reuse tbl0, tbl1 for next layer
                        tbl0[0] = T::from(1);
                        tbl1[0] = T::from(1);
                        for k in 0..12 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl0[i + mask] = tbl0[i] * tau[k];
                                tbl0[i] = tbl0[i] - tbl0[i + mask];
                                tbl1[i + mask] = tbl1[i] * rho[k];
                                tbl1[i] = tbl1[i] - tbl1[i + mask];
                            }
                        }
                        let mut i = 0;
                        for x in 0..=4 {
                            for y in 0..=4 {
                                for z in 0..64 {
                                    let a_idx = (15 + y - 3 * x) * 3 % 5;
                                    let a_idy = x;
                                    let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;
                                    let c0_idx = (a_idx + 4) % 5;
                                    let c0_idz = a_idz;
                                    let c1_idx = (a_idx + 1) % 5;
                                    let c1_idz = (a_idz + 63) % 64;
                                    sum += tbl0[i] * tbl1[(a_idx * 5 + a_idy) * 64 + a_idz];
                                    for c_idy in 0..=4 {
                                        sum += tbl0[i] * tbl1[(c0_idx * 5 + c_idy) * 64 + c0_idz];
                                    }
                                    for c_idy in 0..=4 {
                                        sum += tbl0[i] * tbl1[(c1_idx * 5 + c_idy) * 64 + c1_idz];
                                    }
                                    i += 1;
                                }
                            }
                        }

                        let mut relay = one(tau[11]) * one(rho[11]) *
                            zer(tau[10]) * zer(rho[10]) *
                            zer(tau[9]) * zer(rho[9]);
                        for i in 0..9 {
                            relay *= eq2(tau[i], rho[i]);
                        }
                        let mut block_eq = T::from(1);
                        for i in 12..tau.len() {
                            block_eq *= eq2(tau[i], rho[i]);
                        }
                        [T::from(0), (relay + sum) * block_eq, T::from(0), T::from(0)]
                    }));
            }
            /*
            for k in 0..(2 << log_leaf_size) {
                let ofs = if d != 0 { k << 12 } else { k << 9 };
                for x in 0..=4 {
                    for y in 0..=4 {
                        for z in 0..64 {
                            // # ρ and π steps
                            // B[y,2*x+3*y] = rot(A[x,y], r[x,y])
                            // invert:
                            // let x' = (y-3x)/2, y' = x
                            // B[x, y] = rot[A[x', y'], r[x', y']
                            // B[x, y, z] = A[x', y', z']
                            let a_idx = (15 + y - 3 * x) * 3 % 5;
                            let a_idy = x;
                            // assert_eq!((a_idx * 2 + a_idy * 3) % 5, y);
                            let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;

                            // # θ step
                            // C[x] = A[x,0] xor A[x,1] xor A[x,2] xor A[x,3] xor A[x,4]
                            // A[x,y] = A[x,y] xor C[x-1] xor rot(C[x+1],1)
                            // A[x, y, z] = A[x, y, z] xor
                            //              A[x - 1, 0, z] xor ... xor A[x - 1, 4, z] xor
                            //              A[x + 1, 0, z - 1) xor ... xor A[x + 1, 4, z - 1]
                            let c0_idx = (a_idx + 4) % 5;
                            let c0_idz = a_idz;
                            let c1_idx = (a_idx + 1) % 5;
                            let c1_idz = (a_idz + 63) % 64;

                            let g = &mut layer.gates[(k << 12) + (x * 5 + y) * 64 + z];
                            g.gtype = GateType::SUM;
                            g.inputs = vec![];
                            // A[x, y, z]
                            if d != 0 || a_idx * 5 + a_idy < 8 {
                                g.inputs.push((ofs + (a_idx * 5 + a_idy) * 64 + a_idz) as u32);
                            }
                            // C[x - 1, z] = sum_y A[x - 1, y, z]
                            for c0_idy in 0..=4 {
                                if d != 0 || c0_idx * 5 + c0_idy < 8 {
                                    g.inputs.push((ofs + (c0_idx * 5 + c0_idy) * 64 + c0_idz) as u32);
                                }
                            }
                            // rot(C[x + 1, 1] = rot(sum_y A[x + 1, y], 1)
                            for c1_idy in 0..=4 {
                                if d != 0 || c1_idx * 5 + c1_idy < 8 {
                                    g.inputs.push((ofs + (c1_idx * 5 + c1_idy) * 64 + c1_idz) as u32);
                                }
                            }
                        }
                    }
                }
                for i in 0..512 {
                    let g = &mut layer.gates[(k << 12) + 2048 + i];
                    if d == 0 && k < (1 << log_leaf_size) {
                        let idx = (ofs + i) as u32;
                        g.gtype = GateType::P1D2C0;
                        g.inputs = vec![idx];
                        g.utils = vec![0];
                    } else if d == 0 {
                        let idx = (ofs + i) as u32;
                        g.gtype = GateType::RELAY;
                        g.inputs = vec![idx];
                    } else {
                        g.gtype = GateType::RELAY;
                        g.inputs = vec![((k << 12) + 2048 + i) as u32];
                    }
                }
            }
            */
            depth += 1;
        }

        // non-linear layers
        {
            c.layers.push(Layer::<T>::default());
            let layer = &mut c.layers[depth];
            let log_size =
                if d != 23 {
                    log_leaf_size + 1 + 11 + 1
                } else {
                    log_leaf_size + 1 + 10
                };
            layer.bit_len = log_size;
            // last round
            if d == 23 {
                layer.sumcheck_init[0] = Arc::new(
                    move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                        tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size + 2)]; 2];
                        let tbl_prev_p = &tbl.tbl_prev_point[0];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 10;
                            for y in 0..4 {
                                for z in 0..64 {
                                    let i = ofs +  y * 64 + z;
                                    let w = (ofs << 2) + y * 64 + z;
                                    let u = (ofs << 2) + (y + 5) * 64 + z;
                                    let v = (ofs << 2) + (y + 10) * 64 + z;
                                    tbl.tbl_coeff[1][u] -= tbl_prev_p[i] * tbl.tbl_value[v];
                                    tbl.tbl_coeff[0][u] += tbl_prev_p[i] * tbl.tbl_value[v];
                                    tbl.tbl_coeff[1][w] += tbl_prev_p[i];
                                    if y == 0 && (RC[d] & (1 << z) != 0) {
                                        tbl.tbl_coeff[0][w] += tbl_prev_p[i];
                                    }
                                }
                            }
                            for z in 0..512 {
                                let i = ofs + 512 + z;
                                let w = (ofs << 2) + 2048 + z;
                                tbl.tbl_coeff[1][w] += tbl_prev_p[i];
                            }
                        }
                    });
                layer.sumcheck_init[1] = Arc::new(
                    move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                        tbl.tbl_coeff = vec![vec![]; 4];
                        tbl.tbl_coeff[0] = vec![T::from(0)];
                        tbl.tbl_coeff[1] = vec![T::from(0); 1 << (log_size + 2)];
                        tbl.tbl_coeff[2] = vec![T::from(0); 1 << (log_size + 2)];
                        tbl.tbl_coeff[3] = vec![T::from(0); 1 << (log_size + 2)];
                        let tbl_prev_p = &tbl.tbl_prev_point[0];
                        let tbl_p = &tbl.tbl_point[0];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 10;
                            for y in 0..4 {
                                for z in 0..64 {
                                    let i = ofs +  y * 64 + z;
                                    let w = (ofs << 2) + y * 64 + z;
                                    let u = (ofs << 2) + (y + 5) * 64 + z;
                                    let v = (ofs << 2) + (y + 10) * 64 + z;
                                    let tmp = tbl_prev_p[i] * tbl_p[u];
                                    tbl.tbl_coeff[3][v] -= tmp;
                                    tbl.tbl_coeff[1][v] += tmp;
                                    let tmp = tbl_prev_p[i] * tbl_p[w];
                                    tbl.tbl_coeff[2][0] += tmp;

                                    if y == 0 && (RC[d] & (1 << z) != 0) {
                                        tbl.tbl_coeff[0][0] += tmp;
                                    }
                                }
                            }
                            for z in 0..512 {
                                let i = ofs + 512 + z;
                                let w = (ofs << 2) + 2048 + z;
                                tbl.tbl_coeff[2][0] += tbl_prev_p[i] * tbl_p[w];
                            }
                        }
                    });
                layer.evaluate[0] = Arc::new(
                    move |_: &Layer<T>, inputs: &Vec<T>| {
                        let mut res = vec![T::from(0); 1 << log_size];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 10;
                            for y in 0..4 {
                                for z in 0..64 {
                                    let idx = ofs +  y * 64 + z;
                                    let idx0 = (ofs << 2) + y * 64 + z;
                                    let idx1 = (ofs << 2) + (y + 5) * 64 + z;
                                    let idx2 = (ofs << 2) + (y + 10) * 64 + z;
                                    res[idx] = (T::from(1) - inputs[idx1]) * inputs[idx2] + inputs[idx0];

                                    if y == 0 && (RC[d] & (1 << z) != 0) {
                                        res[idx] += T::from(1);
                                    }
                                }
                            }
                            for i in 0..512 {
                                let idx0 = ofs + 512 + i;
                                let idx1 = (ofs << 2) + 2048 + i;
                                res[idx0] = inputs[idx1];
                            }
                        }
                        res
                    });
                layer.predicate.push(Arc::new(
                        move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>| {
                        let mut tbl0 = vec![T::from(0); 1 << 10];
                        let mut tbl1 = vec![T::from(0); 1 << 12];
                        let mut tbl2 = vec![T::from(0); 1 << 12];
                        // TODO XXX reuse tbl0, tbl1 for next layer
                        tbl0[0] = T::from(1);
                        tbl1[0] = T::from(1);
                        tbl2[0] = T::from(1);
                        for k in 0..10 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl0[i + mask] = tbl0[i] * tau[k];
                                tbl0[i] = tbl0[i] - tbl0[i + mask];
                            }
                        }
                        for k in 0..12 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl1[i + mask] = tbl1[i] * rho[k];
                                tbl1[i] = tbl1[i] - tbl1[i + mask];
                                tbl2[i + mask] = tbl2[i] * sgm[k];
                                tbl2[i] = tbl2[i] - tbl2[i + mask];
                            }
                        }
                        let mut chi_mul = T::from(0);
                        let mut chi_add = T::from(0);
                        let mut chi_cst = T::from(0);
                        let mut i = 0;
                        {
                            for y in 0..4 {
                                for z in 0..64 {
                                    if y == 0 && (RC[d] & (1 << z) != 0) {
                                        chi_cst += tbl0[i] * tbl1[i];
                                    }
                                    chi_mul += tbl0[i] * tbl1[(y + 5) * 64 + z] * tbl2[(y + 10) * 64 + z];
                                    chi_add += tbl0[i] * tbl1[i];
                                    i += 1;
                                }
                            }
                        }

                        let mut sgm_zer = tbl2[0];
                        for i in 12..sgm.len() {
                            sgm_zer *= zer(sgm[i]);
                        }

                        let mut relay = one(tau[9]) * one(rho[11]) *
                            zer(rho[10]) * zer(rho[9]);
                        for i in 0..9 {
                            relay *= eq2(tau[i], rho[i]);
                        }
                        let mut block_eq = T::from(1);
                        for i in 10..tau.len() {
                            block_eq *= eq2(tau[i], rho[i + 2]);
                        }
                        let mut block_eq3 = T::from(1);
                        for i in 10..tau.len() {
                            block_eq3 *= eq3(tau[i], rho[i + 2], sgm[i + 2]);
                        }
                        chi_mul *= block_eq3;
                        [chi_mul, (relay + chi_add) * block_eq * sgm_zer, chi_mul, chi_cst * block_eq * sgm_zer]
                    }));
            } else {
                layer.sumcheck_init[0] = Arc::new(
                    move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                        tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_size]; 2];
                        let tbl_prev_p = &tbl.tbl_prev_point[0];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 12;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let w = ofs + (x * 5 + y) * 64 + z;
                                        let u = ofs + (((x + 1) % 5) * 5 + y) * 64 + z;
                                        let v = ofs + (((x + 2) % 5) * 5 + y) * 64 + z;
                                        tbl.tbl_coeff[1][u] -= tbl_prev_p[w] * tbl.tbl_value[v];
                                        tbl.tbl_coeff[0][u] += tbl_prev_p[w] * tbl.tbl_value[v];
                                        tbl.tbl_coeff[1][w] += tbl_prev_p[w];
                                        if x == 0 && y == 0 && (RC[d] & (1 << z) != 0) {
                                            tbl.tbl_coeff[0][w] += tbl_prev_p[w];
                                        }
                                    }
                                }
                            }
                            for i in 0..512 {
                                let w = ofs + 2048 + i;
                                tbl.tbl_coeff[1][w] += tbl_prev_p[w];
                            }
                        }
                    });
                layer.sumcheck_init[1] = Arc::new(
                    move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                        tbl.tbl_coeff = vec![vec![]; 4];
                        tbl.tbl_coeff[0] = vec![T::from(0)];
                        tbl.tbl_coeff[1] = vec![T::from(0); 1 << log_size];
                        tbl.tbl_coeff[2] = vec![T::from(0); 1 << log_size];
                        tbl.tbl_coeff[3] = vec![T::from(0); 1 << log_size];
                        let tbl_prev_p = &tbl.tbl_prev_point[0];
                        let tbl_p = &tbl.tbl_point[0];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 12;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let w = ofs + (x * 5 + y) * 64 + z;
                                        let u = ofs + (((x + 1) % 5) * 5 + y) * 64 + z;
                                        let v = ofs + (((x + 2) % 5) * 5 + y) * 64 + z;
                                        let tmp = tbl_prev_p[w] * tbl_p[u];
                                        tbl.tbl_coeff[3][v] -= tmp;
                                        tbl.tbl_coeff[1][v] += tmp;
                                        let tmp = tbl_prev_p[w] * tbl_p[w];
                                        tbl.tbl_coeff[2][0] += tmp;

                                        if x == 0 && y == 0 && (RC[d] & (1 << z) != 0) {
                                            tbl.tbl_coeff[0][0] += tmp;
                                        }
                                    }
                                }
                            }
                            for i in 0..512 {
                                let w = ofs + 2048 + i;
                                tbl.tbl_coeff[2][0] += tbl_prev_p[w] * tbl_p[w];
                            }
                        }
                    });
                layer.evaluate[0] = Arc::new(
                    move |_: &Layer<T>, inputs: &Vec<T>| {
                        let mut res = vec![T::from(0); 1 << log_size];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 12;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let idx0 = ofs + (x * 5 + y) * 64 + z;
                                        let idx1 = ofs + (((x + 1) % 5) * 5 + y) * 64 + z;
                                        let idx2 = ofs + (((x + 2) % 5) * 5 + y) * 64 + z;
                                        res[idx0] = (T::from(1) - inputs[idx1]) * inputs[idx2] + inputs[idx0];

                                        if x == 0 && y == 0 && (RC[d] & (1 << z) != 0) {
                                            res[idx0] += T::from(1);
                                        }
                                    }
                                }
                            }
                            for i in 0..512 {
                                let idx = ofs + 2048 + i;
                                res[idx] = inputs[idx];
                            }
                        }
                        res
                    });
                layer.predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>| {
                        let mut tbl0 = vec![T::from(0); 1 << 12];
                        let mut tbl1 = vec![T::from(0); 1 << 12];
                        let mut tbl2 = vec![T::from(0); 1 << 12];
                        // TODO XXX reuse tbl0, tbl1 for next layer
                        tbl0[0] = T::from(1);
                        tbl1[0] = T::from(1);
                        tbl2[0] = T::from(1);
                        for k in 0..12 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl0[i + mask] = tbl0[i] * tau[k];
                                tbl0[i] = tbl0[i] - tbl0[i + mask];
                                tbl1[i + mask] = tbl1[i] * rho[k];
                                tbl1[i] = tbl1[i] - tbl1[i + mask];
                                tbl2[i + mask] = tbl2[i] * sgm[k];
                                tbl2[i] = tbl2[i] - tbl2[i + mask];
                            }
                        }
                        let mut chi_mul = T::from(0);
                        let mut chi_add = T::from(0);
                        let mut chi_cst = T::from(0);
                        let mut i = 0;
                        for x in 0..=4 {
                            for y in 0..=4 {
                                for z in 0..64 {
                                    if x == 0 && y == 0 && (RC[d] & (1 << z) != 0) {
                                        chi_cst += tbl0[i] * tbl1[i];
                                    }

                                    chi_mul += tbl0[i] * tbl1[(((x + 1) % 5) * 5 + y) * 64 + z] * tbl2[(((x + 2) % 5) * 5 + y) * 64 + z];
                                    chi_add += tbl0[i] * tbl1[i];
                                    i += 1;
                                }
                            }
                        }

                        let mut sgm_zer = tbl2[0];
                        for i in 12..sgm.len() {
                            sgm_zer *= zer(sgm[i]);
                        }

                        let mut relay = one(tau[11]) * one(rho[11]) *
                            zer(tau[10]) * zer(rho[10]) *
                            zer(tau[9]) * zer(rho[9]);
                        for i in 0..9 {
                            relay *= eq2(tau[i], rho[i]);
                        }
                        let mut block_eq = T::from(1);
                        for i in 12..tau.len() {
                            block_eq *= eq2(tau[i], rho[i]);
                        }
                        let mut block_eq3 = T::from(1);
                        for i in 12..tau.len() {
                            block_eq3 *= eq3(tau[i], rho[i], sgm[i]);
                        }
                        chi_mul *= block_eq3;
                        [chi_mul, (relay + chi_add) * block_eq * sgm_zer, chi_mul, chi_cst * block_eq * sgm_zer]
                    }));
            }
            /*
            for k in 0..(2 << log_leaf_size) {
                let ofs = k << 12;
                if d != 23 {
                    for x in 0..=4 {
                        for y in 0..=4 {
                            for z in 0..64 {
                                let g = &mut layer.gates[ofs + (x * 5 + y) * 64 + z];

                                // # ι step
                                // A[0,0] = A[0,0] xor RC
                                if x == 0 && y == 0 && (RC[d] & (1 << z) != 0) {
                                    g.gtype = GateType::NKECCAKX;
                                } else {
                                    g.gtype = GateType::KECCAKX;
                                }

                                // # χ step
                                // A[x,y] = B[x,y] xor ((not B[x+1,y]) and B[x+2,y])
                                // A[x, y, z] = B[x, y, z] + (1 - B[x + 1, y, z]) * B[x + 2, y, z]
                                g.inputs = vec![
                                    (ofs + (((x + 1) % 5) * 5 + y) * 64 + z) as u32,
                                    (ofs + (((x + 2) % 5) * 5 + y) * 64 + z) as u32,
                                    (ofs + (x * 5 + y) * 64 + z) as u32,
                                ];
                            }
                        }
                    }
                    for i in 0..512 {
                        let g = &mut layer.gates[ofs + 2048 + i];
                        g.gtype = GateType::RELAY;
                        g.inputs = vec![(ofs + 2048 + i) as u32];
                    }
                // last round
                } else {
                    for y in 0..4 {
                        for z in 0..64 {
                            let g = &mut layer.gates[(k << 10) + y * 64 + z];

                            // # ι step
                            // A[0,0] = A[0,0] xor RC
                            if y == 0 && (RC[d] & (1 << z) != 0) {
                                g.gtype = GateType::NKECCAKX;
                            } else {
                                g.gtype = GateType::KECCAKX;
                            }

                            // # χ step
                            // A[x,y] = B[x,y] xor ((not B[x+1,y]) and B[x+2,y])
                            // A[x, y, z] = B[x, y, z] + (1 - B[x + 1, y, z]) * B[x + 2, y, z]
                            g.inputs = vec![
                                (ofs + (y + 5) * 64 + z) as u32,
                                (ofs + (y + 10) * 64 + z) as u32,
                                (ofs + y * 64 + z) as u32,
                            ];
                        }
                    }

                    for i in 0..512 {
                        let idx = (ofs + 2048 + i) as u32;
                        let g = &mut layer.gates[(k << 10) + 512 + i];
                        g.gtype = GateType::RELAY;
                        g.inputs = vec![idx];
                    }
                }
            }
            */
            depth += 1;
        }
    }

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_leaf_size + 1 + 9;
        layer.bit_len = log_size;
        layer.val_degree = vec![1];
        layer.num_table = vec![1];
        layer.num_phase = 1;
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                let mut relay = one(rho[9]);
                let mut last_block = relay;
                for i in 9..(log_leaf_size + 9 + 1) {
                    relay *= eq2(tau[i], rho[i + 1]);
                    last_block *= one(tau[i]) * one(rho[i + 1]);
                }
                relay = relay * eq2(tau[8], rho[8]) - last_block * one(tau[8]) * one(rho[8]);

                let mut eq_low8 = T::from(1);
                for i in 0..8 {
                    eq_low8 *= eq2(tau[i], rho[i]);
                }

                let mut sum = one(tau[log_leaf_size + 9]) * zer(rho[9]) * zer(rho[8]);
                let mut last_block = sum;
                for i in 8..(log_leaf_size + 9) {
                    sum *= eq2(tau[i], rho[i + 2]);
                    last_block *= one(tau[i]) * one(rho[i + 2]);
                }
                sum -= last_block;

                [T::from(0), (relay + sum) * eq_low8, T::from(0), T::from(0)]
            }));
        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size + 1)]];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in 0..(2 << log_leaf_size) {
                    for i in 0..512 {
                        if k == ((2 << log_leaf_size) - 1) && i == 256 {
                            break;
                        }
                        let idx0 = (k << 10) + 512 + i;
                        let idx1 = (k << 9) + i;
                        tbl.tbl_coeff[0][idx0] += tbl_prev_p[idx1];
                    }
                }
                let ofs = 1 << (log_leaf_size + 9);
                for k in 0..((2 << log_leaf_size) - 1) {
                    for i in 0..256 {
                        let idx0 = (k << 10) + i;
                        let idx1 = ofs + (k << 8) + i;
                        tbl.tbl_coeff[0][idx0] += tbl_prev_p[idx1];
                    }
                }
            });
        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_size];
                for k in 0..(2 << log_leaf_size) {
                    for i in 0..512 {
                        if k == ((2 << log_leaf_size) - 1) && i == 256 {
                            break;
                        }
                        let idx0 = (k << 10) + 512 + i;
                        let idx1 = (k << 9) + i;
                        res[idx1] += inputs[idx0];
                    }
                }
                let ofs = 1 << (log_leaf_size + 9);
                for k in 0..((2 << log_leaf_size) - 1) {
                    for i in 0..256 {
                        let idx0 = (k << 10) + i;
                        let idx1 = ofs + (k << 8) + i;
                        res[idx1] += inputs[idx0];
                    }
                }
                res
            });

        /*
        // mt leaf bittest
        for k in 0..(1 << log_leaf_size) {
            for i in 0..512 {
                let idx = ((k << 10) + 512 + i) as u32;
                let g = &mut layer.gates[(k << 9) + i];
                g.gtype = GateType::RELAY;
                g.inputs = vec![idx];
            }
        }

        // mt inner hash witness
        let mut prev_ofs0 = 0;
        let mut prev_ofs1 = (1 << log_leaf_size) << 10;
        let mut ofs = (1 << log_leaf_size) << 9;
        for d in 0..log_leaf_size {
            for _ in 0..(1 << (log_leaf_size - d - 1)) {
                for i in 0..256 {
                    let idx0 = (prev_ofs0 + i) as u32;
                    let idx1 = (prev_ofs1 + 512 + i) as u32;
                    let g = &mut layer.gates[ofs + i];
                    g.gtype = GateType::SUM;
                    g.inputs = vec![idx0, idx1];
                }
                prev_ofs0 += 1024;
                ofs += 256;
                for i in 0..256 {
                    let idx0 = (prev_ofs0 + i) as u32;
                    let idx1 = (prev_ofs1 + 512 + 256 + i) as u32;
                    let g = &mut layer.gates[ofs + i];
                    g.gtype = GateType::SUM;
                    g.inputs = vec![idx0, idx1];
                }
                prev_ofs0 += 1024;
                prev_ofs1 += 1024;
                ofs += 256;
            }
        }

        // mt roothash
        {
            for i in 0..256 {
                let idx0 = (prev_ofs0 + i) as u32;
                let idx1 = (prev_ofs1 + 512 + i) as u32;
                let g = &mut layer.gates[ofs + i];
                g.gtype = GateType::SUM;
                g.inputs = vec![idx0, idx1];
            }
        }
        */

        depth += 1;
    }

    c.depth = depth;
    c
}

pub fn rot(x: u64, r: usize) -> u64 {
    (x << r) | (x >> (64 - r))
}

pub fn gen_sha3mt_witness<T: Field>(inputs: &Vec<[u64; 8]>, log_leaf_size: usize) -> Vec<T> {
    let mut witness = Vec::with_capacity(512 << (log_leaf_size + 1));
    for i in 0..inputs.len() {
        for w in 0..8 {
            let mut inp = inputs[i][w] as usize;
            for _ in 0..64 {
                witness.push(T::from(inp & 1));
                inp >>= 1;
            }
        }
    }
    let mut inputs = inputs.clone();

    for d in 0..=log_leaf_size {
        let width = 1 << (log_leaf_size - d);
        let mut state = vec![[0u64; 25]; width];
        for i in 0..width {
            for j in 0..8 {
                state[i][j] = inputs[i][j];
            }

            for r in 0..24 {
                // θ step
                let mut c = [0u64; 5];
                for x in 0..=4 {
                    c[x] =
                        state[i][x * 5 + 0] ^
                        state[i][x * 5 + 1] ^
                        state[i][x * 5 + 2] ^
                        state[i][x * 5 + 3] ^
                        state[i][x * 5 + 4];
                }
                let mut a = [0u64; 25];
                for x in 0..=4 {
                    for y in 0..=4 {
                        a[x * 5 + y] = state[i][x * 5 + y] ^ c[(x + 4) % 5] ^ rot(c[(x + 1) % 5], 1);
                    }
                }
                // if d == 0 && r == 0 {
                //     for i in 0..25 {
                //         println!("a {:02}: {:016x}", i, a[i]);
                //     }
                // }

                // ρ and π steps
                let mut b = [0u64; 25];
                for x in 0..=4 {
                    for y in 0..=4 {
                        b[y * 5 + (2 * x + 3 * y) % 5] = rot(a[x * 5 + y], ROTATE[x][y]);
                    }
                }
                // if d == 0 && r == 0 {
                //     for i in 0..25 {
                //         println!("b {:02}: {:016x}", i, b[i]);
                //     }
                // }

                // χ step
                for x in 0..=4 {
                    for y in 0..=4 {
                        state[i][x * 5 + y] =
                            b[x * 5 + y] ^
                            ((!b[((x + 1) % 5) * 5 + y]) & b[((x + 2) % 5) * 5 + y]);
                    }
                }

                // ι step
                state[i][0] ^= RC[r] as u64;

                // if d == 0 && r == 22 {
                //     for j in 0..25 {
                //         println!("b {:02}: {:016x}", j, state[i][j]);
                //     }
                // }
            }
            inputs[i / 2][(i % 2) * 4 + 0] = state[i][0];
            inputs[i / 2][(i % 2) * 4 + 1] = state[i][1];
            inputs[i / 2][(i % 2) * 4 + 2] = state[i][2];
            inputs[i / 2][(i % 2) * 4 + 3] = state[i][3];
            for k in 0..4 {
                let mut inp = state[i][k] as usize;
                for _ in 0..64 {
                    witness.push(T::from(inp & 1));
                    inp >>= 1;
                }
            }
        }
    }
    for _ in 0..256 {
        witness.push(T::from(0));
    }

    witness
}
