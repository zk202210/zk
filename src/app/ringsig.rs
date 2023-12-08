use std::sync::Arc;

use crate::field::*;
use crate::circuit::*;
use crate::gkr::*;
use crate::util::*;

use super::sha3mt::ROTATE;
use super::sha3mt::RC;
use super::sha3mt::rot;

// sha3-256 merkle tree proof without padding
// leaf inputs are 128-bit
// inner node inputs are 256-bit
// path left-right selectors are 1-bit for each layer
pub fn build_circuit<'a, T: Field>(path_size: usize) -> Circuit<'a, T> {
    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);

    let log_num_hash = ceil_log(path_size);

    let log_leaf_size = log_num_hash - 1;

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

    // first block 512 bits: 128 bit ring-key, 1 random, 63 bit selectors, 64 bit zero, 256 bit root-hash;
    // rest blocks 512 bits: 256 bit prev-hash, 256 bit proof hash;
    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_num_hash + 9;
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

    // first block: 128 bit ring-key, 128*3 bits zero, 256 bit root-hash, 256 bit zeros,
    // rest blocks: 256 selected left hash, 256 selected right hash, 256 bit prev-hash, 256 bit proof hash
    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_num_hash + 9 + 1;
        layer.bit_len = log_size;
        /*
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];
        for k in 0..(1 << log_num_hash) {
            if k == 0 {
                for i in 0..128 {
                    // relay of ring-key
                    let g = &mut layer.gates[i];
                    g.gtype = GateType::RELAY;
                    g.inputs = vec![i as u32];
                }
                for i in 512..(512 + 256) {
                    // relay of root-hash
                    let g = &mut layer.gates[i];
                    g.gtype = GateType::RELAY;
                    g.inputs = vec![(i - 256) as u32];
                }
            } else {
                let bk = k as u32;
                let sel_idx = bk + 128;
                let bk = bk << 9;
                for i in 0..256 {
                    // selected of left child
                    let g = &mut layer.gates[(k << 10) + i];
                    g.gtype = GateType::SELECTOR;
                    let i = i as u32;
                    g.inputs = vec![sel_idx, bk + i, bk + i + 256];
                }
                for i in 256..512 {
                    // selected of right child
                    let g = &mut layer.gates[(k << 10) + i];
                    g.gtype = GateType::SELECTOR;
                    let i = (i - 256) as u32;
                    g.inputs = vec![sel_idx, bk + i + 256, bk + i];
                }
                for i in 512..1024 {
                    // relay of two hashes
                    let g = &mut layer.gates[(k << 10) + i];
                    g.gtype = GateType::RELAY;
                    let i = i as u32;
                    g.inputs = vec![bk + i - 512];
                }
            }
        }
        */
        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size - 1)]; 2];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in 0..(1 << log_num_hash) {
                    if k == 0 {
                        for i in 0..128 {
                            // relay of ring-key
                            tbl.tbl_coeff[1][i] += tbl_prev_p[i];
                        }
                        for i in 512..(512 + 256) {
                            // relay of root-hash
                            tbl.tbl_coeff[1][i - 256] += tbl_prev_p[i];
                        }
                    } else {
                        let sel_idx = k + 128;
                        let k0 = k << 9;
                        let k = k << 10;
                        for i in 0..256 {
                            // selected of left child
                            let v = k0 + i;
                            let w = v + 256;
                            tbl.tbl_coeff[1][sel_idx] -= (tbl_prev_p[k + i] + tbl_prev_p[k + 256 + i]) *
                                (tbl.tbl_value[v] - tbl.tbl_value[w]);
                            tbl.tbl_coeff[0][sel_idx] += tbl_prev_p[k + i] * tbl.tbl_value[v];
                            tbl.tbl_coeff[0][sel_idx] += tbl_prev_p[k + 256 + i] * tbl.tbl_value[w];
                        }
                        for i in 512..1024 {
                            tbl.tbl_coeff[1][k0 + i - 512]
                                += tbl_prev_p[k + i];
                        }
                    }
                }
            });

        layer.sumcheck_init[1] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size - 1)]; 4];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                let tbl_p = &tbl.tbl_point[0];
                for k in 0..(1 << log_num_hash) {
                    if k == 0 {
                        for i in 0..128 {
                            // relay of ring-key
                            tbl.tbl_coeff[2][0] += tbl_prev_p[i] * tbl_p[i];
                        }
                        for i in 512..(512 + 256) {
                            // relay of root-hash
                            tbl.tbl_coeff[2][0] += tbl_prev_p[i] * tbl_p[i - 256];
                        }
                    } else {
                        let sel_idx = k + 128;
                        let k0 = k << 9;
                        let k = k << 10;
                        for i in 0..256 {
                            // selected of left child
                            let v = k0 + i;
                            let w = v + 256;
                            let tmp = tbl_prev_p[k + i] * tbl_p[sel_idx];
                            tbl.tbl_coeff[3][v] -= tmp;
                            tbl.tbl_coeff[1][v] += tmp;
                            tbl.tbl_coeff[3][w] += tmp;

                            let tmp = tbl_prev_p[k + 256 + i] * tbl_p[sel_idx];
                            tbl.tbl_coeff[3][w] -= tmp;
                            tbl.tbl_coeff[1][w] += tmp;
                            tbl.tbl_coeff[3][v] += tmp;
                        }
                        for i in 512..1024 {
                            tbl.tbl_coeff[2][0] += tbl_prev_p[k + i] * tbl_p[k0 + i - 512];
                        }
                    }
                }
            });



        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_size];
                for k in 0..(1 << log_num_hash) {
                    if k == 0 {
                        for i in 0..128 {
                            // relay of ring-key
                            res[i] = inputs[i];
                        }
                        for i in 512..(512 + 256) {
                            res[i] = inputs[i - 256];
                        }
                    } else {
                        let sel_idx = k + 128;
                        let k0 = k << 9;
                        let k = k << 10;
                        for i in 0..256 {
                            // selected of left child
                            if inputs[sel_idx] == T::from(0) {
                                res[k + i] = inputs[k0 + i];
                                res[k + i + 256] = inputs[k0 + i + 256];
                            } else {
                                res[k + i] = inputs[k0 + i + 256];
                                res[k + i + 256] = inputs[k0 + i];
                            }
                        }
                        for i in 512..1024 {
                            res[k + i] = inputs[k0 + i - 512];
                        }
                    }
                }
                res
            });

        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>| {
                // selector
                // (1-u)v + uw && (1-u)w + uv
                let mut block_eq = T::from(1);
                let mut block_0 = T::from(1);
                for i in 0..log_num_hash {
                    block_eq *= eq3(tau[10 + i], rho[i], sgm[9 + i]);
                    block_0 *= zer(tau[10 + i]) * zer(rho[i]) * zer(sgm[9 + i]);
                }
                block_eq += block_0;
                for i in 0..log_num_hash {
                    block_eq *= zer(rho[i + 9]);
                }
                for i in log_num_hash..7 {
                    block_eq *= zer(rho[i]);
                }
                block_eq *= one(rho[7]) * zer(rho[8]);

                let mut idx_eq = T::from(1);
                for i in 0..8 {
                    idx_eq *= eq2(tau[i], sgm[i]);
                }

                let selector = block_eq * idx_eq * zer(tau[9]);

                // relay
                let mut zer_sgm = T::from(1);
                for i in 0..sgm.len() {
                    zer_sgm *= zer(sgm[i]);
                }
                let mut block_eq = T::from(1);
                let mut block_0 = T::from(1);
                for i in 0..log_num_hash {
                    block_eq *= eq2(tau[10 + i], rho[9 + i]);
                    block_0 *= zer(tau[10 + i]) * zer(rho[9 + i]);
                }
                let mut idx_eq = T::from(1);
                for i in 0..7 {
                    idx_eq *= eq2(tau[i], rho[i]);
                }
                block_eq += block_0;
                let relay = zer_sgm * idx_eq * (
                    eq2(tau[7], rho[7]) * one(tau[9]) * (
                        block_eq * eq2(tau[8], rho[8]) + block_0 * zer(tau[8]) * one(rho[8])) +
                    block_0 * zer(tau[7]) * zer(rho[7]) * zer(tau[8]) * zer(rho[8]) * zer(tau[9]));
                // let relay = (block_eq * eq2(tau[7], rho[7]) * eq2(tau[8], rho[8]) * one(tau[9]) +
                //              block_0 * eq2(tau[7], rho[7]) * zer(tau[8]) * one(rho[8]) * one(tau[9])) +
                //              block_0 * zer(tau[7]) * zer(rho[7]) * zer(tau[8]) * zer(rho[8]) * zer(tau[9]) +
                //     * idx_eq * zer_sgm ;

                [selector, relay, selector * eq2(tau[8], sgm[8]), T::from(0)]

            }));
        depth += 1;
    }

    for d in 0..24 {
        // linear layers
        {
            c.layers.push(Layer::<T>::default());
            let layer = &mut c.layers[depth];
            // 1600 states => 2048
            let log_size = log_num_hash + 11 + 1;
            layer.bit_len = log_size;
            // layer.gates = vec![Gate::default(); 1 << log_size];
            // layer.values = vec![T::from(0); 1 << log_size];
            // first round
            if d == 0 {
                layer.utils = vec![T::from(1)];
                layer.num_phase = 1;
                layer.num_table = vec![2];
                layer.val_degree = vec![2];
                layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_d2_no_const);
                layer.sumcheck_init[0] = Arc::new(
                    move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                        tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size - 2)]; 2];
                        let tbl_prev_p = &tbl.tbl_prev_point[0];
                        for k in 0..(1 << log_num_hash) {
                            let ofs = k << 10;
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
                            for z in 0..256 {
                                let u = (k << 10) + 512 + z;
                                let i = (k << 12) + 2048 + z;
                                tbl.tbl_coeff[0][u] += tbl_prev_p[i];

                                let u = u + 256;
                                let i = i + 256;
                                tbl.tbl_coeff[1][u] += tbl_prev_p[i];
                                tbl.tbl_coeff[0][u] -= tbl_prev_p[i];
                            }
                        }
                    });
                layer.evaluate[0] = Arc::new(
                    move |_: &Layer<T>, inputs: &Vec<T>| {
                        let mut res = vec![T::from(0); 1 << log_size];
                        for k in 0..(2 << log_leaf_size) {
                            let ofs = k << 10;
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
                            for z in 0..256 {
                                let u = (k << 10) + 512 + z;
                                let i = (k << 12) + 2048 + z;
                                res[i] = inputs[u];

                                let u = u + 256;
                                let i = i + 256;
                                res[i] = inputs[u] * (inputs[u] - T::from(1));
                            }
                        }
                        res
                    });
                // TODO XXX, reuse tbl0 & tbl1
                layer.predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                        let mut sum = T::from(0);
                        let mut tbl0 = vec![T::from(0); 1 << 12];
                        let mut tbl1 = vec![T::from(0); 1 << 10];
                        tbl0[0] = T::from(1);
                        tbl1[0] = T::from(1);
                        for k in 0..12 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl0[i + mask] = tbl0[i] * tau[k];
                                tbl0[i] = tbl0[i] - tbl0[i + mask];
                            }
                        }
                        for k in 0..10 {
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

                        let mut relay = one(tau[11]) * zer(tau[10]) * zer(tau[9]) * one(rho[9]);
                        for i in 0..8 {
                            relay *= eq2(tau[i], rho[i]);
                        }
                        let mut block_eq = T::from(1);
                        for i in 12..tau.len() {
                            block_eq *= eq2(tau[i], rho[i - 2]);
                        }
                        relay *= block_eq;
                        [relay * one(tau[8]) * one(rho[8]),
                         relay * eq2(tau[8], rho[8]) + sum * block_eq,
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
            // layer.gates = vec![Gate::default(); 1 << log_size];
            // layer.values = vec![T::from(0); 1 << log_size];
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
        // TODO XXX
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_num_hash + 9;
        layer.bit_len = log_size;
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];
        // layer.val_degree = vec![1];
        // layer.num_table = vec![1];
        // layer.num_phase = 1;

        for k in 0..path_size {
            let curr_b = (k as u32) << 10;
            let next_b = if k != path_size - 1 { curr_b + 1024 + 512 } else { 512 };
            for i in 0..256 {
                let g = &mut layer.gates[(k << 9) + i];
                g.gtype = GateType::SUM;
                let i = i as u32;
                g.inputs = vec![curr_b + i, next_b + i];
            }
            for i in 256..512 {
                let g = &mut layer.gates[(k << 9) + i];
                g.gtype = GateType::RELAY;
                let i = i as u32;
                g.inputs = vec![curr_b + 512 + i];
            }
        }

        // layer.predicate.push(Arc::new(
        //     move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
        //         // relay
        //         let mut relay = T::from(1);

        //         [T::from(0), (relay + sum) * eq_low8, T::from(0), T::from(0)]
        //     }));

        // layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
        // layer.sumcheck_init[0] = Arc::new(
        //     move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
        //         tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size + 1)]];
        //         let tbl_prev_p = &tbl.tbl_prev_point[0];

        //         for k in 0..path_size {
        //             let curr_b = k << 10;
        //             let next_b = if k != path_size - 1 { curr_b + 1024 + 512 } else { 512 };
        //             let k = k << 9;
        //             for i in 0..256 {
        //                 tbl.tbl_coeff[0][curr_b + i] += tbl_prev_p[k + i];
        //                 tbl.tbl_coeff[0][next_b + i] -= tbl_prev_p[k + i];
        //             }
        //             for i in 256..512 {
        //                 tbl.tbl_coeff[0][curr_b + 512 + i] += tbl_prev_p[k + i];
        //             }
        //         }
        //     });
        // layer.evaluate[0] = Arc::new(
        //     move |_: &Layer<T>, inputs: &Vec<T>| {
        //         let mut res = vec![T::from(0); 1 << log_size];

        //         for k in 0..path_size {
        //             let curr_b = k << 10;
        //             let next_b = if k != path_size - 1 { curr_b + 1024 + 512 } else { 512 };
        //             let k = k << 9;
        //             for i in 0..256 {
        //                 res[k + i] = inputs[curr_b + i] - inputs[next_b + i];
        //             }
        //             for i in 256..512 {
        //                 res[k + i] = inputs[curr_b + 512 + i];
        //             }
        //         }
        //         res
        //     });

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

pub fn gen_ringsig_witness<T: Field>(inputs: &Vec<[u64; 4]>, log_num_hash: usize) -> Vec<T> {
    let mut witness = vec![T::from(0); 512 << log_num_hash];

    {
        for w in 0..2 {
            let mut inp = inputs[0][w] as usize;
            for i in 0..64 {
                witness[w * 64 + i] = T::from(inp & 1);
                inp >>= 1;
            }
        }
        witness[128] = T::random();
        {
            let w = 2;
            let mut inp = inputs[0][3] as usize;
            for i in 0..(inputs.len() - 1) {
                witness[w * 64 + i + 1] = T::from(inp & 1);
                inp >>= 1;
            }
        }
    }

    for b in 1..inputs.len() {
        for w in 0..4 {
            let mut inp = inputs[b][w] as usize;
            for i in 0..64 {
                witness[b * 512 + 256 + w * 64 + i] = T::from(inp & 1);
                inp >>= 1;
            }
        }
    }

    let mut inputs = inputs.clone();

    for i in 0..inputs.len() {
        let mut state = vec![0u64; 25];
        if i == 0 {
            state[0] = inputs[0][0];
            state[1] = inputs[0][1];
        } else {
            let selector = if witness[128 + i] == T::from(0) { 0 } else { 1 };
            for j in 0..2 {
                for k in 0..4 {
                    state[j * 4 + k] = inputs[i - 1 + (j ^ selector)][k];
                }
            }
        }

        for r in 0..24 {
            // θ step
            let mut c = [0u64; 5];
            for x in 0..=4 {
                c[x] =
                    state[x * 5 + 0] ^
                    state[x * 5 + 1] ^
                    state[x * 5 + 2] ^
                    state[x * 5 + 3] ^
                    state[x * 5 + 4];
            }
            let mut a = [0u64; 25];
            for x in 0..=4 {
                for y in 0..=4 {
                    a[x * 5 + y] = state[x * 5 + y] ^ c[(x + 4) % 5] ^ rot(c[(x + 1) % 5], 1);
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
                    state[x * 5 + y] =
                        b[x * 5 + y] ^
                        ((!b[((x + 1) % 5) * 5 + y]) & b[((x + 2) % 5) * 5 + y]);
                }
            }

            // ι step
            state[0] ^= RC[r] as u64;

            // if d == 0 && r == 22 {
            //     for j in 0..25 {
            //         println!("b {:02}: {:016x}", j, state[j]);
            //     }
            // }
        }

        inputs[i][0] = state[0];
        inputs[i][1] = state[1];
        inputs[i][2] = state[2];
        inputs[i][3] = state[3];
        if i != inputs.len() - 1 {
            for w in 0..4 {
                let mut inp = state[w] as usize;
                for j in 0..64 {
                    witness[(i + 1) * 512 + w * 64 + j] = T::from(inp & 1);
                    inp >>= 1;
                }
            }
        } else {
            for w in 0..4 {
                let mut inp = state[w] as usize;
                for j in 0..64 {
                    witness[256 + w * 64 + j] = T::from(inp & 1);
                    inp >>= 1;
                }
            }
        }
    }

    witness
}
