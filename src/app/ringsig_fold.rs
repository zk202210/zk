use std::sync::Arc;

use crate::field::*;
use crate::circuit::*;
use crate::gkr::*;
use crate::util::*;

use super::sha3mt::ROTATE;
use super::sha3mt::RC;
use super::sha3mt::rot;


pub fn build_circuit<'a, T: Field>(path_size: usize, log_fold: usize) -> Circuit<'a, T> {
    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);
    let log_num_hash = ceil_log(path_size);

    let mut path_size_bits: Vec<usize> = vec![];
    let mut tmp = path_size;
    let mut idx = 0;
    while tmp != 0 {
        if (tmp & 1) == 1 {
            path_size_bits.push(idx);
        }
        tmp >>= 1;
        idx += 1;
    }

    let mut path_size_minus_bits: Vec<usize> = vec![];
    let mut tmp = path_size - 1;
    let mut idx = 0;
    while tmp != 0 {
        if (tmp & 1) == 1 {
            path_size_minus_bits.push(idx);
        }
        tmp >>= 1;
        idx += 1;
    }

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_num_hash + 11 + log_fold;
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

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_num_hash + 11 + log_fold;
        layer.bit_len = log_size;
        // layer.gates = vec![Gate::default(); 1 << log_size];
        // layer.values = vec![T::from(0); 1 << log_size];

        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_size]; 2];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in 0..path_size {
                    if k == 0 {
                        for i in 0..128 {
                            tbl.tbl_coeff[1][i] += tbl_prev_p[i];
                        }
                        for i in 512..(512 + 256) {
                            tbl.tbl_coeff[1][i - 256] += tbl_prev_p[i];
                        }
                    } else {
                        let sel_idx = k + 128;
                        let k = k << (11 + (log_fold as u32));
                        for i in 0..256 {
                            let v = k + i;
                            let w = v + 256;
                            tbl.tbl_coeff[1][sel_idx] -= (tbl_prev_p[k + i] + tbl_prev_p[k + 256 + i]) *
                                (tbl.tbl_value[v] - tbl.tbl_value[w]);
                            tbl.tbl_coeff[0][sel_idx] += tbl_prev_p[k + i] * tbl.tbl_value[v];
                            tbl.tbl_coeff[0][sel_idx] += tbl_prev_p[k + 256 + i] * tbl.tbl_value[w];
                        }
                        for i in 512..(512 + 256) {
                            tbl.tbl_coeff[1][k + i - 512] += tbl_prev_p[k + i];
                        }
                        for i in (512 + 256)..(512 + 512) {
                            tbl.tbl_coeff[1][k + i - 512] -= tbl_prev_p[k + i] * tbl.tbl_value[k + i - 512];
                            tbl.tbl_coeff[0][k + i - 512] += tbl_prev_p[k + i] * tbl.tbl_value[k + i - 512];
                        }
                    }
                    for m in 1..(1 << log_fold) {
                        let k = ((k << log_fold) + m) << 11;
                        for i in 0..2048 {
                            tbl.tbl_coeff[1][k + i] += tbl_prev_p[k + i];
                        }
                    }
                }
            });
        layer.sumcheck_init[1] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_size]; 4];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                let tbl_p = &tbl.tbl_point[0];
                for k in 0..path_size {
                    if k == 0 {
                        for i in 0..128 {
                            tbl.tbl_coeff[2][0] += tbl_prev_p[i] * tbl_p[i];
                        }
                        for i in 512..(512 + 256) {
                            tbl.tbl_coeff[2][0] += tbl_prev_p[i] * tbl_p[i - 256];
                        }
                    } else {
                        let sel_idx = k + 128;
                        let k = k << (11 + (log_fold as u32));
                        for i in 0..256 {
                            let v = k + i;
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
                        for i in 512..(512 + 256) {
                            tbl.tbl_coeff[2][0] += tbl_prev_p[k + i] * tbl_p[k + i - 512];
                        }
                        for i in (512 + 256)..(512 + 512) {
                            let tmp = tbl_prev_p[k + i] * tbl_p[k + i - 512];
                            tbl.tbl_coeff[3][k + i - 512] -= tmp;
                            tbl.tbl_coeff[1][k + i - 512] += tmp;
                        }
                    }
                    for m in 1..(1 << log_fold) {
                        let k = ((k << log_fold) + m) << 11;
                        for i in 0..2048 {
                            tbl.tbl_coeff[2][0] += tbl_prev_p[k + i] * tbl_p[k + i];
                        }
                    }
                }
            });

        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_size];
                for k in 0..path_size {
                    if k == 0 {
                        for i in 0..128 {
                            // relay of ring-key
                            res[i] = inputs[i];
                        }
                        for i in 512..(512 + 256) {
                            // relay of root-hash
                            res[i] = inputs[i - 256];
                        }
                    } else {
                        let sel_idx = k + 128;
                        let k = k << (11 + (log_fold as u32));
                        for i in 0..256 {
                            if inputs[sel_idx] == T::from(0) {
                                res[k + i] = inputs[k + i];
                                res[k + i + 256] = inputs[k + i + 256];
                            } else {
                                res[k + i] = inputs[k + i + 256];
                                res[k + i + 256] = inputs[k + i];
                            }
                        }
                        for i in 512..(512 + 256) {
                            res[k + i] = inputs[k + i - 512];
                        }
                        for i in (512 + 256)..(512 + 512) {
                            res[k + i] = inputs[k + i - 512] *
                                (T::from(1) - inputs[k + i - 512]);
                        }
                    }
                    for m in 1..(1 << log_fold) {
                        let k = ((k << log_fold) + m) << 11;
                        for i in 0..2048 {
                            res[k + i] = inputs[k + i];
                        }
                    }
                }
                res
            });

        let c_path_size_bits = path_size_bits.clone();
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>| {
                let mut block_eq2 = T::from(1);
                let mut block_eq3 = T::from(1);
                let mut idx = 0;
                let tau_ofs = log_fold + 11;
                let rho_ofs = log_fold + 11;
                let sgm_ofs = log_fold + 11;
                for i in 0..c_path_size_bits.len() {
                    let idx0 = if i > 0 { idx + 1 } else { 0 };
                    idx = c_path_size_bits[i];
                    if i == 0 {
                        for i in 0..idx {
                            block_eq2 *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                            block_eq3 *= eq3(tau[tau_ofs + i], rho[rho_ofs + i], sgm[sgm_ofs + i]);
                        }
                        if tau_ofs + idx < tau.len() {
                            let tmp = zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                            block_eq2 *= tmp;
                            block_eq3 *= tmp * zer(sgm[sgm_ofs + idx]);
                        }
                    } else {
                        let mut tmp2 = T::from(1);
                        let mut tmp3 = T::from(1);
                        for i in 0..idx {
                            tmp2 *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                            tmp3 *= eq3(tau[tau_ofs + i], rho[rho_ofs + i], sgm[sgm_ofs + i]);
                        }
                        let tmp = zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                        tmp2 *= tmp;
                        tmp3 *= tmp * zer(sgm[sgm_ofs + idx]);

                        for i in idx0..idx {
                            let tmp = zer(tau[tau_ofs + i]) * zer(rho[rho_ofs + i]);
                            block_eq2 *= tmp;
                            block_eq3 *= tmp * zer(sgm[sgm_ofs + i]);
                        }
                        let tmp =  one(tau[tau_ofs + idx]) * one(rho[rho_ofs + idx]);
                        block_eq2 *= tmp;
                        block_eq2 += tmp2;
                        block_eq3 *= tmp * one(sgm[sgm_ofs + idx]);
                        block_eq3 += tmp3;
                    }
                }

                let mut block_first = T::from(1);
                let mut block_first3 = T::from(1);
                for i in tau_ofs..tau.len() {
                    let tmp = zer(tau[i]) * zer(rho[i]);
                    block_first *= tmp;
                    block_first3 *= tmp * zer(sgm[i]);
                }
                let block_rest = block_eq2 - block_first;
                let block_rest3 = block_eq3 - block_first3;

                let mut fold_eq2 = T::from(1);
                let mut fold_eq3 = T::from(1);
                for i in 11..(log_fold + 11) {
                    fold_eq2 *= eq2(tau[i], rho[i]);
                    fold_eq3 *= eq3(tau[i], rho[i], sgm[i]);
                }
                let mut fold_first = T::from(1);
                let mut fold_first3 = T::from(1);
                for i in 11..(log_fold + 11) {
                    let tmp = zer(tau[i]) * zer(rho[i]);
                    fold_first *= tmp;
                    fold_first3 *= tmp * zer(sgm[i]);
                }
                let fold_rest = fold_eq2 - fold_first;

                let mut sgm_zer = T::from(1);
                for i in 0..sgm.len() {
                    sgm_zer *= zer(sgm[i]);
                }

                let mut eq2_128 = T::from(1);
                for i in 0..7 {
                    eq2_128 *= eq2(tau[i], rho[i]);
                }
                let eq2_256 = eq2_128 * eq2(tau[7], rho[7]);
                let relay_first_first =
                    (eq2_128 * zer(tau[7]) * zer(rho[7]) * zer(tau[9]) * zer(rho[8]) +
                     eq2_256 * one(tau[9]) * one(rho[8])) *
                    zer(tau[8]) * zer(rho[9]) * zer(tau[10]) * zer(rho[10]);

                let relay_rest_first =
                    eq2_256 * zer(rho[8]) * zer(rho[9]) * zer(rho[10]) *
                    zer(tau[8]) * one(tau[9]) * zer(tau[10]);

                let mut eq2_2048 = eq2_256;
                for i in 8..11 {
                    eq2_2048 *= eq2(tau[i], rho[i]);
                }
                let relay_all_rest = eq2_2048;

                let mut eq3_256 = T::from(1);
                for i in 0..8 {
                    eq3_256 *= eq3(tau[i], rho[i], sgm[i]);
                }

                let mut naab = eq3_256 *
                    one(rho[8]) * one(sgm[8]) *
                    zer(rho[9]) * zer(sgm[9]) *
                    zer(rho[10]) * zer(sgm[10]) *
                    one(tau[8]) * one(tau[9]) * zer(tau[10]);
                naab *= block_rest3 * fold_first3;

                let mut block_sel_eq3 = T::from(1);
                let mut idx = 0;
                let tau_ofs = log_fold + 11;
                let rho_ofs = 0;
                let sgm_ofs = log_fold + 11;
                for i in 0..c_path_size_bits.len() {
                    let idx0 = if i > 0 { idx + 1 } else { 0 };
                    idx = c_path_size_bits[i];
                    if i == 0 {
                        for i in 0..idx {
                            block_sel_eq3 *= eq3(tau[tau_ofs + i], rho[rho_ofs + i], sgm[sgm_ofs + i]);
                        }
                        if tau_ofs + idx < tau.len() {
                            block_sel_eq3 *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]) * zer(sgm[sgm_ofs + idx]);
                        }
                    } else {
                        let mut tmp3 = T::from(1);
                        for i in 0..idx {
                            tmp3 *= eq3(tau[tau_ofs + i], rho[rho_ofs + i], sgm[sgm_ofs + i]);
                        }
                        tmp3 *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]) * zer(sgm[sgm_ofs + idx]);

                        for i in idx0..idx {
                            block_sel_eq3 *= zer(tau[tau_ofs + i]) * zer(rho[rho_ofs + i]) * zer(sgm[sgm_ofs + i]);
                        }
                        block_sel_eq3 *= one(tau[tau_ofs + idx]) * one(rho[rho_ofs + idx]) * one(sgm[sgm_ofs + idx]);
                        block_sel_eq3 += tmp3;
                    }
                }

                let mut block_sel_first3 = T::from(1);
                for i in tau_ofs..tau.len() {
                    block_sel_first3 *= zer(tau[i]) * zer(rho[i - tau_ofs]) * zer(sgm[i]);
                }
                let block_sel_rest3 = block_sel_eq3 - block_sel_first3;
                let mut selector = block_sel_rest3 * one(rho[7]);
                for i in (tau.len() - tau_ofs)..7 {
                    selector *= zer(rho[i]);
                }
                for i in 8..11 {
                    selector *= zer(rho[i]);
                }
                for i in (11 + log_fold)..rho.len() {
                    selector *= zer(rho[i]);
                }
                for i in 0..8 {
                    selector *= eq2(tau[i], sgm[i]);
                }
                selector *= fold_first3 * zer(tau[9]) * zer(sgm[9]) * zer(tau[10]) * zer(sgm[10]);

                [naab + selector,
                 (relay_all_rest * fold_rest * block_eq2 +
                  relay_rest_first * fold_first * block_rest +
                  relay_first_first * fold_first * block_first) * sgm_zer,
                 naab + selector * eq2(tau[8], sgm[8]),
                 T::from(0)]
            }));

        /* gate description
        for k in 0..path_size {
            if k == 0 {
                for i in 0..128 {
                    // relay of ring-key
                    let g = &mut layer.gates[i];
                    g.gtype = GateType::RELAY;
                    g.inputs = vec![i as u32];
                }
                for i in 256..(256 + 256) {
                    // relay of root-hash
                    let g = &mut layer.gates[i + 256];
                    g.gtype = GateType::RELAY;
                    // g.inputs = vec![(i - 256) as u32];
                    g.inputs = vec![i as u32];
                }
            } else {
                let bk = k as u32;
                let sel_idx = bk + 128;
                let bk = bk << (11 + (log_fold as u32));
                for i in 0..256 {
                    // selected of left child
                    let g = &mut layer.gates[(k << (11 + log_fold)) + i];
                    g.gtype = GateType::SELECTOR;
                    let i = i as u32;
                    g.inputs = vec![sel_idx, bk + i, bk + i + 256];
                }
                for i in 256..512 {
                    // selected of right child
                    let g = &mut layer.gates[(k << (11 + log_fold)) + i];
                    g.gtype = GateType::SELECTOR;
                    let i = (i - 256) as u32;
                    g.inputs = vec![sel_idx, bk + i + 256, bk + i];
                }
                for i in 512..(512 + 256) {
                    // relay of two hashes
                    let g = &mut layer.gates[(k << (11 + log_fold)) + i];
                    g.gtype = GateType::RELAY;
                    let i = i as u32;
                    g.inputs = vec![bk + i - 512];
                }
                for i in (512 + 256)..(512 + 512) {
                    // relay of two hashes
                    let g = &mut layer.gates[(k << (11 + log_fold)) + i];
                    g.gtype = GateType::NAAB;
                    let i = i as u32;
                    g.inputs = vec![bk + i - 512, bk + i - 512];
                }
            }
            for m in 1..(1 << log_fold) {
                let bk = ((k << log_fold) + m) << 11;
                for i in 0..2048 {
                    let g = &mut layer.gates[bk + i];
                    g.gtype = GateType::RELAY;
                    g.inputs = vec![(bk + i) as u32];
                }
            }
        }
        */
        depth += 1;
    }

    for d in 0..(24 / (1 << log_fold)) {
        // linear layers
        {
            c.layers.push(Layer::<T>::default());
            let layer = &mut c.layers[depth];
            let log_size = log_num_hash + log_fold + 12;
            layer.bit_len = log_size;
            // layer.gates = vec![Gate::default(); 1 << log_size];
            // layer.values = vec![T::from(0); 1 << log_size];
            layer.num_phase = 1;
            layer.num_table = vec![1];
            layer.val_degree = vec![1];
            layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
            layer.sumcheck_init[0] = Arc::new(
                move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                    tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size - if d == 0 { 1 } else { 0 })]; 1];
                    let tbl_prev_p = &tbl.tbl_prev_point[0];
                    for k in 0..path_size {
                        for j in 0..(1 << log_fold) {
                            let ofs = if d != 0 { ((k << log_fold) + j) << 12 } else { ((k << log_fold) + j) << 11 };
                            let mut i = ((k << log_fold) + j) << 12;
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
                                        if d != 0 || j != 0 || a_idx * 5 + a_idy < 8 {
                                            tbl.tbl_coeff[0][ofs + (a_idx * 5 + a_idy) * 64 + a_idz] += tbl_prev_p[i];
                                        }
                                        for c_idy in 0..=4 {
                                            if d != 0 || j != 0 || c0_idx * 5 + c_idy < 8 {
                                                tbl.tbl_coeff[0][ofs + (c0_idx * 5 + c_idy) * 64 + c0_idz] += tbl_prev_p[i];
                                            }
                                        }
                                        for c_idy in 0..=4 {
                                            if d != 0 || j != 0 || c1_idx * 5 + c_idy < 8 {
                                                tbl.tbl_coeff[0][ofs + (c1_idx * 5 + c_idy) * 64 + c1_idz] += tbl_prev_p[i];
                                            }
                                        }
                                        i += 1;
                                    }
                                }
                            }
                            if d == 0 && j == 0 {
                                for i in 0..512 {
                                    tbl.tbl_coeff[0][ofs + 512 + i] += tbl_prev_p[(ofs << 1) + 2048 + i];
                                }
                            } else if d == 0 {
                                for i in 0..2048 {
                                    tbl.tbl_coeff[0][ofs + i] += tbl_prev_p[(ofs << 1) + 2048 + i];
                                }
                            } else {
                                for i in 0..2048 {
                                    tbl.tbl_coeff[0][ofs + 2048 + i] += tbl_prev_p[ofs + 2048 + i];
                                }
                            }
                        }
                    }
                });
            layer.evaluate[0] = Arc::new(
                move |_: &Layer<T>, inputs: &Vec<T>| {
                    let mut res = vec![T::from(0); 1 << log_size];
                    for k in 0..path_size {
                        for j in 0..(1 << log_fold) {
                            let ofs = if d != 0 { ((k << log_fold) + j)  << 12 } else { ((k << log_fold) + j) << 11 };
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let a_idx = (15 + y - 3 * x) * 3 % 5;
                                        let a_idy = x;
                                        // assert_eq!((a_idx * 2 + a_idy * 3) % 5, y);
                                        let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;

                                        let c0_idx = (a_idx + 4) % 5;
                                        let c0_idz = a_idz;
                                        let c1_idx = (a_idx + 1) % 5;
                                        let c1_idz = (a_idz + 63) % 64;

                                        if (d != 0 || j != 0) || a_idx * 5 + a_idy < 8 {
                                            res[(((k << log_fold) + j) << 12) + (x * 5 + y) * 64 + z] += inputs[ofs + (a_idx * 5 + a_idy) * 64 + a_idz];
                                        }
                                        for c0_idy in 0..=4 {
                                            if (d != 0 || j != 0) || c0_idx * 5 + c0_idy < 8 {
                                                res[(((k << log_fold) + j) << 12) + (x * 5 + y) * 64 + z] += inputs[ofs + (c0_idx * 5 + c0_idy) * 64 + c0_idz];
                                            }
                                        }
                                        for c1_idy in 0..=4 {
                                            if (d != 0 || j != 0) || c1_idx * 5 + c1_idy < 8 {
                                                res[(((k << log_fold) + j) << 12) + (x * 5 + y) * 64 + z] += inputs[ofs + (c1_idx * 5 + c1_idy) * 64 + c1_idz];
                                            }
                                        }
                                    }
                                }
                            }
                            if d == 0 && j == 0 {
                                for i in 0..512 {
                                    res[(ofs << 1) + 2048 + i] = inputs[ofs + 512 + i];
                                }
                            } else if d == 0 {
                                for i in 0..2048 {
                                    res[(ofs << 1) + 2048 + i] = inputs[ofs + i];
                                }
                            } else {
                                for i in 0..2048 {
                                    res[ofs + 2048 + i] = inputs[ofs + 2048 + i];
                                }
                            }
                        }
                    }
                    res
                });

            if d == 0 {
                let c_path_size_bits = path_size_bits.clone();
                layer.predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                        let mut block_eq2 = T::from(1);
                        let mut idx = 0;
                        let tau_ofs = log_fold + 12;
                        let rho_ofs = log_fold + 11;
                        for i in 0..c_path_size_bits.len() {
                            let idx0 = if i > 0 { idx + 1 } else { 0 };
                            idx = c_path_size_bits[i];
                            if i == 0 {
                                for i in 0..idx {
                                    block_eq2 *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                                }
                                if tau_ofs + idx < tau.len() {
                                    block_eq2 *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                                }
                            } else {
                                let mut tmp2 = T::from(1);
                                for i in 0..idx {
                                    tmp2 *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                                }
                                tmp2 *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                                for i in idx0..idx {
                                    block_eq2 *= zer(tau[tau_ofs + i]) * zer(rho[rho_ofs + i]);
                                }
                                block_eq2 *= one(tau[tau_ofs + idx]) * one(rho[rho_ofs + idx]);
                                block_eq2 += tmp2;
                            }
                        }

                        let mut fold_eq2 = T::from(1);
                        for i in 11..(log_fold + 11) {
                            fold_eq2 *= eq2(tau[i + 1], rho[i]);
                        }
                        let mut fold_first = T::from(1);
                        for i in 11..(log_fold + 11) {
                            fold_first *= zer(tau[i + 1]) * zer(rho[i]);
                        }
                        let fold_rest = fold_eq2 - fold_first;

                        let mut tbl0 = vec![T::from(0); 1 << 12];
                        let mut tbl1 = vec![T::from(0); 1 << 11];
                        // TODO XXX reuse tbl0, tbl1 for next layer
                        tbl0[0] = T::from(1);
                        tbl1[0] = T::from(1);
                        for k in 0..12 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl0[i + mask] = tbl0[i] * tau[k];
                                tbl0[i] = tbl0[i] - tbl0[i + mask];
                            }
                        }
                        for k in 0..11 {
                            let mask = 1 << k;
                            for i in 0..mask {
                                tbl1[i + mask] = tbl1[i] * rho[k];
                                tbl1[i] = tbl1[i] - tbl1[i + mask];
                            }
                        }
                        let mut sum_rest = T::from(0);
                        let mut sum_first = T::from(0);
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
                                    sum_rest += tbl0[i] * tbl1[(a_idx * 5 + a_idy) * 64 + a_idz];
                                    if a_idx * 5 + a_idy < 8 {
                                        sum_first += tbl0[i] * tbl1[(a_idx * 5 + a_idy) * 64 + a_idz];
                                    }
                                    for c_idy in 0..=4 {
                                        sum_rest += tbl0[i] * tbl1[(c0_idx * 5 + c_idy) * 64 + c0_idz];
                                        if c0_idx * 5 + c_idy < 8 {
                                            sum_first += tbl0[i] * tbl1[(c0_idx * 5 + c_idy) * 64 + c0_idz];
                                        }
                                    }
                                    for c_idy in 0..=4 {
                                        sum_rest += tbl0[i] * tbl1[(c1_idx * 5 + c_idy) * 64 + c1_idz];
                                        if c1_idx * 5 + c_idy < 8 {
                                            sum_first += tbl0[i] * tbl1[(c1_idx * 5 + c_idy) * 64 + c1_idz];
                                        }
                                    }
                                    i += 1;
                                }
                            }
                        }

                        let mut eq_512 = T::from(1);
                        for i in 0..9 {
                            eq_512 *= eq2(rho[i], tau[i]);
                        }
                        let mut eq_2048 = eq_512;
                        for i in 9..11 {
                            eq_2048 *= eq2(rho[i], tau[i]);
                        }
                        let relay_first = eq_512 * one(rho[9]) * zer(tau[9]) * zer(tau[10]) * zer(rho[10]) * one(tau[11]);
                        let relay_rest = eq_2048 * one(tau[11]);

                        [T::from(0),
                         ((relay_first + sum_first) * fold_first +
                          (relay_rest + sum_rest) * fold_rest) * block_eq2,
                         T::from(0), T::from(0)]
                    }));
            } else {
                let c_path_size_bits = path_size_bits.clone();
                layer.predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                        let mut block_eq2 = T::from(1);
                        let mut idx = 0;
                        let tau_ofs = log_fold + 12;
                        let rho_ofs = log_fold + 12;
                        for i in 0..c_path_size_bits.len() {
                            let idx0 = if i > 0 { idx + 1 } else { 0 };
                            idx = c_path_size_bits[i];
                            if i == 0 {
                                for i in 0..idx {
                                    block_eq2 *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                                }
                                if tau_ofs + idx < tau.len() {
                                    block_eq2 *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                                }
                            } else {
                                let mut tmp2 = T::from(1);
                                for i in 0..idx {
                                    tmp2 *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                                }
                                tmp2 *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                                for i in idx0..idx {
                                    block_eq2 *= zer(tau[tau_ofs + i]) * zer(rho[rho_ofs + i]);
                                }
                                block_eq2 *= one(tau[tau_ofs + idx]) * one(rho[rho_ofs + idx]);
                                block_eq2 += tmp2;
                            }
                        }

                        let mut fold_eq2 = T::from(1);
                        for i in 12..(log_fold + 12) {
                            fold_eq2 *= eq2(tau[i], rho[i]);
                        }

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
                        let mut sum = T::from(0);
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

                        let mut relay = T::from(1);
                        for i in 0..11 {
                            relay *= eq2(rho[i], tau[i]);
                        }
                        relay *= one(rho[11]) * one(tau[11]);

                        [T::from(0),
                         (relay + sum) * block_eq2 * fold_eq2,
                         T::from(0), T::from(0)]
                    }));
            }

            /* gate description
            for k in 0..path_size {
                for j in 0..(1 << log_fold) {
                    let ofs = if d != 0 { ((k << log_fold) + j)  << 12 } else { ((k << log_fold) + j) << 11 };
                    for x in 0..=4 {
                        for y in 0..=4 {
                            for z in 0..64 {
                                let a_idx = (15 + y - 3 * x) * 3 % 5;
                                let a_idy = x;
                                // assert_eq!((a_idx * 2 + a_idy * 3) % 5, y);
                                let a_idz = (z - ROTATE[a_idx][a_idy]) % 64;

                                let c0_idx = (a_idx + 4) % 5;
                                let c0_idz = a_idz;
                                let c1_idx = (a_idx + 1) % 5;
                                let c1_idz = (a_idz + 63) % 64;

                                let g = &mut layer.gates[(((k << log_fold) + j) << 12) + (x * 5 + y) * 64 + z];
                                g.gtype = GateType::SUM;
                                g.inputs = vec![];

                                if (d != 0 || j != 0) || a_idx * 5 + a_idy < 8 {
                                    g.inputs.push((ofs + (a_idx * 5 + a_idy) * 64 + a_idz) as u32);
                                }
                                for c0_idy in 0..=4 {
                                    if (d != 0 || j != 0) || c0_idx * 5 + c0_idy < 8 {
                                        g.inputs.push((ofs + (c0_idx * 5 + c0_idy) * 64 + c0_idz) as u32);
                                    }
                                }
                                for c1_idy in 0..=4 {
                                    if (d != 0 || j != 0) || c1_idx * 5 + c1_idy < 8 {
                                        g.inputs.push((ofs + (c1_idx * 5 + c1_idy) * 64 + c1_idz) as u32);
                                    }
                                }
                            }
                        }
                    }
                    if d == 0 && j == 0 {
                        for i in 0..512 {
                            let g = &mut layer.gates[(ofs << 1) + 2048 + i];
                            g.gtype = GateType::RELAY;
                            g.inputs = vec![(ofs + 512 + i) as u32];
                        }
                    } else if d == 0 {
                        for i in 0..2048 {
                            let g = &mut layer.gates[(ofs << 1) + 2048 + i];
                            g.gtype = GateType::RELAY;
                            g.inputs = vec![(ofs + i) as u32];
                        }
                    } else {
                        for i in 0..2048 {
                            let g = &mut layer.gates[ofs + 2048 + i];
                            g.gtype = GateType::RELAY;
                            g.inputs = vec![(ofs + 2048 + i) as u32];
                        }
                    }
                }
            }
            */
            depth += 1;
        }

        // nonlinear layers
        {
            c.layers.push(Layer::<T>::default());
            let layer = &mut c.layers[depth];
            let log_size = log_num_hash + log_fold + 12;
            layer.bit_len = log_size;
            // layer.gates = vec![Gate::default(); 1 << log_size];
            layer.sumcheck_init[0] = Arc::new(
                move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                    tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_size]; 2];
                    let tbl_prev_p = &tbl.tbl_prev_point[0];
                    for k in 0..path_size {
                        for j in 0..(1 << log_fold) {
                            let ofs = ((k << log_fold) + j) << 12;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let w = ofs + (x * 5 + y) * 64 + z;
                                        let u = ofs + (((x + 1) % 5) * 5 + y) * 64 + z;
                                        let v = ofs + (((x + 2) % 5) * 5 + y) * 64 + z;
                                        tbl.tbl_coeff[1][u] -= tbl_prev_p[w] * tbl.tbl_value[v];
                                        tbl.tbl_coeff[0][u] += tbl_prev_p[w] * tbl.tbl_value[v];
                                        tbl.tbl_coeff[1][w] += tbl_prev_p[w];
                                        let dd = j * (24 / (1 << log_fold)) + d;
                                        if x == 0 && y == 0 && (RC[dd] & (1 << z) != 0) {
                                            tbl.tbl_coeff[0][w] += tbl_prev_p[w];
                                        }
                                    }
                                }
                            }
                            for i in 0..2048 {
                                let w = ofs + 2048 + i;
                                tbl.tbl_coeff[1][w] += tbl_prev_p[w];
                            }
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
                    for k in 0..path_size {
                        for j in 0..(1 << log_fold) {
                            let ofs = ((k << log_fold) + j) << 12;
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

                                        let dd = j * (24 / (1 << log_fold)) + d;
                                        if x == 0 && y == 0 && (RC[dd] & (1 << z) != 0) {
                                            tbl.tbl_coeff[0][0] += tmp;
                                        }
                                    }
                                }
                            }
                            for i in 0..2048 {
                                let w = ofs + 2048 + i;
                                tbl.tbl_coeff[2][0] += tbl_prev_p[w] * tbl_p[w];
                            }
                        }
                    }
                });
            layer.evaluate[0] = Arc::new(
                move |_: &Layer<T>, inputs: &Vec<T>| {
                    let mut res = vec![T::from(0); 1 << log_size];
                    for k in 0..path_size {
                        for j in 0..(1 << log_fold) {
                            let ofs = ((k << log_fold) + j) << 12;
                            for x in 0..=4 {
                                for y in 0..=4 {
                                    for z in 0..64 {
                                        let idx0 = ofs + (x * 5 + y) * 64 + z;
                                        let idx1 = ofs + (((x + 1) % 5) * 5 + y) * 64 + z;
                                        let idx2 = ofs + (((x + 2) % 5) * 5 + y) * 64 + z;
                                        res[idx0] = (T::from(1) - inputs[idx1]) * inputs[idx2] + inputs[idx0];

                                        let dd = j * (24 / (1 << log_fold)) + d;
                                        if x == 0 && y == 0 && (RC[dd] & (1 << z) != 0) {
                                            res[idx0] += T::from(1);
                                        }
                                    }
                                }
                            }
                            for i in 0..2048 {
                                let idx = ofs + 2048 + i;
                                res[idx] = inputs[idx];
                            }
                        }
                    }
                    res
                });
            let c_path_size_bits = path_size_bits.clone();
            layer.predicate.push(Arc::new(
                move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>| {
                    let mut block_eq2 = T::from(1);
                    let mut block_eq3 = T::from(1);
                    let mut idx = 0;
                    let tau_ofs = log_fold + 12;
                    let rho_ofs = log_fold + 12;
                    let sgm_ofs = log_fold + 12;
                    for i in 0..c_path_size_bits.len() {
                        let idx0 = if i > 0 { idx + 1 } else { 0 };
                        idx = c_path_size_bits[i];
                        if i == 0 {
                            for i in 0..idx {
                                block_eq2 *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                                block_eq3 *= eq3(tau[tau_ofs + i], rho[rho_ofs + i], sgm[sgm_ofs + i]);
                            }
                            if tau_ofs + idx < tau.len() {
                                let tmp = zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                                block_eq2 *= tmp;
                                block_eq3 *= tmp * zer(sgm[sgm_ofs + idx]);
                            }
                        } else {
                            let mut tmp2 = T::from(1);
                            let mut tmp3 = T::from(1);
                            for i in 0..idx {
                                tmp2 *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                                tmp3 *= eq3(tau[tau_ofs + i], rho[rho_ofs + i], sgm[sgm_ofs + i]);
                            }
                            let tmp = zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                            tmp2 *= tmp;
                            tmp3 *= tmp * zer(sgm[sgm_ofs + idx]);

                            for i in idx0..idx {
                                let tmp = zer(tau[tau_ofs + i]) * zer(rho[rho_ofs + i]);
                                block_eq2 *= tmp;
                                block_eq3 *= tmp * zer(sgm[sgm_ofs + i]);
                            }
                            let tmp =  one(tau[tau_ofs + idx]) * one(rho[rho_ofs + idx]);
                            block_eq2 *= tmp;
                            block_eq2 += tmp2;
                            block_eq3 *= tmp * one(sgm[sgm_ofs + idx]);
                            block_eq3 += tmp3;
                        }
                    }

                    let mut tbl0 = vec![T::from(0); 1 << 12];
                    let mut tbl1 = vec![T::from(0); 1 << 12];
                    let mut tbl2 = vec![T::from(0); 1 << 12];
                    // TODO XXX reuse tbl0, tbl1 for next layer
                    tbl0[0] = T::from(1);
                    tbl1[0] = T::from(1);
                    tbl2[0] = T::from(1);
                    for k in 0..(12) {
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
                                chi_mul += tbl0[i] * tbl1[(((x + 1) % 5) * 5 + y) * 64 + z] * tbl2[(((x + 2) % 5) * 5 + y) * 64 + z];
                                chi_add += tbl0[i] * tbl1[i];
                                i += 1;
                            }
                        }
                    }
                    for j in 0..(1 << log_fold) {
                        for z in 0..64 {
                            let dd = j * (24 / (1 << log_fold)) + d;
                            if RC[dd] & (1 << z) != 0 {
                                chi_cst += tbl0[z] * tbl1[z] *
                                    expand_binary_eq2(j, &tau[12..tau_ofs].to_vec(), &rho[12..rho_ofs].to_vec());
                            }
                        }
                    }

                    // relay
                    let mut relay = T::from(1);
                    for i in 0..11 {
                        relay *= eq2(tau[i], rho[i]);
                    }
                    relay *= one(tau[11]) * one(rho[11]);

                    let mut fold_eq2 = T::from(1);
                    for i in 12..(log_fold + 12) {
                        fold_eq2 *= eq2(tau[i], rho[i]);
                    }

                    let mut fold_eq3 = T::from(1);
                    for i in 12..(log_fold + 12) {
                        fold_eq3 *= eq3(tau[i], rho[i], sgm[i]);
                    }

                    let mut sgm_zer = T::from(1);
                    for i in 0..sgm.len() {
                        sgm_zer *= zer(sgm[i]);
                    }

                    chi_mul *= block_eq3 * fold_eq3;
                    [chi_mul, (relay + chi_add) * block_eq2 * sgm_zer * fold_eq2,
                     chi_mul, chi_cst * block_eq2 * sgm_zer]
                }));

            /* gate description
            for k in 0..path_size {
                for j in 0..(1 << log_fold) {
                    let ofs = ((k << log_fold) + j) << 12;
                    for x in 0..=4 {
                        for y in 0..=4 {
                            for z in 0..64 {
                                let g = &mut layer.gates[ofs + (x * 5 + y) * 64 + z];

                                // # ι step
                                // A[0,0] = A[0,0] xor RC
                                let dd = j * (24 / (1 << log_fold)) + d;
                                if x == 0 && y == 0 && (RC[dd] & (1 << z) != 0) {
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
                    for i in 0..2048 {
                        let g = &mut layer.gates[ofs + 2048 + i];
                        g.gtype = GateType::RELAY;
                        g.inputs = vec![(ofs + 2048 + i) as u32];
                    }
                }
            }
            */
            depth += 1;
        }
    }

    // last layer
    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_num_hash + log_fold + 11;
        layer.bit_len = log_size;
        // layer.gates = vec![Gate::default(); 1 << log_size];
        // layer.values = vec![T::from(0); 1 << log_size];

        layer.val_degree = vec![1];
        layer.num_table = vec![1];
        layer.num_phase = 1;
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_size + 1)]; 1];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                // for k in 0..(1 << log_num_hash) {
                for k in 0..path_size {
                    for i in 0..((1 << log_fold) - 1) {
                        for j in 0..2048 {
                            let x = (((k << log_fold) + i) << 11) + j;
                            let u = (((k << log_fold) + i) << 12) + j;
                            let v = (((k << log_fold) + i + 1) << 12) + j + 2048;
                            tbl.tbl_coeff[0][u] += tbl_prev_p[x];
                            tbl.tbl_coeff[0][v] += tbl_prev_p[x];
                        }
                    }
                    if k != path_size - 1 {
                        for j in 0..256 {
                            let i = (1 << log_fold) - 1;
                            let x = (((k << log_fold) + i) << 11) + j;
                            let u = (((k << log_fold) + i) << 12) + j;
                            let v = ((((k + 1) << log_fold) + 0) << 12) + j + 2048;
                            tbl.tbl_coeff[0][u] += tbl_prev_p[x];
                            tbl.tbl_coeff[0][v] += tbl_prev_p[x];
                        }
                    } else {
                        for j in 0..256 {
                            let i = (1 << log_fold) - 1;
                            let x = (((k << log_fold) + i) << 11) + j;
                            let u = (((k << log_fold) + i) << 12) + j;
                            let v = (((0 << log_fold) + 0) << 12) + j + 2048;
                            tbl.tbl_coeff[0][u] += tbl_prev_p[x];
                            tbl.tbl_coeff[0][v] += tbl_prev_p[x];
                        }
                    }

                    for j in 0..256 {
                        let i = (1 << log_fold) - 1;
                        let x = (((k << log_fold) + i) << 11) + j + 256;
                        let u = (((k << log_fold) + 0) << 12) + 2048 + j + 256;
                        tbl.tbl_coeff[0][u] += tbl_prev_p[x];
                    }
                }
            });

        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_size];
                // for k in 0..(1 << log_num_hash) {
                for k in 0..path_size {
                    for i in 0..((1 << log_fold) - 1) {
                        for j in 0..2048 {
                            res[(((k << log_fold) + i) << 11) + j] =
                                inputs[(((k << log_fold) + i) << 12) + j] -
                                inputs[(((k << log_fold) + i + 1) << 12) + j + 2048];
                        }
                    }
                    if k != path_size - 1 {
                        for j in 0..256 {
                            let i = (1 << log_fold) - 1;
                            res[(((k << log_fold) + i) << 11) + j] =
                                inputs[(((k << log_fold) + i) << 12) + j] -
                                inputs[((((k + 1) << log_fold) + 0) << 12) + j + 2048];
                        }
                    } else {
                        for j in 0..256 {
                            let i = (1 << log_fold) - 1;
                            res[(((k << log_fold) + i) << 11) + j] =
                                inputs[(((k << log_fold) + i) << 12) + j] -
                                inputs[(((0 << log_fold) + 0) << 12) + j + 2048];
                        }
                    }

                    for j in 0..256 {
                        let i = (1 << log_fold) - 1;
                        res[(((k << log_fold) + i) << 11) + j + 256] =
                            inputs[(((k << log_fold) + 0) << 12) + 2048 + j + 256];
                    }
                }
                res
            });

        /* TODO fastest 0 <= k < path_size
        // DP[path_size] = min{
        //          0<=k<floor(path_size) and (floor(path_size) + DP[path_size - floor(path_size)],
        //          0<=k<ceil(path_size)  and  (ceil(path_size) - DP[ceil(path_size) - path_size],
        // }
        // (listof(log_range, log_offset), num_mul)
        let mut dp: Vec<Vec<(isize, isize)>> = vec![vec![]; path_size + 1];
        let mut dp_num_mul: Vec<usize> = vec![0; path_size + 1];
        let mut log_idx = 0;
        while (1 << log_idx) <= path_size {
            dp[1 << log_idx].push((log_idx, -1));
            dp_num_mul[1 << log_idx] = log_num_hash - log_idx;
            log_idx += 1;
        }
        for i in 1..(path_size + 1) {
            let log_i = ceil_log(i);
            if (1 << log_i) != i {
                let floor = 1 << (log_i - 1);
                let ceil = 1 << log_i;
                let floor_num = dp_num_mul[floor] + dp_num_mul[i - floor];
                let ceil_num = dp_num_mul[ceil] + dp_num_mul[ceil - i];
                if floor_num < ceil_num {
                    dp_num_mul[i] = floor_num;
                    dp[i] = dp[floor].clone();
                    for d in dp[i - floor].iter() {
                        match d {
                            (log_range, log_offset) => {
                                dp[i].push(log_range, )
                            },
                            _ => unreachable!(),
                        }
                    }
                }
            }
        }
        // end TODO
        */

        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                // tau j: 0~10
                // tau i: 11~log_fold+10
                // tau k: log_fold+11 ~ log_fold+log_num_hash+10
                // rho j: 0~11
                // rho i: 12~log_fold+11
                // rho k: log_fold+12 ~ log_fold+log_num_hash+11
                let mut eq_blockid = T::from(1);
                let mut idx = 0;
                let tau_ofs = log_fold + 11;
                let rho_ofs = log_fold + 12;
                for i in 0..path_size_bits.len() {
                    let idx0 = if i > 0 { idx + 1 } else { 0 };
                    idx = path_size_bits[i];
                    if i == 0 {
                        for i in 0..idx {
                            eq_blockid *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                        }
                        if tau_ofs + idx < tau.len() {
                            eq_blockid *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                        }
                    } else {
                        let mut tmp = T::from(1);
                        for i in 0..idx {
                            tmp *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                        }
                        tmp *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);

                        for i in idx0..idx {
                            eq_blockid *= zer(tau[tau_ofs + i]) * zer(rho[rho_ofs + i]);
                        }
                        eq_blockid *= one(tau[tau_ofs + idx]) * one(rho[rho_ofs + idx]);
                        eq_blockid += tmp;
                    }
                }

                let mut eq_blockid_nomax = T::from(1);
                let mut idx = 0;
                let tau_ofs = log_fold + 11;
                let rho_ofs = log_fold + 12;
                for i in 0..path_size_minus_bits.len() {
                    let idx0 = if i > 0 { idx + 1 } else { 0 };
                    idx = path_size_minus_bits[i];
                    if i == 0 {
                        for i in 0..idx {
                            eq_blockid_nomax *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                        }
                        if tau_ofs + idx < tau.len() {
                            eq_blockid_nomax *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                        }
                    } else {
                        let mut tmp = T::from(1);
                        for i in 0..idx {
                            tmp *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                        }
                        tmp *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);

                        for i in idx0..idx {
                            eq_blockid_nomax *= zer(tau[tau_ofs + i]) * zer(rho[rho_ofs + i]);
                        }
                        eq_blockid_nomax *= one(tau[tau_ofs + idx]) * one(rho[rho_ofs + idx]);
                        eq_blockid_nomax += tmp;
                    }
                }
                if path_size_minus_bits.len() == 0 {
                    eq_blockid_nomax = T::from(0);
                }

                let mut blockid_inc = T::from(0);
                for i in 0..(path_size - 1) {
                    let i1 = i + 1;
                    blockid_inc +=
                        expand_binary_eq(i, &tau[tau_ofs..].to_vec()) *
                        expand_binary_eq(i1, &rho[rho_ofs..].to_vec());
                }
                // TODO faster blockid_inc?
                /*
                let mut blockid_inc = T::from(1);
                let mut idx = 0;
                let tau_ofs = log_fold + 11;
                let rho_ofs = log_fold + 12;
                for i in 0..path_size_bits.len() {
                    let idx0 = if i > 0 { idx + 1 } else { 0 };
                    idx = path_size_bits[i];

                    // let mut eq_foldid_hi = vec![T::from(1); idx - 1];
                    // eq_foldid_hi[idx - 1] = eq2(tau[11 + idx], rho[12 + idx]);
                    // for i in 2..=idx {
                    //     eq_foldid_hi[idx - i] = eq_foldid_hi[idx - i + 1] *
                    //         eq2(tau[11 + idx - i + 1], rho[12 + idx - i + 1]);
                    // }

                    if i == 0 {
                        for i in 0..idx {
                            blockid_inc *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                        }
                        blockid_inc *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);
                    } else {
                        let mut tmp = T::from(1);
                        for i in 0..idx {
                            tmp *= eq2(tau[tau_ofs + i], rho[rho_ofs + i]);
                        }
                        tmp *= zer(tau[tau_ofs + idx]) * zer(rho[rho_ofs + idx]);

                        for i in idx0..idx {
                            blockid_inc *= zer(tau[tau_ofs + i]) * zer(rho[rho_ofs + i]);
                        }
                        blockid_inc *= one(tau[tau_ofs + idx]) * one(rho[rho_ofs + idx]);
                        blockid_inc += tmp;
                    }
                }
                if path_size_minus_bits.len() == 0 {
                    blockid_inc = T::from(0);
                }
                // end TODO
                */

                let mut eq_foldid_hi = vec![T::from(1); log_fold];
                eq_foldid_hi[log_fold - 1] = eq2(tau[tau_ofs - 1], rho[rho_ofs - 1]);
                for i in 2..=log_fold {
                    eq_foldid_hi[log_fold - i] = eq_foldid_hi[log_fold - i + 1] *
                        eq2(tau[tau_ofs - i], rho[rho_ofs - i]);
                }
                let mut foldid_inc = zer(tau[tau_ofs - 1]) * one(rho[rho_ofs - 1]);
                for i in 2..=log_fold {
                    foldid_inc =
                        eq_foldid_hi[log_fold - i + 1] * zer(tau[tau_ofs - i]) * one(rho[rho_ofs - i]) +
                        foldid_inc * one(tau[tau_ofs - i]) * zer(rho[rho_ofs - i]);
                }

                let mut eq_foldid = T::from(1);
                for i in 11..(11 + log_fold) {
                    eq_foldid *= eq2(tau[i], rho[i + 1]);
                }
                let mut foldid_max = T::from(1);
                for i in 11..(11 + log_fold) {
                    foldid_max *= one(tau[i]) * one(rho[i + 1]);
                }
                let foldid_nomax = eq_foldid + foldid_max;

                let mut eq_256 = T::from(1);
                for i in 0..8 {
                    eq_256 *= eq2(tau[i], rho[i]);
                }
                for i in 9..11 {
                    eq_256 *= zer(tau[i]) * zer(rho[i]);
                }

                // fold witness check
                let mut fold_witness = T::from(1);
                for i in 0..11 {
                    fold_witness *= eq2(tau[i], rho[i]);
                }
                fold_witness *=
                    (zer(rho[11]) * foldid_nomax +
                     one(rho[11]) * foldid_inc) * eq_blockid;

                // block witness check
                let mut block_witness =
                    foldid_max * eq_blockid_nomax *
                    zer(rho[8]) * zer(tau[8]) * zer(rho[11]);
                let mut tmp = T::from(1);
                for i in 11..(11 + log_fold) {
                    tmp *= one(tau[i]) * zer(rho[i + 1]);
                }
                block_witness += tmp * blockid_inc * zer(rho[8]) * zer(tau[8]) * one(rho[11]);

                // k == path_size - 1, compare roothash
                let mut tmp = path_size - 1;
                let mut roothash_o = T::from(1);
                let mut roothash_i0 = T::from(1);
                let mut roothash_i1 = T::from(1);
                for i in (log_fold + 11)..tau.len() {
                    roothash_i0 *= zer(rho[i + 1]);
                    if (tmp & 1) == 1 {
                        roothash_o *= one(tau[i]);
                        roothash_i1 *= one(rho[i + 1]);
                    } else {
                        roothash_o *= zer(tau[i]);
                        roothash_i1 *= zer(rho[i + 1]);
                    }
                    tmp >>= 1;
                }
                for i in 11..(log_fold + 11) {
                    roothash_o *= one(tau[i]);
                    roothash_i0 *= zer(rho[i + 1]);
                    roothash_i1 *= one(rho[i + 1]);
                }
                let roothash =
                    (roothash_i0 * one(rho[11]) +
                     roothash_i1 * zer(rho[11])) *
                    roothash_o * zer(tau[8]) * zer(rho[8]);

                // relay
                // 0 <= tau_j < 256;
                // rho_j = 2048 + tau_j;
                // tau_i = (1 << log_fold) - 1;
                // rho_i = 0;
                // 0 <= tau_k < path_size;
                // rho_k = tau_k;
                let mut relay = T::from(1);
                // for i in 0..8 {
                //     relay *= eq2(tau[i], rho[i]);
                // }
                relay *= one(tau[8]) * one(rho[8]);
                // for i in 9..11 {
                //     relay *= zer(tau[i]) * zer(rho[i]);
                // }
                relay *= one(rho[11]);
                for i in 11..(log_fold + 11) {
                    relay *= one(tau[i]) * zer(rho[i + 1]);
                }
                relay *= eq_blockid;

                [T::from(0), (relay + roothash + block_witness) * eq_256 + fold_witness, T::from(0), T::from(0)]
            }));

        /* gate description
        // for k in 0..(1 << log_num_hash) {
        for k in 0..path_size {
            for i in 0..((1 << log_fold) - 1) {
                for j in 0..2048 {
                    let g = &mut layer.gates[(((k << log_fold) + i) << 11) + j];
                    g.gtype = GateType::SUM;
                    g.inputs = vec![
                        ((((k << log_fold) + i) << 12) + j) as u32,
                        ((((k << log_fold) + i + 1) << 12) + j + 2048) as u32,
                    ]
                }
            }
            if k != path_size - 1 {
                for j in 0..256 {
                    let i = (1 << log_fold) - 1;
                    let g = &mut layer.gates[(((k << log_fold) + i) << 11) + j];
                    g.gtype = GateType::SUM;
                    g.inputs = vec![
                        ((((k << log_fold) + i) << 12) + j) as u32,
                        (((((k + 1) << log_fold) + 0) << 12) + j + 2048) as u32,
                    ]
                }
            } else {
                for j in 0..256 {
                    let i = (1 << log_fold) - 1;
                    let g = &mut layer.gates[(((k << log_fold) + i) << 11) + j];
                    g.gtype = GateType::SUM;
                    g.inputs = vec![
                        ((((k << log_fold) + i) << 12) + j) as u32,
                        ((((0 << log_fold) + 0) << 12) + j + 2048) as u32,
                    ]
                }
            }

            for j in 0..256 {
                let i = (1 << log_fold) - 1;
                let g = &mut layer.gates[(((k << log_fold) + i) << 11) + j + 256];
                g.gtype = GateType::RELAY;
                g.inputs = vec![((((k << log_fold) + 0) << 12) + 2048 + j + 256) as u32];
            }
        }
        */
        depth += 1;
    }

    c.depth = depth;
    c
}


pub fn gen_ringsig_witness<T: Field>(inputs: &Vec<[u64; 4]>, log_num_hash: usize, log_fold: usize) -> Vec<T> {
    let mut witness = vec![T::from(0); 1 << (log_num_hash + log_fold + 11)];

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
                witness[(b << (11 + log_fold)) + 256 + w * 64 + i] = T::from(inp & 1);
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

            // ρ and π steps
            let mut b = [0u64; 25];
            for x in 0..=4 {
                for y in 0..=4 {
                    b[y * 5 + (2 * x + 3 * y) % 5] = rot(a[x * 5 + y], ROTATE[x][y]);
                }
            }

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

            if r != 23 && ((r + 1) % (24 >> log_fold)) == 0 {
                let idx = (r + 1) / (24 >> log_fold);
                for w in 0..25 {
                    let mut st = state[w] as usize;
                    for j in 0..64 {
                        witness[(i << (11 + log_fold)) + (idx << 11) + w * 64 + j] = T::from(st & 1);
                        st >>= 1;
                    }
                }
            }

            // if r == 1 && i == 0 {
            //     for w in 0..25 {
            //         println!("{:016x}", state[w]);
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
                    witness[((i + 1) << (11 + log_fold)) + w * 64 + j] = T::from(inp & 1);
                    inp >>= 1;
                }
            }
        } else {
            // roothash
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
