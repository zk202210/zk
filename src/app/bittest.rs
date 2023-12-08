use std::sync::Arc;

use crate::field::*;
use crate::circuit::*;
use crate::gkr::*;
use crate::util::*;

// sha3-256 merkle tree without padding
// leaf inputs are 512-bit each
pub fn build_circuit<'a, T: Field>(log_size: usize) -> Circuit<'a, T> {
    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
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
        // 1600 states => 2048
        layer.bit_len = log_size;
        // first round
        //
        layer.utils = vec![T::from(1)];
        layer.num_phase = 1;
        layer.num_table = vec![2];
        layer.val_degree = vec![2];
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_d2_no_const);
        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_size]; 2];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in 0..(1 << log_size) {
                    tbl.tbl_coeff[1][k] += tbl_prev_p[k];
                    tbl.tbl_coeff[0][k] -= tbl_prev_p[k];
                }
            });
        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_size];
                for k in 0..(1 << log_size) {
                    res[k] = inputs[k] * (inputs[k] - T::from(1));
                }
                res
            });
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                let mut relay = T::from(1);
                for i in 0..log_size {
                    relay *= eq2(tau[i], rho[i]);
                }
                [relay, relay, T::from(0), T::from(0)]
            }));
        depth += 1;
    }
    c.depth = depth;
    c
}
