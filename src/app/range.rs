use std::sync::Arc;

use crate::field::*;
use crate::circuit::*;
use crate::gkr::*;

pub fn build_circuit<'a, T: Field>(log_range: usize, fast_mode: bool) -> Circuit<'a, T> {
    if fast_mode {
        build_circuit_fast::<T>(log_range)
    } else {
        build_circuit_slow::<T>(log_range)
    }
}

fn build_circuit_fast<'a, T: Field>(log_range: usize) -> Circuit<'a, T> {
    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);

    let mut vanishing = Vec::with_capacity(log_range);
    for i in 0..log_range {
        vanishing.push(T::from(1 << i));
    }

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = 0;
        layer.bit_len = log_size;
        // TODO: XXX useless gates at input layer
        layer.gates.push(Gate::default());
        layer.gates[0].gtype = GateType::IN;
        layer.gates[0].inputs.push(0);
        layer.values.push(T::from(0));
        depth += 1;
    }

    for d in 0..log_range {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = 0;
        layer.bit_len = log_size;
        layer.utils.push(vanishing[d]);
        for i in (d + 1)..log_range {
            vanishing[i] = vanishing[i] * (vanishing[i] - vanishing[d]);
        }
        layer.gates.push(Gate::default());
        layer.gates[0].gtype = GateType::P1D2C0;
        layer.gates[0].inputs.push(0);
        layer.gates[0].utils.push(0);
        layer.values.push(T::from(0));
        layer.num_phase = 1;
        layer.num_table = vec![2];
        layer.val_degree = vec![2];
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_d2_no_const);
        depth += 1;
    }
    c.depth = depth;
    c
}

fn build_circuit_slow<'a, T: Field>(log_range: usize) -> Circuit<'a, T> {
    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = 0;
        layer.bit_len = log_size;
        // TODO: XXX useless gates at input layer
        layer.gates.push(Gate::default());
        layer.gates[0].gtype = GateType::IN;
        layer.gates[0].inputs.push(0);
        layer.values.push(T::from(0));
        depth += 1;
    }

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_range - 1;
        layer.bit_len = log_size;
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];
        layer.utils = Vec::with_capacity(1 << log_range);
        for i in 0..(1 << log_range) {
            layer.utils.push(T::from(i));
        }
        for i in 0..(1 << log_size) {
            let g = &mut layer.gates[i];
            g.gtype = GateType::ADDMUL;
            let i = i as u32;
            g.inputs = vec![0, 0];
            g.utils = vec![2 * i, 2 * i + 1];
        }
        depth += 1;
    }
    for d in 0..(log_range - 1) {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        let log_size = log_range - 2 - d;
        layer.bit_len = log_size;
        layer.gates = vec![Gate::default(); 1 << log_size];
        layer.values = vec![T::from(0); 1 << log_size];
        for i in 0..(1 << log_size) {
            let g = &mut layer.gates[i];
            g.gtype = GateType::MUL;
            let i = i as u32;
            g.inputs = vec![2 * i, 2 * i + 1];
        }
        depth += 1;
    }
    c.depth = depth;
    c
}
