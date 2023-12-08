use rand::Rng;
use std::sync::Arc;
use std::time::Instant;

use crate::field::Field;
use crate::field::Coset;
use crate::hash::*;
use crate::channel::Channel;
use crate::channel::ChannelType;
use crate::channel::Tx;
use crate::channel::recv;
use crate::tape::Tape;
use crate::poly::*;
use crate::circuit::*;
use crate::util::*;
use crate::vpd::Prover as VPDProver;
use crate::vpd::Verifier as VPDVerifier;
use crate::statistic::*;

fn pred_not<T: Field>(tau: &Vec<T>, rho: &Vec<T>, _sgm: &Vec<T>) -> [T; 4] {
    let mut s = T::from(1);
    // bit-order: high...low
    // sgm == 0
    // tau[u0] is NOT of rho[u]
    // tau[u1] is RELAY of rho[u]
    //
    // for i in 0..sgm.len() {
    //     s *= T::from(1) - sgm[i];
    // }
    for i in 0..rho.len() {
        // s *= (T::from(1) - tau[i + 1]) * (T::from(1) - rho[i]) + tau[i + 1] * rho[i];
        // char-2 field opt
        // s *= T::from(1) + tau[i + 1] + rho[i];
        s *= T::xor(tau[i + 1], rho[i]);
    }

    // res[1] += tau[0] * s - (T::from(1) - tau[0]) * s;
    let res1 = if T::CHAR == 2 { s } else {
        let mul = tau[0] * s;
        mul + mul - s
    };
    [T::from(0), res1, T::from(0), (T::from(1) - tau[0]) * s]
}

fn pred_comb<T: Field>(grp: usize, dgr: usize, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
    // tau[g | b | a], rho[g | 0 | a], sgm[g | 1 | b]
    let mut s = T::from(1);

    // g in tau, rho, sgm
    for i in 0..grp {
        s *=
            (T::from(1) - rho[i + 1 + dgr]) *
            (T::from(1) - sgm[i + 1 + dgr]) *
            (T::from(1) - tau[i + dgr + dgr]) +
            rho[i + 1 + dgr] *
            sgm[i + 1 + dgr] *
            tau[i + dgr + dgr];
    }

    // 0 in rho, 1 in sgm
    s *= (T::from(1) - rho[dgr]) * sgm[dgr];

    // a, b in tau, rho, sgm
    for i in 0..dgr {
        // s *= (T::from(1) - rho[i]) * (T::from(1) - tau[i])
        //     + rho[i] * tau[i];
        // s *= (T::from(1) - sgm[i]) * (T::from(1) - tau[i + dgr])
        //     + sgm[i] * tau[i + dgr];
        // char-2 field opt
        // s *= (T::from(1) + rho[i] + tau[i]) * (T::from(1) + sgm[i] + tau[i + dgr]);
        s *= T::xor(rho[i], tau[i]) * T::xor(sgm[i], tau[i + dgr]);
    }

    [s, T::from(0), T::from(0), T::from(0)]
}

fn pred_comb_scale<T: Field>(grp: usize, dga: usize, dgb: usize, scalar: T, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
    // tau[g | b | a], rho[g | 0 | a], sgm[g | 1 |0..0 b]
    let mut s =
        (T::from(1) - scalar) *
        (T::from(1) - rho[1 + dga]) *
        (T::from(1) - sgm[1 + dga]) *
        (T::from(1) - tau[dga + dgb]) +
        scalar *
        rho[1 + dga] *
        sgm[1 + dga] *
        tau[dga + dgb];

    // g in tau, rho, sgm
    for i in 1..grp {
        s *=
            (T::from(1) - rho[i + 1 + dga]) *
            (T::from(1) - sgm[i + 1 + dga]) *
            (T::from(1) - tau[i + dga + dgb]) +
            rho[i + 1 + dga] *
            sgm[i + 1 + dga] *
            tau[i + dga + dgb];
    }

    // 0 in rho, 1 in sgm
    s *= (T::from(1) - rho[dga]) * sgm[dga];

    // a, b in tau, rho, sgm
    for i in 0..dga {
        // s *= (T::from(1) - rho[i]) * (T::from(1) - tau[i])
        //     + rho[i] * tau[i];
        // s *= (T::from(1) - sgm[i]) * (T::from(1) - tau[i + dgr])
        //     + sgm[i] * tau[i + dgr];
        // char-2 field opt
        // s *= T::from(1) + rho[i] + tau[i];
        s *= T::xor(rho[i], tau[i]);
    }

    for i in 0..dgb {
        s *= T::xor(sgm[i], tau[i + dga]);
    }
    for i in dgb..dga {
        s *= T::from(1) - sgm[i];
    }

    [s, T::from(0), T::from(0), T::from(0)]
}

// fn pred_scalar<T: Field>(scalar: T, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
//     let h = tau.len() - 1;
//     let h_r = rho.len() - 1;
//     let mut s = T::from(1);
//     for i in 0..sgm.len() {
//         s *= T::from(1) - sgm[i];
//     }
//     for i in 0..h {
//         s *= T::from(1) + rho[i] + tau[i];
//     }
//     // high bits (except the highest) of r are zeros
//     for i in h..h_r {
//         s *= T::from(1) - rho[i];
//     }
//     s *=
//         (T::from(1) - scalar) *
//         (T::from(1) - tau[h]) *
//         (T::from(1) - rho[h_r]) +
//         scalar * tau[h] * rho[h_r];
//     [T::from(0), s, T::from(0), T::from(0)]
// }



// fn pred_dup<T: Field>(tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
//     let mut s = T::from(1);
//     for i in 0..sgm.len() {
//         s *= T::from(1) - sgm[i];
//     }
// 
//     for i in 0..rho.len() {
//         s *= T::from(1) + rho[i] + tau[i];
//     }
//     [T::from(0), s, T::from(0), T::from(0)]
// }

fn _pred_fft<T: Field>(len: usize, factors: Vec<T>, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
    // rho[hi | 0 | lo: len], sgm[hi | 1 | lo: len]
    // -> tau[hi | 0 | lo: len] : LINEAR gate w/ 1, factor(hi)
    // -> tau[hi | 1 | lo: len] : LINEAR gate w/ 1, 1 + factor(hi)

    let mut lo = T::from(1);
    for i in 0..len {
        lo *=
            (T::from(1) - rho[i]) *
            (T::from(1) - sgm[i]) *
            (T::from(1) - tau[i]) +
            rho[i] *
            sgm[i] *
            tau[i];
    }
    // rho=0, sgm=1
    lo *= (T::from(1) - rho[len]) * sgm[len];

    // multiplication tree
    let mut lin0 = T::from(1);
    for i in (len + 1)..tau.len() {
        lin0 *=
            (T::from(1) - rho[i]) *
            (T::from(1) - sgm[i]) *
            (T::from(1) - tau[i]) +
            rho[i] * sgm[i] * tau[i];
    }

    let mut evn_hi1 = T::from(0);

    for k in (len + 1)..tau.len() {
        let mut tmp = T::from(1);
        for i in (len + 1)..tau.len() {
            if i != k {
                tmp *=
                    (T::from(1) - rho[i]) *
                    (T::from(1) - sgm[i]) *
                    (T::from(1) - tau[i]) +
                    rho[i] * sgm[i] * tau[i];
            }
        }
        evn_hi1 += factors[k - len - 1] * rho[k] * sgm[k] * tau[k] * tmp;
    }

    [T::from(0), lin0 * lo,
    (evn_hi1 + lin0 * tau[len]) * lo,
    T::from(0)]
}


// TODO non-full layer, assuming dummy gates appended
pub fn build_circuit<'a, T: Field>(
    // len = n_row * n_col
    // full input length is 2 * len
    log_len: usize,
    log_log_len: usize,
    scalar: T)
-> Circuit<'a, T> {

    let mut depth = 0;
    let mut c = Circuit::<T>::new(depth);

    // input layer
    {
        c.layers.push(Layer::<T>::default());
        // two halves of the circuit for A and B
        c.layers[depth].bit_len = log_log_len + 1;
        for i in 0..(2 << log_log_len) {
            c.layers[depth].gates.push(Gate::default());
            c.layers[depth].gates[i].gtype = GateType::IN;
            c.layers[depth].gates[i].inputs.push(i as u32);
            c.layers[depth].values.push(T::from(0));
        }
        depth += 1;
    }

    // x -> (1-x) x
    {
        c.layers.push(Layer::<T>::default());
        c.layers[depth].bit_len = log_log_len + 2;
        for i in 0..(2 << log_log_len) {
            c.layers[depth].gates.push(Gate::default());
            c.layers[depth].gates.push(Gate::default());
            c.layers[depth].gates[2 * i].gtype = GateType::NOT;
            c.layers[depth].gates[2 * i].inputs.push(i as u32);
            c.layers[depth].gates[2 * i + 1].gtype = GateType::RELAY;
            c.layers[depth].gates[2 * i + 1].inputs.push(i as u32);
            c.layers[depth].values.push(T::from(0));
            c.layers[depth].values.push(T::from(0));
            c.layers[depth].num_phase = 1;
        }
        c.layers[depth].predicate.push(Arc::new(pred_not));
        depth += 1;
    }

    // combinations
    for d in 0..(log_log_len - 1) {
        c.layers.push(Layer::<T>::default());
        // 2 * 2^2^0 * 2^k
        // -> 2 * 2^2^1 * 2^(k-1)   (d = 0)
        // -> 2 * 2^2^2 * 2^(k-2)   (d = 1)
        c.layers[depth].bit_len = 1 + log_log_len - (d + 1) + (1 << (d + 1));
        // merge 2*2^d wires with 2*2^d wires
        // in 2^(log_log_len - d) groups for each half of the circuit
        let mut idx = 0;
        let mut idx0 = 0;
        for _ in 0..(2 << (log_log_len - d - 1)) {
            for k in 0..(1 << (1 << d)) {
                for i in 0..(1 << (1 << d)) {
                    c.layers[depth].gates.push(Gate::default());
                    c.layers[depth].values.push(T::from(0));
                    c.layers[depth].gates[idx].gtype = GateType::MUL;
                    c.layers[depth].gates[idx].inputs.push(idx0 + i);
                    c.layers[depth].gates[idx].inputs.push(idx0 + (1 << (1 << d)) + k);
                    idx += 1;
                }
            }
            idx0 += 2 << (1 << d);
        }
        // assert_eq!(idx, 1 << c.layers[depth].bit_len);
        c.layers[depth].predicate.push(Arc::new(
                move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                pred_comb(log_log_len - d, 1 << d, tau, rho, sgm)));
        depth += 1;
    }

    // combination and scaling halves
    {
        let d = log_log_len - 1;
        c.layers.push(Layer::<T>::default());
        // c.layers[depth].bit_len = 1 + log_log_len - (d + 1) + (1 << (d + 1));
        c.layers[depth].bit_len = 1 + log_len;

        c.layers[depth].utils = vec![scalar, T::from(1) - scalar];

        let mut idx = 0;
        let mut idx0 = 0;
        {
            for k in 0..(1 << (log_len - (1 << d))) {
                for i in 0..(1 << (1 << d)) {
                    c.layers[depth].gates.push(Gate::default());
                    c.layers[depth].values.push(T::from(0));
                    c.layers[depth].gates[idx].gtype = GateType::SCAMUL;
                    c.layers[depth].gates[idx].inputs.push(idx0 + i);
                    c.layers[depth].gates[idx].inputs.push(idx0 + (1 << (1 << d)) + k);
                    c.layers[depth].gates[idx].utils.push(1);
                    // c.layers[depth].gates[idx].utils[0] = T::from(1) - scalar;
                    idx += 1;
                }
            }
            idx0 += 2 << (1 << d);
        }
        {
            for k in 0..(1 << (log_len - (1 << d))) {
                for i in 0..(1 << (1 << d)) {
                    c.layers[depth].gates.push(Gate::default());
                    c.layers[depth].values.push(T::from(0));
                    c.layers[depth].gates[idx].gtype = GateType::SCAMUL;
                    c.layers[depth].gates[idx].inputs.push(idx0 + i);
                    c.layers[depth].gates[idx].inputs.push(idx0 + (1 << (1 << d)) + k);
                    c.layers[depth].gates[idx].utils.push(0);
                    // c.layers[depth].gates[idx].utils[0] = scalar;
                    idx += 1;
                }
            }
        }
        // assert_eq!(idx, 1 << c.layers[depth].bit_len);
        c.layers[depth].predicate.push(Arc::new(
                move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                pred_comb_scale(log_log_len - d, 1 << d, log_len - (1 << d), scalar, tau, rho, sgm)));
        depth += 1;
    }
    c.depth = depth;

    // append_ifft_fft(&mut c, log_l, &coset._fft_factors0, &coset._fft_factors1);
    // append_selector(&mut c, &query_pts);
    c
}

pub struct Prover<'a, T: Field, Ch: Channel<Output = T>> {
    tbl_a: Vec<T>,
    tbl_b: Vec<T>,
    log_len: usize,
    inputs: &'a Vec<T>,
    outputs: &'a Vec<T>,
    q_inputs: Vec<T>,
    ch: &'a mut Ch,
}

impl<'a, T: Field, Ch: Channel<Output = T>> Prover<'a, T, Ch> {
    pub fn new(log_len: usize, inputs: &'a Vec<T>, outputs: &'a Vec<T>, ch: &'a mut Ch) -> Self {
        Self {
            tbl_a: vec![],
            tbl_b: vec![],
            log_len: log_len,
            inputs: inputs,
            outputs: outputs,
            q_inputs: vec![],
            ch: ch,
        }
    }

    fn eval_poly_output(&mut self, tau: &Vec<T>) -> T {
        let mut r0 = vec![T::from(0); 1 << (2 * self.log_len)];
        r0[0] = T::from(1);
        for k in 0..(2 * self.log_len) {
            let mask = 1 << k;
            for i in 0..mask {
                r0[i + mask] = r0[i] * tau[k];
                r0[i] = r0[i] - r0[i + mask];
            }
        }

        let mut res = T::from(0);
        for i in 0..(1 << (2 * self.log_len)) {
            res += r0[i] * self.outputs[i];
        }
        res
    }

    fn sumcheck_init(&mut self, tau: &Vec<T>) {
        let size = 1 << (2 * self.log_len);
        self.tbl_a = vec![T::from(0); size / 2];
        self.tbl_b = vec![T::from(0); size / 2];

        for u in 0..(size / 2) {
            let u0 = u << 1;
            let v0 = u0 + size;
            self.tbl_a[u] = self.inputs[u0] + (self.inputs[u0 + 1] - self.inputs[u0]) * tau[0];
            self.tbl_b[u] = self.inputs[v0] + (self.inputs[v0 + 1] - self.inputs[v0]) * tau[self.log_len];
        }

        for i in 1..self.log_len {
            let mask = 1 << (2 * self.log_len - i - 1);
            for u in 0..mask {
                let u0 = u << 1;
                self.tbl_a[u0 + 1] = self.tbl_a[u0 + 1] - self.tbl_a[u0];
                self.tbl_b[u0 + 1] = self.tbl_a[u0 + 1] - self.tbl_b[u0];
                self.tbl_a[u] = self.tbl_a[u0] + self.tbl_a[u0 + 1] * tau[i];
                self.tbl_b[u] = self.tbl_b[u0] + self.tbl_b[u0 + 1] * tau[i + self.log_len];
            }
        }
    }

    #[cfg(not(debug_assertions))]
    fn sumcheck_eval(&mut self, i: usize) -> [T; 2] {
        let mut coef = [T::from(0); 2];
        for u in (0..(1 << (self.log_len - i))).step_by(2) {
            let v = u + 1;
            self.tbl_a[v] = self.tbl_a[v] - self.tbl_a[u];
            self.tbl_b[v] = self.tbl_b[v] - self.tbl_b[u];
            coef[0] += self.tbl_a[u] * self.tbl_b[u];
            coef[1] += self.tbl_a[v] * self.tbl_b[v];
        }
        coef
    }

    #[cfg(debug_assertions)]
    fn sumcheck_eval(&mut self, i: usize) -> [T; 3] {
        let mut coef = [T::from(0); 3];
        for u in (0..(1 << (self.log_len - i))).step_by(2) {
            let v = u + 1;
            self.tbl_a[v] = self.tbl_a[v] - self.tbl_a[u];
            self.tbl_b[v] = self.tbl_b[v] - self.tbl_b[u];
            coef[0] += self.tbl_a[u] * self.tbl_b[u];
            // swap coef[1] with coef[2]
            coef[1] += self.tbl_a[v] * self.tbl_b[v];
            coef[2] += self.tbl_a[u] * self.tbl_b[v] + self.tbl_a[v] * self.tbl_b[u];
        }
        coef
    }

    fn sumcheck_update(&mut self, i: usize, r: T) {
        let mask = 1 << (self.log_len - i - 1);
        for u in 0..mask {
            let u0 = u << 1;
            self.tbl_a[u] = self.tbl_a[u0] + self.tbl_a[u0 + 1] * r;
            self.tbl_b[u] = self.tbl_b[u0] + self.tbl_b[u0 + 1] * r;
        }
    }

    fn sumcheck_final(&mut self) -> (T, T) {
        (self.tbl_a[0], self.tbl_b[0])
    }

    pub fn prove_w_flpcp(&mut self) {
        self.ch.send(0, Tx::<T>::Input(self.inputs.to_vec()));
        let tau: Vec<T> = recv!(self.ch, MultiRandom, 1);

        let mut r = vec![T::from(0); self.log_len];

        self.sumcheck_init(&tau);

        let fr = self.eval_poly_output(&tau);
        self.ch.send(0, Tx::<T>::OLayer(fr));
        let _ = recv!(self.ch, @Syn, 1);

        for i in 0..self.log_len {
            let coef = self.sumcheck_eval(i);
            self.ch.send(0, Tx::<T>::Sumcheck(coef.to_vec()));
            r[i] = recv!(self.ch, Random, 1);
            self.sumcheck_update(i, r[i]);

            // let tmp = fr - coef[0] - coef[0] - coef[1];
            // fr = r[i] * ((r[i] * coef[1]) + tmp) + coef[0];
        }

        let (v_a, v_b) = self.sumcheck_final();

        self.ch.send(0, Tx::<T>::Layer([v_a, v_b].to_vec()));
        let _ = recv!(self.ch, @Syn, 1);

        self.ch.send(0, Tx::<T>::LayerPred([T::from(1), T::from(0), T::from(0), T::from(0)]));
        let _ = recv!(self.ch, @Syn, 1);

        self.ch.send(0, Tx::<T>::FLPCP);
        let _ = recv!(self.ch, @Syn, 1);

        self.ch.send(0, Tx::<T>::MPCitH);
        let _ = recv!(self.ch, @Syn, 1);
    }

    pub fn prove(&mut self) {
        self.ch.send(0, Tx::<T>::Input(self.inputs.to_vec()));
        self.main_prove();
    }

    pub fn main_prove(&mut self) -> [T; 2] {
        let mpc0 = mpc_prv_get();
        let start = Instant::now();

        let tau: Vec<T> = recv!(self.ch, MultiRandom, 1);

        let mut r = vec![T::from(0); self.log_len];

        self.sumcheck_init(&tau);

        let fr = self.eval_poly_output(&tau);
        self.ch.send(0, Tx::<T>::OLayer(fr));
        let _ = recv!(self.ch, @Syn, 1);

        for i in 0..self.log_len {
            let coef = self.sumcheck_eval(i);
            self.ch.send(0, Tx::<T>::Sumcheck(coef.to_vec()));
            r[i] = recv!(self.ch, Random, 1);
            self.sumcheck_update(i, r[i]);

            // let tmp = fr - coef[0] - coef[0] - coef[1];
            // fr = r[i] * ((r[i] * coef[1]) + tmp) + coef[0];
        }

        let log_log_len = ceil_log(self.log_len * 2);
        self.q_inputs = vec![T::from(0); 2 << log_log_len];
        for i in 0..self.log_len {
            self.q_inputs[i] = tau[i];
            self.q_inputs[i + (1 << log_log_len)] = tau[i + self.log_len];
        }
        for i in 0..self.log_len {
            self.q_inputs[i + self.log_len] = r[i];
            self.q_inputs[i + self.log_len + (1 << log_log_len)] = r[i];
        }

        let (v_a, v_b) = self.sumcheck_final();

        self.ch.send(0, Tx::<T>::FinalLayer([v_a, v_b].to_vec()));
        let _ = recv!(self.ch, @Syn, 1);

        self.ch.send(0, Tx::<T>::FinalLayerPred([T::from(1), T::from(0), T::from(0), T::from(0)]));
        let _ = recv!(self.ch, @Syn, 1);

        self.ch.send(0, Tx::<T>::MPCitH0);
        let _ = recv!(self.ch, @Syn, 1);

        // self.ch.send(0, Tx::<T>::LayerPred([T::from(1), T::from(0), T::from(0), T::from(0)]));
        // let _ = recv!(self.ch, @Syn, 1);

        // self.ch.send(0, Tx::<T>::FLPCP);
        // let _ = recv!(self.ch, @Syn, 1);

        self.ch.send(0, Tx::<T>::MPCitH);
        let _ = recv!(self.ch, @Syn, 1);

        let elapsed = start.elapsed();
        let in_ms = elapsed.as_micros() as usize;
        let mpc1 = mpc_prv_get();
        let flp_ms = flp_prv_get();
        gkr_prv_inc(in_ms - (mpc1 - mpc0) - flp_ms);
        [v_a, v_b]
    }

    pub fn prove_w_vpd(&'a mut self, n_ldt_iter: usize, coset: &'a Coset<T>) {
        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);

        let mut prover = VPDProver::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let l_commit = prover.commit(self.inputs);

        self.ch = prover.ch;
        self.ch.send(0, Tx::<T>::LDTCommitHash(l_commit));

        // println!("Prover input commit done.");

        let v = self.main_prove();

        // println!("Prover main proof done.");

        self.ch.send(0, Tx::<T>::Syn);
        let scalar = recv!(self.ch, Random, 1);
        let mu = v[0] * (T::from(1) - scalar) + v[1] * scalar;

        let log_log_len = ceil_log(self.log_len * 2);
        let mut q_circ = build_circuit(self.log_len * 2, log_log_len, scalar);
        prover.ch = self.ch;
        let _ = prover.open(mu, &mut q_circ, &self.q_inputs);
    }

    pub fn prove_w_vpdv2(&'a mut self, n_ldt_iter: usize, coset: &'a Coset<T>) {
        use crate::vpdv2::Prover as VPDProver;
        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);

        let mut prover = VPDProver::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let l_commit = prover.commit(self.inputs);

        self.ch = prover.ch;
        self.ch.send(0, Tx::<T>::LDTCommitHash(l_commit));

        // println!("Prover input commit done.");

        let v = self.main_prove();

        // println!("Prover main proof done.");

        self.ch.send(0, Tx::<T>::Syn);
        let scalar = recv!(self.ch, Random, 1);
        let mu = v[0] * (T::from(1) - scalar) + v[1] * scalar;

        let log_log_len = ceil_log(self.log_len * 2);
        let mut q_circ = build_circuit(self.log_len * 2, log_log_len, scalar);
        prover.ch = self.ch;
        let _ = prover.open(mu, &mut q_circ, &self.q_inputs);
    }



    pub fn prove_w_mpcvpd(&'a mut self, coset: &'a Coset<T>, n_ldt_iter: usize, n_party: usize, n_mpc_iter: usize) {
        use crate::mpcvpd::Prover as VPDProver;
        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);

        let mut prover = VPDProver::<T, Ch>::new(self.ch, coset, n_ldt_iter, n_party, n_mpc_iter);
        let l_commit = prover.commit(self.inputs);

        self.ch = prover.ch;
        self.ch.send(0, Tx::<T>::LDTCommitHash(l_commit));

        // println!("Prover input commit done.");

        let v = self.main_prove();

        // println!("Prover main proof done.");

        self.ch.send(0, Tx::<T>::Syn);
        let scalar = recv!(self.ch, Random, 1);
        let mu = v[0] * (T::from(1) - scalar) + v[1] * scalar;

        let log_log_len = ceil_log(self.log_len * 2);
        let mut q_circ = build_circuit(self.log_len * 2, log_log_len, scalar);
        prover.ch = self.ch;
        let _ = prover.open(mu, &mut q_circ, &self.q_inputs);
    }
}

pub struct Verifier<'a, T: Field, Ch: Channel<Output = T>> {
    tbl_a: Vec<T>,
    tbl_b: Vec<T>,
    log_len: usize,
    q_inputs: Vec<T>,
    ch: &'a mut Ch,
}

impl<'a, T: Field, Ch: Channel<Output = T>> Verifier<'a, T, Ch> {
    pub fn new(log_len: usize, ch: &'a mut Ch) -> Self {
        Self {
            tbl_a: vec![],
            tbl_b: vec![],
            log_len: log_len,
            q_inputs: vec![],
            ch: ch,
        }
    }

    // fn eval_poly_output(&mut self, output: &'a Vec<T>, tau: &Vec<T>) -> T {
    //     let mut r0 = vec![T::from(0); 1 << (2 * self.log_len)];
    //     r0[0] = T::from(1);
    //     for k in 0..(2 * self.log_len) {
    //         let mask = 1 << k;
    //         for i in 0..mask {
    //             r0[i + mask] = r0[i] * tau[k];
    //             r0[i] = r0[i] - r0[i + mask];
    //         }
    //     }

    //     let mut res = T::from(0);
    //     for i in 0..(1 << (2 * self.log_len)) {
    //         res += r0[i] * output[i];
    //     }
    //     res
    // }

    fn eval_poly_output(&mut self, outputs: &'a Vec<T>, tau: &Vec<T>) -> T {
        self.sumcheck_init(&tau);

        let mask = (1 << self.log_len) - 1;
        let mut res = T::from(0);
        for i in 0..(1 << (2 * self.log_len)) {
            res += outputs[i] * self.tbl_a[i & mask] * self.tbl_b[i >> self.log_len];
        }
        res
    }

    fn sumcheck_init(&mut self, tau: &Vec<T>) {
        self.tbl_a = vec![T::from(0); 1 << (2 * self.log_len)];
        self.tbl_b = vec![T::from(0); 1 << (2 * self.log_len)];
        self.tbl_a[0] = T::from(1);
        self.tbl_b[0] = T::from(1);

        for k in 0..self.log_len {
            let mask = 1 << k;
            for i in 0..mask {
                self.tbl_a[i + mask] = self.tbl_a[i] * tau[k];
                self.tbl_a[i] = self.tbl_a[i] - self.tbl_a[i + mask];
                self.tbl_b[i + mask] = self.tbl_b[i] * tau[k + self.log_len];
                self.tbl_b[i] = self.tbl_b[i] - self.tbl_b[i + mask];
            }
        }
    }

    fn sumcheck_final(&mut self, r: &Vec<T>) {
        for k in 0..self.log_len {
            let mask = 1 << (k + self.log_len);
            for i in 0..mask {
                self.tbl_a[i + mask] = self.tbl_a[i] * r[k];
                self.tbl_a[i] = self.tbl_a[i] - self.tbl_a[i + mask];
                self.tbl_b[i + mask] = self.tbl_b[i] * r[k];
                self.tbl_b[i] = self.tbl_b[i] - self.tbl_b[i + mask];
            }
        }
    }

    fn input_test(&mut self, v_ab: &Vec<T>, input: &Vec<T>, r: &Vec<T>) -> bool {
        for k in 0..self.log_len {
            let mask = 1 << (k + self.log_len);
            for i in 0..mask {
                self.tbl_a[i + mask] = self.tbl_a[i] * r[k];
                self.tbl_a[i] = self.tbl_a[i] - self.tbl_a[i + mask];
                self.tbl_b[i + mask] = self.tbl_b[i] * r[k];
                self.tbl_b[i] = self.tbl_b[i] - self.tbl_b[i + mask];
            }
        }
        let mut v_a_test = v_ab[0];
        let mut v_b_test = v_ab[1];
        for i in 0..(1 << (2 * self.log_len)) {
            v_a_test -= self.tbl_a[i] * input[i];
            v_b_test -= self.tbl_b[i] * input[i + (1 << (2 * self.log_len))];
        }
        v_a_test == T::from(0) && v_b_test == T::from(0)
    }

    pub fn verify(&mut self, outputs: &'a Vec<T>) -> bool {
        let input: Vec<T> = recv!(self.ch, Input, 0);
        let mut tau = vec![T::from(0); 2 * self.log_len];
        for i in 0..self.log_len { tau[i] = T::random(); }
        self.ch.send(1, Tx::<T>::MultiRandom(tau.clone()));

        let mut fr = self.eval_poly_output(&outputs, &tau);
        let _fr: T = recv!(self.ch, OLayer, 0);
        assert_eq!(fr, _fr);
        self.ch.send(1, Tx::<T>::Syn);

        let mut r = vec![T::from(0); self.log_len];

        for i in 0..self.log_len {
            let coef: Vec<T> = recv!(self.ch, Sumcheck, 0);
            r[i] = T::random();
            self.ch.send(1, Tx::<T>::Random(r[i]));
            let tmp = fr - coef[0] - coef[0] - coef[1];
            #[cfg(debug_assertions)] {
                assert_eq!(coef[2], tmp);
                fr = r[i] * ((r[i] * coef[1]) + tmp) + coef[0];
            }
            #[cfg(not(debug_assertions))] {
                fr = r[i] * ((r[i] * coef[1]) + tmp) + coef[0];
            }
        }
        let v_ab: Vec<T> = recv!(self.ch, FinalLayer, 0);
        self.ch.send(1, Tx::<T>::Syn);

        // let _pred: [T; 4] = recv!(self.ch, LayerPred, 0);
        // self.ch.send(1, Tx::<T>::Syn);

        if fr != v_ab[0] * v_ab[1] {
            println!("fail layer mul test");
            return false;
        }

        let _ = recv!(self.ch, FinalLayerPred, 0);
        self.ch.send(1, Tx::<T>::Syn);

        let _ = recv!(self.ch, @MPCitH0, 0);
        self.ch.send(1, Tx::<T>::Syn);

        // let _ = recv!(self.ch, @FLPCP, 0);
        // self.ch.send(1, Tx::<T>::Syn);

        let _ = recv!(self.ch, @MPCitH, 0);
        self.ch.send(1, Tx::<T>::Syn);

        self.input_test(&v_ab, &input, &r)
    }

    pub fn zk_verify_w_flpcp(&mut self, outputs: &'a Vec<T>, n_party: usize, tape_len: usize) -> bool {
        let mut coef_icr: Vec<Vec<T>> = vec![];

        let _input_icr: Vec<T> = recv!(self.ch, Input, 0);
        let mut tau = vec![T::from(0); 2 * self.log_len];
        for i in 0..self.log_len { tau[i] = T::random(); }
        self.ch.send(1, Tx::<T>::MultiRandom(tau.clone()));

        let fr = self.eval_poly_output(&outputs, &tau);
        // let _fr: T = recv!(self.ch, OLayer, 0);
        // assert_eq!(fr, _fr);
        // self.ch.send(1, Tx::<T>::Syn);

        let mut r = vec![T::from(0); self.log_len];

        for i in 0..self.log_len {
            let _coef_icr = recv!(self.ch, Sumcheck, 0);
            coef_icr.push(_coef_icr);
            r[i] = T::random();
            self.ch.send(1, Tx::<T>::Random(r[i]));
            // let tmp = fr - coef[0] - coef[0] - coef[1];
            // #[cfg(debug_assertions)] {
            //     assert_eq!(coef[2], tmp);
            //     fr = r[i] * ((r[i] * coef[2]) + coef[1]) + coef[0];
            // }
            // #[cfg(not(debug_assertions))] {
            //     fr = r[i] * ((r[i] * coef[1]) + tmp) + coef[0];
            // }
        }
        let v_ab_icr = recv!(self.ch, Layer, 0);
        self.ch.send(1, Tx::<T>::Syn);

        let poly_icr = recv!(self.ch, FLPCPCommit, 0);

        let pr = T::random();
        self.ch.send(1, Tx::<T>::Random(pr));


        let (h, mut a, mut b, mut c) = recv!(self.ch, FLPCPOpen, 0);

        let missing_party: u8 = rand::thread_rng().gen::<u8>() % (n_party as u8);
        self.ch.send(1, Tx::<T>::Party(vec![missing_party as usize]));

        let (stc, mut sto) = recv!(self.ch, Seed, 0);

        // begin verify
        assert!(sto.verify(&stc));
        let mut tapes: Vec<Tape> = vec![];
        for p in 0..n_party {
            tapes.push(Tape::new(sto.get_seed(0, p), tape_len * T::SIZE));
        }


        // TODO: transpose shares for SIMD
        // TODO: prove inputs
        let mut input_shares: Vec<Vec<T>> = vec![vec![]; n_party];

        for p in 0..n_party {
            for _ in 0..(1 << (2 * self.log_len + 1)) {
                let share = tapes[p].get::<T>();
                input_shares[p].push(share);
                // if p == 0 {
                //     inputs_shares[0][i] += input_icr[i];
                // }
            }
        }

        self.sumcheck_final(&r);

        let mut v_a_shares: Vec<Vec<T>> = vec![vec![T::from(0); 2]; n_party];
        let mut v_b_shares: Vec<Vec<T>> = vec![vec![T::from(0); 2]; n_party];
        let mut fr_shares:  Vec<Vec<T>> = vec![vec![T::from(0); 2]; n_party];
        let mut v_a_poly_shares: Vec<Vec<T>> = vec![];
        let mut v_b_poly_shares: Vec<Vec<T>> = vec![];
        let mut fr_poly_shares: Vec<Vec<T>> = vec![];

        fr_shares[0][0] = fr;

        for i in 0..self.log_len {
            for p in 0..n_party {
                let mut coef0_share = tapes[p].get::<T>();
                let mut coef2_share = tapes[p].get::<T>();
                if p == 0 {
                    coef0_share += coef_icr[i][0];
                    coef2_share += coef_icr[i][1];
                }
                let coef1_share = fr_shares[p][0] - coef0_share - coef0_share - coef2_share;
                fr_shares[p][0] =
                    (coef2_share
                     * r[i] + coef1_share)
                    * r[i] + coef0_share;
            }
        }

        for p in 0..n_party {
            v_a_shares[p][0] = tapes[p].get::<T>();
            v_b_shares[p][0] = tapes[p].get::<T>();
            if p == 0 {
                v_a_shares[p][0] += v_ab_icr[0];
                v_b_shares[p][0] += v_ab_icr[1];
            }
        }


        let d = 1;
        let x_values_to_m = first_n_field_elements::<T>(d + 1);
        let precompute_to_m = precompute_lagrange_polynomials::<T>(&x_values_to_m);
        let x_values_to_2m = first_n_field_elements::<T>(2 * d + 1);
        let precompute_to_2m = precompute_lagrange_polynomials::<T>(&x_values_to_2m);

        let mut p_shares: Vec<Vec<T>> = vec![vec![T::from(0); 2 * d + 1]; n_party];

        for p in 0..n_party {
            v_a_shares[p][1] = tapes[p].get::<T>();
            v_b_shares[p][1] = tapes[p].get::<T>();
            fr_shares[p][1] = tapes[p].get::<T>();

            v_a_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &v_a_shares[p]));
            v_b_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &v_b_shares[p]));
            fr_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &fr_shares[p]));
            p_shares[p][1] = tapes[p].get::<T>();
            p_shares[p][2] = tapes[p].get::<T>();
            if p == 0 {
                p_shares[p][1] += poly_icr[0];
                p_shares[p][2] += poly_icr[1];
            }
        }


        let mut lagrange_poly_m_eval_r = vec![T::from(0); d + 1];
        let mut lagrange_poly_2m_eval_r = vec![T::from(0); 2 * d + 1];
        for k in 0..(d + 1) {
            lagrange_poly_m_eval_r[k] = poly_eval(&precompute_to_m[k], pr);
        }
        for k in 0..(2 * d + 1) {
            lagrange_poly_2m_eval_r[k] = poly_eval(&precompute_to_2m[k], pr);
        }

        let mut a_shares: Vec<T> = vec![T::from(0); n_party];
        let mut b_shares: Vec<T> = vec![T::from(0); n_party];
        let mut c_shares: Vec<T> = vec![T::from(0); n_party];
        let mut d_shares: Vec<T> = vec![T::from(0); n_party];
        let mut d = a * b - c;
        for p in 0..n_party {
            if p != missing_party as usize {
                a_shares[p] = dot_product(&lagrange_poly_m_eval_r, &v_a_shares[p]);
                b_shares[p] = dot_product(&lagrange_poly_m_eval_r, &v_b_shares[p]);
                c_shares[p] = dot_product(&lagrange_poly_m_eval_r, &fr_shares[p]);
                d_shares[p] = dot_product(&lagrange_poly_2m_eval_r, &p_shares[p]);
            }
            a -= a_shares[p];
            b -= b_shares[p];
            c -= c_shares[p];
            d -= d_shares[p];
        }
        a_shares[missing_party as usize] = a;
        b_shares[missing_party as usize] = b;
        c_shares[missing_party as usize] = c;
        d_shares[missing_party as usize] = d;

        a_shares.append(&mut b_shares);
        a_shares.append(&mut c_shares);
        a_shares.append(&mut d_shares);

        let mut h_ = SHA256::from(0);
        h_.commit(&a_shares);

        h_ == h
    }

    pub fn zk_verify_wo_inputs(&mut self, output: &'a Vec<T>, n_party: usize, tape_len: usize, input_commit: bool) -> (bool, [T; 2]) {
        if self.ch.get_type() == ChannelType::MPCitH {
            self.zk_verify_wo_inputs_mith(output, n_party, tape_len, input_commit)
        } else if self.ch.get_type() == ChannelType::VOLEitH {
            self.zk_verify_wo_inputs_mith(output, n_party, tape_len, input_commit)
            // self.zk_verify_wo_inputs_vith(output, n_party, tape_len, input_commit)
        } else {
            unreachable!()
        }
    }

    /*
    pub fn zk_verify_wo_inputs_vith(&mut self, output: &'a Vec<T>, n_party: usize, tape_len: usize, input_commit: bool) -> (bool, [T; 2]) {
        assert!(!input_commit);

        let n = if n_party == 256 { 32 } else if n_party == 16 { 64 } else { unreachable!() };
        let m = n / 2;
        let log_n = ceil_log(n);
        let log_m = ceil_log(m);

        let mut res = true;

        let mut coef_icr: Vec<Vec<T>> = vec![];

        let mut tau = vec![T::from(0); 2 * self.log_len];
        for i in 0..self.log_len { tau[i] = T::random(); }
        self.ch.send(1, Tx::<T>::MultiRandom(tau.clone()));

        let fr = self.eval_poly_output(&output, &tau);

        let mut r = vec![T::from(0); self.log_len];

        for i in 0..self.log_len {
            let _coef_icr = recv!(self.ch, Sumcheck, 0);
            mpc_bw_inc(T::SIZE * _coef_icr.len());
            coef_icr.push(_coef_icr);
            r[i] = T::random();
            self.ch.send(1, Tx::<T>::Random(r[i]));
        }

        let log_log_len = ceil_log(self.log_len * 2);
        self.q_inputs = vec![T::from(0); 2 << log_log_len];
        for i in 0..self.log_len {
            self.q_inputs[i] = tau[i];
            self.q_inputs[i + (1 << log_log_len)] = tau[i + self.log_len];
        }
        for i in 0..self.log_len {
            self.q_inputs[i + self.log_len] = r[i];
            self.q_inputs[i + self.log_len + (1 << log_log_len)] = r[i];
        }

        let v_ab = recv!(self.ch, FinalLayer, 0);
        mpc_bw_inc(T::SIZE * v_ab.len());
        self.ch.send(1, Tx::<T>::Syn);

        let h = recv!(self.ch, FinalLayerOpen, 0);
        mpc_bw_inc(SHA256::SIZE);
        self.ch.send(1, Tx::<T>::Syn);

        let _ = recv!(self.ch, @MPCitH, 0);
        let alpha = T::random();
        self.ch.send(1, Tx::<T>::Random(alpha));

        let (vith_c, vith_s, hv) = recv!(self.ch, VOLEOpen, 0);
        mpc_bw_inc(SHA256::SIZE + T::SIZE * (vith_c.len() + vith_s.len()));

        let mut missings: Vec<usize> = vec![0; n];
        for i in 0..n {
            missings[i] = rand::thread_rng().gen::<usize>() % n_party;
        }
        self.ch.send(1, Tx::<T>::Party(missings.clone()));

        let (stc, mut sto) = recv!(self.ch, Seed, 0);
        mpc_bw_inc(sto.proof_size);
        self.ch.send(1, Tx::<T>::Syn);

        // begin verify
        // TODO XXX
        let start = Instant::now();

    }
    */

    pub fn zk_verify_wo_inputs_mith(&mut self, output: &'a Vec<T>, n_party: usize, tape_len: usize, input_commit: bool) -> (bool, [T; 2]) {
        let mut coef_icr: Vec<Vec<T>> = vec![];

        let mut tau = vec![T::from(0); 2 * self.log_len];
        for i in 0..self.log_len { tau[i] = T::random(); }
        self.ch.send(1, Tx::<T>::MultiRandom(tau.clone()));

        let fr = self.eval_poly_output(&output, &tau);
        // let _fr: T = recv!(self.ch, OLayer, 0);
        // assert_eq!(fr, _fr);
        // self.ch.send(1, Tx::<T>::Syn);

        let mut r = vec![T::from(0); self.log_len];

        for i in 0..self.log_len {
            let _coef_icr = recv!(self.ch, Sumcheck, 0);
            mpc_bw_inc(T::SIZE * _coef_icr.len());
            coef_icr.push(_coef_icr);
            r[i] = T::random();
            self.ch.send(1, Tx::<T>::Random(r[i]));
            // let tmp = fr - coef[0] - coef[0] - coef[1];
            // #[cfg(debug_assertions)] {
            //     assert_eq!(coef[2], tmp);
            //     fr = r[i] * ((r[i] * coef[2]) + coef[1]) + coef[0];
            // }
            // #[cfg(not(debug_assertions))] {
            //     fr = r[i] * ((r[i] * coef[1]) + tmp) + coef[0];
            // }
        }

        let log_log_len = ceil_log(self.log_len * 2);
        self.q_inputs = vec![T::from(0); 2 << log_log_len];
        for i in 0..self.log_len {
            self.q_inputs[i] = tau[i];
            self.q_inputs[i + (1 << log_log_len)] = tau[i + self.log_len];
        }
        for i in 0..self.log_len {
            self.q_inputs[i + self.log_len] = r[i];
            self.q_inputs[i + self.log_len + (1 << log_log_len)] = r[i];
        }

        let v_ab = recv!(self.ch, FinalLayer, 0);
        mpc_bw_inc(T::SIZE * 2);
        self.ch.send(1, Tx::<T>::Syn);

        let h = recv!(self.ch, FinalLayerOpen, 0);
        mpc_bw_inc(SHA256::SIZE);
        mpc_bw_hash_inc(1);
        self.ch.send(1, Tx::<T>::Syn);

        // let poly_icr = recv!(self.ch, FLPCPCommit, 0);

        // let mut pr = T::random();
        // self.ch.send(1, Tx::<T>::Random(pr));


        // let (h, mut a, mut b, mut c) = recv!(self.ch, FLPCPOpen, 0);


        let missing_party: u8 = rand::thread_rng().gen::<u8>() % (n_party as u8);
        let _ = recv!(self.ch, @MPCitH, 0);
        self.ch.send(1, Tx::<T>::Party(vec![missing_party as usize]));

        let (stc, mut sto) = recv!(self.ch, Seed, 0);
        // TODO: seed without missing party and seed tree
        mpc_bw_inc(sto.proof_size);
        self.ch.send(1, Tx::<T>::Syn);

        // begin verify
        let start = Instant::now();

        assert!(sto.verify(&stc));
        let mut tapes: Vec<Tape> = vec![];
        for p in 0..n_party {
            tapes.push(Tape::new(sto.get_seed(0, p), tape_len * T::SIZE));
        }

        // TODO: transpose shares for SIMD
        // TODO: prove inputs
        if input_commit {
            let mut input_shares: Vec<Vec<T>> = vec![vec![]; n_party];
            for p in 0..n_party {
                for _ in 0..(1 << (2 * self.log_len + 1)) {
                    let share = tapes[p].get::<T>();
                    input_shares[p].push(share);
                    // if p == 0 {
                    //     inputs_shares[0][i] += input_icr[i];
                    // }
                }
            }
            self.sumcheck_final(&r);
        }


        // let mut v_a_shares: Vec<Vec<T>> = vec![T::from(0); 2]; n_party];
        // let mut v_b_shares: Vec<Vec<T>> = vec![T::from(0); 2]; n_party];
        let mut fr_shares:  Vec<T> = vec![T::from(0); n_party];
        // let mut v_a_poly_shares: Vec<Vec<T>> = vec![];
        // let mut v_b_poly_shares: Vec<Vec<T>> = vec![];
        // let mut fr_poly_shares: Vec<Vec<T>> = vec![];

        fr_shares[0] = fr;

        for i in 0..self.log_len {
            for p in 0..n_party {
                let mut coef0_share = tapes[p].get::<T>();
                let mut coef2_share = tapes[p].get::<T>();
                if p == 0 {
                    coef0_share += coef_icr[i][0];
                    coef2_share += coef_icr[i][1];
                }
                let coef1_share = fr_shares[p] - coef0_share - coef0_share - coef2_share;
                fr_shares[p] =
                    (coef2_share
                     * r[i] + coef1_share)
                    * r[i] + coef0_share;
            }
        }

        fr_shares[0] -= v_ab[0] * v_ab[1];

        let mut tmp = T::from(0);
        for p in 0..n_party {
            if p != missing_party as usize {
                tmp += fr_shares[p];
            }
        }

        fr_shares[missing_party as usize] = T::from(0) - tmp;


        // for p in 0..n_party {
        //     v_a_shares[p][0] = tapes[p].get::<T>();
        //     v_b_shares[p][0] = tapes[p].get::<T>();
        //     if p == 0 {
        //         v_a_shares[p][0] += v_ab_icr[0];
        //         v_b_shares[p][0] += v_ab_icr[1];
        //     }
        // }

        /*
        let d = 1;
        let x_values_to_m = first_n_field_elements::<T>(d + 1);
        let precompute_to_m = precompute_lagrange_polynomials::<T>(&x_values_to_m);
        let x_values_to_2m = first_n_field_elements::<T>(2 * d + 1);
        let precompute_to_2m = precompute_lagrange_polynomials::<T>(&x_values_to_2m);

        let mut p_shares: Vec<Vec<T>> = vec![vec![T::from(0); 2 * d + 1]; n_party];

        for p in 0..n_party {
            v_a_shares[p][1] = tapes[p].get::<T>();
            v_b_shares[p][1] = tapes[p].get::<T>();
            fr_shares[p][1] = tapes[p].get::<T>();

            v_a_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &v_a_shares[p]));
            v_b_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &v_b_shares[p]));
            fr_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &fr_shares[p]));
            p_shares[p][1] = tapes[p].get::<T>();
            p_shares[p][2] = tapes[p].get::<T>();
            if p == 0 {
                p_shares[p][1] += poly_icr[0];
                p_shares[p][2] += poly_icr[1];
            }
        }


        let mut lagrange_poly_m_eval_r = vec![T::from(0); d + 1];
        let mut lagrange_poly_2m_eval_r = vec![T::from(0); 2 * d + 1];
        for k in 0..(d + 1) {
            lagrange_poly_m_eval_r[k] = poly_eval(&precompute_to_m[k], pr);
        }
        for k in 0..(2 * d + 1) {
            lagrange_poly_2m_eval_r[k] = poly_eval(&precompute_to_2m[k], pr);
        }

        let mut a_shares: Vec<T> = vec![T::from(0); n_party];
        let mut b_shares: Vec<T> = vec![T::from(0); n_party];
        let mut c_shares: Vec<T> = vec![T::from(0); n_party];
        let mut d_shares: Vec<T> = vec![T::from(0); n_party];
        let mut d = a * b - c;
        for p in 0..n_party {
            if p != missing_party as usize {
                a_shares[p] = dot_product(&lagrange_poly_m_eval_r, &v_a_shares[p]);
                b_shares[p] = dot_product(&lagrange_poly_m_eval_r, &v_b_shares[p]);
                c_shares[p] = dot_product(&lagrange_poly_m_eval_r, &fr_shares[p]);
                d_shares[p] = dot_product(&lagrange_poly_2m_eval_r, &p_shares[p]);
            }
            a -= a_shares[p];
            b -= b_shares[p];
            c -= c_shares[p];
            d -= d_shares[p];
        }
        a_shares[missing_party as usize] = a;
        b_shares[missing_party as usize] = b;
        c_shares[missing_party as usize] = c;
        d_shares[missing_party as usize] = d;

        a_shares.append(&mut b_shares);
        a_shares.append(&mut c_shares);
        a_shares.append(&mut d_shares);

        */
        let mut h_ = SHA256::from(0);
        h_.commit(&fr_shares);

        assert!(h_ == h);

        let elapsed = start.elapsed();
        let in_ms = elapsed.as_micros() as usize;
        mpc_vrf_inc(in_ms);

        (true, [v_ab[0], v_ab[1]])
    }

    pub fn zk_verify(&mut self, outputs: &'a Vec<T>, n_party: usize, tape_len: usize) -> bool {
        let _input_icr: Vec<T> = recv!(self.ch, Input, 0);
        let (res, _) = self.zk_verify_wo_inputs(outputs, n_party, tape_len, true);
        res
    }

    pub fn zk_verify_w_vpdv2(&mut self, outputs: &'a Vec<T>, n_party: usize, tape_len: usize, n_ldt_iter: usize, coset: &'a Coset::<T>) -> bool {
        use crate::vpdv2::Verifier as VPDVerifier;

        let l_commit = recv!(self.ch, LDTCommitHash, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);

        let (res, v_ab) = self.zk_verify_wo_inputs(outputs, n_party, tape_len, false);
        assert!(res);

        let _ = recv!(self.ch, @Syn, 0);
        let scalar = T::random();
        self.ch.send(1, Tx::<T>::Random(scalar.clone()));
        let mu = v_ab[0] * (T::from(1) - scalar) + v_ab[1] * scalar;

        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);
        let log_log_len = ceil_log(self.log_len * 2);

        let mut q_circ = build_circuit(self.log_len * 2, log_log_len, scalar);

        let mut verifier = VPDVerifier::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let (res, _) = verifier.verify(mu, l_commit, &mut q_circ, &self.q_inputs);
        assert!(res);
        // assert_eq!(mu, v_ab[0] * (T::from(1) - scalar) + v_ab[1] * scalar);
        true
    }

    pub fn zk_verify_w_vpd(&mut self, outputs: &'a Vec<T>, n_party: usize, tape_len: usize, n_ldt_iter: usize, coset: &'a Coset::<T>) -> bool {
        let l_commit = recv!(self.ch, LDTCommitHash, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);

        let (res, v_ab) = self.zk_verify_wo_inputs(outputs, n_party, tape_len, false);
        assert!(res);

        let _ = recv!(self.ch, @Syn, 0);
        let scalar = T::random();
        self.ch.send(1, Tx::<T>::Random(scalar.clone()));
        let mu = v_ab[0] * (T::from(1) - scalar) + v_ab[1] * scalar;

        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);
        let log_log_len = ceil_log(self.log_len * 2);

        let mut q_circ = build_circuit(self.log_len * 2, log_log_len, scalar);

        let mut verifier = VPDVerifier::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let (res, _) = verifier.verify(mu, l_commit, &mut q_circ, &self.q_inputs);
        assert!(res);
        // assert_eq!(mu, v_ab[0] * (T::from(1) - scalar) + v_ab[1] * scalar);
        true
    }

    pub fn zk_verify_w_mpcvpd(&mut self, outputs: &'a Vec<T>, coset: &'a Coset::<T>, n_party: usize, tape_len: usize, n_ldt_iter: usize, n_mpc_iter: usize) -> bool {
        use crate::mpcvpd::Verifier as VPDVerifier;

        let stc = recv!(self.ch, SeedTree, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);
        self.ch.send(1, Tx::<T>::Syn);

        let lr_coeff_icr = recv!(self.ch, InputZK, 0);
        vpd_bw_inc(lr_coeff_icr.len() * lr_coeff_icr[0].len() * T::SIZE);
        self.ch.send(1, Tx::<T>::Syn);

        let l_commit = recv!(self.ch, LDTCommitHash, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);

        let (res, v_ab) = self.zk_verify_wo_inputs(outputs, n_party, tape_len, false);
        assert!(res);

        let _ = recv!(self.ch, @Syn, 0);
        let scalar = T::random();
        self.ch.send(1, Tx::<T>::Random(scalar.clone()));
        let mu = v_ab[0] * (T::from(1) - scalar) + v_ab[1] * scalar;

        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);
        let log_log_len = ceil_log(self.log_len * 2);

        let mut q_circ = build_circuit(self.log_len * 2, log_log_len, scalar);

        let mut verifier = VPDVerifier::<T, Ch>::new(self.ch, coset, n_ldt_iter, n_party, n_mpc_iter);
        let (res, _) = verifier.verify(mu, l_commit, stc, &lr_coeff_icr, &mut q_circ, &self.q_inputs);
        assert!(res);
        // assert_eq!(mu, v_ab[0] * (T::from(1) - scalar) + v_ab[1] * scalar);
        true
    }
}
