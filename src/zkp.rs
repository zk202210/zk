use std::sync::Arc;
// use std::time::Instant;

use crate::field::Field;
use crate::field::Coset;
use crate::hash::*;
use crate::channel::Channel;
use crate::channel::Tx;
use crate::channel::recv;
use crate::circuit::*;
use crate::util::*;
use crate::statistic::*;
use crate::gkr::Prover as GKRProver;
use crate::gkr::Verifier as GKRVerifier;


pub struct Prover<'a, 'b, T: Field, Ch: Channel<Output = T>> {
    log_len: usize,
    circ: &'a mut Circuit<'b, T>,
    inputs: &'a Vec<T>,
    // outputs: &'a Vec<T>,
    ch: &'a mut Ch,
}

impl<'a, 'b, T: Field, Ch: Channel<Output = T>> Prover<'a, 'b, T, Ch> {
    pub fn new(log_len: usize, circ: &'a mut Circuit<'b, T>, inputs: &'a Vec<T>, ch: &'a mut Ch) -> Self {
        Self {
            log_len: log_len,
            circ: circ,
            inputs: inputs,
            // outputs: outputs,
            ch: ch,
        }
    }

    pub fn prove_plain(&'a mut self, n_ldt_iter: usize, coset: &'a Coset<T>) {
        use crate::plainvpd::Prover as VPDProver;

        let mut prover = VPDProver::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let l_commit = prover.commit(self.inputs);

        self.ch = prover.ch;
        self.ch.send(0, Tx::<T>::LDTCommitHash(l_commit));

        let log_log_len = ceil_log(self.log_len);

        let mu: T;
        let q_scalar: T;
        let mut q_inputs: Vec<T>;
        {
            let mut prover = GKRProver::<T, Ch>::new(self.circ, self.ch);
            prover.main_prove();
            q_scalar = prover.q_scalar;
            if prover.v_ab.len() > 1 {
                mu = prover.v_ab[0] * (T::from(1) - prover.q_scalar) + prover.v_ab[1] * prover.q_scalar;
                q_inputs = Vec::with_capacity(2 << log_log_len);

                for i in 0..(prover.q_inputs.len() / 2) {
                    q_inputs.push(prover.q_inputs[i]);
                }
                for _ in (prover.q_inputs.len() / 2)..(1 << log_log_len) {
                    q_inputs.push(T::from(0));
                }
                let offset = prover.q_inputs.len() / 2;
                for i in 0..(prover.q_inputs.len() / 2) {
                    q_inputs.push(prover.q_inputs[i + offset]);
                }
                for _ in (prover.q_inputs.len() / 2)..(1 << log_log_len) {
                    q_inputs.push(T::from(0));
                }
            } else {
                mu = prover.v_ab[0];
                q_inputs = Vec::with_capacity(2 << log_log_len);
                // q_inputs = Vec::with_capacity(1 << log_log_len);
                q_inputs.append(&mut prover.q_inputs);
                // TODO XXX: one phase q_circ
                // q_inputs.resize(1 << log_log_len, T::from(0));
                q_inputs.resize(2 << log_log_len, T::from(0));
            }
        }
        self.circ.clean();

        // self.ch.send(0, Tx::<T>::Syn);
        // let scalar = recv!(self.ch, Random, 1);

        let mut q_circ = build_circuit(self.log_len, log_log_len, q_scalar);
        prover.ch = self.ch;
        let _ = prover.open(mu, &mut q_circ, &q_inputs);
    }

    pub fn prove_v2(&'a mut self, n_ldt_iter: usize, coset: &'a Coset<T>) {
        use crate::vpdv2::Prover as VPDProver;

        let mut prover = VPDProver::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let l_commit = prover.commit(self.inputs);

        self.ch = prover.ch;
        self.ch.send(0, Tx::<T>::LDTCommitHash(l_commit));

        let log_log_len = ceil_log(self.log_len);

        let mu: T;
        let q_scalar: T;
        let mut q_inputs: Vec<T>;
        {
            let mut prover = GKRProver::<T, Ch>::new(self.circ, self.ch);
            prover.main_prove();
            q_scalar = prover.q_scalar;
            if prover.v_ab.len() > 1 {
                mu = prover.v_ab[0] * (T::from(1) - prover.q_scalar) + prover.v_ab[1] * prover.q_scalar;
                q_inputs = Vec::with_capacity(2 << log_log_len);

                for i in 0..(prover.q_inputs.len() / 2) {
                    q_inputs.push(prover.q_inputs[i]);
                }
                for _ in (prover.q_inputs.len() / 2)..(1 << log_log_len) {
                    q_inputs.push(T::from(0));
                }
                let offset = prover.q_inputs.len() / 2;
                for i in 0..(prover.q_inputs.len() / 2) {
                    q_inputs.push(prover.q_inputs[i + offset]);
                }
                for _ in (prover.q_inputs.len() / 2)..(1 << log_log_len) {
                    q_inputs.push(T::from(0));
                }
            } else {
                mu = prover.v_ab[0];
                q_inputs = Vec::with_capacity(2 << log_log_len);
                q_inputs.append(&mut prover.q_inputs);
                // TODO XXX: one phase q_circ
                // q_inputs.resize(1 << log_log_len, T::from(0));
                q_inputs.resize(2 << log_log_len, T::from(0));
            }
        }
        self.circ.clean();

        let mut q_circ = build_circuit(self.log_len, log_log_len, q_scalar);
        prover.ch = self.ch;
        let _ = prover.open(mu, &mut q_circ, &q_inputs);
    }

    pub fn prove(&'a mut self, n_ldt_iter: usize, coset: &'a Coset<T>) {
        use crate::vpd::Prover as VPDProver;
        // use crate::plainvpd_slow::Prover as VPDProver;

        let mut prover = VPDProver::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let l_commit = prover.commit(self.inputs);

        self.ch = prover.ch;
        self.ch.send(0, Tx::<T>::LDTCommitHash(l_commit));

        let log_log_len = ceil_log(self.log_len);

        let mu: T;
        let q_scalar: T;
        let mut q_inputs: Vec<T>;
        {
            let mut prover = GKRProver::<T, Ch>::new(self.circ, self.ch);
            prover.main_prove();
            q_scalar = prover.q_scalar;
            if prover.v_ab.len() > 1 {
                mu = prover.v_ab[0] * (T::from(1) - prover.q_scalar) + prover.v_ab[1] * prover.q_scalar;
                q_inputs = Vec::with_capacity(2 << log_log_len);

                for i in 0..(prover.q_inputs.len() / 2) {
                    q_inputs.push(prover.q_inputs[i]);
                }
                for _ in (prover.q_inputs.len() / 2)..(1 << log_log_len) {
                    q_inputs.push(T::from(0));
                }
                let offset = prover.q_inputs.len() / 2;
                for i in 0..(prover.q_inputs.len() / 2) {
                    q_inputs.push(prover.q_inputs[i + offset]);
                }
                for _ in (prover.q_inputs.len() / 2)..(1 << log_log_len) {
                    q_inputs.push(T::from(0));
                }
            } else {
                mu = prover.v_ab[0];
                q_inputs = Vec::with_capacity(2 << log_log_len);
                // q_inputs = Vec::with_capacity(1 << log_log_len);
                q_inputs.append(&mut prover.q_inputs);
                // TODO XXX: one phase q_circ
                // q_inputs.resize(1 << log_log_len, T::from(0));
                q_inputs.resize(2 << log_log_len, T::from(0));
            }
        }
        self.circ.clean();

        // self.ch.send(0, Tx::<T>::Syn);
        // let scalar = recv!(self.ch, Random, 1);

        let mut q_circ = build_circuit(self.log_len, log_log_len, q_scalar);
        prover.ch = self.ch;
        let _ = prover.open(mu, &mut q_circ, &q_inputs);
    }

    pub fn prove_mpcvpd(&'a mut self, coset: &'a Coset<T>, n_ldt_iter: usize, n_party: usize, n_mpc_iter: usize) {
        use crate::mpcvpd::Prover as VPDProver;

        let mut prover = VPDProver::<T, Ch>::new(self.ch, coset, n_ldt_iter, n_party, n_mpc_iter);
        let l_commit = prover.commit(self.inputs);

        self.ch = prover.ch;
        self.ch.send(0, Tx::<T>::LDTCommitHash(l_commit));

        let log_log_len = ceil_log(self.log_len);

        let mu: T;
        let q_scalar: T;
        let mut q_inputs: Vec<T>;
        {
            let mut prover = GKRProver::<T, Ch>::new(self.circ, self.ch);
            prover.main_prove();
            q_scalar = prover.q_scalar;
            if prover.v_ab.len() > 1 {
                mu = prover.v_ab[0] * (T::from(1) - prover.q_scalar) + prover.v_ab[1] * prover.q_scalar;
                q_inputs = Vec::with_capacity(2 << log_log_len);

                for i in 0..(prover.q_inputs.len() / 2) {
                    q_inputs.push(prover.q_inputs[i]);
                }
                for _ in (prover.q_inputs.len() / 2)..(1 << log_log_len) {
                    q_inputs.push(T::from(0));
                }
                let offset = prover.q_inputs.len() / 2;
                for i in 0..(prover.q_inputs.len() / 2) {
                    q_inputs.push(prover.q_inputs[i + offset]);
                }
                for _ in (prover.q_inputs.len() / 2)..(1 << log_log_len) {
                    q_inputs.push(T::from(0));
                }
            } else {
                mu = prover.v_ab[0];
                q_inputs = Vec::with_capacity(2 << log_log_len);
                // q_inputs = Vec::with_capacity(1 << log_log_len);
                q_inputs.append(&mut prover.q_inputs);
                // TODO XXX: one phase q_circ
                // q_inputs.resize(1 << log_log_len, T::from(0));
                q_inputs.resize(2 << log_log_len, T::from(0));
            }
        }
        self.circ.clean();

        // self.ch.send(0, Tx::<T>::Syn);
        // let scalar = recv!(self.ch, Random, 1);

        // TODO XXX: one phase q_circ
        let mut q_circ = build_circuit(self.log_len, log_log_len, q_scalar);
        prover.ch = self.ch;
        let _ = prover.open(mu, &mut q_circ, &q_inputs);
    }
}


pub struct Verifier<'a, 'b, T: Field, Ch: Channel<Output = T>> {
    // tbl_a: Vec<T>,
    // tbl_b: Vec<T>,
    circ: &'a mut Circuit<'b, T>,
    log_len: usize,
    ch: &'a mut Ch,
}

impl<'a, 'b, T: Field, Ch: Channel<Output = T>> Verifier<'a, 'b, T, Ch> {
    pub fn new(log_len: usize, circ: &'a mut Circuit<'b, T>, ch: &'a mut Ch) -> Self {
        Self {
            // tbl_a: vec![],
            // tbl_b: vec![],
            circ: circ,
            log_len: log_len,
            ch: ch,
        }
    }
    pub fn verify_plain(&mut self, outputs: &'a Vec<T>, n_party: usize, tape_len: usize, n_ldt_iter: usize, coset: &'a Coset::<T>) -> bool {
        use crate::plainvpd::Verifier as VPDVerifier;

        let l_commit = recv!(self.ch, LDTCommitHash, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);

        let log_log_len = ceil_log(self.log_len);

        let mu: T;
        let q_scalar: T;
        let q_inputs = Vec::with_capacity(2 << log_log_len);
        {
            // TODO XXX q_inputs
            let mut verifier = GKRVerifier::<T, Ch>::new(self.circ, self.ch);
            // let res = verifier.zk_verify(outputs, n_party, tape_len, false);
            let res = verifier.zk_verify(outputs, n_party, tape_len);
            assert!(res);

            if verifier.v_ab.len() > 1 {
                mu = verifier.v_ab[0] * (T::from(1) - verifier.q_scalar) + verifier.v_ab[1] * verifier.q_scalar;
            } else {
                mu = verifier.v_ab[0];
            }
            q_scalar = verifier.q_scalar;
        }
        self.circ.clean();

        // let _ = recv!(self.ch, @Syn, 0);
        // let scalar = T::random();
        // self.ch.send(1, Tx::<T>::Random(scalar.clone()));

        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);

        let mut q_circ = build_circuit(self.log_len, log_log_len, q_scalar);

        let mut verifier = VPDVerifier::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let (res, _) = verifier.verify(mu, l_commit, &mut q_circ, &q_inputs);
        assert!(res);
        true
    }

    pub fn verify_v2(&mut self, outputs: &'a Vec<T>, n_party: usize, tape_len: usize, n_ldt_iter: usize, coset: &'a Coset::<T>) -> bool {
        use crate::vpdv2::Verifier as VPDVerifier;

        let l_commit = recv!(self.ch, LDTCommitHash, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);

        let log_log_len = ceil_log(self.log_len);

        let mu: T;
        let q_scalar: T;
        let q_inputs = Vec::with_capacity(2 << log_log_len);
        {
            // TODO XXX q_inputs
            let mut verifier = GKRVerifier::<T, Ch>::new(self.circ, self.ch);
            // let res = verifier.zk_verify(outputs, n_party, tape_len, false);
            let res = verifier.zk_verify(outputs, n_party, tape_len);
            assert!(res);

            if verifier.v_ab.len() > 1 {
                mu = verifier.v_ab[0] * (T::from(1) - verifier.q_scalar) + verifier.v_ab[1] * verifier.q_scalar;
            } else {
                mu = verifier.v_ab[0];
            }
            q_scalar = verifier.q_scalar;
        }
        self.circ.clean();

        let mut q_circ = build_circuit(self.log_len, log_log_len, q_scalar);

        let mut verifier = VPDVerifier::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let (res, _) = verifier.verify(mu, l_commit, &mut q_circ, &q_inputs);
        assert!(res);
        true
    }

    pub fn verify(&mut self, outputs: &'a Vec<T>, n_party: usize, tape_len: usize, n_ldt_iter: usize, coset: &'a Coset::<T>) -> bool {
        use crate::vpd::Verifier as VPDVerifier;
        // use crate::plainvpd_slow::Verifier as VPDVerifier;

        let l_commit = recv!(self.ch, LDTCommitHash, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);

        let log_log_len = ceil_log(self.log_len);

        let mu: T;
        let q_scalar: T;
        let q_inputs = Vec::with_capacity(2 << log_log_len);
        {
            // TODO XXX q_inputs
            let mut verifier = GKRVerifier::<T, Ch>::new(self.circ, self.ch);
            // let res = verifier.zk_verify(outputs, n_party, tape_len, false);
            let res = verifier.zk_verify(outputs, n_party, tape_len);
            assert!(res);

            if verifier.v_ab.len() > 1 {
                mu = verifier.v_ab[0] * (T::from(1) - verifier.q_scalar) + verifier.v_ab[1] * verifier.q_scalar;
            } else {
                mu = verifier.v_ab[0];
            }
            q_scalar = verifier.q_scalar;
        }
        self.circ.clean();

        // let _ = recv!(self.ch, @Syn, 0);
        // let scalar = T::random();
        // self.ch.send(1, Tx::<T>::Random(scalar.clone()));

        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);

        let mut q_circ = build_circuit(self.log_len, log_log_len, q_scalar);

        let mut verifier = VPDVerifier::<T, Ch>::new(self.ch, coset, n_ldt_iter);
        let (res, _) = verifier.verify(mu, l_commit, &mut q_circ, &q_inputs);
        assert!(res);
        true
    }

    pub fn verify_mpcvpd(&mut self, outputs: &'a Vec<T>, coset: &'a Coset::<T>, n_party: usize, tape_len: usize, n_ldt_iter: usize, n_mpc_iter: usize) -> bool {
        use crate::mpcvpd::Verifier as VPDVerifier;

        let stc = recv!(self.ch, SeedTree, 0);
        self.ch.send(1, Tx::<T>::Syn);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);

        let lr_coeff_icr = recv!(self.ch, InputZK, 0);
        self.ch.send(1, Tx::<T>::Syn);
        vpd_bw_inc(lr_coeff_icr.len() * lr_coeff_icr[0].len() * T::SIZE);

        let l_commit = recv!(self.ch, LDTCommitHash, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);

        let log_log_len = ceil_log(self.log_len);

        let mu: T;
        let q_scalar: T;
        let q_inputs = Vec::with_capacity(2 << log_log_len);
        {
            // TODO XXX q_inputs
            let mut verifier = GKRVerifier::<T, Ch>::new(self.circ, self.ch);
            // let res = verifier.zk_verify(outputs, n_party, tape_len, false);
            let res = verifier.zk_verify(outputs, n_party, tape_len);
            assert!(res);

            if verifier.v_ab.len() > 1 {
                mu = verifier.v_ab[0] * (T::from(1) - verifier.q_scalar) + verifier.v_ab[1] * verifier.q_scalar;
            } else {
                mu = verifier.v_ab[0];
            }
            q_scalar = verifier.q_scalar;
        }
        self.circ.clean();

        // let _ = recv!(self.ch, @Syn, 0);
        // let scalar = T::random();
        // self.ch.send(1, Tx::<T>::Random(scalar.clone()));

        // let log_l = self.log_len * 2 + 1 + log_rho;
        // let coset = Coset::<T>::init(log_l, log_rho);

        let mut q_circ = build_circuit(self.log_len, log_log_len, q_scalar);

        let mut verifier = VPDVerifier::<T, Ch>::new(self.ch, coset, n_ldt_iter, n_party, n_mpc_iter);
        let (res, _) = verifier.verify(mu, l_commit, stc, &lr_coeff_icr, &mut q_circ, &q_inputs);
        assert!(res);
        true
    }
}


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
        s *= T::from(1) + tau[i + 1] + rho[i];
    }
    // res[1] += tau[0] * s - (T::from(1) - tau[0]) * s;
    // char-2 field opt
    [T::from(0), s, T::from(0), (T::from(1) - tau[0]) * s]
    // [s, (T::from(1) - tau[0]) * s, T::from(0), T::from(0)]
}

fn pred_not_fprime<T: Field>(tau: &Vec<T>, rho: &Vec<T>, _sgm: &Vec<T>) -> [T; 4] {
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
        let tmp = tau[i + 1] * rho[i];
        s *= T::from(1) - tau[i + 1] - rho[i] + tmp + tmp;
    }
    // res[1] += tau[0] * s - (T::from(1) - tau[0]) * s;
    let tmp = tau[0] * s;
    [T::from(0), tmp + tmp - s, T::from(0), s - tmp]
    // [s, (T::from(1) - tau[0]) * s, T::from(0), T::from(0)]
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
        s *= (T::from(1) - rho[i] - tau[i]) * (T::from(1) - sgm[i] - tau[i + dgr]);
    }

    [s, T::from(0), T::from(0), T::from(0)]
}

fn pred_comb_fprime<T: Field>(grp: usize, dgr: usize, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
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
        let tmp = rho[i] * tau[i];
        s *= T::from(1) - rho[i] - tau[i] + tmp + tmp;
        let tmp = sgm[i] * tau[i + dgr];
        s *= T::from(1) - sgm[i] - tau[i + dgr] + tmp + tmp;
    }

    [s, T::from(0), T::from(0), T::from(0)]
}

fn pred_comb_scale_fprime<T: Field>(grp: usize, dga: usize, dgb: usize, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
    // tau[g | b | a], rho[g | 0 | a], sgm[g | 1 |0..0 b]
    let mut s = T::from(1);

    // g in tau, rho, sgm
    for i in 0..grp {
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
        let tmp = rho[i] * tau[i];
        s *= T::from(1) - rho[i] - tau[i] + tmp + tmp;
    }

    for i in 0..dgb {
        let tmp = sgm[i] * tau[i + dga];
        s *= T::from(1) - sgm[i] - tau[i + dga] + tmp + tmp;
    }
    for i in dgb..dga {
        s *= T::from(1) - sgm[i];
    }

    [s, T::from(0), T::from(0), T::from(0)]
}

fn pred_comb_scale<T: Field>(grp: usize, dga: usize, dgb: usize, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
    // tau[g | b | a], rho[g | 0 | a], sgm[g | 1 |0..0 b]
    let mut s = T::from(1);

    // g in tau, rho, sgm
    for i in 0..grp {
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
        s *= T::from(1) + rho[i] + tau[i];
    }

    for i in 0..dgb {
        s *= T::from(1) + sgm[i] + tau[i + dga];
    }
    for i in dgb..dga {
        s *= T::from(1) - sgm[i];
    }

    [s, T::from(0), T::from(0), T::from(0)]
}

fn pred_linear_comb<T: Field>(log_len: usize, scalar: T, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
    let mut s = T::from(1);
    s *= sgm[log_len] * (T::from(1) - rho[log_len]);
    for i in 0..log_len {
        s *= tau[i] * rho[i] * sgm[i] +
            (T::from(1) - tau[i]) *
            (T::from(1) - rho[i]) *
            (T::from(1) - sgm[i]);
    }
    [T::from(0), (T::from(1) - scalar) * s, scalar * s, T::from(0)]
}

// XXX mu, nu
pub fn build_circuit<'a, T: Field>(log_len: usize, log_log_len: usize, scalar: T) -> Circuit<'a, T> {
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
        if T::CHAR == 2 {
            c.layers[depth].predicate.push(Arc::new(pred_not));
        } else {
            c.layers[depth].predicate.push(Arc::new(pred_not_fprime));
        }
        depth += 1;
    }

    // combinations
    for d in 0..(log_log_len - 1) {
        c.layers.push(Layer::<T>::default());
        // 2^2^0 * 2^k
        // -> 2^2^1 * 2^(k-1)   (d = 0)
        // -> 2^2^2 * 2^(k-2)   (d = 1)
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
        if T::CHAR == 2 {
            c.layers[depth].predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                    pred_comb(log_log_len - d, 1 << d, tau, rho, sgm)));
        } else {
            c.layers[depth].predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                    pred_comb_fprime(log_log_len - d, 1 << d, tau, rho, sgm)));
        }
        depth += 1;
    }

    // combination and scaling halves
    // XXX
    {
        let d = log_log_len - 1;
        c.layers.push(Layer::<T>::default());
        // c.layers[depth].bit_len = log_log_len - (d + 1) + (1 << (d + 1));
        c.layers[depth].bit_len = 1 + log_len;

        let mut idx = 0;
        let mut idx0 = 0;
        {
            for k in 0..(1 << (log_len - (1 << d))) {
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
        {
            for k in 0..(1 << (log_len - (1 << d))) {
                for i in 0..(1 << (1 << d)) {
                    c.layers[depth].gates.push(Gate::default());
                    c.layers[depth].values.push(T::from(0));
                    c.layers[depth].gates[idx].gtype = GateType::MUL;
                    c.layers[depth].gates[idx].inputs.push(idx0 + i);
                    c.layers[depth].gates[idx].inputs.push(idx0 + (1 << (1 << d)) + k);
                    c.layers[depth].gates[idx].utils.push(0);
                    // c.layers[depth].gates[idx].utils[0] = scalar;
                    idx += 1;
                }
            }
        }
        // assert_eq!(idx, 1 << c.layers[depth].bit_len);
        if T::CHAR == 2 {
            c.layers[depth].predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                    pred_comb_scale(log_log_len - d, 1 << d, log_len - (1 << d), tau, rho, sgm)));
        } else {
            c.layers[depth].predicate.push(Arc::new(
                    move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                    pred_comb_scale_fprime(log_log_len - d, 1 << d, log_len - (1 << d), tau, rho, sgm)));
        }
        depth += 1;
    }

    // Linear combination
    {
        c.layers.push(Layer::<T>::default());
        // c.layers[depth].bit_len = log_log_len - (d + 1) + (1 << (d + 1));
        c.layers[depth].bit_len = log_len;
        c.layers[depth].utils = vec![T::from(1) - scalar, scalar];

        {
            for i in 0..(1 << log_len) {
                c.layers[depth].gates.push(Gate::default());
                c.layers[depth].values.push(T::from(0));
                c.layers[depth].gates[i].gtype = GateType::LINEAR;
                c.layers[depth].gates[i].inputs.push(i as u32);
                c.layers[depth].gates[i].inputs.push((i + (1 << log_len)) as u32);
                c.layers[depth].gates[i].utils = vec![0, 1];
            }
        }

        c.layers[depth].predicate.push(Arc::new(
                move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                pred_linear_comb(log_len, scalar, tau, rho, sgm)));
        depth += 1;
    }

    c.depth = depth;
    c
}
