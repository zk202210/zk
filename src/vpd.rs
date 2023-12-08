use std::time::Instant;

use crate::field::Field;
use crate::field::Coset;
use crate::hash::*;
use crate::channel::Channel;
use crate::channel::Tx;
use crate::channel::recv;
use crate::merkle::MT;
use crate::merkle::MTO;
use crate::circuit::*;
use crate::ldt::Prover as LDTProver;
use crate::ldt::Verifier as LDTVerifier;
use crate::gkr::Prover as GKRProver;
use crate::gkr::Verifier as GKRVerifier;
use crate::statistic::*;

// The zkVPD protocol as proposed in Virgo's paper
// with masking and LDT over a linear combination of polynomials.
// Since the degree will be doubled compared with the input polynomial,
// RS code rate will be doubled as well: from 1/r to 1/(r/2)

pub struct Prover<'a, T: Field, Ch: Channel<Output = T>> {
    pub n_ldt_iter: usize,
    pub ch: &'a mut Ch,
    log_len: usize,
    log_l: usize,
    coset: &'a Coset<T>,
    mt_l: MT<T>,
    mt_s: MT<T>,
    mt_h: MT<T>,
    input_ex: Vec<T>,
}

impl<'a, T: Field, Ch: Channel<Output = T>> Prover<'a, T, Ch> {
    pub fn new(ch: &'a mut Ch, coset: &'a Coset<T>, n_ldt_iter: usize) -> Self {
        // let mut log_kappa: usize = 63 - ((kappa as u64).leading_zeros() as usize);
        // if log_kappa != ((kappa as u64).trailing_zeros() as usize) {
        //     log_kappa += 1;
        // }
        Self {
            n_ldt_iter: n_ldt_iter,
            ch: ch,
            log_len: coset.log_len,
            log_l: coset.log_l,
            // log_kappa: log_kappa,
            coset: coset,
            mt_l: MT::<T>::new(coset.log_l - coset.log_layer_ratio, coset.log_layer_ratio),
            mt_s: MT::<T>::new(coset.log_l - coset.log_layer_ratio, coset.log_layer_ratio),
            mt_h: MT::<T>::new(coset.log_l - coset.log_layer_ratio, coset.log_layer_ratio),
            input_ex: vec![],
        }
    }

    pub fn commit(&mut self, inputs: &Vec<T>) -> SHA256 {

let start = Instant::now();
        // each fri iter open a subcoset number of points
        let k = self.n_ldt_iter << self.coset.log_layer_ratio;
        let mut r_coeff: Vec<T> = Vec::with_capacity(1 << self.log_len);
        for _ in 0..k {
            r_coeff.push(T::random());
        }
        r_coeff.resize(1 << self.log_len, T::from(0));
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let r = (self.coset.fft)(self.coset, &r_coeff, self.log_len, self.log_l);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("    ┌──────#fft of mask r(x) {:>2}->{:>2}:        {:9.3} ms", self.log_len, self.log_l, in_ms as f64 / 1000.0);

    #[cfg(feature = "microbench")] let start0 = Instant::now();
        self.input_ex = (self.coset.interpolate)(self.coset, inputs, self.log_len, self.log_l);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("    ├──────#interpolation of l(x) {:>2}->{:>2}:   {:9.3} ms", self.log_len, self.log_l, in_ms as f64 / 1000.0);

        for i in 0..self.input_ex.len() {
            self.input_ex[i] = self.input_ex[i] + r[i] * self.coset.vanishing_eval[i];
        }

    #[cfg(feature = "microbench")] let start0 = Instant::now();
        (self.coset.rscode_commit)(self.coset, &mut self.mt_l, &self.input_ex, 0);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("    ├──────#merkle tree commit of l'(x) {:>2}: {:9.3} ms", self.log_l, in_ms as f64 / 1000.0);

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
#[cfg(feature = "microbench")]
println!("  ┏━┷━*vpd commit in {} ms", in_ms as f64 / 1000.0);
#[cfg(not(feature = "microbench"))]
println!("  ┏━━━*vpd commit in {} ms", in_ms as f64 / 1000.0);
vpd_prv_inc(in_ms);

        self.mt_l.nodes[0]
    }

    pub fn open(&mut self, mu: T, q_circ: &mut Circuit<'a, T>, q_inputs: &Vec<T>) {
        // send all msg in plaintext
        self.ch.force_plain();

        // step 1: commit mask for l*q

        channel0_rst();
let start = Instant::now();

        let n = 1 << self.log_len;
        let k = self.n_ldt_iter << self.coset.log_layer_ratio;

        q_circ.eval(q_inputs);
        let q_depth = q_circ.depth;

    let start0 = Instant::now();
        if T::CHAR == 2 {
            append_ifft(q_circ, self.log_l, &self.coset.fft_coeff);
            // append_fft(q_circ, self.log_l, &self.coset.fft_coeff);
        } else {
            append_ifft_fprime(q_circ, self.log_l, self.coset);
            // append_fft_fprime(q_circ, self.log_l, self.coset);
        }

    let elapsed = start0.elapsed();
    let exclude_time = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ┌──────#(excluded) q circ append ifft/fft layers: {} ms", exclude_time as f64 / 1000.0);

    #[cfg(feature = "microbench")] let start0 = Instant::now();
        #[cfg(not(feature = "memsave"))]
        q_circ.continue_eval(q_depth);

        #[cfg(feature = "memsave")]
        q_circ.continue_eval_only(q_depth);

        let q_depth = q_circ.depth;

    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#q circ evaluation:              {:9.3} ms", in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();

        let q_ex = (self.coset.fft)(
            self.coset,
            &q_circ.layers[q_depth - 1].values,
            self.log_len, self.log_l);

    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#fft of q(x) {:>2}->{:>2}:             {:9.3} ms", self.log_len, self.log_l, in_ms as f64 / 1000.0);

        let mut s_coeff: Vec<T> = Vec::with_capacity(4 << self.log_len);
        for _ in 0..(2 * n + k - 1) {
            s_coeff.push(T::random());
        }

        // lq_s degree: 2n - 2 + k
        // lq_s = h * vanishing + g
        // g degree: n - 1
        // h degree: n - 2 + k
        let (m, log_m) = if 2 * n + k - 1 <= 4 * n {
                (4, 2)
            } else if 2 * n + k - 1 <= 8 * n {
                (8, 3)
            } else if 2 * n + k - 1 <= 16 * n {
                (16, 4)
            } else if 2 * n + k - 1 <= 32 * n {
                (32, 5)
            } else {
                // TODO
                unreachable!()
            };

        let log_m2 = if n + k - 1 <= 2 * n {
                1
            } else if n + k - 1 <= 4 * n {
                2
            } else if n + k - 1 <= 8 * n {
                3
            } else if n + k - 1 <= 16 * n {
                4
            } else if n + k - 1 <= 32 * n {
                5
            } else {
                // TODO
                unreachable!()
            };

        s_coeff.resize(m << self.log_len, T::from(0));
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let s_ex = (self.coset.fft)(self.coset, &s_coeff, self.log_len + log_m, self.log_l);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#fft of mask s(x) {:>2}->{:>2}:        {:9.3} ms", self.log_len + log_m, self.log_l, in_ms as f64 / 1000.0);
        let mut sum_s = T::from(0);
        if T::CHAR == 2 {
            for i in 0..n {
                sum_s += s_ex[i];
            }
        } else {
            for i in 0..n {
                sum_s += s_ex[i << (self.coset.log_rho)];
            }
        }

    #[cfg(feature = "microbench")] let start0 = Instant::now();
        (self.coset.rscode_commit)(self.coset, &mut self.mt_s, &s_ex, 0);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#merkle tree commit of s(x) {:>2}:  {:9.3} ms", self.log_l, in_ms as f64 / 1000.0);

        self.ch.send(0, Tx::<T>::VPDOpen0((mu, sum_s, self.mt_s.nodes[0])));
        let alpha = recv!(self.ch, Random, 1);

        // step 2: division and commit each part

        let mut lq: Vec<T> = Vec::with_capacity(1 << self.log_l);
        let mut lq_s: Vec<T> = Vec::with_capacity(1 << self.log_l);
        for i in 0..(1 << self.log_l) {
            // let q = q_circ.layers[q_depth - 1].values[i];
            let tmp = self.input_ex[i] * q_ex[i];
            lq.push(tmp);
            lq_s.push(tmp * alpha + s_ex[i]);
        }
        // clear mem
        self.input_ex = vec![];

    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let mut lq_s_coeff =
            if T::CHAR == 2 {
                (self.coset.ifft)(self.coset, &lq_s, self.log_len + log_m)
            } else {
                (self.coset.ifft)(self.coset, &lq_s, self.log_l)
            };
        lq_s_coeff.truncate(2 * n + k - 1);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#ifft of (αl'q+s)(x) {:>2}:         {:9.3} ms", self.log_len + log_m, in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let mut h_coeff = (self.coset.poly_over_vanishing)(self.coset, &lq_s_coeff);
        h_coeff.resize(1 << (self.log_len + log_m2), T::from(0));
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#polynomial division:            {:9.3} ms", in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let mut h_ex = (self.coset.fft)(self.coset, &h_coeff, self.log_len + log_m2, self.log_l);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#fft of quotient h(x) {:>2}->{:>2}:    {:9.3} ms", self.log_len + log_m2, self.log_l, in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        (self.coset.rscode_commit)(self.coset, &mut self.mt_h, &h_ex, 0);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#merkle tree commit of h(x) {:>2}:  {:9.3} ms", self.log_l, in_ms as f64 / 1000.0);

        self.ch.send(0, Tx::<T>::VPDOpen1(self.mt_h.nodes[0]));
        let randoms = recv!(self.ch, MultiRandom, 1);

        // step 3: random combination of lifted polynomials
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        // prove that lq, h, s, p has degree less than 2n+k-1, n+k-1, 2n+k-1, n-1
        // offset: 0, n, 0, n+k

        let p_n_1 = (mu * alpha + sum_s) * self.coset.ksi;

        // let mut p_ex = lq_s.to_vec();
        let mut p_ex = Vec::with_capacity(1 << self.log_l);
        if T::CHAR == 2 {
            for i in 0..(1 << self.log_l) {
                // p_ex[i] = (p_ex[i] - self.coset.vanishing_eval[i] * h_ex[i] - p_n_1 * self.coset.x_n_1[i]) * self.coset.x_n_k[i];
                p_ex.push((lq_s[i] - self.coset.vanishing_eval[i] * h_ex[i] - p_n_1 * self.coset.x_n_1[i]) * self.coset.x_n_k[i]);
            }
        } else {
            for i in 0..(1 << self.log_l) {
                // p_ex[i] = (p_ex[i] - self.coset.vanishing_eval[i] * h_ex[i] - p_n_1) * self.coset.x_n_1[i] * self.coset.x_n_k[i];
                p_ex.push((lq_s[i] - self.coset.vanishing_eval[i] * h_ex[i] - p_n_1) * self.coset.x_n_1[i] * self.coset.x_n_k[i]);
            }
        }

        for i in 0..(1 << self.log_l) {
            h_ex[i] *= self.coset.x_n[i];
        }

        // random combinations
        let mut c_ex: Vec<T> = Vec::with_capacity(1 << self.log_l);
        for i in 0..(1 << self.log_l) {
            c_ex.push(randoms[0] * lq[i] + randoms[1] * h_ex[i] + randoms[2] * s_ex[i] + randoms[3] * p_ex[i]);
        }
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#poly degree lift and combine:   {:9.3} ms", in_ms as f64 / 1000.0);

// let elapsed = start1.elapsed();
// let us = elapsed.as_micros() as usize as f64 / 1000.0;
// println!("  *vpd frac 2: {} ms", us);

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
let wait_ms = channel0_get();
channel0_rst();
#[cfg(feature = "microbench")]
println!("  ┣━┷━*vpd open time w/o ldt/gkr: {} ms", (in_ms - wait_ms - exclude_time) as f64 / 1000.0);
#[cfg(not(feature = "microbench"))]
println!("  ┣━━━*vpd open time w/o ldt/gkr: {} ms", (in_ms - wait_ms - exclude_time) as f64 / 1000.0);
vpd_prv_inc(in_ms - wait_ms - exclude_time);

let start = Instant::now();

        let log_degree = log_m + self.log_len;
        let mut ldt_prover = LDTProver::<T, Ch>::new(self.ch, &self.coset, log_degree, 2 * n + k - 1);
        // let mut ldt_prover = LDTProver::<T, Ch>::new(self.ch, &self.coset, log_degree, (1 << log_degree));

        let query_pts = ldt_prover.prove(&mut c_ex, self.n_ldt_iter);

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
let wait_ms = channel0_get();
channel0_rst();
// #[cfg(feature = "microbench")]
// println!("  ┣━┷━*vpd ldt in {} ms", (in_ms - wait_ms) as f64 / 1000.0);
// #[cfg(not(feature = "microbench"))]
println!("  ┣━━━*vpd ldt in {} ms", (in_ms - wait_ms) as f64 / 1000.0);
vpd_prv_inc(in_ms - wait_ms);
// let start = Instant::now();

        let mut query_cids = Vec::with_capacity(self.n_ldt_iter);
        for p in query_pts.iter() {
            let (cid, _) = (self.coset.subcoset_id)(self.coset, *p);
            query_cids.push(cid);
        }

        let mut mto_l = MTO::<T>::new(&self.mt_l, self.n_ldt_iter);
        mto_l.ids = query_cids.to_vec();
        self.mt_l.open(&mut mto_l);
        let mut mto_s = MTO::<T>::new(&self.mt_s, self.n_ldt_iter);
        mto_s.ids = query_cids.to_vec();
        self.mt_s.open(&mut mto_s);
        let mut mto_h = MTO::<T>::new(&self.mt_h, self.n_ldt_iter);
        mto_h.ids = query_cids.to_vec();
        self.mt_h.open(&mut mto_h);

        self.mt_l.clean();
        self.mt_s.clean();
        self.mt_h.clean();

        // assert no statistics in LDT Prover/Verifier
        self.ch.send(0, Tx::<T>::LDTQuery(vec![mto_l, mto_s, mto_h]));
        let _ = recv!(self.ch, @Syn, 1);

let start = Instant::now();

        if T::CHAR == 2 {
            let q_power_series = self.coset.additive_power_series(&query_pts, self.coset.log_len);
            append_poly_eval(q_circ, &q_power_series);
        } else {
            let q_power_series = self.coset.multiplicative_power_series(&query_pts, self.coset.log_len);
            append_poly_eval_fprime(q_circ, &q_power_series);
        };

        q_circ.continue_eval(q_depth);

let elapsed = start.elapsed();
let eval_ms = elapsed.as_micros() as usize;

        let q_depth = q_circ.depth;

        self.ch.send(0, Tx::<T>::Output(q_circ.layers[q_depth - 1].values.to_vec()));
        let _ = recv!(self.ch, @Syn, 1);

// let elapsed = start.elapsed();
// let in_ms = elapsed.as_micros() as usize;
// let wait_ms = channel0_get();
channel0_rst();
// println!("  *vpd frac {} ms", (in_ms - wait_ms) as f64 / 1000.0);
// vpd_prv_inc(in_ms - wait_ms);

let start = Instant::now();
        // assert no statistics in plain GKR
        let mut prover = GKRProver::<T, Ch>::new(q_circ, self.ch);
        prover.prove();

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
let wait_ms = channel0_get();
channel0_rst();
// #[cfg(feature = "microbench")]
// println!("  ┣━┷━*vpd gkr in {} ms", (in_ms - wait_ms) as f64 / 1000.0);
// #[cfg(not(feature = "microbench"))]
println!("  ┣━━━*vpd gkr in {} ms", (in_ms + eval_ms - wait_ms) as f64 / 1000.0);
vpd_prv_inc(in_ms + eval_ms - wait_ms);
    }
}


pub struct Verifier<'a, T: Field, Ch: Channel<Output = T>> {
    pub n_ldt_iter: usize,
    pub ch: &'a mut Ch,
    log_len: usize,
    log_l: usize,
    coset: &'a Coset<T>,
}

impl<'a, T: Field, Ch: Channel<Output = T>> Verifier<'a, T, Ch> {
    pub fn new(ch: &'a mut Ch, coset: &'a Coset<T>, n_ldt_iter: usize) -> Self {
        Self {
            n_ldt_iter: n_ldt_iter,
            ch: ch,
            log_len: coset.log_len,
            log_l: coset.log_l,
            coset: coset,
        }
    }

    // TODO make q_inputs useful
    pub fn verify(&mut self, mu: T, l_commit: SHA256, q_circ: &mut Circuit<'a, T>, _q_inputs: &Vec<T>) -> (bool, T) {
        let (_, sum_s, s_commit) = recv!(self.ch, VPDOpen0, 0);

        // without bw for mu
        vpd_bw_inc(T::SIZE + SHA256::SIZE);
        let alpha = T::random();
        self.ch.send(1, Tx::<T>::Random(alpha));

        let h_commit = recv!(self.ch, VPDOpen1, 0);
        vpd_bw_inc(SHA256::SIZE);
        let randoms = vec![T::random(), T::random(), T::random(), T::random()];
        self.ch.send(1, Tx::<T>::MultiRandom(randoms.to_vec()));

        let n = 1 << self.log_len;
        let k = self.n_ldt_iter << self.coset.log_layer_ratio;
        let log_m = if 2 * n + k - 1 <= 4 * n {
                2
            } else if 2 * n + k - 1 <= 8 * n {
                3
            } else if 2 * n + k - 1 <= 16 * n {
                4
            } else if 2 * n + k - 1 <= 32 * n {
                5
            } else {
                // TODO
                unreachable!()
            };
        let log_degree = log_m + self.log_len;
        let mut ldt_verifier = LDTVerifier::<T, Ch>::new(self.ch, &self.coset, log_degree, 2 * n + k - 1);
        // let mut ldt_verifier = LDTVerifier::<T, Ch>::new(self.ch, &self.coset, log_degree, (1 << log_degree));
        let (query_pts, mtos) = ldt_verifier.verify(self.n_ldt_iter);

        let mut query_cids = Vec::with_capacity(self.n_ldt_iter);
        let mut query_eids = Vec::with_capacity(self.n_ldt_iter);
        for p in query_pts.iter() {
            let (cid, eid) = (self.coset.subcoset_id)(self.coset, *p);
            query_cids.push(cid);
            query_eids.push(eid);
        }

        let new_mtos = recv!(self.ch, LDTQuery, 0);
        // vpd_bw_inc(new_mtos[0].proof_size + new_mtos[1].proof_size + new_mtos[2].proof_size);
        vpd_bw_inc(new_mtos[0].get_avg_proof_size(self.coset.log_l) * 3);
        self.ch.send(1, Tx::<T>::Syn);

        let q = recv!(self.ch, Output, 0);
        vpd_bw_inc(T::SIZE * q.len());
        self.ch.send(1, Tx::<T>::Syn);

let start = Instant::now();

        assert!(new_mtos[0].verify(&l_commit));
        assert!(new_mtos[1].verify(&s_commit));
        assert!(new_mtos[2].verify(&h_commit));

        for i in 0..self.n_ldt_iter {
            let x = query_pts[i];
            let eid = query_eids[i];

            let x_n_1  = self.coset.x_n_1[x];
            let x_n_k = self.coset.x_n_k[x];
            let x_n = self.coset.x_n[x];

            let c_x = mtos[0].get_value(i)[eid];
            let l_x = new_mtos[0].get_value(i)[eid];
            let s_x = new_mtos[1].get_value(i)[eid];
            let h_x = new_mtos[2].get_value(i)[eid];

            let q_x = q[i];
            let lq_x = l_x * q_x;

            // p * ksi = alq + s - Wx^{n-1} - h Z
            let p_x =
                if T::CHAR == 2 {
                    alpha * lq_x + s_x - h_x * self.coset.vanishing_eval[x] -
                        (mu * alpha + sum_s) * x_n_1 * self.coset.ksi
                } else {
                    (alpha * lq_x + s_x - h_x * self.coset.vanishing_eval[x] -
                     (mu * alpha + sum_s) * self.coset.ksi) * x_n_1
                };

            // c_ex = randoms[0] * lq + randoms[1] * h[i] + randoms[2] * s[i] + randoms[3] * p_ex[i]);
            assert_eq!(c_x,
                       randoms[0] * lq_x +
                       randoms[1] * h_x * x_n +
                       randoms[2] * s_x +
                       randoms[3] * p_x * x_n_k);
        }

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
vpd_vrf_inc(in_ms);

        if T::CHAR == 2 {
            append_ifft(q_circ, self.log_l, &self.coset.fft_coeff);
            // append_fft(q_circ, self.log_l, &self.coset.fft_coeff);
        } else {
            append_ifft_fprime(q_circ, self.log_l, self.coset);
            // append_fft_fprime(q_circ, self.log_l, self.coset);
        }

let start = Instant::now();
        let q_power_series =
            if T::CHAR == 2 {
                self.coset.additive_power_series(&query_pts, self.log_len)
            } else {
                self.coset.multiplicative_power_series(&query_pts, self.log_len)
            };
let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
vpd_vrf_inc(in_ms);

        // append_selector(q_circ, &query_pts);
        if T::CHAR == 2 {
            append_poly_eval(q_circ, &q_power_series);
        } else {
            append_poly_eval_fprime(q_circ, &q_power_series);
        };

        let mut verifier = GKRVerifier::<T, Ch>::new(&q_circ, self.ch);
        let res = verifier.verify(&q);

        assert!(res);
        (true, mu)
    }
}
