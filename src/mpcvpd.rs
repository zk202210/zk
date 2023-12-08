use std::time::Instant;

use rand::Rng;

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
use crate::tape::*;
use crate::statistic::*;
use crate::util::*;

pub struct Prover<'a, T: Field, Ch: Channel<Output = T>> {
    pub n_ldt_iter: usize,
    pub n_mpc_iter: usize,
    pub n_party: usize,
    pub ch: &'a mut Ch,
    pub seedtree: SeedTree,
    pub tapes: Vec<Vec<Tape>>,
    log_len: usize,
    log_l: usize,
    coset: &'a Coset<T>,
    mt_l: MT<T>,
    mt_h: MT<T>,
    lp_ex: Vec<T>,
    lr_coeff_shares: Vec<Vec<Vec<T>>>,
    l_ex: Vec<T>,
}

impl<'a, T: Field, Ch: Channel<Output = T>> Prover<'a, T, Ch> {
    pub fn new(ch: &'a mut Ch, coset: &'a Coset<T>, n_ldt_iter: usize, n_party: usize, n_mpc_iter: usize) -> Self {

let start = Instant::now();
        let mseed = Seed::random();
        let st = SeedTree::new(&mseed, n_party, n_mpc_iter);

        let mut tapes = vec![Vec::with_capacity(n_party); n_mpc_iter];
        let tape_len = n_ldt_iter * (2 + (1 << coset.log_layer_ratio));
        for i in 0..n_mpc_iter {
            for p in 0..n_party {
                let seed = st.get_seed(i, p);
                tapes[i].push(Tape::new(seed, tape_len * T::SIZE));
            }
        }

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
vpd_prv_inc(in_ms);
#[cfg(feature = "microbench")] println!("    ┌──────#gen seeds random tape {:>4}:     {:9.3} ms", tape_len, in_ms as f64 / 1000.0);

        Self {
            n_ldt_iter,
            n_mpc_iter,
            n_party,
            ch: ch,
            seedtree: st,
            tapes,
            log_len: coset.log_len,
            log_l: coset.log_l,
            // log_kappa: log_kappa,
            coset: coset,
            mt_l: MT::<T>::new(coset.log_l - coset.log_layer_ratio, coset.log_layer_ratio),
            mt_h: MT::<T>::new(coset.log_l - coset.log_layer_ratio, coset.log_layer_ratio),
            lp_ex: vec![],
            lr_coeff_shares: vec![],
            l_ex: vec![],
        }
    }

    pub fn commit(&mut self, inputs: &Vec<T>) -> SHA256 {

channel0_rst();
let start = Instant::now();
        // seed tree commitment
        let stc = self.seedtree.commit();

        self.ch.send(0, Tx::<T>::SeedTree(stc));
        let _ = recv!(self.ch, @Syn, 1);

    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let mut l_coeff = (self.coset.ifft)(self.coset, inputs, self.log_len);
        self.l_ex = (self.coset.fft)(self.coset, &l_coeff, self.log_len, self.log_l);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("    ├──────#interpolation of l(x) {:>2}->{:>2}:   {:9.3} ms", self.log_len, self.log_l, in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let lp_coeff = &mut l_coeff;
        let mut lr_coeff = Vec::with_capacity(self.n_ldt_iter);
        for i in 0..self.n_ldt_iter {
            let r = T::random();
            lp_coeff[i] -= r;
            lr_coeff.push(r);
        }
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("    ├──────#decompose l(x) = lᵣ(x) + l'(x): {:9.3} ms", in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let mut lr_coeff_icr = vec![lr_coeff.clone(); self.n_mpc_iter];
        let mut lr_coeff_shares = vec![vec![Vec::with_capacity(self.n_party); self.n_mpc_iter]; lr_coeff.len()];
        for i in 0..self.n_mpc_iter {
            for j in 0..lr_coeff.len() {
                for p in 0..self.n_party {
                    let share = self.tapes[i][p].get();
                    lr_coeff_icr[i][j] -= share;
                    lr_coeff_shares[j][i].push(share);
                }
                lr_coeff_shares[j][i][0] += lr_coeff_icr[i][j];
            }
        }
        self.lr_coeff_shares = lr_coeff_shares;
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("    ├──────#mpc commit of lᵣ(x) coeffs:     {:9.3} ms", in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let lp_ex = (self.coset.fft)(self.coset, lp_coeff, self.log_len, self.log_l);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("    ├──────#fft of l'(x) {:>2}->{:>2}:            {:9.3} ms", self.log_len, self.log_l, in_ms as f64 / 1000.0);
        self.ch.send(0, Tx::<T>::InputZK(lr_coeff_icr.clone()));
        let _ = recv!(self.ch, @Syn, 1);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        (self.coset.rscode_commit)(self.coset, &mut self.mt_l, &lp_ex, 0);
        self.lp_ex = lp_ex;
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("    ├──────#merkle tree commit of l'(x) {:>2}: {:9.3} ms", self.log_l, in_ms as f64 / 1000.0);

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
let wait_ms = channel0_get();
let in_ms = in_ms - wait_ms;
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

        // // zk LDT requre another mask TODO XXX TODO
        // let mut s_coeff: Vec<T> = Vec::with_capacity(2 << self.log_len);
        // for _ in 0..(2 * n + k - 1) {
        //     s_coeff.push(T::random());
        // }

        // step 2: division and commit each part

        let mut lq: Vec<T> = Vec::with_capacity(1 << self.log_l);
        // let mut lq_s: Vec<T> = Vec::with_capacity(1 << self.log_l);
        for i in 0..(1 << self.log_l) {
            // let q = q_circ.layers[q_depth - 1].values[i];
            let tmp = self.l_ex[i] * q_ex[i];
            lq.push(tmp);
            // lq_s.push(tmp * alpha + s_ex[i]);
        }
        // clear mem
        self.l_ex = vec![];

    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let mut lq_coeff =
            if T::CHAR == 2 {
                (self.coset.ifft)(self.coset, &lq, self.log_len + 1)
            } else {
                (self.coset.ifft)(self.coset, &lq, self.log_l)
            };
        lq_coeff.truncate(2 * n);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#ifft of lq(x) {:>2}:               {:9.3} ms", self.log_len + 1, in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let (mut h_coeff, mut p_coeff) = (self.coset.poly_over_vanishing_w_remainder)(self.coset, &lq_coeff);
        if T::CHAR == 2 {
            // assert_eq!(p_coeff[n - 1], mu * self.coset.ksi);
            p_coeff[n - 1] = T::from(0);
        } else {
            p_coeff = p_coeff[1..].to_vec();
            p_coeff.resize(n, T::from(0));
        }
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#polynomial division lq = hZ + g:{:9.3} ms", in_ms as f64 / 1000.0);

    // working on quotient h(x)
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let hp_coeff = &mut h_coeff;
        let mut hr_coeff = Vec::with_capacity(self.n_ldt_iter);
        for i in 0..self.n_ldt_iter {
            let r = T::random();
            hp_coeff[i] -= r;
            hr_coeff.push(r);
        }
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#decompose quotient h = hᵣ + h': {:9.3} ms", in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let mut hr_coeff_icr = vec![hr_coeff.clone(); self.n_mpc_iter];
        let mut hr_coeff_shares = vec![vec![Vec::with_capacity(self.n_party); self.n_mpc_iter]; hr_coeff.len()];
        for i in 0..self.n_mpc_iter {
            for j in 0..hr_coeff.len() {
                for p in 0..self.n_party {
                    let share = self.tapes[i][p].get();
                    hr_coeff_icr[i][j] -= share;
                    hr_coeff_shares[j][i].push(share);
                }
                hr_coeff_shares[j][i][0] += hr_coeff_icr[i][j];
            }
        }
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#mpc commit of hᵣ(x) coeffs:     {:9.3} ms", in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let hp_ex = (self.coset.fft)(self.coset, &hp_coeff, self.log_len, self.log_l);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#fft of quotient h'(x) {:>2}->{:>2}:   {:9.3} ms", self.log_len, self.log_l, in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        (self.coset.rscode_commit)(self.coset, &mut self.mt_h, &hp_ex, 0);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#merkle tree commit of h'(x) {:>2}: {:9.3} ms", self.log_l, in_ms as f64 / 1000.0);

    // working on modified remainder p(x)
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let pp_coeff = &mut p_coeff;
        let mut pr_coeff = Vec::with_capacity(self.n_ldt_iter << self.coset.log_layer_ratio);
        for i in 0..(self.n_ldt_iter << self.coset.log_layer_ratio) {
            let r = T::random();
            pp_coeff[i] -= r;
            pr_coeff.push(r);
        }
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#decompose p = pᵣ + p':          {:9.3} ms", in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let mut pr_coeff_icr = vec![pr_coeff.clone(); self.n_mpc_iter];
        let mut pr_coeff_shares = vec![vec![Vec::with_capacity(self.n_party); self.n_mpc_iter]; pr_coeff.len()];
        for i in 0..self.n_mpc_iter {
            for j in 0..pr_coeff.len() {
                for p in 0..self.n_party {
                    let share = self.tapes[i][p].get();
                    pr_coeff_icr[i][j] -= share;
                    pr_coeff_shares[j][i].push(share);
                }
                pr_coeff_shares[j][i][0] += pr_coeff_icr[i][j];
            }
        }
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#mpc commit of pᵣ(x) coeffs:     {:9.3} ms", in_ms as f64 / 1000.0);
    #[cfg(feature = "microbench")] let start0 = Instant::now();
        let pp_ex = (self.coset.fft)(self.coset, &pp_coeff, self.log_len, self.log_l);
    #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "microbench")] println!("  ┃ ├──────#fft of p'(x) {:>2}->{:>2}:            {:9.3} ms", self.log_len, self.log_l, in_ms as f64 / 1000.0);

        self.ch.send(0, Tx::<T>::InputZK(hr_coeff_icr.clone()));
        let _ = recv!(self.ch, @Syn, 1);

        self.ch.send(0, Tx::<T>::InputZK(pr_coeff_icr.clone()));
        let _ = recv!(self.ch, @Syn, 1);

        self.ch.send(0, Tx::<T>::VPDOpen1(self.mt_h.nodes[0]));
        let _ = recv!(self.ch, @Syn, 1);
    //     let randoms = recv!(self.ch, MultiRandom, 1);

    //     // step 3: random combination of lifted polynomials
    // #[cfg(feature = "microbench")] let start0 = Instant::now();
    //     // prove that lq, h, p has degree less than 2n-1, n, n-1
    //     // random combinations
    //     let mut c_ex: Vec<T> = Vec::with_capacity(1 << self.log_l);
    //     for i in 0..(1 << self.log_l) {
    //         c_ex.push(randoms[0] * lq[i] +
    //                   randoms[1] * hp_ex[i] * self.coset.x_n[i] +
    //                   randoms[2] * pp_ex[i] * self.coset.x_n_k[i]);
    //     }
    // #[cfg(feature = "microbench")] let elapsed = start0.elapsed();
    // #[cfg(feature = "microbench")] let in_ms = elapsed.as_micros() as usize;
    // #[cfg(feature = "microbench")] println!("  ┃ ├──────#poly degree lift and combine:   {:9.3} ms", in_ms as f64 / 1000.0);

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

        let log_degree = 1 + self.log_len;
        let mut ldt_prover = LDTProver::<T, Ch>::new(self.ch, &self.coset, log_degree, n - 1);
        // let mut ldt_prover = LDTProver::<T, Ch>::new(self.ch, &self.coset, log_degree, (1 << log_degree));

        let query_pts = ldt_prover.prove(&mut pp_ex.clone(), self.n_ldt_iter);

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

let start = Instant::now();
        let mut query_cids = Vec::with_capacity(self.n_ldt_iter);
        for p in query_pts.iter() {
            let (cid, _) = (self.coset.subcoset_id)(self.coset, *p);
            query_cids.push(cid);
        }

        let mut mto_l = MTO::<T>::new(&self.mt_l, self.n_ldt_iter);
        mto_l.ids = query_cids.to_vec();
        self.mt_l.open(&mut mto_l);
        let mut mto_h = MTO::<T>::new(&self.mt_h, self.n_ldt_iter);
        mto_h.ids = query_cids.to_vec();
        self.mt_h.open(&mut mto_h);

        let q_power_series =
            if T::CHAR == 2 {
                // assert n_ldt_iter * layer_ratio < (1 << log_len)
                self.coset.additive_power_series(&query_pts, self.coset.log_len)
            } else {
                self.coset.multiplicative_power_series(&query_pts, self.coset.log_len)
            };

        // begin mpc open
        let mut output_shares = Vec::with_capacity(self.n_party * self.n_mpc_iter * self.n_ldt_iter);
        for i in 0..self.n_ldt_iter {
            let x = query_pts[i];

            let mut power_ex = vec![T::from(1); self.n_ldt_iter << self.coset.log_layer_ratio];
            for k in 0..(self.n_ldt_iter << self.coset.log_layer_ratio) {
                let mut x = k;
                for j in 0..ceil_log(self.n_ldt_iter << self.coset.log_layer_ratio) {
                    if (x & 1) == 1 {
                        power_ex[k] *= q_power_series[i][j];
                    }
                    x >>= 1;
                }
            }
            let x_n_1  = self.coset.x_n_1[x];
            // let x_n_k = self.coset.x_n_k[x];
            // let x_n = self.coset.x_n[x];

            let p_x = pp_ex[x];
            let q_x = q_ex[x];
            let lp_x = self.lp_ex[x];
            let hp_x = hp_ex[x];
            let z_x = self.coset.vanishing_eval[x];

            // tmp = pp_x * x_n_k * randoms[2] + lr_x * q_x * randoms[0]
            // let tmp = c_x - randoms[0] * lp_x * q_x - randoms[1] * hp_x * x_n;
            let tmp = p_x;

            for ii in 0..self.n_mpc_iter {
                // let mut output = T::from(0);
                for p in 0..self.n_party {
                    let mut hr_x_share = T::from(0);
                    let mut lr_x_share = T::from(0);
                    let mut pr_x_share = T::from(0);
                    for j in 0..self.n_ldt_iter {
                        lr_x_share += self.lr_coeff_shares[j][ii][p] * power_ex[j];
                        hr_x_share += hr_coeff_shares[j][ii][p] * power_ex[j];
                    }
                    for j in 0..(self.n_ldt_iter << self.coset.log_layer_ratio) {
                        pr_x_share += pr_coeff_shares[j][ii][p] * power_ex[j];
                    }
                    let l_x_share = if p != 0 { lr_x_share } else { lr_x_share + lp_x };
                    let lq_x_share = l_x_share * q_x;
                    let h_x_share = if p != 0 { hr_x_share } else { hr_x_share + hp_x };
                    let p_x_share =
                        if T::CHAR == 2 {
                            if p == 0 {
                                lq_x_share - h_x_share * z_x - mu * self.coset.ksi * x_n_1
                            } else {
                                lq_x_share - h_x_share * z_x
                            }
                        } else {
                            if p == 0 {
                                (lq_x_share - h_x_share * z_x - mu * self.coset.ksi) * x_n_1
                            } else {
                                (lq_x_share - h_x_share * z_x) * x_n_1
                            }
                        };
                    let pp_x_share = p_x_share - pr_x_share;

                    let mut output_share = pp_x_share; // * x_n_k * randoms[2] + lr_x_share * q_x * randoms[0];
                    if p == 0 {
                        output_share -= tmp;
                    }
                    // output += output_share;
                    output_shares.push(output_share);
                }
                // assert_eq!(output, T::from(0));
            }
        }
        let mut h = SHA256::from(0);
        h.commit(&output_shares);

        self.mt_l.clean();
        self.mt_h.clean();
let elapsed = start.elapsed();
let mpc_ms = elapsed.as_micros() as usize;
vpd_prv_inc(mpc_ms);

        self.ch.send(0, Tx::<T>::Hash(h));
        let _ = recv!(self.ch, @Syn, 1);

        // assert no statistics in LDT Prover/Verifier
        self.ch.send(0, Tx::<T>::LDTQuery(vec![mto_l, mto_h]));
        let _ = recv!(self.ch, @Syn, 1);

let start = Instant::now();

        if T::CHAR == 2 {
            append_poly_eval(q_circ, &q_power_series);
        } else {
            append_poly_eval_fprime(q_circ, &q_power_series);
        };

        q_circ.continue_eval(q_depth);

let elapsed = start.elapsed();
let eval_ms = elapsed.as_micros() as usize;

        let q_depth = q_circ.depth;

        self.ch.send(0, Tx::<T>::Output(q_circ.layers[q_depth - 1].values.to_vec()));
        let missings = recv!(self.ch, Party, 1);

let start = Instant::now();
        let mut sto = SeedTreeOpen::new(&self.seedtree);
        sto.missings = missings;
        self.seedtree.open(&mut sto);
let elapsed = start.elapsed();
let mpc_seed_ms = elapsed.as_micros() as usize;
vpd_prv_inc(mpc_ms);
println!("  ┣━━━*vpd mpc open in {} ms", (mpc_ms + mpc_seed_ms) as f64 / 1000.0);

        self.ch.send(0, Tx::<T>::Seed((SHA256::from(0), sto)));
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
    pub n_mpc_iter: usize,
    pub n_party: usize,
    pub ch: &'a mut Ch,
    log_len: usize,
    log_l: usize,
    coset: &'a Coset<T>,
}

impl<'a, T: Field, Ch: Channel<Output = T>> Verifier<'a, T, Ch> {
    pub fn new(ch: &'a mut Ch, coset: &'a Coset<T>, n_ldt_iter: usize, n_party: usize, n_mpc_iter: usize) -> Self {
        Self {
            n_ldt_iter,
            n_mpc_iter,
            n_party,
            ch,
            log_len: coset.log_len,
            log_l: coset.log_l,
            coset: coset,
        }
    }

    // TODO make q_inputs useful
    pub fn verify(&mut self, mu: T, l_commit: SHA256, stc: SHA256, lr_coeff_icr: &Vec<Vec<T>>, q_circ: &mut Circuit<'a, T>, _q_inputs: &Vec<T>) -> (bool, T) {
        // let (_, sum_s, s_commit) = recv!(self.ch, VPDOpen0, 0);

        // without bw for mu
        // vpd_bw_inc(T::SIZE + SHA256::SIZE);
        // let alpha = T::random();
        // self.ch.send(1, Tx::<T>::Random(alpha));

        let hr_coeff_icr = recv!(self.ch, InputZK, 0);
        self.ch.send(1, Tx::<T>::Syn);
        vpd_bw_inc(hr_coeff_icr.len() * hr_coeff_icr[0].len() * T::SIZE);

        let pr_coeff_icr = recv!(self.ch, InputZK, 0);
        self.ch.send(1, Tx::<T>::Syn);
        vpd_bw_inc(pr_coeff_icr.len() * pr_coeff_icr[0].len() * T::SIZE);

        let h_commit = recv!(self.ch, VPDOpen1, 0);
        self.ch.send(1, Tx::<T>::Syn);
        vpd_bw_inc(SHA256::SIZE);
        // let randoms = vec![T::random(), T::random(), T::random()];
        // self.ch.send(1, Tx::<T>::MultiRandom(randoms.to_vec()));

        let n = 1 << self.log_len;

        let log_degree = 1 + self.log_len;
        let mut ldt_verifier = LDTVerifier::<T, Ch>::new(self.ch, &self.coset, log_degree, 2 * n - 1);
        // let mut ldt_verifier = LDTVerifier::<T, Ch>::new(self.ch, &self.coset, log_degree, (1 << log_degree));
        let (query_pts, mtos) = ldt_verifier.verify(self.n_ldt_iter);

        let mut query_cids = Vec::with_capacity(self.n_ldt_iter);
        let mut query_eids = Vec::with_capacity(self.n_ldt_iter);
        for p in query_pts.iter() {
            let (cid, eid) = (self.coset.subcoset_id)(self.coset, *p);
            query_cids.push(cid);
            query_eids.push(eid);
        }

        let h = recv!(self.ch, Hash, 0);
        self.ch.send(1, Tx::<T>::Syn);

        let new_mtos = recv!(self.ch, LDTQuery, 0);
        // vpd_bw_inc(new_mtos[0].proof_size + new_mtos[1].proof_size + new_mtos[2].proof_size);
        vpd_bw_inc(new_mtos[0].get_avg_proof_size(self.coset.log_l) * 3);
        self.ch.send(1, Tx::<T>::Syn);

        let q = recv!(self.ch, Output, 0);
        vpd_bw_inc(T::SIZE * q.len());

        let mut missings = Vec::with_capacity(self.n_mpc_iter);
        for _ in 0..self.n_mpc_iter {
            let missing = (rand::thread_rng().gen::<u16>() as usize) % self.n_party;
            missings.push(missing);
        }
        self.ch.send(1, Tx::<T>::Party(missings.clone()));

        let (_, mut sto) = recv!(self.ch, Seed, 0);
        vpd_bw_inc(sto.proof_size);
        self.ch.send(1, Tx::<T>::Syn);

let start = Instant::now();

        assert!(new_mtos[0].verify(&l_commit));
        assert!(new_mtos[1].verify(&h_commit));
        assert!(sto.verify(&stc));

        // begin vpd mpc verify
        // gen tapes
        let mut tapes = vec![Vec::with_capacity(self.n_party); self.n_mpc_iter];
        let tape_len = self.n_ldt_iter * (2 + (1 << self.coset.log_layer_ratio));
        for i in 0..self.n_mpc_iter {
            for p in 0..self.n_party {
                let seed = sto.get_seed(i, p);
                tapes[i].push(Tape::new(seed, tape_len * T::SIZE));
            }
        }
        // shares
        let mut lr_coeff_shares = vec![vec![Vec::with_capacity(self.n_party); self.n_mpc_iter]; self.n_ldt_iter];
        for i in 0..self.n_mpc_iter {
            for j in 0..self.n_ldt_iter {
                for p in 0..self.n_party {
                    let share = tapes[i][p].get::<T>();
                    lr_coeff_shares[j][i].push(share);
                }
                lr_coeff_shares[j][i][0] += lr_coeff_icr[i][j];
            }
        }
        let mut hr_coeff_shares = vec![vec![Vec::with_capacity(self.n_party); self.n_mpc_iter]; self.n_ldt_iter];
        for i in 0..self.n_mpc_iter {
            for j in 0..self.n_ldt_iter {
                for p in 0..self.n_party {
                    let share = tapes[i][p].get::<T>();
                    hr_coeff_shares[j][i].push(share);
                }
                hr_coeff_shares[j][i][0] += hr_coeff_icr[i][j];
            }
        }
        let mut pr_coeff_shares = vec![vec![Vec::with_capacity(self.n_party); self.n_mpc_iter]; self.n_ldt_iter << self.coset.log_layer_ratio];
        for i in 0..self.n_mpc_iter {
            for j in 0..(self.n_ldt_iter << self.coset.log_layer_ratio) {
                for p in 0..self.n_party {
                    let share = tapes[i][p].get::<T>();
                    pr_coeff_shares[j][i].push(share);
                }
                pr_coeff_shares[j][i][0] += pr_coeff_icr[i][j];
            }
        }
        let mut output_shares = Vec::with_capacity(self.n_party * self.n_mpc_iter * self.n_ldt_iter);

        let q_power_series =
            if T::CHAR == 2 {
                self.coset.additive_power_series(&query_pts, self.log_len)
            } else {
                self.coset.multiplicative_power_series(&query_pts, self.log_len)
            };

        for i in 0..self.n_ldt_iter {
            let x = query_pts[i];
            let eid = query_eids[i];

            let mut power_ex = vec![T::from(1); self.n_ldt_iter << self.coset.log_layer_ratio];
            for k in 0..(self.n_ldt_iter << self.coset.log_layer_ratio) {
                let mut x = k;
                for j in 0..ceil_log(self.n_ldt_iter << self.coset.log_layer_ratio) {
                    if (x & 1) == 1 {
                        power_ex[k] *= q_power_series[i][j];
                    }
                    x >>= 1;
                }
            }

            let x_n_1  = self.coset.x_n_1[x];
            // let x_n_k = self.coset.x_n_k[x];
            // let x_n = self.coset.x_n[x];

            let c_x = mtos[0].get_value(i)[eid];
            let lp_x = new_mtos[0].get_value(i)[eid];
            let hp_x = new_mtos[1].get_value(i)[eid];
            let q_x = q[i];
            let z_x = self.coset.vanishing_eval[x];

            let tmp = c_x; // - randoms[0] * lp_x * q_x - randoms[1] * hp_x * x_n;

            for ii in 0..self.n_mpc_iter {
                let mut output = T::from(0);
                let mut missing_idx = 0;
                for p in 0..self.n_party {
                    if p == missings[ii] {
                        missing_idx = output_shares.len();
                        output_shares.push(T::from(0));
                        continue;
                    }
                    let mut hr_x_share = T::from(0);
                    let mut lr_x_share = T::from(0);
                    let mut pr_x_share = T::from(0);
                    for j in 0..self.n_ldt_iter {
                        lr_x_share += lr_coeff_shares[j][ii][p] * power_ex[j];
                        hr_x_share += hr_coeff_shares[j][ii][p] * power_ex[j];
                    }
                    for j in 0..(self.n_ldt_iter << self.coset.log_layer_ratio) {
                        pr_x_share += pr_coeff_shares[j][ii][p] * power_ex[j];
                    }

                    let l_x_share = if p != 0 { lr_x_share } else { lr_x_share + lp_x };
                    let lq_x_share = l_x_share * q_x;
                    let h_x_share = if p != 0 { hr_x_share } else { hr_x_share + hp_x };
                    let p_x_share =
                        if T::CHAR == 2 {
                            if p == 0 {
                                lq_x_share - h_x_share * z_x - mu * self.coset.ksi * x_n_1
                            } else {
                                lq_x_share - h_x_share * z_x
                            }
                        } else {
                            if p == 0 {
                                (lq_x_share - h_x_share * z_x - mu * self.coset.ksi) * x_n_1
                            } else {
                                (lq_x_share - h_x_share * z_x) * x_n_1
                            }
                        };
                    let pp_x_share = p_x_share - pr_x_share;

                    let mut output_share = pp_x_share; // * x_n_k * randoms[2] + lr_x_share * q_x * randoms[0];
                    if p == 0 {
                        output_share -= tmp;
                    }
                    output += output_share;
                    output_shares.push(output_share);
                }
                output_shares[missing_idx] = T::from(0) - output;
            }
        }
        let mut _h = SHA256::from(0);
        _h.commit(&output_shares);
        assert_eq!(h, _h);

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
