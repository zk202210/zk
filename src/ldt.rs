use std::time::Instant;

use rand::Rng;

use crate::field::Field;
use crate::field::Coset;
use crate::hash::*;
use crate::merkle::MT;
use crate::merkle::MTO;
use crate::channel::*;
use crate::util::*;
use crate::statistic::*;


pub struct Prover<'a, T: Field, Ch: Channel<Output = T>> {
    ch: &'a mut Ch,
    coset: &'a Coset<T>,
    log_degree: usize,
    degree: usize,
}

impl<'a, T: Field, Ch: Channel<Output = T>> Prover<'a, T, Ch> {
    pub fn new(ch: &'a mut Ch, coset: &'a Coset<T>, log_degree: usize, degree: usize) -> Self {
        Self {
            ch: ch,
            coset: coset,
            log_degree: log_degree,
            degree: degree,
        }
    }

    // rscode will be destroyed
    pub fn prove(&mut self, rscode: &mut Vec<T>, n_ldt_iter: usize) -> Vec<usize> {
        // commit phase
        let mut roots: Vec<SHA256> = Vec::new();
        let mut r: Vec<T> = Vec::with_capacity(self.log_degree);
        let mut mts: Vec<MT<T>> = Vec::new();

        {
            let mut mt = MT::<T>::new(self.coset.log_l - self.coset.log_layer_ratio, self.coset.log_layer_ratio);
            (self.coset.rscode_commit)(self.coset, &mut mt, &rscode, 0);

            let root = mt.nodes[0];
            mts.push(mt);
            roots.push(root);
            self.ch.send(0, Tx::<T>::LDTCommitHash(root));
        }

        let mut codelen = rscode.len();
        let mut degree = self.degree;
        let mut j = 0;
        let mut mt_size = 1 << (self.coset.log_l - self.coset.log_layer_ratio);
        loop {
            let degree_parity = degree & ((1 << self.coset.log_layer_ratio) - 1);

            // assert!(self.coset.log_l - j >= self.coset.log_layer_ratio)
            // BUGGED if we have a extremely short inputs with high layer_ratio:
            // the degree will not decrease into the range of (1, log_layer_ratio)
            // before log_l-j decreased to be less than log_layer_raio.
            // Currently just ignore this case.
            let new_r = recv!(self.ch, Random, 1);
            r.push(new_r);
            let mut mt = MT::<T>::new(self.coset.log_l - j - self.coset.log_layer_ratio - self.coset.log_layer_ratio, self.coset.log_layer_ratio);
            mt_size >>= self.coset.log_layer_ratio;

            assert_eq!(mt_size, mt.size);

            (self.coset.rscode_commit_next_layer)(self.coset, &mut mt, rscode, r[j / self.coset.log_layer_ratio], j, degree_parity);
            codelen >>= self.coset.log_layer_ratio;
            rscode.truncate(codelen);

            let root = mt.nodes[0];
            mts.push(mt);
            roots.push(root);

            degree >>= self.coset.log_layer_ratio;
            degree = if degree_parity == 0 { degree } else { degree + 1 };
            j += self.coset.log_layer_ratio;

            if n_ldt_iter >= mt_size {
                break;
            }

            self.ch.send(0, Tx::<T>::LDTCommitHash(root));
        }

        self.ch.send(0, Tx::<T>::LDTCommitValue(rscode.clone()));

        // clear mem
        *rscode = vec![];

        // query phase
        let query_pts = recv!(self.ch, LDTQueryPoint, 1);

        let mut mtos: Vec<MTO<T>> = Vec::new();

        let mut cids = Vec::with_capacity(n_ldt_iter);
        for i in 0..n_ldt_iter {
            let (cid, _) = (self.coset.subcoset_id)(self.coset, query_pts[i]);
            cids.push(cid);
        }
        for j in 0..(mts.len() - 1) {
            let mut mto = MTO::<T>::new(&mts[j], n_ldt_iter);
            mto.ids = cids.clone();
            mts[j].open(&mut mto);
            mtos.push(mto);
            for i in 0..n_ldt_iter {
                let (cid, _) = (self.coset.subcoset_id_next_layer)(self.coset, cids[i], j * self.coset.log_layer_ratio);
                cids[i] = cid;
            }
        }
        println!("  ┃\tldt peak mem {}", procinfo());
        for j in 0..mts.len() {
            mts[j].clean();
        }
        self.ch.send(0, Tx::<T>::LDTQuery(mtos));
        let _ = recv!(self.ch, @Syn, 1);

        query_pts
    }
}

pub struct Verifier<'a, T: Field, Ch: Channel<Output = T>> {
    ch: &'a mut Ch,
    coset: &'a Coset<T>,
    log_degree: usize,
    degree: usize,
}

impl<'a, T: Field, Ch: Channel<Output = T>> Verifier<'a, T, Ch> {
    pub fn new(ch: &'a mut Ch, coset: &'a Coset<T>, log_degree: usize, degree: usize) -> Self {
        Self {
            ch: ch,
            coset: coset,
            log_degree: log_degree,
            degree: degree,
        }
    }

    pub fn verify(&mut self, n_ldt_iter: usize) -> (Vec<usize>, Vec<MTO<T>>) {
        let mut roots: Vec<SHA256> = Vec::new();
        let mut r: Vec<T> = Vec::with_capacity(self.log_degree);
        let mut query_pts: Vec<usize> = Vec::with_capacity(n_ldt_iter);

        let root = recv!(self.ch, LDTCommitHash, 0);
        vpd_bw_inc(SHA256::SIZE);
        vpd_bw_hash_inc(1);
        roots.push(root);
        let mut mt_size = 1 << (self.coset.log_l - self.coset.log_layer_ratio);

        let mut degree = self.degree;
        let mut n_mt = 1;
        loop {
            let degree_parity = degree & ((1 << self.coset.log_layer_ratio) - 1);

            mt_size >>= self.coset.log_layer_ratio;

            degree >>= self.coset.log_layer_ratio;
            degree = if degree_parity == 0 { degree } else { degree + 1 };

            if n_ldt_iter >= mt_size { // || (n_ldt_iter << self.coset.log_layer_ratio) >= degree {
                break;
            }

            let new_r = T::random();
            r.push(new_r);
            self.ch.send(1, Tx::<T>::Random(new_r));
            let root = recv!(self.ch, LDTCommitHash, 0);
            vpd_bw_inc(SHA256::SIZE);
            vpd_bw_hash_inc(1);
            n_mt += 1;

            roots.push(root);
        }

        let new_r = T::random();
        r.push(new_r);
        self.ch.send(1, Tx::<T>::Random(new_r));
        let final_val = recv!(self.ch, LDTCommitValue, 0);
        vpd_bw_inc(T::SIZE * final_val.len());

        // generate sorted query pts, no duplication
        let mut i = 0;
        while i < n_ldt_iter {
            let raw: u64 = rand::thread_rng().gen();
            let pt = (raw as usize) % (1 << self.coset.log_l);
            match query_pts.binary_search(&pt) {
                Ok(_) => {},
                Err(pos) => {
                    query_pts.insert(pos, pt);
                    i += 1;
                },
            }
        }

        self.ch.send(1, Tx::<T>::LDTQueryPoint(query_pts.to_vec()));
        let mtos = recv!(self.ch, LDTQuery, 0);
        for i in 0..mtos.len() {
            // vpd_bw_inc(mtos[i].proof_size);
            vpd_bw_inc(mtos[i].get_avg_proof_size(self.coset.log_l));
        }

        self.ch.send(1, Tx::<T>::Syn);

let start = Instant::now();

        // verify
        let mut cids = Vec::with_capacity(n_ldt_iter);
        for i in 0..n_ldt_iter {
            let (cid, _) = (self.coset.subcoset_id)(self.coset, query_pts[i]);
            cids.push(cid);
        }
        let mut degree = self.degree;
        let mut j = 0;
        for k in 0..(n_mt - 1) {
    #[cfg(feature = "nanobench")] let start0 = Instant::now();
            assert!(mtos[k].verify(&roots[k]));
            let degree_parity = degree & ((1 << self.coset.log_layer_ratio) - 1);

            for i in 0..n_ldt_iter {
                let (cid, eid) = (self.coset.subcoset_id_next_layer)(self.coset, cids[i], j);
                assert_eq!(
                    (self.coset.next_layer)(
                        self.coset, j, cids[i], &mtos[k].get_value(i), r[k], degree_parity),
                    mtos[k + 1].get_value(i)[eid]);
                cids[i] = cid;
            }
            j += self.coset.log_layer_ratio;
    #[cfg(feature = "nanobench")] let elapsed = start0.elapsed();
    #[cfg(feature = "nanobench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "nanobench")] println!("  ┃         vrf ldt dgr {:>8} mt sz {:>2}:  {:9.3} ms @────", degree, mtos[k].leaf_degree, in_ms as f64 / 1000.0);
            degree >>= self.coset.log_layer_ratio;
            degree = if degree_parity == 0 { degree } else { degree + 1 };
        }
        {
    #[cfg(feature = "nanobench")] let start0 = Instant::now();
            let k = n_mt - 1;
            let degree_parity = degree & ((1 << self.coset.log_layer_ratio) - 1);
            assert!(mtos[k].verify(&roots[k]));

            for i in 0..n_ldt_iter {
                let (cid, eid) = (self.coset.subcoset_id_next_layer)(self.coset, cids[i], j);
                assert_eq!(
                    (self.coset.next_layer)(
                        self.coset, j, cids[i], &mtos[k].get_value(i), r[k], degree_parity),
                    final_val[(self.coset.normal_id)(self.coset, cid, eid, j)]);
            }
    #[cfg(feature = "nanobench")] let elapsed = start0.elapsed();
    #[cfg(feature = "nanobench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "nanobench")] println!("  ┃         vrf ldt dgr {:>8} mt sz {:>2}:  {:9.3} ms @────", degree, mtos[k].leaf_degree, in_ms as f64 / 1000.0);
            degree >>= self.coset.log_layer_ratio;
            degree = if degree_parity == 0 { degree } else { degree + 1 };
        }
    #[cfg(feature = "nanobench")] let start0 = Instant::now();

        let coeffs = (self.coset.ifft_small)(self.coset, &final_val, ceil_log(final_val.len()));
        for i in degree..coeffs.len() {
            assert_eq!(coeffs[i], T::from(0));
        }
    #[cfg(feature = "nanobench")] let elapsed = start0.elapsed();
    #[cfg(feature = "nanobench")] let in_ms = elapsed.as_micros() as usize;
    #[cfg(feature = "nanobench")] println!("  ┃         vrf ifft of degree({:>3}) poly {:>2}:{:9.3} ms @────", degree, ceil_log(final_val.len()), in_ms as f64 / 1000.0);

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
vpd_vrf_inc(in_ms);

        (query_pts, mtos)
    }
}
