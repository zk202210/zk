use std::sync::Arc;

use super::*;
use crate::util::*;
use crate::merkle::MT;

// a group of additive cosets (log_len layers)
// layer 0   basis: 1, 2, 4, 8, ..., (1 << (log_len - 1))
// ...
// layer i   basis: β₁, β₂, β₃, ...
// layer i+1 basis: β₂(β₁+β₂), β₃(β₁+β₃), ...
// basis0 stores the first basis of each layer
// element[i] stores all of the elements of the layer-i coset
#[derive(Clone)]
pub struct Coset<T: Field> {
    pub log_len: usize,
    pub log_rho: usize,
    pub log_l: usize,
    pub log_layer_ratio: usize,
    basis0: Vec<T>,
    basis0_inv: Vec<T>,
    p: Vec<T>,
    p_inv: Vec<T>,
    pub pow_2_inv: Vec<T>,
    pub element: Vec<Vec<T>>,
    pub pack_element: Vec<Vec<T::P4>>,
    pub pack_element_inv: Vec<Vec<T::P4>>,
    pub fft_coeff: Vec<Vec<T>>,
    pub fft_coeff_small: Vec<Vec<Vec<T>>>,
    // scalar multiplier needed in polynomial degree shift
    // in the basis defined in https://arxiv.org/pdf/1503.05761.pdf
    // we only need thress shift:
    //   n~2n -> 0~n
    //  2n~3n -> 0~n
    //  2n~3n -> n~2n
    //  the last two use the same scalar
    // where n := log_len and assuming kappa < n so that 2n+kappa < 3n
    fft_basis_shift_scalar: [T; 5],
    // pub vanishing_coeff: Vec<T>,
    pub vanishing_eval: Vec<T>,
    // pub vanishing_rev_inv: Vec<(usize, T)>,
    // pts representation of monomials
    pub x_n_k: Vec<T>,
    pub x_n: Vec<T>,
    pub x_n_1: Vec<T>,
    pub ksi: T,
    pub fft: Arc<dyn Fn(&Coset<T>, &Vec<T>, usize, usize) -> Vec<T> + Sync + Send>,
    pub ifft: Arc<dyn Fn(&Coset<T>, &Vec<T>, usize) -> Vec<T> + Sync + Send>,
    pub ifft_small: Arc<dyn Fn(&Coset<T>, &Vec<T>, usize) -> Vec<T> + Sync + Send>,
    pub interpolate: Arc<dyn Fn(&Coset<T>, &Vec<T>, usize, usize) -> Vec<T> + Sync + Send>,
    pub poly_over_vanishing: Arc<dyn Fn(&Coset<T>, &Vec<T>) -> Vec<T> + Sync + Send>,
    pub poly_over_vanishing_w_remainder: Arc<dyn Fn(&Coset<T>, &Vec<T>) -> (Vec<T>, Vec<T>) + Sync + Send>,
    pub next_layer: Arc<dyn Fn(&Coset<T>, usize, usize, &Vec<T>, T, usize) -> T + Sync + Send>,
    pub next_layer_final: Arc<dyn Fn(&Coset<T>, usize, usize, &Vec<T>, T, usize) -> T + Sync + Send>,
    pub rscode_commit: Arc<dyn Fn(&Coset<T>, &mut MT<T>, &Vec<T>, usize) + Sync + Send>,
    pub rscode_commit_next_layer: Arc<dyn Fn(&Coset<T>, &mut MT<T>, &mut Vec<T>, T, usize, usize) + Sync + Send>,
    pub subcoset_id: Arc<dyn Fn(&Coset<T>, usize) -> (usize, usize) + Sync + Send>,
    pub subcoset_id_next_layer: Arc<dyn Fn(&Coset<T>, usize, usize) -> (usize, usize) + Sync + Send>,
    pub normal_id: Arc<dyn Fn(&Coset<T>, usize, usize, usize) -> usize + Sync + Send>,
}

fn fft_vanishing<T: Field>(basis0: &Vec<T>, basis_len: usize, x: &T) -> T {
    let mut res = *x;
    for i in 0..basis_len {
        res *= res + basis0[i];
    }
    res
}

// fn fft_slow_coeff<T: Field>(log_l: usize) -> (Vec<Vec<T>>, Vec<Vec<T>>) {
//     let m = log_l;
//     let n = 1 << m;
// 
//     let mut recursed_betas: Vec<T> = Vec::with_capacity((m + 1) * m / 2);
//     let mut betas2: Vec<T> = Vec::with_capacity(m);
//     let mut betas: Vec<Vec<T>> = Vec::with_capacity(m);
//     let mut idx = 0usize;
// 
//     for i in 0..m {
//         betas2.push(T::from(1 << i));
//     }
// 
//     for j in 0..m {
//         betas.push(Vec::with_capacity(n >> j));
//         let beta = betas2[m - 1 - j];
//         let betainv = beta.inv();
//         let mut betai = T::from(1);
//         for _ in (0..n).step_by(1 << j) {
//             betas[j].push(betai);
//             betai *= beta;
//         }
//         assert_eq!(betas[j].len(), n >> j);
// 
//         for i in 0..(m - 1 - j) {
//             let newbeta = betas2[i] * betainv;
//             recursed_betas.push(newbeta);
//             betas2[i] = newbeta * (newbeta - T::from(1));
//             idx += 1;
//         }
//     }
// 
//     let mut subset_sums: Vec<Vec<T>> = Vec::with_capacity(m);
//     for j in 0..m {
//         idx -= j;
//         let mut sums: Vec<T> = vec![T::from(0); 1 << j];
//         for i in 0..j {
//             let m = 1 << i;
//             for k in 0..m {
//                 sums[m + k] = sums[k] + recursed_betas[idx + i];
//             }
//         }
//         subset_sums.push(sums);
//     }
// 
//     (subset_sums, betas)
// }


impl<T: Field> Coset<T> {
    // (1 << log_rho) will be the RS code rate
    pub fn init(log_l: usize, log_rho: usize, log_layer_ratio: usize) -> Self {
        if T::CHAR == 2 {
            Self::additive_coset_init(log_l, log_rho, log_layer_ratio)
        } else {
            Self::multiplicative_coset_init(log_l, log_rho, log_layer_ratio)
        }
    }
    pub fn additive_coset_init(log_l: usize, log_rho: usize, log_layer_ratio: usize) -> Self {
        let log_len = log_l - log_rho;
        let mut basis0 = vec![T::from(0); log_l];
        let mut basis0_inv = vec![T::from(0); log_l];
        let mut element = vec![vec![]; log_l];
        let mut basis = vec![T::from(0); log_l];

        let mut fft_coeff: Vec<Vec<T>> = Vec::with_capacity(log_l);

        for i in 0..(1usize << log_l) {
            element[0].push(T::from(i));
        }
        for i in 0..log_l {
            basis[i] = T::from(1usize << i);
        }
        basis0[0] = basis[0];
        basis0_inv[0] = basis[0].inv();

        for j in 0..(log_l - 1) {
            element[j + 1] = vec![T::from(0); 1usize << (log_l - j - 1)];
            for i in 0..(log_l - j - 1) {
                basis[i] = (basis0[j] + basis[i + 1]) * basis[i + 1];
            }
            for i in 0..(1usize << (log_l - j - 1)) {
                element[j + 1][i] = element[j][2 * i] * (
                    element[j][2 * i] + basis0[j]);
            }
            basis0[j + 1] = basis[0];
            basis0_inv[j + 1] = basis[0].inv();
        }

        // Remark: fft_coeff is xor homomophism
        let mut p = Vec::with_capacity(log_l);
        let mut p_inv = Vec::with_capacity(log_l);
        for i in 0..log_l {
            let denom = fft_vanishing(&basis0, i, &T::from(1usize << i));
            let denom_inv = denom.inv();
            p.push(denom);
            p_inv.push(denom_inv);

            let log_min = log_len.min(i + 1);
            let log_coeff_len = log_l - log_min;
            fft_coeff.push(Vec::with_capacity(1 << log_coeff_len));
            let mut fft_coeff_base: Vec<T> = Vec::with_capacity(log_coeff_len);
            for j in log_min..log_l {
                fft_coeff_base.push(fft_vanishing(&basis0, i, &T::from(1 << j)) * denom_inv);
            }
            fft_coeff[i].push(T::from(0));
            for j in log_min..log_l {
                for k in 0..(1 << (j - log_min)) {
                    let tmp = fft_coeff[i][k] + fft_coeff_base[j - log_min];
                    fft_coeff[i].push(tmp);
                }
            }
        }
        let mut p_n_1 = T::from(1);
        for i in 0..log_len {
            p_n_1 *= p[i];
        }

        // let mut fft_basis_conv_p = Vec::with_capacity(1 << log_l);
        // let mut fft_basis_conv_p_inv = Vec::with_capacity(1 << log_l);
        // fft_basis_conv_p.push(T::from(1));
        // fft_basis_conv_p_inv.push(T::from(1));
        // for i in 0..log_l {
        //     for k in 0..(1 << i) {
        //         let tmp = fft_basis_conv_p[k] * p[i];
        //         fft_basis_conv_p.push(tmp);
        //         let tmp = fft_basis_conv_p_inv[k] * p_inv[i];
        //         fft_basis_conv_p_inv.push(tmp);
        //     }
        // }
        let fft_basis_shift_scalar =
            // TODO: weird for small log_rho
            if log_rho < 2 {
                [p_inv[log_len],
                 T::from(0), T::from(0), T::from(0), T::from(0)]
            } else if log_rho < 4 {
                [p_inv[log_len],
                 p_inv[log_len + 1] * p[log_len],
                 p_inv[log_len + 2] * p[log_len] * p[log_len + 1],
                 p_inv[log_len + 3] * p[log_len] * p[log_len + 1] * p[log_len + 2],
                 p_inv[log_len + 4] * p[log_len] * p[log_len + 1] * p[log_len + 2] * p[log_len + 3]]
            } else {
                [p_inv[log_len],
                 p_inv[log_len + 1] * p[log_len],
                 p_inv[log_len + 2] * p[log_len] * p[log_len + 1],
                 p_inv[log_len + 3] * p[log_len] * p[log_len + 1] * p[log_len + 2],
                 T::from(0)]
            };

        let mut res = Self {
            log_len: log_len,
            log_rho: log_rho,
            log_l: log_l,
            log_layer_ratio: log_layer_ratio,
            basis0,
            basis0_inv,
            element,
            pack_element: vec![],
            pack_element_inv: vec![],
            pow_2_inv: vec![],
            p: vec![],
            p_inv: vec![],
            fft_coeff: fft_coeff,
            fft_coeff_small: vec![],
            fft_basis_shift_scalar: fft_basis_shift_scalar,
            // vanishing_coeff: vec![],
            vanishing_eval: vec![T::from(0); 1 << log_l],
            // vanishing_rev_inv: vec![],
            x_n_k: vec![],
            x_n: vec![],
            x_n_1: vec![],
            ksi: p_n_1,
            fft: Arc::new(Self::additive_fft),
            ifft: Arc::new(Self::additive_ifft),
            ifft_small: Arc::new(Self::additive_ifft_small),
            interpolate: Arc::new(Self::additive_interpolate),
            poly_over_vanishing: Arc::new(Self::poly_over_additive_vanishing),
            poly_over_vanishing_w_remainder: Arc::new(Self::poly_over_additive_vanishing_w_remainder),
            next_layer: Arc::new(move |a, b, c, d, e, f|
                Self::additive_next_layer(a, b, c, d, e, f, log_layer_ratio)),
            next_layer_final: Arc::new(Self::additive_next_layer_final),
            rscode_commit_next_layer: Arc::new(Self::additive_rscode_commit_next_layer),
            rscode_commit: Arc::new(Self::additive_rscode_commit),
            subcoset_id: Arc::new(Self::additive_subcoset_id),
            subcoset_id_next_layer: Arc::new(Self::additive_subcoset_id_next_layer),
            normal_id: Arc::new(Self::additive_normal_id),
        };

        let mut vanishing_coeff_novel = vec![T::from(0); 1 << (log_len + 1)];
        vanishing_coeff_novel[1 << log_len] = p[log_len];
        res.vanishing_eval = res.additive_fft(&vanishing_coeff_novel, log_len + 1, log_l);
        for i in 0..log_len {
            res.ksi *= res.basis0_inv[i];
        }

        res.p = p;
        res.p_inv = p_inv;


        // vanishing_eval[x ^ y] = vanishing_eval[x] + vanishing[y]
        // for i in 0..log_rho {
        //     let mut tmp = T::from(1);
        //     let u = 1 << (log_l - log_rho + i);
        //     let x = T::from(u);
        //     for k in 0..(1 << (log_l - log_rho)) {
        //         tmp *= x - T::from(k);
        //     }
        //     for j in ((1 << i)..(1 << log_rho)).step_by(2 << i) {
        //         for k in j..(j + (1 << i)) {
        //             res.vanishing_eval[k << (log_l - log_rho)] += tmp;
        //         }
        //     }
        // }
        // for i in (1 << (log_l - log_rho))..(1 << log_l) {
        //     res.vanishing_eval[i] = res.vanishing_eval[(i >> (log_l - log_rho)) << (log_l - log_rho)];
        // }

        // res.vanishing_coeff = res.ifft_slow(&res.vanishing_eval, log_l);
        // res.vanishing_coeff.truncate(1 + (1 << (log_l - log_rho)));

        // res.ksi = res.vanishing_coeff[1].inv();
        // res.zeta = res.vanishing_coeff[res.vanishing_coeff.len() - 1].inv();

        res
    }

    pub fn multiplicative_coset_init(log_l: usize, log_rho: usize, log_layer_ratio: usize) -> Self {
        let log_len = log_l - log_rho;

        let mut basis0 = vec![T::from(0); log_l + 1];
        let mut basis0_inv = vec![T::from(0); log_l + 1];
        let mut element = vec![vec![]; log_l + 1];
        // let mut fft_coeff: Vec<Vec<T>> = Vec::with_capacity(log_l);

        basis0[log_l] = T::from(1);
        basis0_inv[log_l] = T::from(1);
        element[log_l] = vec![T::from(1)];
        let mut root = - T::from(1);
        let mut root_inv = - T::from(1);
        basis0[log_l - 1] = root;
        basis0_inv[log_l - 1] = root_inv;
        element[log_l - 1] = vec![T::from(1), root];
        for i in 2..=log_l {
            root = T::sqrt(root);
            root_inv = T::sqrt(root_inv);
            if root_inv * root != T::from(1) {
                root_inv = - root_inv;
            }
            // assert_eq!(root * root_inv, T::from(1));
            basis0[log_l - i] = root;
            basis0_inv[log_l - i] = root_inv;
        }

        for i in 1..log_l {
            for j in 0..(1 << i) {
                let ele = element[log_l - i][j];
                element[log_l - 1 - i].push(ele);
                element[log_l - 1 - i].push(ele * basis0[log_l - 1 - i]);
            }
        }

        // pow_2_inv[i] = 2^{-i}
        let mut pow_2_inv = Vec::with_capacity(log_l + 1);
        pow_2_inv.push(T::from(1));
        pow_2_inv.push(T::from(2).inv());
        for i in 2..=log_l {
            pow_2_inv.push(pow_2_inv[1] * pow_2_inv[i - 1]);
        }

        // for i in 0..log_l {
        //     assert_eq!(element[i][1], basis0[i]);
        //     assert_eq!(element[i].len(), 1 << (log_l - i));
        //     assert_eq!(element[i][(1 << (log_l - i)) - 1], basis0_inv[i]);
        // }

        let mut res = Self {
            log_len: log_len,
            log_rho: log_rho,
            log_l: log_l,
            log_layer_ratio: log_layer_ratio,
            basis0: basis0,
            basis0_inv: basis0_inv,
            element: element,
            pack_element: vec![],
            pack_element_inv: vec![],
            pow_2_inv: pow_2_inv,
            p: vec![],
            p_inv: vec![],
            fft_coeff: vec![],
            fft_coeff_small: vec![],
            fft_basis_shift_scalar: [T::from(0); 5],
            vanishing_eval: vec![T::from(0); 1 << log_l],
            x_n_k: vec![],
            x_n: vec![],
            x_n_1: vec![],
            ksi: T::from(0),
            fft: Arc::new(Self::multiplicative_fft),
            ifft: Arc::new(Self::multiplicative_ifft),
            ifft_small: Arc::new(Self::multiplicative_ifft),
            interpolate: Arc::new(Self::multiplicative_interpolate),
            poly_over_vanishing: Arc::new(Self::poly_over_multiplicative_vanishing),
            poly_over_vanishing_w_remainder: Arc::new(Self::poly_over_multiplicative_vanishing_w_remainder),
            next_layer: Arc::new(move |a, b, c, d, e, f|
                Self::multiplicative_next_layer(a, b, c, d, e, f, log_layer_ratio)),
            next_layer_final: Arc::new(Self::multiplicative_next_layer_final),
            rscode_commit_next_layer: Arc::new(Self::multiplicative_rscode_commit_next_layer),
            rscode_commit: Arc::new(Self::multiplicative_rscode_commit),
            subcoset_id: Arc::new(Self::multiplicative_subcoset_id),
            subcoset_id_next_layer: Arc::new(Self::multiplicative_subcoset_id_next_layer),
            normal_id: Arc::new(Self::multiplicative_normal_id),
        };

        // precompute pack elements for fft/ifft
        if T::CHAR == (1 << 61) - 1 {
            res.fft = Arc::new(Self::multiplicative_fft_pack4);
            res.ifft = Arc::new(Self::multiplicative_ifft_pack4);
            res.pack_element = Vec::with_capacity(log_l - 1);
            res.pack_element_inv = Vec::with_capacity(log_l - 1);
            for i in 0..(log_l - 1) {
                res.pack_element.push(Vec::with_capacity(1 << (log_l - i - 2)));
                res.pack_element_inv.push(Vec::with_capacity(1 << (log_l - i - 2)));
                for j in 0..(1 << (log_l - i - 2)) {
                    res.pack_element[i].push(T::P4::pack(unsafe {
                        std::mem::transmute(&res.element[i][4 * j])
                    }));
                    let m = 1 << (log_l - i);
                    res.pack_element_inv[i].push(T::P4::pack(&[
                        res.element[i][(m - 4 * j - 0) % m],
                        res.element[i][(m - 4 * j - 1) % m],
                        res.element[i][(m - 4 * j - 2) % m],
                        res.element[i][(m - 4 * j - 3) % m]
                    ]));
                }
            }
        }

        // x^n - 1
        let mut vanishing_coeff = vec![T::from(0); 1 << (log_len + 1)];
        vanishing_coeff[1 << log_len] = T::from(1);
        vanishing_coeff[0] = - T::from(1);
        res.vanishing_eval = res.multiplicative_fft(&vanishing_coeff, log_len + 1, log_l);

        // ksi = 1 / n
        res.ksi = res.pow_2_inv[log_len];


        res
    }

    pub fn precompute(&mut self, n_ldt_iter: usize, mpc_vpd_flag: bool) {
        // self.fft = Arc::new(move |inputs, log_input_size, log_output_size| {
        //     self.add_fft(inputs, log_input_size, log_output_size)
        // });
        // compute pt repr of monomials
        // x^{2n-k}, x^{3n-k}, x^{3n+1}, x^{n - 1}
        // has changed to the corresponding novel basis
        let k = if mpc_vpd_flag { 0 } else { n_ldt_iter << self.log_layer_ratio };
        let n = 1 << self.log_len;

        assert!(2 * n + k - 1 <= (1 << (self.log_len + self.log_rho)));

        let (log_m, m) =
            if 2 * n + k - 1 <= 2 * n {
                (1, 2)
            } else if 2 * n + k - 1 <= 4 * n {
                (2, 4)
            } else if 2 * n + k - 1 <= 8 * n {
                (3, 8)
            } else if 2 * n + k - 1 <= 16 * n {
                (4, 16)
            } else if 2 * n + k - 1 <= 32 * n {
                (5, 32)
            } else {
                // TODO
                unreachable!()
            };

        // only init these monomial when using Virgo's proposed VPD protocol
        if !mpc_vpd_flag {
            let mut monomial = vec![T::from(0); m * n];
            monomial[n] = T::from(1);
            self.x_n = (self.fft)(self, &monomial, self.log_len + log_m, self.log_l);

            monomial[n] = T::from(0);
            monomial[n + k] = T::from(1);
            if T::CHAR == 2 {
                self.x_n_k = (self.fft)(self, &monomial, self.log_len + log_m, self.log_l);
            } else {
                // add a non-zero const term to lifting polynomial x_n_k
                // so that g(x)/x * x_n_k is low-degree iff const term of g(x) is zero
                monomial[0] = T::from(1);
                self.x_n_k = (self.fft)(self, &monomial, self.log_len + log_m, self.log_l);
            }
        }

        if T::CHAR == 2 {
            // self.x_n_1 = x^{n-1}
            let mut monomial = vec![T::from(0); n];
            monomial[n - 1] = T::from(1);
            self.x_n_1 = (self.fft)(self, &monomial, self.log_len, self.log_l);
        } else {
            // self.x_n_1 = x^{-1}
            let n = 1 << self.log_l;
            for i in 0..n {
                self.x_n_1.push(self.element[0][(n - i) % n]);
            }
        }

        // compute pt repr of Rev(Z)^{-1} mod x^{2n+k-1 - n + 1}
        // let mut f = self.vanishing_coeff.to_vec();
        // f.reverse();
        // self.vanishing_rev_inv = poly_mod_inv_sparse(&f, n + k);

        // let inv = poly_mod_inv(&f, n + k, self);
        // for i in 1..inv.len() {
        //     if inv[i] != T::from(0) {
        //         self.vanishing_rev_inv.push((i, inv[i]));
        //     }
        // }
        //
        // let log_len = ceil_log(2 * n + k - 1 + n + k - 1 + 1);
        // self.vanishing_rev_inv = self.fft_slow(&inv, log_len);

        // let mut b = vec![T::from(0); n + k - 1];
        // b[0] = T::from(1);

        // println!("f len {}", f.len());
        // println!("inv len {}", inv.len());
        // let mut p = poly_mul(&f, &inv);
        // p.truncate(n + k - 1);
        // assert_eq!(p, b);


        // let f_pts = self.fft_slow(&f, self.log_len + 2);
        // let mut a_pts = self.vanishing_rev_inv.to_vec();
        // for i in 0..(4 << self.log_len) {
        //     a_pts[i] *= f_pts[i];
        // }
        // let mut a = self.ifft_slow(&a_pts, self.log_len + 2);
        // a.truncate(n + k - 1);
        // assert_eq!(a, b);

        if T::CHAR == 2 {
            // each ifft/fft is for a perticular coset
            // below generates coeffs for small size (~kappa * ratio) ifft
            // assert max_log_len < self.log_len
            let max_log_len = ceil_log(n_ldt_iter) + self.log_layer_ratio;
            for log_len in self.log_layer_ratio..=max_log_len {
                let mut fft_coeff: Vec<Vec<T>> = Vec::with_capacity(log_len);
                let basis0 = self.basis0[(self.log_l - log_len)..].to_vec();
                let element = &self.element[self.log_l - log_len];
                for i in 0..log_len {
                    let denom = fft_vanishing(&basis0, i, &element[1 << i]);
                    let denom_inv = denom.inv();
                    let log_coeff_len = log_len - i - 1;
                    fft_coeff.push(Vec::with_capacity(1 << log_coeff_len));
                    let mut fft_coeff_base: Vec<T> = Vec::with_capacity(log_coeff_len);
                    for j in (i + 1)..log_len {
                        fft_coeff_base.push(fft_vanishing(&basis0, i, &element[1 << j]) * denom_inv);
                    }
                    fft_coeff[i].push(T::from(0));
                    for j in (i + 1)..log_len {
                        for k in 0..(1 << (j - (i + 1))) {
                            let tmp = fft_coeff[i][k] + fft_coeff_base[j - (i + 1)];
                            fft_coeff[i].push(tmp);
                        }
                    }
                }
                self.fft_coeff_small.push(fft_coeff);
            }
        }
    }

    pub fn multiplicative_normal_id(&self, cid: usize, eid: usize, lvl: usize) -> usize {
        cid + (eid << (self.log_l - lvl - 2 * self.log_layer_ratio))
    }

    pub fn multiplicative_subcoset_id(&self, pts: usize) -> (usize, usize) {
        let subcoset_mask = (1 << (self.log_l - self.log_layer_ratio)) - 1;
        (pts & subcoset_mask, pts >> (self.log_l - self.log_layer_ratio))
    }

    pub fn multiplicative_subcoset_id_next_layer(&self, pts: usize, lvl: usize) -> (usize, usize) {
        let subcoset_mask = (1 << (self.log_l - lvl - 2 * self.log_layer_ratio)) - 1;
        (pts & subcoset_mask, pts >> (self.log_l - lvl - 2 * self.log_layer_ratio))
    }

    pub fn multiplicative_next_layer_final(&self, lvl: usize, cid: usize, subcoset_values: &Vec<T>, r: T, degree: usize) -> T {
        let log_ratio = ceil_log(degree);
        let gap = self.log_layer_ratio - log_ratio;
        let mut res = vec![T::from(0); 1 << gap];
        let ratio = 1 << log_ratio;
        for i in 0..res.len() {
            let mut tmp = Vec::with_capacity(ratio);
            for j in (i..(1 << self.log_layer_ratio)).step_by(1 << gap) {
                tmp.push(subcoset_values[j]);
            }
            res[i] = self.multiplicative_next_layer(
                lvl, cid + (i << (self.log_l - lvl - self.log_layer_ratio)),
                &tmp, r, degree, log_ratio);

            // additional verifier checks, much better to output another bool
            if i >= 1 {
                assert_eq!(res[i], res[i - 1]);
            }
        }
        res[0]
    }

    pub fn multiplicative_next_layer(&self, lvl: usize, cid: usize, subcoset_values: &Vec<T>, alpha: T, parity: usize, log_ratio: usize) -> T {
        let mut rscode = subcoset_values.to_vec();
        let mut parity = parity;
        let mut len = 1 << log_ratio;
        let log_stride = self.log_l - lvl - log_ratio;
        for k in 0..log_ratio {
            for i in 0..(len / 2) {
                let fu = rscode[i];
                let fv = rscode[i + len / 2];
                let y = self.element[lvl + k + 1][cid + (i << log_stride)];
                let alpha = if (parity & 1) == 0 { alpha } else { alpha * y };
                // self.element[lvl + k][cid + (i << log_stride)].inv())
                let u_inv = self.element[lvl + k][((1 << (self.log_l - lvl - k)) - cid - (i << log_stride)) % (1 << (self.log_l - lvl - k))];
                let fy = ((fu + fv) + (fu - fv) * alpha * u_inv) * self.pow_2_inv[1];
                rscode[i] = fy;
            }
            len >>= 1;
            parity = (parity >> 1) + (parity & 1);
        }
        return rscode[0];
    }

    pub fn multiplicative_rscode_commit(&self, mt: &mut MT<T>, rscode: &Vec<T>, _lvl: usize) {
        let stride = 1 << (self.log_l - self.log_layer_ratio);
        for i in 0..mt.size {
            mt.set_leaf_from_slice(i, &rscode[i..], stride, 1 << self.log_layer_ratio);
        }
        mt.commit();
    }

    pub fn multiplicative_rscode_commit_next_layer(&self, mt: &mut MT<T>, rscode: &mut Vec<T>, alpha: T, lvl: usize, parity: usize) {
        let mut parity = parity;
        // assert_eq!(rscode.len(), mt.size << (2 * self.log_layer_ratio));
        for k in 0..(self.log_layer_ratio - 1) {
            let codelen = mt.size << (2 * self.log_layer_ratio - k - 1);
            for i in 0..codelen {
                let fu = rscode[i];
                let fv = rscode[i + codelen];
                // assert_eq!(self.element[lvl + k][i] + self.element[lvl + k][i + codelen], T::from(0));
                let y = self.element[lvl + k + 1][i];
                let alpha = if (parity & 1) == 0 { alpha } else { alpha * y };
                let u_inv = self.element[lvl + k][((1 << (self.log_l - lvl - k)) - i) % (1 << (self.log_l - lvl - k))];
                // self.element[lvl + k][i].inv()
                let fy = ((fu + fv) + (fu - fv) * alpha * u_inv) * self.pow_2_inv[1];
                rscode[i] = fy;
            }
            parity = (parity >> 1) + (parity & 1);
        }
        let k = self.log_layer_ratio - 1;

        let ratio = 1 << self.log_layer_ratio;
        let stride = mt.size;
        let codelen = mt.size << (2 * self.log_layer_ratio - k - 1);
        for i in 0..mt.size {
            let mut leaf = vec![T::from(0); ratio];
            for j in 0..ratio {
                let fu = rscode[i + stride * j];
                let fv = rscode[i + stride * j + codelen];
                // assert_eq!(self.element[lvl + k][i + stride * j] +
                //            self.element[lvl + k][i + stride * j + codelen],
                //            T::from(0));
                // assert_eq!(self.element[lvl + k][i + stride * j] * self.element[lvl + k][i + stride * j],
                //            self.element[lvl + k + 1][i + stride * j]);
                let y = self.element[lvl + k + 1][i + stride * j];
                let alpha = if (parity & 1) == 0 { alpha } else { alpha * y };
                // let u_inv = self.element[lvl + k][i + stride * j].inv();
                let u_inv = self.element[lvl + k][((1 << (self.log_l - lvl - k)) - i - stride * j) % (1 << (self.log_l - lvl - k))];
                let fy = ((fu + fv) + (fu - fv) * alpha * u_inv) * self.pow_2_inv[1];
                rscode[i + stride * j] = fy;
                leaf[j] = fy;
            }
            mt.set_leaf(i, &leaf);
        }
        mt.commit();
    }


    pub fn additive_normal_id(&self, cid: usize, eid: usize, _lvl: usize) -> usize {
        (cid << self.log_layer_ratio) + eid
    }

    pub fn additive_subcoset_id(&self, pts: usize) -> (usize, usize) {
        (pts >> self.log_layer_ratio, pts & ((1 << self.log_layer_ratio) - 1))
    }

    pub fn additive_subcoset_id_next_layer(&self, cid: usize, _lvl: usize) -> (usize, usize) {
        (cid >> self.log_layer_ratio, cid & ((1 << self.log_layer_ratio) - 1))
    }

    pub fn additive_next_layer_final(&self, lid: usize, cid: usize, subcoset_values: &Vec<T>, r: T, degree: usize) -> T {
        let log_ratio = ceil_log(degree);
        let gap = self.log_layer_ratio - log_ratio;
        let mut res = vec![T::from(0); 1 << gap];
        let ratio = 1 << log_ratio;
        for i in 0..res.len() {
            res[i] = self.additive_next_layer(lid, (cid << gap) + i, &subcoset_values[(i * ratio)..((i + 1) * ratio)].to_vec(), r, degree, log_ratio);

            // additional verifier checks, much better to output another bool
            if i >= 1 {
                assert_eq!(res[i], res[i - 1]);
            }
        }
        res[0]
    }

    pub fn additive_next_layer(&self, lvl: usize, cid: usize, subcoset_values: &Vec<T>, alpha: T, parity: usize, log_ratio: usize) -> T {
        let mut rscode = subcoset_values.to_vec();
        let mut parity = parity;
        let mut cid = cid << log_ratio;
        let mut len = 1 << log_ratio;
        for k in 0..log_ratio {
            let mut ofs = 0;
            for i in 0..(len / 2) {
                let u = rscode[ofs];
                let v = rscode[ofs + 1];
                let alpha = if (parity & 1) == 0 { alpha } else { alpha * self.element[lvl + k + 1][(cid >> 1) + i] };
                let y = (u + v) * self.basis0_inv[lvl + k] * (alpha + self.element[lvl + k][cid + (i << 1)]) + u;
                rscode[i] = y;
                ofs += 2;
            }
            len >>= 1;
            cid >>= 1;
            parity = (parity >> 1) + (parity & 1);
        }
        return rscode[0];
    }

    pub fn additive_rscode_commit(&self, mt: &mut MT<T>, rscode: &Vec<T>, _lvl: usize) {
        let ratio = 1 << self.log_layer_ratio;
        for i in 0..mt.size {
            mt.set_leaf_from_slice(i, &rscode[(ratio * i)..], 1, ratio);
        }
        mt.commit();
    }

    pub fn additive_rscode_commit_next_layer(&self, mt: &mut MT<T>, rscode: &mut Vec<T>, alpha: T, lvl: usize, parity: usize) {
        let mut parity = parity;
        // assert_eq!(rscode.len(), mt.size << (2 * self.log_layer_ratio));
        for k in 0..(self.log_layer_ratio - 1) {
            let ratio = 2;
            let mut ofs = 0;
            for i in 0..(mt.size << (2 * self.log_layer_ratio - k - 2)) {
                for j in 0..ratio {
                    let u = rscode[ofs];
                    let v = rscode[ofs + 1];
                    let alpha = if (parity & 1) == 0 { alpha } else { alpha * self.element[lvl + k + 1][ratio * i + j] };
                    let y = (u + v) * self.basis0_inv[lvl + k] * (alpha + self.element[lvl + k][ofs]) + u;
                    // assert_eq!(ofs/2, ratio * i + j);
                    rscode[ratio * i + j] = y;
                    ofs += ratio;
                }
            }
            parity = (parity >> 1) + (parity & 1);
        }
        let k = self.log_layer_ratio - 1;

        let ratio = 1 << self.log_layer_ratio;
        let mut ofs = 0;
        for i in 0..mt.size {
            let mut leaf = vec![T::from(0); ratio];
            for j in 0..ratio {
                let u = rscode[ofs];
                let v = rscode[ofs + 1];
                let alpha = if (parity & 1) == 0 { alpha } else { alpha * self.element[lvl + k + 1][ratio * i + j] };
                let y = (u + v) * self.basis0_inv[lvl + k] * (alpha + self.element[lvl + k][ofs]) + u;
                rscode[ratio * i + j] = y;
                leaf[j] = y;
                ofs += 2;
            }
            // TODO change to set_leaf_from_slice and omit vect leaf
            mt.set_leaf(i, &leaf);
        }
        mt.commit();
    }

    // O(n^2) division
    // pub fn poly_over_vanish(&self, inputs: &Vec<T>) -> (Vec<T>, Vec<T>) {
    //     let deg_z = self.vanishing_coeff.len() - 1;
    //     let mut q = inputs[deg_z..].to_vec();
    //     let mut r = inputs[0..deg_z].to_vec();
    //     let linv = self.zeta;
    //     for i in (0..q.len()).rev() {
    //         let twist = q[i] * linv;
    //         q[i] = twist;
    //         for j in (0..(deg_z)).rev() {
    //             if i + j >= deg_z {
    //                 q[i + j - deg_z] = q[i + j - deg_z] - twist * self.vanishing_coeff[j];
    //             } else {
    //                 r[i + j] = r[i + j] - twist * self.vanishing_coeff[j];
    //             }
    //         }
    //     }
    //     // assert_eq!(*inputs, poly_add(&poly_mul(&q, &self.vanishing_coeff), &r));
    //     (q, r)
    // }

    // fast division : div and conq
    // ~ 10^3 ms
    // pub fn poly_over_vanish(&self, inputs: &Vec<T>) -> Vec<T> {
    //     poly_div_with_leading_inv(inputs, &self.vanishing_coeff, &self.zeta, self)
    // }

    // faster division: using inverse of reverse of Z
    // ~ 10^2 ms
    // ref: ch2 http://people.csail.mit.edu/madhu/ST12/scribe/lect06.pdf
    // one mistake in the lecture notes: deg(f) - deg(g) should be ... + 1
    // input: polynomial of degree 2n + k - 1 over vanishing polynomial of degree n
    // pub fn poly_over_vanish(&self, inputs: &Vec<T>) -> Vec<T> {
    //     let n = 1 << self.log_len;
    //     let k = inputs.len() - 2 * n;

    //     let log_len = ceil_log(2 * n + k - 1 + n + k - 1 + 1);

    //     let mut f = inputs.to_vec();
    //     f.reverse();

    //     // let mut q = poly_mul(&f, &self.vanishing_rev_inv);
    //     let mut pts = self.fft_slow(&f, log_len);
    //     for i in 0..(1 << log_len) {
    //         pts[i] *= self.vanishing_rev_inv[i];
    //     }
    //     let mut q = self.ifft_slow(&pts, log_len);

    //     q.truncate(n + k);
    //     q.reverse();
    //     q
    // }

    // fastest division: since the inv rev Z is very sparse
    // ~0.x ms
    // with log|H| + 2 non-zero coeff for 7 < log|H| (can be proved?)
    // a direct partial multiplication can be super fast
    // pub fn poly_monomial_over_vanish(&self, inputs: &Vec<T>) -> Vec<T> {
    //     let n = 1 << self.log_len;
    //     let k = inputs.len() - 2 * n;

    //     let mut q = inputs[(2 * n + k - (n + k))..(2 * n + k)].to_vec();

    //     for i in 0..self.vanishing_rev_inv.len() {
    //         let (i, coeff) = self.vanishing_rev_inv[i];
    //         if i > n + k {
    //             break;
    //         }
    //         for j in i..(n + k) {
    //             // q[j] += coeff * f[j - i]
    //             // but reversed
    //             q[n + k - 1 - j] += coeff * inputs[2 * n + k - 1 - (j - i)];
    //         }
    //     }
    //     q
    // }

    // ref: http://www.math.clemson.edu/~sgao/papers/GM10.pdf
    // https://github.com/scipr-lab/libiop/blob/master/libiop/algebra/fft.tcc
    /*
    pub fn fft_slow(&self, inputs: &Vec<T>, log_output_size: usize) -> Vec<T> {
        self.fft_slow_offset(inputs, log_output_size, inputs.len(), 0)
    }

    pub fn fft_slow_offset(&self, inputs: &Vec<T>, log_output_size: usize, input_len: usize, offset: usize) -> Vec<T> {
        let m = log_output_size;
        let n = 1 << log_output_size;
        let mut res: Vec<T> = vec![T::from(0); n];

        let mut recursed_betas: Vec<T> = vec![T::from(0); (m + 1) * m / 2];
        let mut recursed_betas_ptr = 0usize;
        // std::vector<FieldT> recursed_shifts(m, FieldT(0));

        let mut betas2: Vec<T> = vec![T::from(0); log_output_size];
        // FieldT shift2 = domain.shift();

        for i in 0..input_len {
            res[i + offset] = inputs[i];
        }

        for i in 0..log_output_size {
            betas2[i] = T::from(1 << i);
        }

        for j in 0..m {
            let beta = betas2[m - 1 - j];
            let mut betai = T::from(1);

            /* twist by beta. TODO: this can often be elided by a careful choice of betas */
            for ofs in (0..n).step_by(1 << j) {
                for p in 0..(1 << j) {
                    res[ofs + p] = res[ofs + p] * betai;
                    // res[ofs + p] = res[ofs + p] * self.fft_betas[j][ofs >> j];
                }
                betai *= beta;
            }

            /* perform radix conversion */
            let mut stride = n / 4;
            while stride >= (1 << j) {
                for ofs in (0..n).step_by(stride * 4) {
                    for i in 0..stride {
                        res[ofs + 2 * stride + i] = res[ofs + 2 * stride + i] + res[ofs + 3 * stride + i];
                        res[ofs + 1 * stride + i] = res[ofs + 1 * stride + i] + res[ofs + 2 * stride + i];
                    }
                }
                stride >>= 1;
            }

            /* compute deltas used in the reverse process */
            let betainv = beta.inv();
            for i in 0..(m - 1 - j) {
                let newbeta = betas2[i] * betainv;
                recursed_betas[recursed_betas_ptr] = newbeta;
                betas2[i] = newbeta * newbeta - newbeta;
                recursed_betas_ptr += 1;
            }

            // FieldT newshift = shift2 * betainv;
            // recursed_shifts[j] = newshift;
            // shift2 = newshift.squared() - newshift;
        }

        // bitreverse_vector<FieldT>(S);
        for k in 0..n {
            // use bit_reverse::ParallelReverse;
            let rk = swap_bits(k, m);
            if k < rk {
                res.swap(k, rk);
            }
        }

        /* unwind the recursion */
        for j in 0..m {
            // recursed_betas_ptr -= j;
            /* note that this devolves to empty range for the first loop iteration */
            // let popped_betas = &recursed_betas[recursed_betas_ptr..(recursed_betas_ptr + j)];
            // std::vector<FieldT> popped_betas = std::vector<FieldT>(recursed_betas.begin()+recursed_betas_ptr, recursed_betas.begin()+recursed_betas_ptr+j);
            // const FieldT popped_shift = recursed_shifts[m-1-j];
            // std::vector<FieldT> sums = all_subset_sums<FieldT>(popped_betas, popped_shift);
            let mut sums: Vec<T> = vec![T::from(0); 1 << j];
            for i in 0..j {
                let m = 1 << i;
                for k in 0..m {
                    sums[m + k] = sums[k] + recursed_betas[recursed_betas_ptr + i];
                }
            }

            let stride = 1 << j;
            for ofs in (0..n).step_by(2 * stride) {
                for i in 0..stride {
                    res[ofs + i] = res[ofs + i] + res[ofs + stride + i] * sums[i];
                    // res[ofs + i] = res[ofs + i] + res[ofs + stride + i] * self.fft_subset_sums[j][i];
                    res[ofs + stride + i] = res[ofs + stride + i] + res[ofs + i];
                }
            }
        }
        // assert_eq!(recursed_betas_ptr, 0);

        res
    }

    pub fn ifft_slow(&self, inputs: &Vec<T>, log_size: usize) -> Vec<T> {
        let m = log_size;
        let n = 1 << m;

        let mut res: Vec<T> = inputs.to_vec();
        res.truncate(1 << log_size);
        let mut recursed_twists: Vec<T> = vec![T::from(0); m];
        let mut betas2: Vec<T> = vec![T::from(0); m];
        for i in 0..m {
            betas2[i] = T::from(1 << i);
        }

        for j in 0..m {
            let beta = betas2[m - 1 - j];
            let betai = beta.inv();
            recursed_twists[j] = betai;

            let mut newbetas: Vec<T> = vec![T::from(0); m - 1 - j];

            for i in 0..(m - 1 - j) {
                let newbeta = betas2[i] * betai;
                newbetas[i] = newbeta;
                betas2[i] = newbeta * newbeta - newbeta;
            }

            // const std::vector<FieldT> sums = all_subset_sums<FieldT>(newbetas, newshift);
            let mut sums: Vec<T> = vec![T::from(0); 1 << (m - 1 - j)];
            for i in 0..(m - 1 - j) {
                let m = 1 << i;
                for k in 0..m {
                    sums[m + k] = sums[k] + newbetas[i];
                }
            }

            let half = 1 << (m - 1 - j);
            for ofs in (0..n).step_by(2 * half) {
                for p in 0..half {
                    res[ofs + half + p] = res[ofs + half + p] + res[ofs + p];
                    res[ofs + p] = res[ofs + p] + res[ofs + half + p] * sums[p];
                }
            }
        }

        // bitreverse_vector<FieldT>(S);
        for k in 0..n {
            let rk = swap_bits(k, m);
            if k < rk {
                res.swap(k, rk);
            }
        }

        for j in 0..m {
            let mut nn = 4 << (m - 1 - j);
            /* perform radix combinations */
            while nn <= n {
                let quarter = nn / 4;
                for ofs in (0..n).step_by(nn) {
                    for i in 0..quarter {
                        res[ofs + 1 * quarter + i] = res[ofs + 1 * quarter + i] + res[ofs + 2 * quarter + i];
                        res[ofs + 2 * quarter + i] = res[ofs + 2 * quarter + i] + res[ofs + 3 * quarter + i];
                    }
                }
                nn *= 2;
            }

            /* twist by \beta^{-1} */
            let betainv = recursed_twists[m - 1 - j];
            let mut betainvi = T::from(1);
            for ofs in (0..n).step_by(1 << (m - 1 - j)) {
                for p in 0..(1 << (m - 1 - j)) {
                    res[ofs + p] = res[ofs + p] * betainvi;
                }
                betainvi *= betainv;
            }
        }

        res
    }
    */

    // pub fn poly_over_vanishing(&self, inputs: &Vec<T>) -> (Vec<T>, Vec<T>) {
    pub fn poly_over_additive_vanishing(&self, inputs: &Vec<T>) -> Vec<T> {
        // let reminder = inputs[0..(1 << self.log_len)].to_vec();
        let mut quotient: Vec<T> = Vec::with_capacity(inputs.len() - (1 << self.log_len));

        let mut k = 1;
        loop {
            // inspired by ref https://arxiv.org/pdf/1503.05761.pdf
            for i in ((2 * k - 1) << self.log_len)..inputs.len().min((2 * k) << self.log_len) {
                quotient.push(inputs[i] * self.fft_basis_shift_scalar[0]);
            }

            if inputs.len() <= ((2 * k) << self.log_len) {
                return quotient;
            }

            let mut tmp = 2 * k;
            let mut shift = 0;
            while (tmp & 1) != 1 {
                shift += 1;
                tmp >>= 1;
            }

            for i in ((2 * k) << self.log_len)..inputs.len().min((2 * k + 1) << self.log_len) {
                let coeff = inputs[i] * self.fft_basis_shift_scalar[shift];
                quotient.push(coeff);
                for j in 0..((1 << shift) - 1) {
                    // move j block backward,
                    // current have (2 * k - 1) blocks
                    // quotient[i - ((2 * k) << self.log_len) + ((2 * k - 2 - j) << self.log_len)] += coeff;
                    quotient[i - ((j + 2) << self.log_len)] += coeff;
                }
            }

            if inputs.len() <= ((2 * k + 1) << self.log_len) {
                return quotient;
            }
            k += 1;
        }
    }

    pub fn poly_over_additive_vanishing_w_remainder(&self, inputs: &Vec<T>) -> (Vec<T>, Vec<T>) {
        // let reminder = inputs[0..(1 << self.log_len)].to_vec();
        let mut quotient: Vec<T> = Vec::with_capacity(inputs.len() - (1 << self.log_len));
        let remainder = inputs[0..(1 << self.log_len)].to_vec();

        let mut k = 1;
        loop {
            // inspired by ref https://arxiv.org/pdf/1503.05761.pdf
            for i in ((2 * k - 1) << self.log_len)..inputs.len().min((2 * k) << self.log_len) {
                quotient.push(inputs[i] * self.fft_basis_shift_scalar[0]);
            }

            if inputs.len() <= ((2 * k) << self.log_len) {
                return (quotient, remainder);
            }

            let mut tmp = 2 * k;
            let mut shift = 0;
            while (tmp & 1) != 1 {
                shift += 1;
                tmp >>= 1;
            }

            for i in ((2 * k) << self.log_len)..inputs.len().min((2 * k + 1) << self.log_len) {
                let coeff = inputs[i] * self.fft_basis_shift_scalar[shift];
                quotient.push(coeff);
                for j in 0..((1 << shift) - 1) {
                    // move j block backward,
                    // current have (2 * k - 1) blocks
                    // quotient[i - ((2 * k) << self.log_len) + ((2 * k - 2 - j) << self.log_len)] += coeff;
                    quotient[i - ((j + 2) << self.log_len)] += coeff;
                }
            }

            if inputs.len() <= ((2 * k + 1) << self.log_len) {
                return (quotient, remainder);
            }
            k += 1;
        }
    }


    pub fn poly_over_multiplicative_vanishing(&self, inputs: &Vec<T>) -> Vec<T> {
        // vanishing polynomial is x^n - 1
        let mut quotient = inputs[(1 << self.log_len)..].to_vec();
        let k = quotient.len() >> self.log_len;
        for i in (k << self.log_len)..quotient.len() {
            let tmp = quotient[i];
            quotient[i - (1 << self.log_len)] += tmp;
        }
        for j in 2..=k {
            for i in 0..(1 << self.log_len) {
                let tmp = quotient[i + ((k - j + 1) << self.log_len)];
                quotient[i + ((k - j) << self.log_len)] += tmp;
            }
        }
        quotient
    }

    pub fn poly_over_multiplicative_vanishing_w_remainder(&self, inputs: &Vec<T>) -> (Vec<T>, Vec<T>) {
        // vanishing polynomial is x^n - 1
        let mut quotient = inputs[(1 << self.log_len)..].to_vec();
        let mut remainder = inputs[0..(1 << self.log_len)].to_vec();
        let k = quotient.len() >> self.log_len;
        for i in (k << self.log_len)..quotient.len() {
            let tmp = quotient[i];
            quotient[i - (1 << self.log_len)] += tmp;
        }
        for j in 2..=k {
            for i in 0..(1 << self.log_len) {
                let tmp = quotient[i + ((k - j + 1) << self.log_len)];
                quotient[i + ((k - j) << self.log_len)] += tmp;
            }
        }
        for i in 0..quotient.len().min(1 << self.log_len) {
            remainder[i] += quotient[i];
        }
        (quotient, remainder)
    }


    // ref: https://arxiv.org/pdf/1404.3458.pdf
    pub fn additive_ifft_small(&self, inputs: &Vec<T>, log_size: usize) -> Vec<T> {
        let fft_coeff = &self.fft_coeff_small[log_size - self.log_layer_ratio];
        let mut res = inputs.to_vec();
        for i in 0..log_size {
            // for j in 0..(inputs.len() / 2) {
            for j in 0..(1 << (log_size - 1)) {
                let m = j % (1usize << i);
                let c = (j - m) * 2;
                let idx0 = c + m;
                let idx1 = idx0 + (1usize << i);
                let tmp = res[idx0];
                res[idx1] += tmp;
                let tmp = fft_coeff[i][c >> (i + 1)] * res[idx1];
                res[idx0] += tmp;
            }
        }
        res
    }

    pub fn additive_ifft(&self, inputs: &Vec<T>, log_size: usize) -> Vec<T> {
        // let mut res = vec![T::from(0); inputs.len()];
        // for i in 0..inputs.len() {
        //     res[i] = inputs[i];
        // }
        let mut res = inputs.to_vec();
        res.truncate(1 << log_size);
        // assert!(log_size <= self.log_l - self.log_rho);
        let log_min = log_size.min(self.log_len);
        // for i in 0..log_size {
        for i in 0..log_min {
            // for j in 0..(inputs.len() / 2) {
            for j in 0..(1 << (log_size - 1)) {
                let m = j % (1usize << i);
                let c = (j - m) * 2;
                let idx0 = c + m;
                let idx1 = idx0 + (1usize << i);
                let tmp = res[idx0];
                res[idx1] += tmp;
                // let tmp = self._fft_factors0[i][c >> (i + 1)] * res[idx1];
                let tmp = self.fft_coeff[i][c >> (i + 1)] * res[idx1];
                res[idx0] += tmp;
            }
        }
        for i in log_min..log_size {
            for j in 0..(1 << (log_size - 1)) {
                let m = j % (1usize << i);
                let c = (j - m) * 2;
                let idx0 = c + m;
                let idx1 = idx0 + (1usize << i);
                let tmp = res[idx0];
                res[idx1] += tmp;
                // let tmp = self._fft_factors0[i][c >> (i + 1)] * res[idx1];
                let tmp = self.fft_coeff[i][c >> self.log_len] * res[idx1];
                res[idx0] += tmp;
            }
        }
        res
    }

    // pub fn fft(&self, inputs: &Vec<T>, log_input_size: usize, log_output_size: usize) -> Vec<T> {
    //     let mut res = vec![T::from(0); 1usize << log_output_size];
    //     for k in 0..(1usize << (log_output_size - log_input_size)) {
    //         let offset = k << log_input_size;
    //         for i in 0..inputs.len() {
    //             res[offset + i] = inputs[i];
    //         }
    //     }

    //     for k in 0..(1usize << (log_output_size - log_input_size)) {
    //         let offset = k << log_input_size;
    //         for i in (self.log_len..log_input_size).rev() {
    //             // let i = log_input_size - 1 - ii;
    //             for j in 0..(1usize << (log_input_size - 1)) {
    //                 let m = j % (1usize << i);
    //                 let c = (j - m) * 2;
    //                 let idx0 = c + m + offset;
    //                 let idx1 = idx0 + (1usize << i);

    //                 let coeff = self.fft_coeff[i][(c ^ offset) >> (i + 1)];

    //                 res[idx0] = res[idx0] + coeff * res[idx1];
    //                 res[idx1] = res[idx1] + res[idx0];
    //             }
    //         }
    //     }
    //     res
    // }

    pub fn additive_fft(&self, inputs: &Vec<T>, log_input_size: usize, log_output_size: usize) -> Vec<T> {
        let k_shift = log_input_size - self.log_len;

        let mut res = vec![T::from(0); 1usize << log_output_size];

        for k in 0..(1usize << (log_output_size - log_input_size)) {
            let offset = k << log_input_size;

            if self.log_len < log_input_size {
                {
                    let i = log_input_size - 1;
                    // let i = log_input_size - 1 - ii;
                    for j in 0..(1usize << (log_input_size - 1)) {
                        let m = j % (1usize << i);
                        let c = (j - m) * 2;
                        let idx0 = c + m;
                        let idx1 = idx0 + (1usize << i);

                        let coeff = self.fft_coeff[i][(c >> self.log_len) ^ (k << (k_shift))];

                        res[idx0 + offset] = inputs[idx0] + coeff * inputs[idx1];
                        res[idx1 + offset] = inputs[idx1] + res[idx0 + offset];
                    }
                }

                for i in (self.log_len..(log_input_size - 1)).rev() {
                    // let i = log_input_size - 1 - ii;
                    for j in 0..(1usize << (log_input_size - 1)) {
                        let m = j % (1usize << i);
                        let c = (j - m) * 2;
                        let idx0 = c + m + offset;
                        let idx1 = idx0 + (1usize << i);

                        let coeff = self.fft_coeff[i][(c >> self.log_len) ^ (k << (k_shift))];

                        res[idx0] = res[idx0] + coeff * res[idx1];
                        res[idx1] = res[idx1] + res[idx0];
                    }
                }

                // self.log_len = self.log_l - self.log_rho
                for i in (0..self.log_len).rev() {
                    for j in 0..(1usize << (log_input_size - 1)) {
                        let m = j % (1usize << i);
                        let c = (j - m) * 2;
                        let idx0 = c + m + offset;
                        let idx1 = idx0 + (1usize << i);

                        let shift = self.log_len - (i + 1);
                        let coeff = self.fft_coeff[i][(c >> (i + 1)) ^ (k << (shift + k_shift))];

                        res[idx0] = res[idx0] + coeff * res[idx1];
                        res[idx1] = res[idx1] + res[idx0];
                    }
                }
            } else {
                // assert_eq!(self.log_len, log_input_size);
                {
                    let i = self.log_len - 1;
                    for j in 0..(1usize << (log_input_size - 1)) {
                        let m = j % (1usize << i);
                        let c = (j - m) * 2;
                        let idx0 = c + m;
                        let idx1 = idx0 + (1usize << i);

                        let shift = self.log_len - (i + 1);
                        let coeff = self.fft_coeff[i][(c >> (i + 1)) ^ (k << (shift + k_shift))];

                        res[idx0 + offset] = inputs[idx0] + coeff * inputs[idx1];
                        res[idx1 + offset] = inputs[idx1] + res[idx0 + offset];
                    }
                }

                // self.log_len = self.log_l - self.log_rho
                for i in (0..(self.log_len - 1)).rev() {
                    for j in 0..(1usize << (log_input_size - 1)) {
                        let m = j % (1usize << i);
                        let c = (j - m) * 2;
                        let idx0 = c + m + offset;
                        let idx1 = idx0 + (1usize << i);

                        let shift = self.log_len - (i + 1);
                        let coeff = self.fft_coeff[i][(c >> (i + 1)) ^ (k << (shift + k_shift))];

                        res[idx0] = res[idx0] + coeff * res[idx1];
                        res[idx1] = res[idx1] + res[idx0];
                    }
                }
            }
        }
        res
    }

    pub fn additive_interpolate(&self, inputs: &Vec<T>, log_input_size: usize, log_output_size: usize) -> Vec<T> {
        self.additive_fft(&self.additive_ifft(inputs, log_input_size), log_input_size, log_output_size)
    }

    pub fn multiplicative_fft_pack4(&self, inputs: &Vec<T>, log_input_size: usize, log_output_size: usize) -> Vec<T> {
    // pub fn multiplicative_fft(&self, inputs: &Vec<T>, log_input_size: usize, log_output_size: usize) -> Vec<T> {
        let mut res = vec![T::from(0); 1 << log_output_size];
        let mut pack_res = Vec::with_capacity(1 << (log_output_size - 2));

        {
            let s = log_output_size - log_input_size;
            let m = 2 << s;
            let lvl = self.log_l - s - 1;
            // let root = self.basis0[lvl];
            // assert_eq!(root, self.element[lvl][m - 1]);
            for k in (0..(1 << log_output_size)).step_by(m) {
                // let mut w = T::from(1);
                for j in 0..(m / 2) {
                    // assert_eq!(w, self.element[lvl][j]);
                    let w = self.element[lvl][j];
                    let idx0 = swap_bits(k >> s, log_input_size);
                    let idx1 = idx0 + (1 << (log_input_size - 1));
                    let u = inputs[idx0];
                    let v = w * inputs[idx1];
                    res[k + j] = u + v;
                    res[k + j + m / 2] = u - v;
                    // w *= root;
                }
            }
        }

        for i in 0..(1 << (log_output_size - 2)) {
            unsafe {
                pack_res.push(T::P4::pack(
                    std::mem::transmute(&res[4 * i] as *const T)
                ));
            }
        }

        for s in (log_output_size - log_input_size + 1)..log_output_size {
            let m = 2 << s;
            let lvl = self.log_l - s - 1;

            for k in (0..(1 << log_output_size)).step_by(m) {
                for j in 0..(m / 8) {
                    let pack_w = self.pack_element[lvl][j];
                    let tmp = pack_w * pack_res[(k >> 2) + j + (m >> 3)];
                    pack_res[(k >> 2) + j + (m >> 3)] = pack_res[(k >> 2) + j] - tmp;
                    pack_res[(k >> 2) + j] += tmp;
                }
            }
        }

        for i in 0..(1 << (log_output_size - 2)) {
            unsafe {
                pack_res[i].unpack(
                    std::mem::transmute(&res[4 * i] as *const T)
                    )
            }
        }
        res
    }

    pub fn multiplicative_fft(&self, inputs: &Vec<T>, log_input_size: usize, log_output_size: usize) -> Vec<T> {
        // assert_eq!(inputs.len(), 1 << log_input_size);
        let mut res = vec![T::from(0); 1 << log_output_size];

        {
            let s = log_output_size - log_input_size;
            let m = 2 << s;
            let lvl = self.log_l - s - 1;
            // let root = self.basis0[lvl];
            // assert_eq!(root, self.element[lvl][m - 1]);
            for k in (0..(1 << log_output_size)).step_by(m) {
                // let mut w = T::from(1);
                for j in 0..(m / 2) {
                    // assert_eq!(w, self.element[lvl][j]);
                    let w = self.element[lvl][j];
                    let idx0 = swap_bits(k >> s, log_input_size);
                    let idx1 = idx0 + (1 << (log_input_size - 1));
                    let u = inputs[idx0];
                    let v = w * inputs[idx1];
                    res[k + j] = u + v;
                    res[k + j + m / 2] = u - v;
                    // w *= root;
                }
            }
        }

        for s in (log_output_size - log_input_size + 1)..log_output_size {
            let m = 2 << s;
            let lvl = self.log_l - s - 1;
            // let root = self.basis0[lvl];
            // assert_eq!(root, self.element[lvl][m - 1]);
            for k in (0..(1 << log_output_size)).step_by(m) {
                // let mut w = T::from(1);
                for j in 0..(m / 2) {
                    // assert_eq!(w, self.element[lvl][j]);
                    let w = self.element[lvl][j];
                    let u = res[k + j];
                    let v = w * res[k + j + m / 2];
                    res[k + j] = u + v;
                    res[k + j + m / 2] = u - v;
                    // w *= root;
                }
            }
        }
        res
    }

    pub fn multiplicative_ifft_pack4(&self, inputs: &Vec<T>, log_size: usize) -> Vec<T> {
        // assert_eq!(inputs.len(), 1 << log_size)
        let mut res = vec![T::from(0); 1 << log_size];

        // let size_inv = T::from(1 << (61 - log_size));
        // let size_inv = T::from(1 << log_size).inv();
        let size_inv = self.pow_2_inv[log_size];

        let size = 1 << log_size;

        {
            for k in (0..size).step_by(2) {
                let u = inputs[swap_bits(k, log_size)] * size_inv;
                let v = inputs[swap_bits(k + 1, log_size)] * size_inv;
                res[k + 0] = u + v;
                res[k + 1] = u - v;
            }
        }

        {
            let lvl = self.log_l - 2;
            for k in (0..size).step_by(4) {
                let u = res[k + 0];
                let v = res[k + 2];
                res[k + 0] = u + v;
                res[k + 2] = u - v;

                let w = self.element[lvl][3];
                let u = res[k + 1];
                let v = w * res[k + 3];
                res[k + 1] = u + v;
                res[k + 3] = u - v;
            }
        }

        let mut pack_res = Vec::with_capacity(1 << (log_size - 2));
        for i in 0..(1 << (log_size - 2)) {
            unsafe {
                pack_res.push(T::P4::pack(
                    std::mem::transmute(&res[4 * i] as *const T)
                ));
            }
        }

        for s in 2..log_size {
            let m = 2 << s;
            let lvl = self.log_l - s - 1;
            for k in (0..size).step_by(m) {
                for j in 0..(m / 8) {
                    let pack_w = self.pack_element_inv[lvl][j];
                    let tmp = pack_w * pack_res[(k >> 2) + j + (m >> 3)];
                    pack_res[(k >> 2) + j + (m >> 3)] = pack_res[(k >> 2) + j] - tmp;
                    pack_res[(k >> 2) + j] += tmp;
                }
            }
        }
        for i in 0..(1 << (log_size - 2)) {
            unsafe {
                pack_res[i].unpack(
                    std::mem::transmute(&res[4 * i] as *const T)
                    )
            }
        }
        res
    }

    pub fn multiplicative_ifft(&self, inputs: &Vec<T>, log_size: usize) -> Vec<T> {
        let mut res = vec![T::from(0); 1 << log_size];

        // let size_inv = T::from(1 << (61 - log_size));
        // let size_inv = T::from(1 << log_size).inv();
        let size_inv = self.pow_2_inv[log_size];

        let size = 1 << log_size;

        for i in 0..inputs.len() {
            res[swap_bits(i, log_size)] = inputs[i] * size_inv;
        }

        for s in 0..log_size {
            let m = 2 << s;
            let lvl = self.log_l - s - 1;
            // let root = self.basis0_inv[lvl];
            // assert_eq!(root, self.element[lvl][m - 1]);
            for k in (0..size).step_by(m) {
                // let mut w = T::from(1);
                for j in 0..(m / 2) {
                    // assert_eq!(w, self.element[lvl][(m - j) % m]);
                    let w = self.element[lvl][(m - j) % m];
                    let u = res[k + j];
                    let v = w * res[k + j + m / 2];
                    res[k + j] = u + v;
                    res[k + j + m / 2] = u - v;
                    // w *= root;
                }
            }
        }
        res
    }

    pub fn multiplicative_interpolate(&self, inputs: &Vec<T>, log_input_size: usize, log_output_size: usize) -> Vec<T> {
        // let size_inv = T::from(1 << (61 - log_size));
        // let size_inv = T::from(1 << log_input_size).inv();
        let size_inv = self.pow_2_inv[log_input_size];

        let inputs =
        {
            let log_size = log_input_size;
            let mut res = inputs.to_vec();

            {
                for k in 0..(1 << (log_size - 1)) {
                    let u = size_inv * res[k];
                    let v = size_inv * res[k + (1 << (log_size - 1))];
                    res[k] = u + v;
                    res[k + (1 << (log_size - 1))] = u - v;
                }
            }
            for s in 1..log_size {
                let m = 2 << s;
                let lvl = self.log_l - s - 1;
                for k in 0..(1 << (log_size - s - 1)) {
                    for j in 0..(1 << s) {
                        let w = self.element[lvl][(m - swap_bits(j, s)) % m];
                        let u = res[k + (j << (log_size - s))];
                        let v = w * res[k + (j << (log_size - s)) + (1 << (log_size - s - 1))];
                        res[k + (j << (log_size - s))] = u + v;
                        res[k + (j << (log_size - s)) + (1 << (log_size - s - 1))] = u - v;
                    }
                }
            }
            res
        };

        let mut res = vec![T::from(0); 1 << log_output_size];

        {
            let s = log_output_size - log_input_size;
            let m = 2 << s;
            let lvl = self.log_l - s - 1;
            for k in (0..(1 << log_output_size)).step_by(m) {
                // let mut w = T::from(1);
                for j in 0..(m / 2) {
                    // assert_eq!(w, self.element[lvl][j]);
                    let w = self.element[lvl][j];
                    let u = inputs[k >> s];
                    let v = w * inputs[(k >> s) + 1];
                    res[k + j] = u + v;
                    res[k + j + m / 2] = u - v;
                    // w *= root;
                }
            }
        }

        for s in (log_output_size - log_input_size + 1)..log_output_size {
            let m = 2 << s;
            let lvl = self.log_l - s - 1;
            for k in (0..(1 << log_output_size)).step_by(m) {
                // let mut w = T::from(1);
                for j in 0..(m / 2) {
                    // assert_eq!(w, self.element[lvl][j]);
                    let w = self.element[lvl][j];
                    let u = res[k + j];
                    let v = w * res[k + j + m / 2];
                    res[k + j] = u + v;
                    res[k + j + m / 2] = u - v;
                    // w *= root;
                }
            }
        }
        res
    }

    // ret [[p^1, p^2, p^4, ... p^{2^len}] for p in pts]
    pub fn multiplicative_power_series(&self, pts: &Vec<usize>, len: usize) -> Vec<Vec<T>> {
        let mut res = vec![Vec::with_capacity(len); pts.len()];
        for k in 0..res.len() {
            res[k].push(self.element[0][pts[k]]);
        }
        for d in 1..len {
            for k in 0..res.len() {
                let tmp = res[k][d - 1] * res[k][d - 1];
                res[k].push(tmp);
            }
        }
        res
    }

    // ret [[X1(p), X2(p), X4(p), ... X{2^len}(p)] for p in pts]
    // X{2^i}(p) = p_inv[i] * sᵢ(p)
    //           = p_inv[i] * sᵢ₋₁(p)(sᵢ₋₁(p) + sᵢ₋₁(vᵢ))
    //           = p_inv[i] * sᵢ₋₁(p)(sᵢ₋₁(p) + p[i - 1])
    pub fn additive_power_series(&self, pts: &Vec<usize>, len: usize) -> Vec<Vec<T>> {
        let mut res = vec![Vec::with_capacity(len); pts.len()];
        for k in 0..res.len() {
            res[k].push(self.element[0][pts[k]]);
        }
        for d in 1..len {
            for k in 0..res.len() {
                let tmp = res[k][d - 1] * (res[k][d - 1] + self.p[d - 1]);
                res[k].push(tmp);
                res[k][d - 1] *= self.p_inv[d - 1];
            }
        }
        for k in 0..res.len() {
            res[k][len - 1] *= self.p_inv[len - 1];
        }
        res
    }
}
