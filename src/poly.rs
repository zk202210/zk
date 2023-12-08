use crate::field::Field;


pub fn first_n_field_elements<T: Field>(n: usize) -> Vec<T> {
    let mut res = vec![T::from(0); n];
    let mut gen = T::from(2);
    for i in 0..n {
        res[i] = gen;
        gen *= T::from(2);
    }
    res
}

fn build_from_roots<T: Field>(roots: &Vec<T>) -> Vec<T> {
    let len = roots.len();
    let mut poly = roots.to_vec();
    poly.push(T::from(0));

    for k in 1..len {
        let tmp = poly[k];
        poly[k] = tmp + poly[k - 1];
        for i in 1..k {
            poly[k - i] = poly[k - i] * tmp + poly[k - i - 1];
        }
        poly[0] *= tmp;
    }
    poly[len] = T::from(1);
    poly
}

pub fn precompute_lagrange_polynomials<T: Field>(x_values: &Vec<T>) -> Vec<Vec<T>> {
    let m = x_values.len();
    let mut res = vec![];

    for k in 0..m {
        let mut x_except_k: Vec<T> = Vec::with_capacity(m);
        let mut denominator = T::from(1);
        for j in 0..m {
            if k != j {
                denominator *= x_values[k] - x_values[j];
                x_except_k.push(- x_values[j]);
            }
        }
        let mut numerator = build_from_roots(&x_except_k);
        for n in numerator.iter_mut() {
            *n /= denominator;
        }
        res.push(numerator);
    }
    res
}

pub fn interpolate_with_precomputation<T: Field>(precompute: &Vec<Vec<T>>, y_values: &Vec<T>) -> Vec<T> {
    let m = y_values.len();
    let n = precompute[0].len();
    let mut res = vec![T::from(0); n];
    for k in 0..m {
        for i in 0..n {
            res[i] +=
                precompute[k][i] *
                y_values[k];
        }
    }
    res
}

pub fn poly_eval<T: Field>(poly: &Vec<T>, point: T) -> T {
    let mut res = T::from(0);
    for i in 0..poly.len() {
        res *= point;
        res += poly[poly.len() - 1 - i];
    }
    res
}


pub fn poly_eval_all<T: Field>(poly: &Vec<T>, log_out_len: usize) -> Vec<T> {
    let mut res: Vec<T> = vec![T::from(0); 1 << log_out_len];
    let poly_len = poly.len();
    for i in 0..(1 << log_out_len) {
        let mut tmp = poly[poly_len - 1];
        for j in 1..poly_len {
            tmp = tmp * T::from(i) + poly[poly_len - j - 1];
        }
        res[i] = tmp;
    }
    res
}


pub fn dot_product<T: Field>(a: &Vec<T>, b: &Vec<T>) -> T {
    let mut res = T::from(0);
    for i in 0..a.len() {
        res += a[i] * b[i];
    }
    res
}

pub fn poly_mul_sub<T: Field>(a: &Vec<T>, b: &Vec<T>, c: &Vec<T>) -> Vec<T> {
    let mut p: Vec<T> = vec![T::from(0); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            p[i + j] += a[i] * b[j];
        }
    }

    for k in 0..c.len() {
        p[k] -= c[k];
    }
    p
}

pub fn poly_add<T: Field>(a: &Vec<T>, b: &Vec<T>) -> Vec<T> {
    // assert(a.len() >= b.len())
    let mut res = a.to_vec();
    for i in 0..b.len() {
        res[i] = res[i] + b[i];
    }
    res
}

pub fn poly_scale<T: Field>(a: &Vec<T>, b: &T) -> Vec<T> {
    let mut res = a.to_vec();
    for i in 0..a.len() {
        res[i] *= *b;
    }
    res
}

pub fn poly_mul<T: Field>(a: &Vec<T>, b: &Vec<T>) -> Vec<T> {
    let mut res: Vec<T> = vec![T::from(0); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            res[i + j] += a[i] * b[j];
        }
    }
    res
}

/*
fn _poly_mul_inplace<T: Field>(res: &mut Vec<T>, r_ofs: usize, a: &Vec<T>, a_ofs: usize, a_len: usize, b: &Vec<T>, b_ofs: usize, b_len: usize, coset: &Coset<T>) {
    // reset the first b_len elements in res
    // and add res to the product of a and b
    let threshold = 1 << 7;
    if a_len <= threshold || b_len <= threshold {
        // TODO: Karatsuba for O(n^{log_2 3})
        for i in 0..a_len {
            for j in 0..b_len {
                if i == 0 {
                    res[r_ofs + i + j] = a[a_ofs + i] * b[b_ofs + j];
                } else {
                    res[r_ofs + i + j] += a[a_ofs + i] * b[b_ofs + j];
                }
            }
        }
    } else {
        // fft ifft
        let log_len = ceil_log(a_len + b_len - 1);
        let len = 1 << log_len;

        let mut a = a[a_ofs..(a_ofs + a_len)].to_vec();
        a.resize(len, T::from(0));
        let mut b = b[b_ofs..(b_ofs + b_len)].to_vec();
        b.resize(len, T::from(0));

        let mut a_pts = coset.fft_slow(&a, log_len);
        let b_pts = coset.fft_slow(&b, log_len);

        for i in 0..len {
            a_pts[i] *= b_pts[i];
        }

        let c = coset.ifft_slow(&a_pts, log_len);

        for i in 0..b_len {
            res[r_ofs + i] = c[i];
        }
        for i in b_len..(a_len + b_len - 1) {
            res[r_ofs + i] += c[i];
        }
    }
}
*/

/*
fn _poly_div_rec_inplace<T: Field>(
    q: &mut Vec<T>, q_ofs: usize,
    p: &mut Vec<T>, p_ofs: usize,
    f: &Vec<T>, f_ofs: usize,
    g: &Vec<T>, g_ofs: usize,
    n: usize,
    g_inv: &T, coset: &Coset<T>)
{

    // TODO: maybe increase the cutting point other than n <= 2
    if n == 1 {
        q[q_ofs] = *g_inv * f[f_ofs];
        p[p_ofs] = f[f_ofs];
        return;
    } else if n == 2 {
        q[q_ofs + 1] = *g_inv * f[f_ofs + 2];
        q[q_ofs + 0] = *g_inv * (f[f_ofs + 1] - q[q_ofs + 1] * g[g_ofs]);
        p[p_ofs + 0] = q[q_ofs] * g[g_ofs];
        p[p_ofs + 1] = f[f_ofs + 1];
        p[p_ofs + 2] = f[f_ofs + 2];
        return;
    }

    let n2 = n / 2;
    let n1 = n - n2;

    _poly_div_rec_inplace(q, n2 + q_ofs, p, p_ofs + 2 * n2, f, f_ofs + 2 * n2, g, g_ofs + n2, n1, g_inv, coset);

    _poly_mul_inplace(p, p_ofs + n2, q, n2 + q_ofs, n1, g, g_ofs, n2, coset);

    let mut t = vec![T::from(0); 2 * n2 - 1];
    for i in 0..n2 {
        p[p_ofs + i] = f[f_ofs + n - 1 + i] - p[p_ofs + n - 1 + i];
        t[n2 - 1 + i] = p[p_ofs + i];
    }

    let mut new_p = vec![T::from(0); 2 * n2 - 1];
    _poly_div_rec_inplace(q, q_ofs, &mut new_p, 0, &t, 0, g, g_ofs + n1, n2, g_inv, coset);

    _poly_mul_inplace(p, p_ofs, g, g_ofs, n1, q, q_ofs, n2, coset);

    for i in 0..(2 * n2 - 1) {
        p[p_ofs + n1 + i] = p[p_ofs + n1 + i] + new_p[i];
    }
}
*/

/*
fn _poly_div_short<T: Field>(
    q: &mut Vec<T>, q_ofs: usize,
    p: &mut Vec<T>, p_ofs: usize,
    f: &Vec<T>, f_ofs: usize,
    g: &Vec<T>, g_ofs: usize,
    m: usize, n: usize,
    g_inv: &T, coset: &Coset<T>)
{
    if m == 2 * n - 1 {
        _poly_div_rec_inplace(
            q, q_ofs, p, p_ofs,
            f, f_ofs, g, g_ofs, n, g_inv, coset);
    } else  {
        let n1 = m - n + 1;
        let n2 = n - n1;
        _poly_div_rec_inplace(
            q, q_ofs, p, p_ofs + n2,
            f, f_ofs + n2, g, g_ofs + n2, n1, g_inv, coset);
        _poly_mul_inplace(p, p_ofs, q, q_ofs, n1, g, g_ofs, n2, coset);
    }
}
*/

// ref: https://doi.org/10.1007/978-3-642-23568-9_16
// https://github.com/flintlib/flint2/blob/trunk/src/n_poly/n_fq_poly_divrem.c
/*
pub fn poly_div_with_leading_inv<T: Field>(f: &Vec<T>, g: &Vec<T>, g_inv: &T, coset: &Coset<T>) -> Vec<T> {
    let mut m = f.len();
    let n = g.len();

    if m < n {
        return vec![T::from(0)];
    } else if m <= 2 * n - 1 {
        let mut q = vec![T::from(0); m - n + 1];
        let mut p = vec![T::from(0); m];
        _poly_div_short(&mut q, 0, &mut p, 0, f, 0, g, 0, m, n, &g_inv, coset);

        return q;
    } else {
        let mut q = vec![T::from(0); m - n + 1];
        let mut p = f.to_vec();
        let mut r = vec![T::from(0); m];
        let k = 2 * n - 1;
        while m >= k {
            let shift = m - k;
            _poly_div_rec_inplace(
                &mut q, shift, &mut r, 0,
                &p, shift, g, 0, n, &g_inv, coset);
            for i in 0..k {
                p[shift + i] = p[shift + i] - r[i];
            }
            m -= n;
        }

        if m >= n {
            _poly_div_short(
                &mut q, 0, &mut r, 0,
                &p, 0, g, 0, m, n, &g_inv, coset);
            for i in 0..(n - 1) {
                p[i] = p[i] - r[i];
            }
        }
        return q;
    }
}

pub fn poly_div_rem_with_leading_inv<T: Field>(f: &Vec<T>, g: &Vec<T>, g_inv: &T, coset: &Coset<T>) -> (Vec<T>, Vec<T>) {
    let mut m = f.len();
    let n = g.len();

    if m < n {
        return (vec![T::from(0)], f.to_vec());
    } else if m <= 2 * n - 1 {
        let mut q = vec![T::from(0); m - n + 1];
        let mut p = vec![T::from(0); m];
        _poly_div_short(&mut q, 0, &mut p, 0, f, 0, g, 0, m, n, &g_inv, coset);

        for i in 0..(n - 1) {
            p[i] = f[i] - p[i];
        }
        p.truncate(n - 1);
        return (q, p);
    } else {
        let mut q = vec![T::from(0); m - n + 1];
        let mut p = f.to_vec();
        let mut r = vec![T::from(0); m];
        let k = 2 * n - 1;
        while m >= k {
            let shift = m - k;
            _poly_div_rec_inplace(
                &mut q, shift, &mut r, 0,
                &p, shift, g, 0, n, &g_inv, coset);
            for i in 0..k {
                p[shift + i] = p[shift + i] - r[i];
            }
            m -= n;
        }

        if m >= n {
            _poly_div_short(
                &mut q, 0, &mut r, 0,
                &p, 0, g, 0, m, n, &g_inv, coset);
            for i in 0..(n - 1) {
                p[i] = p[i] - r[i];
            }
        }
        p.truncate(n - 1);
        return (q, p);
    }
}

pub fn poly_div<T: Field>(f: &Vec<T>, g: &Vec<T>, coset: &Coset<T>) -> (Vec<T>, Vec<T>) {
    let g_inv = g[g.len() - 1].inv();
    poly_div_rem_with_leading_inv(f, g, &g_inv, coset)
}
*/

// compute f^{-1}(x) mod x^degree
// assert f[0] = 1
// O(n^2) algorithm
// pub fn poly_mod_inv<T: Field>(f: &Vec<T>, degree: usize) -> Vec<T> {
//     // assert_eq!(f[0], T::from(1));
//     let mut inv = f.to_vec();
//     inv.resize(degree, T::from(0));
// 
//     for i in 1..inv.len() {
//         inv[i] = T::from(0) - inv[i];
//         for j in (i + 1)..inv.len().min(i + f.len()) {
//             inv[j] = inv[j] + inv[i] * f[j - i];
//         }
//     }
//     inv
// }
pub fn poly_mod_inv_sparse<T: Field>(f: &Vec<T>, degree: usize) -> Vec<(usize, T)> {
    // assert_eq!(f[0], T::from(1));
    let mut non_zero: Vec<(usize, T)> = vec![];
    // let mut inv = f.to_vec();
    // inv.resize(degree, T::from(0));

    for i in 1..degree.min(f.len()) {
        if f[i] != T::from(0) {
            non_zero.push((i, f[i]));
        }
    }
    let mut inv = non_zero.clone();

    let mut i = 0;
    while i < inv.len() {
        let (idx, coeff) = inv[i];
        // for arith field
        // inv[i] = (idx, -coeff)
        i += 1;
        for j in 0..non_zero.len() {
            let (idx0, coeff0) = non_zero[j];
            let new_idx = idx + idx0;
            if new_idx >= degree {
                break;
            }
            let new_coeff = T::from(0) - coeff0 * coeff;

            let mut k = i;
            let mut inserted = false;
            while k < inv.len() {
                let (idx, coeff) = inv[k];
                if idx > new_idx {
                    break;
                } else if idx == new_idx {
                    inserted = true;
                    if coeff + new_coeff == T::from(0) {
                        inv.remove(k);
                    } else {
                        inv[k] = (idx, coeff + new_coeff);
                    }
                    break;
                }
                k += 1
            }
            if !inserted {
                inv.insert(k, (new_idx, new_coeff));
            }
        }
    }

    inv
}


// Hensel lifting
// ref: ch2.1 http://people.csail.mit.edu/madhu/ST12/scribe/lect06.pdf
/*
pub fn poly_mod_inv<T: Field>(f: &Vec<T>, degree: usize, coset: &Coset<T>) -> Vec<T> {
    let log_degree = ceil_log(degree);

    let mut f = f.to_vec();
    f.resize(1 << log_degree, T::from(0));

    let mut inv = vec![T::from(0); 1 << log_degree];
    inv[0] = T::from(1);

    // let mut p = vec![T::from(0); 1 << log_degree];

    for i in 0..log_degree {
        // let d = 1 << i;
        // a  := inv[0..d] = inv f mod x^d
        // c  := inv[d..2d] s.t. a * f = 1 + c x^d
        // h1 := f[d..2d]
        // output
        // b  := inv[d..2d]  s.t. b = -a(a * h1 + c) mod x^d
        // c  := inv[2d..4d] s.t. c = (a || b) * f[0..2d] from x^2d to x^4d
        _poly_lo_mul_add(&mut inv, &f, i, coset);
        _poly_lo_mul(&mut inv, i, coset);

        if i < log_degree - 1 {
            _poly_hi_mul(&mut inv, &f, i + 1, coset);
        }
    }
    inv.truncate(degree);
    inv
}

fn _poly_lo_mul_add<T: Field>(a: &mut Vec<T>, b: &Vec<T>, log_len: usize, coset: &Coset<T>) {
    // a[d..2d] += lo(a[0..d] * b[d..2d])
    let d = 1 << log_len;
    // for i in 0..d {
    //     a[d + i] = T::from(0);
    // }
    if d <= 512 {
        for i in 0..d {
            for j in i..d {
                a[d + j] = a[d + j] + a[i] * b[d + j - i];
            }
        }
    } else {
        let f = a[0..d].to_vec();
        let g = b[d..(2 * d)].to_vec();
        let mut f_pts = coset.fft_slow(&f, log_len + 1);
        let g_pts = coset.fft_slow(&g, log_len + 1);
        for i in 0..(2 * d) {
            f_pts[i] *= g_pts[i];
        }
        let f = coset.ifft_slow(&f_pts, log_len + 1);
        for i in 0..d {
            a[d + i] += f[i];
        }
    }
}
*/

/*
fn _poly_lo_mul<T: Field>(a: &mut Vec<T>, log_len: usize, coset: &Coset<T>) {
    // a[d..2d] *= a[0..d], take lower part
    let d = 1 << log_len;
    if d <= 512 {
        let b = a[d..(2 * d)].to_vec();
        for i in 0..d {
            for j in i..d {
                if i == 0 {
                    a[d + j] = a[i] * b[j - i];
                } else {
                    a[d + j] = a[d + j] + a[i] * b[j - i];
                }
            }
        }
    } else {
        let f = a[0..d].to_vec();
        let g = a[d..(2 * d)].to_vec();
        let mut f_pts = coset.fft_slow(&f, log_len + 1);
        let g_pts = coset.fft_slow(&g, log_len + 1);
        for i in 0..(2 * d) {
            f_pts[i] *= g_pts[i];
        }
        let f = coset.ifft_slow(&f_pts, log_len + 1);
        for i in 0..d {
            a[d + i] = f[i];
        }
    }
}
*/

/*
fn _poly_hi_mul<T: Field>(a: &mut Vec<T>, b: &Vec<T>, log_len: usize, coset: &Coset<T>) {
    // a[d..2d] *= a[0..d] * f[0..d], take high part
    let d = 1 << log_len;
    if d <= 512 {
        for i in 0..d {
            for j in i..d {
                a[2 * d - 2 - j] = a[2 * d - 2 - j] + a[d - 1 - i] * b[d - 1 - (j - i)];
            }
        }
    } else {
        let f = a[0..d].to_vec();
        let g = b[0..d].to_vec();
        let mut f_pts = coset.fft_slow(&f, log_len + 1);
        let g_pts = coset.fft_slow(&g, log_len + 1);
        for i in 0..(2 * d) {
            f_pts[i] *= g_pts[i];
        }
        let f = coset.ifft_slow(&f_pts, log_len + 1);
        for i in 0..d {
            a[d + i] = f[d + i];
        }
    }
}
*/
