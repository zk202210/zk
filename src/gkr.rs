use std::time::Instant;

use rand::Rng;

use super::field::Field;
use super::field::Coset;
use super::circuit::Circuit;
use super::circuit::Layer;
use super::circuit::GateType;
use super::hash::*;
use super::channel::ChannelType;
use super::channel::Channel;
use super::channel::Tx;
use super::channel::recv;
use super::statistic::*;
use super::tape::*;
use super::poly::*;
use super::util::*;


// Bookkeeping tables
pub struct BKT<T> {
    pub tbl_coeff: Vec<Vec<T>>,
    pub tbl_point: Vec<Vec<T>>,
    pub tbl_prev_point: Vec<Vec<T>>,
    // pub tbl_point_vec: Vec<Vec<T>>,
    // pub tbl_prev_point_vec: Vec<Vec<T>>,
    pub tbl_value: Vec<T>,
    pub final_value: Vec<T>,
}

impl<T: Field> BKT<T> {
    pub fn new() -> Self {
        Self {
            tbl_coeff: vec![],
            tbl_point: vec![],
            tbl_prev_point: vec![],
            // tbl_point_vec: vec![],
            // tbl_prev_point_vec: vec![],
            tbl_value: vec![],
            final_value: vec![],
        }
    }

    pub fn swap_tbl_point(&mut self) {
        std::mem::swap(&mut self.tbl_prev_point, &mut self.tbl_point);
        // std::mem::swap(&mut self.tbl_prev_point_vec, &mut self.tbl_point_vec);
        self.tbl_point = vec![];
        // self.tbl_point_vec = vec![];
    }

    pub fn get_final_value(&mut self) {
        self.final_value.push(self.tbl_value[0]);
    }

    pub fn generate_point_table(&mut self, rho: &Vec<T>) {
        let mut tbl: Vec<T> = Vec::with_capacity(1 << rho.len());
        tbl.push(T::from(1));
        for k in 0..rho.len() {
            let mask = 1 << k;
            for i in 0..mask {
                let tbl_i_mask = tbl[i] * rho[k];
                // tbl[i + mask] = tbl[i] * rho[k];
                tbl.push(tbl_i_mask);
                tbl[i] = tbl[i] - tbl_i_mask;
            }
        }

        self.tbl_point.push(tbl);
        // self.tbl_point_vec.push(rho.to_vec());
    }

    pub fn generate_dummy_point_table(&mut self, rho: &Vec<T>) {
        self.tbl_point.push(vec![T::from(1); 1 << rho.len()]);
    }

    pub fn merge_point_table(&mut self, scalar: T) {
        let scalar = [T::from(1) - scalar, scalar];
        for i in 0..self.tbl_point[0].len() {
            self.tbl_point[0][i] =
                scalar[0] * self.tbl_point[0][i] + scalar[1] * self.tbl_point[1][i];
        }
        // for i in 0..self.tbl_point_vec[0].len() {
        //     self.tbl_point_vec[0][i] =
        //         scalar[0] * self.tbl_point_vec[0][i] + scalar[1] * self.tbl_point_vec[1][i];
        // }
        self.tbl_point.truncate(1);
        // self.tbl_point_vec.truncate(1);
    }

    pub fn sumcheck_phase1_init(&mut self, layer: &Layer<T>, back_bit_len: usize) {
        let phase = 0;

        let back_size = 1 << back_bit_len;
        self.tbl_coeff = vec![];
        for _ in 0..layer.num_table[phase] {
            self.tbl_coeff.push(vec![T::from(0); back_size]);
        }
        let ofs = if layer.val_degree[0] == layer.num_table[0] { 0 } else { 1 };
        // ofs == 0 iff no const

        let tbl_prev_p = &self.tbl_prev_point[0];
        let size = 1 << layer.bit_len;
        for i in 0..size {
            let gate = &layer.gates[i];
            match gate.gtype {
                GateType::ADD => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    self.tbl_coeff[1][u] += tbl_prev_p[i];
                    self.tbl_coeff[0][u] += tbl_prev_p[i] * self.tbl_value[v];
                }
                GateType::SUB => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    self.tbl_coeff[1][u] += tbl_prev_p[i];
                    self.tbl_coeff[0][u] -= tbl_prev_p[i] * self.tbl_value[v];
                }
                GateType::XOR => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    self.tbl_coeff[1][u] += tbl_prev_p[i] * (T::from(1) - T::from(2) * self.tbl_value[v]);
                    self.tbl_coeff[0][u] += tbl_prev_p[i] * self.tbl_value[v];
                }
                GateType::LINEAR => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let a = gate.utils[0] as usize;
                    let b = gate.utils[1] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    if ofs == 0 {
                        self.tbl_coeff[ofs][u] += tbl_prev_p[i] * a;
                        self.tbl_coeff[ofs][v] += tbl_prev_p[i] * b;
                    } else {
                        self.tbl_coeff[1][u] += tbl_prev_p[i] * a;
                        self.tbl_coeff[0][u] += tbl_prev_p[i] * self.tbl_value[v] * b;
                    }
                }
                GateType::FLINEAR => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let a = gate.utils[0] as usize;
                    let b = gate.utils[1] as usize;
                    let c = gate.utils[2] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    let c = layer.utils[c];
                    self.tbl_coeff[1][u] += tbl_prev_p[i] * a;
                    self.tbl_coeff[0][u] += tbl_prev_p[i] * (self.tbl_value[v] * b);
                    self.tbl_coeff[0][0] += tbl_prev_p[i] * c;
                }
                GateType::SCALAR => {
                    let u = gate.inputs[0] as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    self.tbl_coeff[ofs][u] += tbl_prev_p[i] * a;
                }
                GateType::MUL => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    self.tbl_coeff[1][u] += tbl_prev_p[i] * self.tbl_value[v];
                }
                GateType::ADDMUL => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let a = gate.utils[0] as usize;
                    let b = gate.utils[1] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    let tmp = tbl_prev_p[i] * (self.tbl_value[v] + b);
                    self.tbl_coeff[1][u] += tmp;
                    // self.tbl_coeff[1][u] += tbl_prev_p[i] * b;
                    self.tbl_coeff[0][u] += tmp * a;

                    // self.tbl_coeff[1][u] += tbl_prev_p[i] * (self.tbl_value[v] + b);
                    // self.tbl_coeff[0][u] += tbl_prev_p[i] * (self.tbl_value[v] + b) * a;
                }
                GateType::SCAMUL => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    self.tbl_coeff[1][u] += tbl_prev_p[i] * self.tbl_value[v] * a;
                }
                GateType::NOT => {
                    let u = gate.inputs[0] as usize;
                    self.tbl_coeff[1][u] -= tbl_prev_p[i];
                    self.tbl_coeff[0][u] += tbl_prev_p[i];
                }
                GateType::ADDCONS => {
                    let u = gate.inputs[0] as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    self.tbl_coeff[1][u] += tbl_prev_p[i];
                    self.tbl_coeff[0][u] += tbl_prev_p[i] * a;
                }
                GateType::SUMCONS => {
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    for j in 0..gate.inputs.len() {
                        let u = gate.inputs[j] as usize;
                        self.tbl_coeff[1][u] += tbl_prev_p[i];
                    }
                    self.tbl_coeff[0][0] += tbl_prev_p[i] * a;
                }
                GateType::SUM => {
                    for j in 0..gate.inputs.len() {
                        let u = gate.inputs[j] as usize;
                        self.tbl_coeff[ofs][u] += tbl_prev_p[i];
                    }
                }
                GateType::EXPSUM => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let mut pow = T::from(1);
                    for j in u..=v {
                        self.tbl_coeff[ofs][j] += pow * tbl_prev_p[i];
                        pow = pow + pow;
                    }
                }
                GateType::SUMMUL => {
                    let v = *gate.inputs.last().unwrap() as usize;
                    for j in 0..(gate.inputs.len() - 1) {
                        let u = gate.inputs[j] as usize;
                        self.tbl_coeff[1][u] += tbl_prev_p[i] * self.tbl_value[v];
                    }
                }
                GateType::SUMLIN => {
                    let v = *gate.inputs.last().unwrap() as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = gate.utils[1] as usize;
                    let b = layer.utils[b];
                    let c = gate.utils[2] as usize;
                    let c = layer.utils[c];
                    let tmp = tbl_prev_p[i] * self.tbl_value[v] * a;
                    for j in 0..(gate.inputs.len() - 1) {
                        let u = gate.inputs[j] as usize;
                        self.tbl_coeff[1][u] += tmp;
                    }
                    self.tbl_coeff[0][0] += tmp * c;
                    self.tbl_coeff[0][0] += tbl_prev_p[i] * b;
                }
                GateType::NAAB => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    self.tbl_coeff[1][u] -= tbl_prev_p[i] * self.tbl_value[v];
                    self.tbl_coeff[0][u] += tbl_prev_p[i] * self.tbl_value[v];
                }
                GateType::RELAY => {
                    let u = gate.inputs[0] as usize;
                    self.tbl_coeff[ofs][u] += tbl_prev_p[i];
                }
                GateType::P1D2C0 => {
                    let u = gate.inputs[0] as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    self.tbl_coeff[ofs + 1][u] += tbl_prev_p[i];
                    self.tbl_coeff[ofs][u] += tbl_prev_p[i] * a;
                }
                GateType::SELECTOR => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let w = gate.inputs[2] as usize;
                    self.tbl_coeff[1][u] -= tbl_prev_p[i] * (self.tbl_value[v] - self.tbl_value[w]);
                    self.tbl_coeff[0][u] += tbl_prev_p[i] * self.tbl_value[v];
                }
                GateType::KECCAKX => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let w = gate.inputs[2] as usize;
                    self.tbl_coeff[1][u] -= tbl_prev_p[i] * self.tbl_value[v];
                    self.tbl_coeff[0][u] += tbl_prev_p[i] * self.tbl_value[v];
                    self.tbl_coeff[1][w] += tbl_prev_p[i];
                }
                GateType::NKECCAKX => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let w = gate.inputs[2] as usize;
                    self.tbl_coeff[1][u] += tbl_prev_p[i] * self.tbl_value[v];
                    self.tbl_coeff[0][u] -= tbl_prev_p[i] * self.tbl_value[v];
                    self.tbl_coeff[1][w] -= tbl_prev_p[i];
                    self.tbl_coeff[0][w] += tbl_prev_p[i];
                }
                GateType::DUMMY => {}
                GateType::IN => {}
            }
        }

    }

    pub fn sumcheck_phase1_eval(&mut self, bit_len: usize, idx: usize) -> Vec<T> {
        let mut coef = vec![T::from(0); 2];
        for u in (0..(1 << (bit_len - idx))).step_by(2) {
            let v = u + 1;
            self.tbl_value[v] = self.tbl_value[v] - self.tbl_value[u];
            self.tbl_coeff[1][v] = self.tbl_coeff[1][v] - self.tbl_coeff[1][u];
            self.tbl_coeff[0][v] = self.tbl_coeff[0][v] - self.tbl_coeff[0][u];
            coef[0] += self.tbl_value[u] * self.tbl_coeff[1][u] + self.tbl_coeff[0][u];
            // coef[1] += self.tbl_value[u] * self.tbl_coeff[1][v] + self.tbl_coeff[1][u] * self.tbl_value[v] + self.tbl_coeff[0][v];
            coef[1] += self.tbl_coeff[1][v] * self.tbl_value[v];
        }
        coef
    }

    pub fn sumcheck_phase1_eval_d2_no_const(&mut self, bit_len: usize, idx: usize) -> Vec<T> {
        let mut coef = vec![T::from(0); 3];
        for u in (0..(1 << (bit_len - idx))).step_by(2) {
            // (c1u + c1v x) (vu + vv x)^2 + (c0u + c0v x) (vu + vv x)
            // = (c1u + c1v x) (vu^2 + vv vu x + vv vu x + vv^2  x^2) + (c0u vu + c0v vu x + c0u vv x + c0v vv x^2)
            // = (c1v vv^2) x^3 + (c1u vv^2 + c1v(vvvu + vvvu)) x^2 + (c1v vu^2 + c1u(vvvu + vvvu)) x + c1u vu^2 + c0v vv x^2 + (c0v vu + c0u vv) x + c0u vu
            // => coef3 = c1v vv^2                                      =  c1v vv vv
            //    coef2 = c1u vv^2 + c1v(vvvu + vvvu) + c0v vv          = (c1u vv + c0v) vv
            //    coef1 = c1v vu^2 + c1u(vvvu + vvvu) + c0v vu + c0u vv
            //    coef0 = c1u vu^2 + c0u vu                             = (c1u vu + c0u) vu
            let v = u + 1;
            self.tbl_value[v] = self.tbl_value[v] - self.tbl_value[u];
            self.tbl_coeff[1][v] = self.tbl_coeff[1][v] - self.tbl_coeff[1][u];
            self.tbl_coeff[0][v] = self.tbl_coeff[0][v] - self.tbl_coeff[0][u];
            coef[0] += self.tbl_value[u] * (self.tbl_coeff[1][u] * self.tbl_value[u] + self.tbl_coeff[0][u]);
            coef[1] += self.tbl_value[v] * (self.tbl_coeff[1][u] * self.tbl_value[v] + self.tbl_coeff[0][v]);
            // TODO non char-2 field
            // let tmp = self.tbl_value[u] * self.tbl_value[v];
            // coef[1] += self.tbl_coeff[1][v] * (tmp + tmp);
            coef[2] += self.tbl_coeff[1][v] * self.tbl_value[v] * self.tbl_value[v];
        }
        coef
    }

    pub fn sumcheck_phase1_eval_no_const(&mut self, bit_len: usize, idx: usize) -> Vec<T> {
        let mut coef = vec![T::from(0); 2];
        for u in (0..(1 << (bit_len - idx))).step_by(2) {
            let v = u + 1;
            self.tbl_value[v] = self.tbl_value[v] - self.tbl_value[u];
            self.tbl_coeff[0][v] = self.tbl_coeff[0][v] - self.tbl_coeff[0][u];
            // self.tbl_coeff[0][v] = self.tbl_coeff[0][v] - self.tbl_coeff[0][u];
            coef[0] += self.tbl_value[u] * self.tbl_coeff[0][u]; // + self.tbl_coeff[0][u];
            // coef[2] += self.tbl_value[u] * self.tbl_coeff[0][v] + self.tbl_coeff[0][u] * self.tbl_value[v];// + self.tbl_coeff[0][v];
            coef[1] += self.tbl_coeff[0][v] * self.tbl_value[v];
        }
        coef
    }

    pub fn sumcheck_phase1_update(&mut self, bit_len: usize, idx: usize, r: T) {
        let mask = 1 << (bit_len - idx - 1);
        for u in 0..mask {
            let u0 = u << 1;
            self.tbl_value[u] = self.tbl_value[u0] + self.tbl_value[u0 + 1] * r;
            for i in 0..self.tbl_coeff.len() {
                self.tbl_coeff[i][u] = self.tbl_coeff[i][u0] + self.tbl_coeff[i][u0 + 1] * r;
            }
        }
    }

    pub fn sumcheck_phase2_init(&mut self, layer: &Layer<T>, back_bit_len: usize) {
        let phase = 1;
        let back_size = 1 << back_bit_len;

        self.tbl_coeff = vec![vec![T::from(0)]];
        for _ in 1..layer.num_table[phase] {
            self.tbl_coeff.push(vec![T::from(0); back_size]);
        }

        let tbl_prev_p = &self.tbl_prev_point[0];
        let tbl_p = &self.tbl_point[0];
        let size = 1 << layer.bit_len;
        for i in 0..size {
            let gate = &layer.gates[i];
            match gate.gtype {
                GateType::ADD => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[1][v] += tmp;
                    self.tbl_coeff[2][v] += tmp;
                }
                GateType::SUB => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[1][v] -= tmp;
                    self.tbl_coeff[2][v] += tmp;
                }
                GateType::XOR => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[1][v] += tmp;
                    self.tbl_coeff[2][v] += tmp;
                    self.tbl_coeff[3][v] -= T::from(2) * tbl_prev_p[i] * tbl_p[u];
                }
                GateType::LINEAR => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    let a = gate.utils[0] as usize;
                    let b = gate.utils[1] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    self.tbl_coeff[2][v] += tmp * a;
                    self.tbl_coeff[1][v] += tmp * b;
                }
                GateType::FLINEAR => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    let a = gate.utils[0] as usize;
                    let b = gate.utils[1] as usize;
                    let c = gate.utils[2] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    let c = layer.utils[c];
                    self.tbl_coeff[2][v] += tmp * a;
                    self.tbl_coeff[1][v] += tmp * b;
                    self.tbl_coeff[0][0] += tbl_prev_p[i] * tbl_p[0] * c;
                }
                GateType::MUL => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    self.tbl_coeff[3][v] += tbl_prev_p[i] * tbl_p[u];
                }
                GateType::ADDMUL => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let a = gate.utils[0] as usize;
                    let b = gate.utils[1] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[3][v] += tmp;
                    self.tbl_coeff[2][v] += tmp * b;
                    self.tbl_coeff[1][v] += tmp * a;
                    self.tbl_coeff[0][0] += tmp * a * b;
                    // self.tbl_coeff[3][v] += tmp;
                    // self.tbl_coeff[2][v] += tmp * b;
                    // self.tbl_coeff[1][v] += tmp * a;
                    // self.tbl_coeff[0][0] += tmp * a * b;
                }
                GateType::SCAMUL => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[3][v] += tmp * a;
                }
                GateType::NAAB => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[3][v] -= tmp;
                    self.tbl_coeff[1][v] += tmp;
                }
                GateType::NOT => {
                    let u = gate.inputs[0] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[2][0] -= tmp;
                    self.tbl_coeff[0][0] += tmp;
                }
                GateType::ADDCONS => {
                    let u = gate.inputs[0] as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[2][0] += tmp;
                    self.tbl_coeff[0][0] += tmp * a;
                }
                GateType::SUMCONS => {
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    for j in 0..gate.inputs.len() {
                        let u = gate.inputs[j] as usize;
                        let tmp = tbl_prev_p[i] * tbl_p[u];
                        self.tbl_coeff[2][0] += tmp;
                    }
                    self.tbl_coeff[0][0] += tbl_prev_p[i] * tbl_p[0] * a;
                }
                GateType::SUM => {
                    for j in 0..gate.inputs.len() {
                        let u = gate.inputs[j] as usize;
                        let tmp = tbl_prev_p[i] * tbl_p[u];
                        self.tbl_coeff[2][0] += tmp;
                    }
                }
                GateType::EXPSUM => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let mut pow = T::from(1);
                    for j in u..=v {
                        let tmp = tbl_prev_p[i] * tbl_p[j];
                        self.tbl_coeff[2][0] += pow * tmp;
                        pow = pow + pow;
                    }
                }
                GateType::SUMMUL => {
                    let v = *gate.inputs.last().unwrap() as usize;
                    for j in 0..(gate.inputs.len() - 1) {
                        let u = gate.inputs[j] as usize;
                        self.tbl_coeff[3][v] += tbl_prev_p[i] * tbl_p[u];
                    }
                }
                GateType::SUMLIN => {
                    let v = *gate.inputs.last().unwrap() as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = gate.utils[1] as usize;
                    let b = layer.utils[b];
                    let c = gate.utils[2] as usize;
                    let c = layer.utils[c];
                    for j in 0..(gate.inputs.len() - 1) {
                        let u = gate.inputs[j] as usize;
                        self.tbl_coeff[3][v] += tbl_prev_p[i] * tbl_p[u] * a;
                    }
                    self.tbl_coeff[1][v] += tbl_prev_p[i] * tbl_p[0] * c * a;
                    self.tbl_coeff[0][0] += tbl_prev_p[i] * tbl_p[0] * b;
                }
                GateType::RELAY => {
                    let u = gate.inputs[0] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    // self.tbl_coeff[2][size - 1] += tmp;
                    self.tbl_coeff[2][0] += tmp;
                }
                GateType::P1D2C0 => {
                    let u = gate.inputs[0] as usize;
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    let tmp = tbl_prev_p[i] * (tbl_prev_p[i] + a) * tbl_p[u];
                    // self.tbl_coeff[2][size - 1] += tmp;
                    self.tbl_coeff[2][0] += tmp;
                }
                GateType::SELECTOR => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let w = gate.inputs[2] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[3][v] -= tmp;
                    self.tbl_coeff[1][v] += tmp;
                    self.tbl_coeff[3][w] += tmp;
                }
                GateType::KECCAKX => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let w = gate.inputs[2] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[3][v] -= tmp;
                    self.tbl_coeff[1][v] += tmp;
                    let tmp = tbl_prev_p[i] * tbl_p[w];
                    self.tbl_coeff[2][0] += tmp;
                }
                GateType::NKECCAKX => {
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let w = gate.inputs[2] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    self.tbl_coeff[3][v] += tmp;
                    self.tbl_coeff[1][v] -= tmp;
                    let tmp = tbl_prev_p[i] * tbl_p[w];
                    self.tbl_coeff[2][0] -= tmp;
                    self.tbl_coeff[0][0] += tmp;
                }
                GateType::SCALAR => {
                    let u = gate.inputs[0] as usize;
                    let tmp = tbl_prev_p[i] * tbl_p[u];
                    let a = gate.utils[0] as usize;
                    let a = layer.utils[a];
                    self.tbl_coeff[2][0] += tmp * a;
                }
                GateType::DUMMY => {}
                GateType::IN => {}
            }
        }
    }

    pub fn sumcheck_phase2_eval(&mut self, bit_len: usize, idx: usize) -> Vec<T> {
        let mut coef = vec![T::from(0); 2];
        for u in (0..(1 << (bit_len - idx))).step_by(2) {
            let v = u + 1;
            self.tbl_value[v] = self.tbl_value[v] - self.tbl_value[u];
            self.tbl_coeff[3][v] = self.tbl_coeff[3][v] - self.tbl_coeff[3][u];
            self.tbl_coeff[1][v] = self.tbl_coeff[1][v] - self.tbl_coeff[1][u];
            self.tbl_coeff[2][v] = self.tbl_coeff[2][v] - self.tbl_coeff[2][u];
            let tmp0 = self.tbl_coeff[3][u] * self.final_value[0] + self.tbl_coeff[1][u];
            let tmp1 = self.tbl_coeff[3][v] * self.final_value[0] + self.tbl_coeff[1][v];
            coef[0] += self.tbl_value[u] * tmp0 + self.final_value[0] * self.tbl_coeff[2][u];
            // coef[1] += self.tbl_value[u] * tmp1 + tmp0 * self.tbl_value[v] + self.Vrho * self.tbl_coeff[2][v];
            coef[1] += self.tbl_value[v] * tmp1;
        }
        // AH0[1] = - AH0[0]
        // coefs + AH0[1] x + AH0[0]
        coef[0] += self.tbl_coeff[0][0];
        // coef[1] -= self.AH0;
        coef
    }

    pub fn sumcheck_phase2_eval_no_mul(&mut self, bit_len: usize, idx: usize) -> Vec<T> {
        let mut coef = vec![T::from(0); 2];
        for u in (0..(1 << (bit_len - idx))).step_by(2) {
            let v = u + 1;
            self.tbl_value[v] = self.tbl_value[v] - self.tbl_value[u];
            self.tbl_coeff[1][v] = self.tbl_coeff[1][v] - self.tbl_coeff[1][u];
            self.tbl_coeff[2][v] = self.tbl_coeff[2][v] - self.tbl_coeff[2][u];
            let tmp0 = self.tbl_coeff[1][u];
            let tmp1 = self.tbl_coeff[1][v];
            coef[0] += self.tbl_value[u] * tmp0 + self.final_value[0] * self.tbl_coeff[2][u];
            // coef[1] += self.tbl_value[u] * tmp1 + tmp0 * self.tbl_value[v] + self.Vrho * self.tbl_coeff[2][v];
            coef[1] += self.tbl_value[v] * tmp1;
        }
        // AH0[1] = - AH0[0]
        // coefs + AH0[1] x + AH0[0]
        coef[0] += self.tbl_coeff[0][0];
        // coef[1] -= self.AH0;
        coef
    }

    pub fn sumcheck_phase2_update(&mut self, bit_len: usize, idx: usize, r: T) {
        let mask = 1 << (bit_len - idx - 1);
        // for u in 0..mask {
        //     let u0 = u << 1;
        //     self.tbl_value[u] = self.tbl_value[u0] + self.tbl_value[u0 + 1] * r;
        //     self.AG[u] = self.AG[u0] + self.AG[u0 + 1] * r;
        //     self.AG0[u] = self.AG0[u0] + self.AG0[u0 + 1] * r;
        //     self.AH[u] = self.AH[u0] + self.AH[u0 + 1] * r;
        // }
        for u in 0..mask {
            let u0 = u << 1;
            self.tbl_value[u] = self.tbl_value[u0] + self.tbl_value[u0 + 1] * r;
            for i in 1..self.tbl_coeff.len() {
                self.tbl_coeff[i][u] = self.tbl_coeff[i][u0] + self.tbl_coeff[i][u0 + 1] * r;
            }
        }
        // AH0[0] = AH0[0] + AH0[1] * r
        //        = AH0[0] - AH0[0] * r
        self.tbl_coeff[0][0] *= T::from(1) - r;
    }
}

pub struct Prover<'a, 'b, 'c, T: Field, Ch: Channel<Output = T>> {
    pub tbl: BKT<T>,
    pub ch: &'a mut Ch,
    pub circ: &'c Circuit<'b, T>,
    pub q_inputs: Vec<T>,
    pub q_scalar: T,
    pub v_ab: Vec<T>,
}


impl<'a, 'b, 'c, T: Field, Ch: Channel<Output = T>> std::ops::Deref for Prover<'a, 'b, 'c, T, Ch> {
    type Target = BKT<T>;

    fn deref(&self) -> &BKT<T> {
        &self.tbl
    }
}

impl<'a, 'b, 'c, T: Field, Ch: Channel<Output = T>> std::ops::DerefMut for Prover<'a, 'b, 'c, T, Ch> {
    // type Target = BKT<T>;

    fn deref_mut(&mut self) -> &mut BKT<T> {
        &mut self.tbl
    }
}



impl<'a, 'b, 'c, T: Field, Ch: Channel<Output = T>> Prover<'a, 'b, 'c, T, Ch> {
    pub fn new(c: &'c Circuit<'b, T>, ch: &'a mut Ch) -> Self {
        Self {
            tbl: BKT::new(),
            ch: ch,
            circ: c,
            q_inputs: vec![],
            q_scalar: T::from(0),
            v_ab: vec![],
        }
    }

    fn eval_poly_output(&self) -> T {
        let bit_len = self.circ.layers[self.circ.depth - 1].bit_len;
        let mut res = T::from(0);
        for i in 0..(1 << bit_len) {
            res += self.tbl_point[0][i] * self.circ.layers[self.circ.depth - 1].values[i];
        }
        res
    }

    pub fn prove(&mut self) {
        // let layer = &self.circ.layers[self.circ.depth - 1];
        // let outputs = layer.values;
        // let layer = &self.circ.layers[self.circ.depth - 1];
        // let back_layer = &self.circ.layers[self.circ.depth - 2];

        self.ch.send(0, Tx::<T>::Input(self.circ.layers[0].values.to_vec()));

        let mut r: Vec<T> = recv!(self.ch, MultiRandom, 1);
        // let mut scalar: Vec<T> = vec![];

#[cfg(feature = "nanobench")] let mut start0 = Instant::now();
#[cfg(feature = "nanobench")] let mut wait_ms0 = channel0_get();

        self.generate_point_table(&r);

        #[cfg(feature = "memsave")]
        let mut values: Vec<T> = vec![];

        for j in 0..(self.circ.depth - 1) {

#[cfg(feature = "nanobench")] {
    let elapsed = start0.elapsed();
    let mut in_ms = elapsed.as_micros() as usize;
    let wait_ms_new = channel0_get();
    in_ms -= wait_ms_new - wait_ms0;
    println!("\t\t#sumcheck layer {}: {} ms", j, in_ms as f64 / 1000.0);
    start0 = Instant::now();
    wait_ms0 = wait_ms_new;
}
            let forword_id = self.circ.depth - 1 - j;
            let layer = &self.circ.layers[forword_id];
            let back_layer = &self.circ.layers[forword_id - 1];
            let bit_len = back_layer.bit_len;

            self.swap_tbl_point();

            #[cfg(feature = "memsave")]
            if layer.values.len() > 0 && back_layer.values.len() == 0 {
                values = layer.values.to_vec();
            }

            for p in 0..layer.num_phase {
                #[cfg(not(feature = "memsave"))] {
                    self.tbl_value = back_layer.values.to_vec();
                }
                #[cfg(feature = "memsave")]
                if back_layer.values.len() == 0 {
                    values = layer.evaluate[1](&back_layer, &values);
                    self.tbl_value = values.to_vec();
                } else {
                    self.tbl_value = back_layer.values.to_vec();
                }
                layer.sumcheck_init[p](self, &layer, bit_len);
                if j == 0 { println!("  ┃\tgkr peak mem {}", procinfo()); }
                r = Vec::with_capacity(bit_len);
                for i in 0..bit_len {
                    let coef = layer.sumcheck_eval[p](self, bit_len, i);
                    self.ch.send(0, Tx::<T>::Sumcheck(coef.to_vec()));
                    r.push(recv!(self.ch, Random, 1));
                    layer.sumcheck_update[p](self, bit_len, i, r[i]);
                }
                self.get_final_value();
                self.generate_point_table(&r);
            }
            self.ch.send(0, Tx::<T>::Layer(self.final_value.clone()));
            self.final_value = vec![];


            if layer.num_phase > 1 && j != self.circ.depth - 2 {
                // println!("{} prover pred {:?}", j,
                //          [self.tbl_coeff[3][0], self.tbl_coeff[2][0],
                //           self.tbl_coeff[1][0], self.tbl_coeff[0][0]]);
                let scalar = recv!(self.ch, Random, 1);
                self.merge_point_table(scalar);
            } else {
                // if layer.num_table[0] == 3 {
                //     println!("{} prover pred {:?}", j,
                //          [self.tbl_coeff[2][0], self.tbl_coeff[1][0],
                //           T::from(0), self.tbl_coeff[0][0]]);
                // } else if layer.num_table[0] == 1 {
                //     println!("{} prover pred {:?}", j,
                //          [T::from(0), self.tbl_coeff[0][0],
                //           T::from(0), T::from(0)]);
                // } else if layer.val_degree[0] == 2 {
                //     println!("{} prover pred {:?}", j,
                //          [self.tbl_coeff[1][0], self.tbl_coeff[0][0],
                //           T::from(0), T::from(0)]);
                // } else {
                //     println!("{} prover pred {:?}", j,
                //          [T::from(0), self.tbl_coeff[1][0],
                //           T::from(0), self.tbl_coeff[0][0]]);
                // }
                let _ = recv!(self.ch, @Syn, 1);
            }
        }
#[cfg(feature = "nanobench")] let elapsed = start0.elapsed();
#[cfg(feature = "nanobench")] let in_ms = elapsed.as_micros() as usize;
#[cfg(feature = "nanobench")] let wait_ms_new = channel0_get();
#[cfg(feature = "nanobench")] println!("\t\t#sumcheck layer last: {} ms", (in_ms + wait_ms0 -  wait_ms_new) as f64 / 1000.0);
    }


    pub fn main_prove(&mut self) {
        let mpc0 = mpc_prv_get();
        let start = Instant::now();
        // let layer = &self.circ.layers[self.circ.depth - 1];
        // let outputs = layer.values;
        // let layer = &self.circ.layers[self.circ.depth - 1];
        // let back_layer = &self.circ.layers[self.circ.depth - 2];

        // self.ch.send(0, Tx::<T>::Input(self.circ.layers[0].values.to_vec()));

        let mut r: Vec<T> = recv!(self.ch, MultiRandom, 1);
        // let mut scalar: Vec<T> = vec![];

#[cfg(feature = "nanobench")] let mut start0 = Instant::now();
#[cfg(feature = "nanobench")] let mut wait_ms0 = channel0_get();

        self.generate_point_table(&r);

        let fr = self.eval_poly_output();
        self.ch.send(0, Tx::<T>::OLayer(fr));
        let _ = recv!(self.ch, @Syn, 1);

        #[cfg(feature = "memsave")]
        let mut values: Vec<T> = vec![];

        for j in 0..(self.circ.depth - 1) {

#[cfg(feature = "nanobench")] {
    let elapsed = start0.elapsed();
    let mut in_ms = elapsed.as_micros() as usize;
    let wait_ms_new = channel0_get();
    in_ms -= wait_ms_new - wait_ms0;
    println!("\t\t#sumcheck layer {}: {} ms", j, in_ms as f64 / 1000.0);
    start0 = Instant::now();
    wait_ms0 = wait_ms_new;
}
            let layer = &self.circ.layers[self.circ.depth - 1 - j];
            let back_layer = &self.circ.layers[self.circ.depth - 2 - j];
            let bit_len = back_layer.bit_len;

            self.swap_tbl_point();

            #[cfg(feature = "memsave")]
            if layer.values.len() > 0 && back_layer.values.len() == 0 {
                values = layer.values.to_vec();
            }

            for p in 0..layer.num_phase {
                #[cfg(not(feature = "memsave"))] {
                    self.tbl_value = back_layer.values.to_vec();
                }
                #[cfg(feature = "memsave")]
                if back_layer.values.len() == 0 {
                    values = layer.evaluate[1](&back_layer, &values);
                    self.tbl_value = values.to_vec();
                } else {
                    self.tbl_value = back_layer.values.to_vec();
                }
    // let start1 = Instant::now();
                layer.sumcheck_init[p](self, &layer, bit_len);
                if j == 0 { println!("  ┃\tzkgkr peak mem {}", procinfo()); }
    // let elapsed = start1.elapsed();
    // let init_time = elapsed.as_micros() as usize;
    // println!("\t\t#sumcheck layer {} phase {} init time: {} ms", j, p, init_time as f64 / 1000.0);
                r = Vec::with_capacity(bit_len);
                for i in 0..bit_len {
                    let coef = layer.sumcheck_eval[p](self, bit_len, i);
                    self.ch.send(0, Tx::<T>::Sumcheck(coef.to_vec()));
                    r.push(recv!(self.ch, Random, 1));
                    layer.sumcheck_update[p](self, bit_len, i, r[i]);
                }
                self.get_final_value();
                self.generate_point_table(&r);

                if j == self.circ.depth - 2 && p == 0 {
                    self.q_inputs = r.clone();
                } else if j == self.circ.depth - 2 {
                    self.q_inputs.append(&mut r);
                }
            }

            // TODO XXX if no secmul, no flpcp for that layer
            if layer.num_phase > 1 && j != self.circ.depth - 2 {

                self.ch.send(0, Tx::<T>::Layer(self.final_value.clone()));
                self.final_value = vec![];

                let scalar = recv!(self.ch, Random, 1);

                self.ch.send(0, Tx::<T>::LayerPred([
                    self.tbl_coeff[3][0], self.tbl_coeff[2][0],
                    self.tbl_coeff[1][0], self.tbl_coeff[0][0]]));
                let _ = recv!(self.ch, @Syn, 1);
                self.ch.send(0, Tx::<T>::LayerScalar(scalar));
                let _ = recv!(self.ch, @Syn, 1);
                self.merge_point_table(scalar);
            } else if j != self.circ.depth - 2 {
                self.ch.send(0, Tx::<T>::Layer(self.final_value.clone()));
                self.final_value = vec![];

                let _ = recv!(self.ch, @Syn, 1);
                // TODO XXX add a function for predicate getter
                if layer.num_table[0] == 3 {
                    // d2 with const
                    self.ch.send(0, Tx::<T>::LayerPred([
                        self.tbl_coeff[2][0], self.tbl_coeff[1][0],
                        T::from(0), self.tbl_coeff[0][0]]));
                } else if layer.num_table[0] == 1 {
                    // d1 without const
                    self.ch.send(0, Tx::<T>::LayerPred([
                        T::from(0), self.tbl_coeff[0][0],
                        T::from(0), T::from(0)]));
                } else if layer.val_degree[0] == 2 {
                    // d2 no const
                    self.ch.send(0, Tx::<T>::LayerPred([
                        self.tbl_coeff[1][0], self.tbl_coeff[0][0],
                        T::from(0), T::from(0)]));
                } else {
                    // d1 with const
                    self.ch.send(0, Tx::<T>::LayerPred([
                        T::from(0), self.tbl_coeff[1][0],
                        T::from(0), self.tbl_coeff[0][0]]));
                }
                let _ = recv!(self.ch, @Syn, 1);

                self.ch.send(0, Tx::<T>::LayerScalar(T::from(0)));
                let _ = recv!(self.ch, @Syn, 1);

            // final layer
            } else if layer.num_phase > 1 {
                self.ch.send(0, Tx::<T>::FinalLayer(self.final_value.clone()));
                self.v_ab = self.final_value.clone();
                self.final_value = vec![];

                let scalar = recv!(self.ch, Random, 1);
                self.q_scalar = scalar;

                self.ch.send(0, Tx::<T>::FinalLayerPred([
                    self.tbl_coeff[3][0], self.tbl_coeff[2][0],
                    self.tbl_coeff[1][0], self.tbl_coeff[0][0]]));
                let _ = recv!(self.ch, @Syn, 1);
                // self.merge_point_table(&scalar);
            } else {
                self.ch.send(0, Tx::<T>::FinalLayer(self.final_value.clone()));
                self.v_ab = self.final_value.clone();
                self.final_value = vec![];
                let _ = recv!(self.ch, @Syn, 1);
                if layer.num_table[0] == 3 {
                    // d2 with const
                    self.ch.send(0, Tx::<T>::FinalLayerPred([
                        self.tbl_coeff[2][0], self.tbl_coeff[1][0],
                        T::from(0), self.tbl_coeff[0][0]]));
                } else if layer.num_table[0] == 1 {
                    // d1 without const
                    self.ch.send(0, Tx::<T>::FinalLayerPred([
                        T::from(0), self.tbl_coeff[0][0],
                        T::from(0), T::from(0)]));
                } else if layer.val_degree[0] == 2 {
                    // d2 no const
                    self.ch.send(0, Tx::<T>::FinalLayerPred([
                        self.tbl_coeff[1][0], self.tbl_coeff[0][0],
                        T::from(0), T::from(0)]));
                } else {
                    // d1 with const
                    self.ch.send(0, Tx::<T>::FinalLayerPred([
                        T::from(0), self.tbl_coeff[1][0],
                        T::from(0), self.tbl_coeff[0][0]]));
                }
                let _ = recv!(self.ch, @Syn, 1);
            }

        }
#[cfg(feature = "nanobench")] let elapsed = start0.elapsed();
#[cfg(feature = "nanobench")] let in_ms = elapsed.as_micros() as usize;
#[cfg(feature = "nanobench")] let wait_ms_new = channel0_get();
#[cfg(feature = "nanobench")] println!("\t\t#sumcheck layer last: {} ms", (in_ms + wait_ms0 -  wait_ms_new) as f64 / 1000.0);


        self.ch.send(0, Tx::<T>::MPCitH0);
        let _ = recv!(self.ch, @Syn, 1);

        if self.circ.depth > 2 {
            // sync
            self.ch.send(0, Tx::<T>::FLPCP);
            let _ = recv!(self.ch, @Syn, 1);
        }

        self.ch.send(0, Tx::<T>::MPCitH);
        let _ = recv!(self.ch, @Syn, 1);

        let elapsed = start.elapsed();
        let in_ms = elapsed.as_micros() as usize;
        let mpc1 = mpc_prv_get();
        let flp_ms = flp_prv_get();
        gkr_prv_inc(in_ms - (mpc1 - mpc0) - flp_ms);
    }
}

pub struct Verifier<'a, 'b, 'c, T: Field, Ch: Channel<Output = T>> {
    pub tbl: BKT<T>,
    ch: &'a mut Ch,
    circ: &'c Circuit<'b, T>,
    pub q_scalar: T,
    pub v_ab: Vec<T>,
}

impl<'a, 'b, 'c, T: Field, Ch: Channel<Output = T>> std::ops::Deref for Verifier<'a, 'b, 'c, T, Ch> {
    type Target = BKT<T>;

    fn deref(&self) -> &BKT<T> {
        &self.tbl
    }
}

impl<'a, 'b, 'c, T: Field, Ch: Channel<Output = T>> std::ops::DerefMut for Verifier<'a, 'b, 'c, T, Ch> {
    // type Target = BKT<T>;

    fn deref_mut(&mut self) -> &mut BKT<T> {
        &mut self.tbl
    }
}

impl<'a, 'b, 'c, T: Field, Ch: Channel<Output = T>> Verifier<'a, 'b, 'c, T, Ch> {
    pub fn new(circ: &'c Circuit<'b, T>, ch: &'a mut Ch) -> Self {
        Self {
            tbl: BKT::new(),
            ch: ch,
            circ: circ,
            q_scalar: T::from(0),
            v_ab: vec![],
        }
    }

    // fn eval_poly_output_old(&mut self, outputs: &Vec<T>, tau: &Vec<T>, bit_len: usize, log_grp_len: usize) -> T {
    //     self.generate_point_table(&tau);

    //     let mut res = T::from(0);
    //     for i in 0..(1 << bit_len) {
    //         res += self.tbl_point[0][i] * outputs[i];
    //     }
    //     res
    // }

    fn eval_poly_output(&mut self, outputs: &Vec<T>, tau: &Vec<T>, bit_len: usize, log_grp_len: usize) -> T {
        self.generate_point_table(&tau[0..(bit_len - log_grp_len)].to_vec());
        self.generate_point_table(&tau[(bit_len - log_grp_len)..].to_vec());

        let mut res = T::from(0);
        for g in 0..(1 << log_grp_len) {
            for i in 0..(1 << (bit_len - log_grp_len)) {
                res += self.tbl_point[1][g] * self.tbl_point[0][i] * outputs[i + (g << (bit_len - log_grp_len))];
            }
        }

        // self.swap_tbl_point();
        // self.generate_point_table(&tau);

        res
    }

    fn predicate_out(&mut self, layer_id: usize, r0: &Vec<Vec<T>>, r: &Vec<Vec<T>>) -> [T; 4] {
        let r0 = &r0[0];

        if self.circ.layers[layer_id].predicate.len() > 0 {

            if layer_id >= 2 && self.circ.layers[layer_id - 1].predicate.len() == 0 {
                self.tbl_point = vec![];
                self.generate_point_table(&r[0]);
                if r.len() > 1 {
                    self.generate_point_table(&r[1]);
                } else {
                    self.generate_dummy_point_table(&r[0]);
                }
            }

            let pred = &self.circ.layers[layer_id].predicate[0];
            if r.len() > 1 {
                return pred(&r0, &r[0], &r[1]);
            } else {
                return pred(&r0, &r[0], &vec![]);
            }
        }

        // default O(C) predicate
        let mut res = [T::from(0); 4];
        let len0 = self.circ.layers[layer_id].bit_len;
        let len1 = self.circ.layers[layer_id - 1].bit_len;
        // assert some total grp num
        let grp_len = self.circ.layers[layer_id - 1].log_grp_len;
        let mut pred = [T::from(0); GateType::LEN];

        self.swap_tbl_point();
        self.generate_point_table(&r[0][0..(len1 - grp_len)].to_vec());
        if r.len() > 1 {
            self.generate_point_table(&r[1][0..(len1 - grp_len)].to_vec());
        } else {
            self.generate_dummy_point_table(&r[0][0..(len1 - grp_len)].to_vec());
        }
        // self.generate_point_table(&r[0][(len1 - grp_len)..].to_vec());
        // self.generate_point_table(&r[1][(len1 - grp_len)..].to_vec());

        let mut grp_pred3 = T::from(1);
        for i in 0..grp_len {
            grp_pred3 *= eq3(r0[len0 - grp_len + i], r[0][len1 - grp_len + i], r[1][len1 - grp_len + i]);
        }

        let mut grp_pred2 = T::from(1);
        for i in 0..grp_len {
            grp_pred2 *= eq2(r0[len0 - grp_len + i], r[0][len1 - grp_len + i]) * (T::from(1) - r[1][len1 - grp_len + i]);
        }

        let layer = &self.circ.layers[layer_id];
        for i in 0..(1 << (len0 - grp_len)) {
            let g = &layer.gates[i];
            let gt = g.gtype;
            match gt {
                GateType::DUMMY => {},
                GateType::IN => {},
                GateType::SCALAR => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                },
                GateType::SCAMUL => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                },
                GateType::ADDMUL => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0] * a * b;
                },
                GateType::LINEAR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    if layer.val_degree[0] == layer.num_table[0] {
                        res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                        res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][v] * self.tbl_point[1][0] * b;
                    } else {
                        res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                        res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    }
                },
                GateType::FLINEAR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let b = g.utils[1] as usize;
                    let c = g.utils[2] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    let c = layer.utils[c];
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][0] * self.tbl_point[1][0] * c;
                },
                GateType::RELAY | GateType::NOT => {
                    let u = g.inputs[0] as usize;
                    pred[gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                },
                GateType::P1D2C0 => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let tmp = self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    res[0] += tmp;
                    res[1] += tmp * a;
                },
                GateType::ADDCONS => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                },
                GateType::SELECTOR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    let nab = self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    let mul = self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][w];
                    res[0] += mul - nab;
                    res[2] += nab;
                },
                GateType::KECCAKX => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][w] * self.tbl_point[1][0];
                    res[0] -= self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                },
                GateType::NKECCAKX => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    res[1] -= self.tbl_prev_point[0][i] * self.tbl_point[0][w] * self.tbl_point[1][0];
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][w] * self.tbl_point[1][0];
                    res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[2] -= self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                },
                GateType::SUMCONS => {
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    for k in 0..g.inputs.len() {
                        let u = g.inputs[k] as usize;
                        pred[gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                        res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    }
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][0] * self.tbl_point[1][0] * a;
                },
                GateType::SUM => {
                    for k in 0..g.inputs.len() {
                        let u = g.inputs[k] as usize;
                        pred[gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    }
                },
                GateType::EXPSUM => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let mut pow = T::from(1);
                    for k in u..=v {
                        pred[gt as usize] += pow * self.tbl_prev_point[0][i] * self.tbl_point[0][k] * self.tbl_point[1][0];
                        pow = pow + pow;
                    }
                },
                GateType::SUMMUL => {
                    let v = *g.inputs.last().unwrap() as usize;
                    for k in 0..(g.inputs.len() - 1) {
                        let u = g.inputs[k] as usize;
                        res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    }
                },
                GateType::SUMLIN => {
                    let v = *g.inputs.last().unwrap() as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    let c = g.utils[2] as usize;
                    let c = layer.utils[c];
                    for k in 0..(g.inputs.len() - 1) {
                        let u = g.inputs[k] as usize;
                        res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    }
                    res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][0] * self.tbl_point[1][v] * c * a;
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][0] * self.tbl_point[1][0] * b;
                },
                _ => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    pred[gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                },
            }
        }
        res[0] += pred[GateType::MUL as usize]
                - T::from(2) * pred[GateType::XOR as usize]
                - pred[GateType::NAAB as usize];
        res[1] += pred[GateType::ADD as usize]
                + pred[GateType::SUB as usize]
                + pred[GateType::XOR as usize];
        res[2] += pred[GateType::ADD as usize]
                - pred[GateType::SUB as usize]
                + pred[GateType::XOR as usize]
                + pred[GateType::NAAB as usize];
        res[3] += pred[GateType::NOT as usize];

        res[0] *= grp_pred3;
        res[1] = res[1] * grp_pred3 + (
            pred[GateType::RELAY as usize] +
            pred[GateType::SUM as usize] +
            pred[GateType::EXPSUM as usize] -
            pred[GateType::NOT as usize]) * grp_pred2;
        res[2] *= grp_pred3;
        res[3] *= grp_pred2;

        res
    }

    // default O(C) predicate
    fn predicate(&mut self, layer_id: usize, scalar: T, r0: &Vec<Vec<T>>, r: &Vec<Vec<T>>) -> [T; 4] {
        let scalar = [T::from(1) - scalar, scalar];

        if self.circ.layers[layer_id].predicate.len() > 0 {


            if layer_id >= 2 && self.circ.layers[layer_id - 1].predicate.len() == 0 {
                self.tbl_point = vec![];
                self.generate_point_table(&r[0]);
                if r.len() > 1 {
                    self.generate_point_table(&r[1]);
                } else {
                    self.generate_dummy_point_table(&r[0]);
                }
            }

            let pred = &self.circ.layers[layer_id].predicate[0];

            let mut res = if r.len() > 1 {
                pred(&r0[0], &r[0], &r[1])
            } else {
                pred(&r0[0], &r[0], &vec![])
            };

            let res1 = if r.len() > 1 {
                pred(&r0[1], &r[0], &r[1])
            } else {
                pred(&r0[1], &r[0], &vec![])
            };

            for i in 0..4 {
                res[i] = scalar[0] * res[i] + scalar[1] * res1[i];
            }
            return res;
        }

        let mut res = [T::from(0); 4];
        let len0 = self.circ.layers[layer_id].bit_len;
        let len1 = self.circ.layers[layer_id - 1].bit_len;
        let grp_len0 = self.circ.layers[layer_id].log_grp_len;
        let grp_len = self.circ.layers[layer_id - 1].log_grp_len;
        let mut pred = [[T::from(0); GateType::LEN]; 2];

        if grp_len0 != grp_len {
            self.swap_tbl_point();
            self.generate_point_table(&r0[0]);
            self.generate_point_table(&r0[1]);
        }

        // self.merge_point_table(scalar);

        self.swap_tbl_point();
        self.generate_point_table(&r[0][0..(len1 - grp_len)].to_vec());
        // self.generate_point_table(&r[1][0..(len1 - grp_len)].to_vec());
        if r.len() > 1 {
            self.generate_point_table(&r[1][0..(len1 - grp_len)].to_vec());
        } else {
            self.generate_dummy_point_table(&r[0][0..(len1 - grp_len)].to_vec());
        }

        let mut grp_pred3 = [T::from(1), T::from(1)];
        for j in 0..2 {
            for i in 0..grp_len {
                grp_pred3[j] *= eq3(r0[j][len0 - grp_len + i], r[0][len1 - grp_len + i], r[1][len1 - grp_len + i]);
            }
        }

        let mut grp_pred2 = [T::from(1), T::from(1)];
        for j in 0..2 {
            for i in 0..grp_len {
                grp_pred2[j] *= eq2(r0[j][len0 - grp_len + i], r[0][len1 - grp_len + i]) * (T::from(1) - r[1][len1 - grp_len + i]);
            }
        }

        let layer = &self.circ.layers[layer_id];
        for i in 0..(1 << (len0 - grp_len)) {
            let g = &layer.gates[i];
            let gt = g.gtype;

            match gt {
                GateType::DUMMY => {},
                GateType::IN => {},
                GateType::SCALAR => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                },
                GateType::SCAMUL => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[0] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                },
                GateType::ADDMUL => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    res[0] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    res[2] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0] * a * b;
                },
                GateType::LINEAR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    if layer.val_degree[0] == layer.num_table[0] {
                        res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                        res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][v] * self.tbl_point[1][0] * b;
                    } else {
                        res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                        res[2] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    }
                },
                GateType::FLINEAR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let b = g.utils[1] as usize;
                    let c = g.utils[2] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    let c = layer.utils[c];
                    res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    res[2] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][0] * self.tbl_point[1][0] * c;
                },
                GateType::RELAY | GateType::NOT => {
                    let u = g.inputs[0] as usize;
                    pred[0][gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    pred[1][gt as usize] += self.tbl_prev_point[1][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                },
                GateType::P1D2C0 => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let tmp = (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0];
                    res[0] += tmp;
                    res[1] += tmp * a;
                },
                GateType::SELECTOR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    let nab = (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v];
                    let mul = (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][w];
                    res[0] += mul - nab;
                    res[2] += nab;
                },
                GateType::KECCAKX => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    let tmp = scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i];
                    res[1] += tmp * self.tbl_point[0][w] * self.tbl_point[1][0];
                    let tmp = tmp * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[0] -= tmp;
                    res[2] += tmp;
                },
                GateType::NKECCAKX => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    let tmp0 = scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i];
                    let tmp = tmp0 * self.tbl_point[0][w] * self.tbl_point[1][0];
                    res[1] -= tmp;
                    res[3] += tmp;
                    let tmp = tmp0 * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[0] += tmp;
                    res[2] -= tmp;
                },
                GateType::ADDCONS => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0];
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                },
                GateType::SUMCONS => {
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    for k in 0..g.inputs.len() {
                        let u = g.inputs[k] as usize;
                        res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0];
                    }
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][0] * self.tbl_point[1][0] * a;
                },
                GateType::SUM => {
                    for k in 0..g.inputs.len() {
                        let u = g.inputs[k] as usize;
                        pred[0][gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                        pred[1][gt as usize] += self.tbl_prev_point[1][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    }
                },
                GateType::EXPSUM => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let mut pow = T::from(1);
                    for k in u..=v {
                        pred[0][gt as usize] += pow * self.tbl_prev_point[0][i] * self.tbl_point[0][k] * self.tbl_point[1][0];
                        pred[1][gt as usize] += pow * self.tbl_prev_point[1][i] * self.tbl_point[0][k] * self.tbl_point[1][0];
                        pow = pow + pow;
                    }
                },
                GateType::SUMMUL => {
                    let v = *g.inputs.last().unwrap() as usize;
                    for k in 0..(g.inputs.len() - 1) {
                        let u = g.inputs[k] as usize;
                        res[0] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v];
                    }
                },
                GateType::SUMLIN => {
                    let v = *g.inputs.last().unwrap() as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    let c = g.utils[2] as usize;
                    let c = layer.utils[c];
                    for k in 0..(g.inputs.len() - 1) {
                        let u = g.inputs[k] as usize;
                        res[0] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    }
                    res[2] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][0] * self.tbl_point[1][v] * c * a;
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][0] * self.tbl_point[1][0] * b;
                },
                _ => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    pred[0][gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    pred[1][gt as usize] += self.tbl_prev_point[1][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                },
            }
        }

        // assert_eq!(res[0], T::from(0));
        // assert_eq!(res[1], T::from(0));
        // assert_eq!(res[2], T::from(0));
        // assert_eq!(res[3], T::from(0));
        for i in 0..2 {
            res[0] += scalar[i] * (
                pred[i][GateType::MUL as usize] -
                T::from(2) * pred[i][GateType::XOR as usize] -
                pred[i][GateType::NAAB as usize]
                ) * grp_pred3[i];
            res[1] += scalar[i] * (
                (
                    pred[i][GateType::ADD as usize] +
                    pred[i][GateType::SUB as usize] +
                    pred[i][GateType::XOR as usize]
                ) * grp_pred3[i] + (
                    pred[i][GateType::SUM as usize] +
                    pred[i][GateType::EXPSUM as usize] +
                    pred[i][GateType::RELAY as usize] -
                    pred[i][GateType::NOT as usize]
                ) * grp_pred2[i]);
            res[2] += scalar[i] * (
                pred[i][GateType::ADD as usize] -
                pred[i][GateType::SUB as usize] +
                pred[i][GateType::XOR as usize] +
                pred[i][GateType::NAAB as usize]
                ) * grp_pred3[i];
            res[3] += scalar[i] * (
                pred[i][GateType::NOT as usize]
                ) * grp_pred2[i];
        }
        res
    }


    // TODO reorganize these two default predicate function
    /*
    fn predicate_out_old(&mut self, layer_id: usize, r0: &Vec<Vec<T>>, r: &Vec<Vec<T>>) -> [T; 4] {
        let r0 = &r0[0];

        if self.circ.layers[layer_id].predicate.len() > 0 {

            if layer_id >= 2 && self.circ.layers[layer_id - 1].predicate.len() == 0 {
                self.tbl_point = vec![];
                self.generate_point_table(&r[0]);
                self.generate_point_table(&r[1]);
            }

            let pred = &self.circ.layers[layer_id].predicate[0];
            if r.len() > 1 {
                return pred(&r0, &r[0], &r[1]);
            } else {
                return pred(&r0, &r[0], &vec![]);
            }
        }

        // default O(C) predicate
        let mut res = [T::from(0); 4];
        let len0 = self.circ.layers[layer_id].bit_len;
        let mut pred = [T::from(0); GateType::LEN];

        self.swap_tbl_point();
        self.generate_point_table(&r[0]);
        self.generate_point_table(&r[1]);

        let layer = &self.circ.layers[layer_id];
        for i in 0..(1 << len0) {
            let g = &layer.gates[i];
            let gt = g.gtype;
            match gt {
                GateType::DUMMY => {},
                GateType::IN => {},
                GateType::SCALAR => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                },
                GateType::SCAMUL => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                },
                GateType::ADDMUL => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0] * a * b;
                },
                GateType::LINEAR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    if layer.val_degree[0] == layer.num_table[0] {
                        res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                        res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][v] * self.tbl_point[1][0] * b;
                    } else {
                        res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                        res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    }
                },
                GateType::FLINEAR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let b = g.utils[1] as usize;
                    let c = g.utils[2] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    let c = layer.utils[c];
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][0] * self.tbl_point[1][0] * c;
                },
                GateType::RELAY | GateType::NOT => {
                    let u = g.inputs[0] as usize;
                    pred[gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                },
                GateType::P1D2C0 => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let tmp = self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    res[0] += tmp;
                    res[1] += tmp * a;
                },
                GateType::ADDCONS => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                },
                GateType::KECCAKX => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][w] * self.tbl_point[1][0];
                    res[0] -= self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                },
                GateType::NKECCAKX => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    res[1] -= self.tbl_prev_point[0][i] * self.tbl_point[0][w] * self.tbl_point[1][0];
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][w] * self.tbl_point[1][0];
                    res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[2] -= self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                },
                GateType::SUMCONS => {
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    for k in 0..g.inputs.len() {
                        let u = g.inputs[k] as usize;
                        pred[gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                        res[1] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    }
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][0] * self.tbl_point[1][0] * a;
                },
                GateType::SUM => {
                    for k in 0..g.inputs.len() {
                        let u = g.inputs[k] as usize;
                        pred[gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    }
                },
                GateType::EXPSUM => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let mut pow = T::from(1);
                    for k in u..=v {
                        pred[gt as usize] += pow * self.tbl_prev_point[0][i] * self.tbl_point[0][k] * self.tbl_point[1][0];
                        pow = pow + pow;
                    }
                },
                GateType::SUMMUL => {
                    let v = *g.inputs.last().unwrap() as usize;
                    for k in 0..(g.inputs.len() - 1) {
                        let u = g.inputs[k] as usize;
                        res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    }
                },
                GateType::SUMLIN => {
                    let v = *g.inputs.last().unwrap() as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    let c = g.utils[2] as usize;
                    let c = layer.utils[c];
                    for k in 0..(g.inputs.len() - 1) {
                        let u = g.inputs[k] as usize;
                        res[0] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    }
                    res[2] += self.tbl_prev_point[0][i] * self.tbl_point[0][0] * self.tbl_point[1][v] * c * a;
                    res[3] += self.tbl_prev_point[0][i] * self.tbl_point[0][0] * self.tbl_point[1][0] * b;
                },
                _ => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    pred[gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                },
            }
        }
        res[0] += pred[GateType::MUL as usize]
                - T::from(2) * pred[GateType::XOR as usize]
                - pred[GateType::NAAB as usize];
        res[1] += pred[GateType::ADD as usize]
                + pred[GateType::SUB as usize]
                + pred[GateType::XOR as usize]
                + pred[GateType::RELAY as usize]
                + pred[GateType::SUM as usize]
                + pred[GateType::EXPSUM as usize]
                - pred[GateType::NOT as usize];
        res[2] += pred[GateType::ADD as usize]
                - pred[GateType::SUB as usize]
                + pred[GateType::XOR as usize]
                + pred[GateType::NAAB as usize];
        res[3] += pred[GateType::NOT as usize];
        res
    }

    // default O(C) predicate
    fn predicate_old(&mut self, layer_id: usize, scalar: T, r0: &Vec<Vec<T>>, r: &Vec<Vec<T>>) -> [T; 4] {
        let scalar = [T::from(1) - scalar, scalar];

        if self.circ.layers[layer_id].predicate.len() > 0 {


            if layer_id >= 2 && self.circ.layers[layer_id - 1].predicate.len() == 0 {
                self.tbl_point = vec![];
                self.generate_point_table(&r[0]);
                self.generate_point_table(&r[1]);
            }

            let pred = &self.circ.layers[layer_id].predicate[0];

            let mut res = if r.len() > 1 {
                pred(&r0[0], &r[0], &r[1])
            } else {
                pred(&r0[0], &r[0], &vec![])
            };

            let res1 = if r.len() > 1 {
                pred(&r0[1], &r[0], &r[1])
            } else {
                pred(&r0[1], &r[0], &vec![])
            };

            for i in 0..4 {
                res[i] = scalar[0] * res[i] + scalar[1] * res1[i];
            }
            return res;
        }

        let mut res = [T::from(0); 4];
        let len0 = self.circ.layers[layer_id].bit_len;
        let mut pred = [[T::from(0); GateType::LEN]; 2];

        // self.merge_point_table(scalar);
        self.swap_tbl_point();
        self.generate_point_table(&r[0]);
        self.generate_point_table(&r[1]);

        let layer = &self.circ.layers[layer_id];
        for i in 0..(1 << len0) {
            let g = &layer.gates[i];
            let gt = g.gtype;

            match gt {
                GateType::DUMMY => {},
                GateType::IN => {},
                GateType::SCALAR => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                },
                GateType::SCAMUL => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[0] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                },
                GateType::ADDMUL => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    res[0] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    res[2] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0] * a * b;
                },
                GateType::LINEAR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    if layer.val_degree[0] == layer.num_table[0] {
                        res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                        res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][v] * self.tbl_point[1][0] * b;
                    } else {
                        res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                        res[2] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    }
                },
                GateType::FLINEAR => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let a = g.utils[0] as usize;
                    let b = g.utils[1] as usize;
                    let c = g.utils[2] as usize;
                    let a = layer.utils[a];
                    let b = layer.utils[b];
                    let c = layer.utils[c];
                    res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    res[2] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * b;
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][0] * self.tbl_point[1][0] * c;
                },
                GateType::RELAY | GateType::NOT => {
                    let u = g.inputs[0] as usize;
                    pred[0][gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    pred[1][gt as usize] += self.tbl_prev_point[1][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                },
                GateType::P1D2C0 => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let tmp = (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0];
                    res[0] += tmp;
                    res[1] += tmp * a;
                },
                GateType::KECCAKX => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    let tmp = scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i];
                    res[1] += tmp * self.tbl_point[0][w] * self.tbl_point[1][0];
                    let tmp = tmp * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[0] -= tmp;
                    res[2] += tmp;
                },
                GateType::NKECCAKX => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let w = g.inputs[2] as usize;
                    let tmp0 = scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i];
                    let tmp = tmp0 * self.tbl_point[0][w] * self.tbl_point[1][0];
                    res[1] -= tmp;
                    res[3] += tmp;
                    let tmp = tmp0 * self.tbl_point[0][u] * self.tbl_point[1][v];
                    res[0] += tmp;
                    res[2] -= tmp;
                },
                GateType::ADDCONS => {
                    let u = g.inputs[0] as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0];
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0] * a;
                },
                GateType::SUMCONS => {
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    for k in 0..g.inputs.len() {
                        let u = g.inputs[k] as usize;
                        res[1] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][0];
                    }
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][0] * self.tbl_point[1][0] * a;
                },
                GateType::SUM => {
                    for k in 0..g.inputs.len() {
                        let u = g.inputs[k] as usize;
                        pred[0][gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                        pred[1][gt as usize] += self.tbl_prev_point[1][i] * self.tbl_point[0][u] * self.tbl_point[1][0];
                    }
                },
                GateType::EXPSUM => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    let mut pow = T::from(1);
                    for k in u..=v {
                        pred[0][gt as usize] += pow * self.tbl_prev_point[0][i] * self.tbl_point[0][k] * self.tbl_point[1][0];
                        pred[1][gt as usize] += pow * self.tbl_prev_point[1][i] * self.tbl_point[0][k] * self.tbl_point[1][0];
                        pow = pow + pow;
                    }
                },
                GateType::SUMMUL => {
                    let v = *g.inputs.last().unwrap() as usize;
                    for k in 0..(g.inputs.len() - 1) {
                        let u = g.inputs[k] as usize;
                        res[0] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v];
                    }
                },
                GateType::SUMLIN => {
                    let v = *g.inputs.last().unwrap() as usize;
                    let a = g.utils[0] as usize;
                    let a = layer.utils[a];
                    let b = g.utils[1] as usize;
                    let b = layer.utils[b];
                    let c = g.utils[2] as usize;
                    let c = layer.utils[c];
                    for k in 0..(g.inputs.len() - 1) {
                        let u = g.inputs[k] as usize;
                        res[0] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][u] * self.tbl_point[1][v] * a;
                    }
                    res[2] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][0] * self.tbl_point[1][v] * c * a;
                    res[3] += (scalar[0] * self.tbl_prev_point[0][i] + scalar[1] * self.tbl_prev_point[1][i]) * self.tbl_point[0][0] * self.tbl_point[1][0] * b;
                },
                _ => {
                    let u = g.inputs[0] as usize;
                    let v = g.inputs[1] as usize;
                    pred[0][gt as usize] += self.tbl_prev_point[0][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                    pred[1][gt as usize] += self.tbl_prev_point[1][i] * self.tbl_point[0][u] * self.tbl_point[1][v];
                },
            }
        }

        for i in 0..2 {
            res[0] += scalar[i] * (
                pred[i][GateType::MUL as usize] -
                T::from(2) * pred[i][GateType::XOR as usize] -
                pred[i][GateType::NAAB as usize]
                );
            res[1] += scalar[i] * (
                pred[i][GateType::ADD as usize] +
                pred[i][GateType::SUB as usize] +
                pred[i][GateType::XOR as usize] +
                pred[i][GateType::SUM as usize] +
                pred[i][GateType::EXPSUM as usize] +
                pred[i][GateType::RELAY as usize] -
                pred[i][GateType::NOT as usize]
                );
            res[2] += scalar[i] * (
                pred[i][GateType::ADD as usize] -
                pred[i][GateType::SUB as usize] +
                pred[i][GateType::XOR as usize] +
                pred[i][GateType::NAAB as usize]
                );
            res[3] += scalar[i] * (
                pred[i][GateType::NOT as usize]
                );
        }
        res
    }
    */

    pub fn verify(&mut self, outputs: &Vec<T>) -> bool {

        let mut res = true;

        let layer = &self.circ.layers[self.circ.depth - 1];
        let mut bit_len = layer.bit_len;

        // TODO: do not send inputs again in VPD setting
        //  and take input from argument
        let inputs: Vec<T> = recv!(self.ch, Input, 0);

        let mut tau = vec![T::from(0); bit_len];
        // TODO: assert no more than two phases
        let mut rs: Vec<Vec<Vec<T>>> = vec![];
        let mut scalars: Vec<T> = vec![];
        let mut final_vals: Vec<Vec<T>> = vec![];
        let mut coefs: Vec<Vec<T>> = vec![];

        for i in 0..bit_len {
            tau[i] = T::random();
        }
        self.ch.send(1, Tx::<T>::MultiRandom(tau.clone()));

        for j in 0..(self.circ.depth - 1) {
            let layer = &self.circ.layers[self.circ.depth - 1 - j];
            let back_layer = &self.circ.layers[self.circ.depth - 2 - j];
            bit_len = back_layer.bit_len;

            let mut r_p: Vec<Vec<T>> = Vec::with_capacity(layer.num_phase);
            for _ in 0..layer.num_phase {
                let mut r = vec![T::from(0); bit_len];
                for i in 0..bit_len {
                    let coef: Vec<T> = recv!(self.ch, Sumcheck, 0);
                    gkr_bw_inc(T::SIZE * coef.len());
                    coefs.push(coef);
                    r[i] = T::random();
                    self.ch.send(1, Tx::<T>::Random(r[i]));
                }
                r_p.push(r);
            }
            rs.push(r_p);

            let final_val = recv!(self.ch, Layer, 0);
            gkr_bw_inc(T::SIZE * final_val.len());
            final_vals.push(final_val);

            if layer.num_phase > 1 && j != self.circ.depth - 2 {
                let scalar =T::random();
                self.ch.send(1, Tx::<T>::Random(scalar));
                scalars.push(scalar);
            } else {
                self.ch.send(1, Tx::<T>::Syn);
            }
        }


        // begin verify
let start = Instant::now();

        let bit_len = self.circ.layers[self.circ.depth - 1].bit_len;
        let log_grp_len = self.circ.layers[self.circ.depth - 1].log_grp_len;
        let mut fr = self.eval_poly_output(&outputs, &tau, bit_len, log_grp_len);

        // let mut coef_iter = coefs.iter();
        let mut idx = 0;
        let mut r_iter = rs.iter();
        let mut scalar_iter = scalars.iter();

        let mut scalar = T::from(0);
        let mut rs = &vec![tau];

        for j in 0..(self.circ.depth - 1) {
            let forword_id = self.circ.depth - 1 - j;
            let layer = &self.circ.layers[forword_id];
            let back_layer = &self.circ.layers[forword_id - 1];
            let bit_len = back_layer.bit_len;

            let old_rs = rs;
            rs = r_iter.next().unwrap();
            for p in 0..layer.num_phase {
                for i in 0..bit_len {
                    if coefs[idx].len() == 2 {
                        let tmp = fr - coefs[idx][0] - coefs[idx][0] - coefs[idx][1];
                        fr = rs[p][i] * ((rs[p][i] * coefs[idx][1]) + tmp) + coefs[idx][0];
                    } else {
                        let tmp = fr - coefs[idx][0] - coefs[idx][0] - coefs[idx][1] - coefs[idx][2];
                        fr = ((coefs[idx][2] * rs[p][i] + coefs[idx][1]) * rs[p][i] + tmp) * rs[p][i] + coefs[idx][0];
                    }
                    idx += 1;
                }
            }

            let pred =
                if old_rs.len() == 1 {
                    self.predicate_out(forword_id, old_rs, rs)
                } else {
                    self.predicate(forword_id, scalar, old_rs, rs)
                };
            // println!("{} verifier pred {:?}", j, pred);

            if layer.num_phase == 1 {
                if fr != final_vals[j][0] * (final_vals[j][0] * pred[0] + pred[1]) + pred[3] {
                    println!("fail at {} {}", j, self.circ.depth - 1 - j);
                    res = false;
                }
                if j != self.circ.depth - 2 {
                    fr = final_vals[j][0];
                }
            } else {
                if fr != (final_vals[j][0] * final_vals[j][1]) * pred[0] + final_vals[j][0] * pred[1] + final_vals[j][1] * pred[2] + pred[3] {
                    println!("fail at {} {}", j, self.circ.depth - 1 - j);
                    res = false;
                }

                if j != self.circ.depth - 2 {
                    scalar = *scalar_iter.next().unwrap();
                    fr = (T::from(1) - scalar) * final_vals[j][0] + scalar * final_vals[j][1];
                }
            }

        }

        let bit_len = self.circ.layers[0].bit_len;

        for p in 0..rs.len() {
            let mut tbl: Vec<T> = Vec::with_capacity(1 << bit_len);
            tbl.push(T::from(1));
            for k in 0..bit_len {
                let mask = 1 << k;
                for i in 0..mask {
                    let tbl_i_mask = tbl[i] * rs[p][k];
                    tbl.push(tbl_i_mask);
                    tbl[i] -= tbl_i_mask;
                }
            }
            for i in 0..(1 << bit_len) {
                final_vals[self.circ.depth - 2][p] -= tbl[i] * inputs[i];
            }

            if final_vals[self.circ.depth - 2][p] != T::from(0) {
                println!("fail at input");
                res = false;
            }
        }

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
gkr_vrf_inc(in_ms);

        res
    }

    pub fn zk_verify(&mut self, outputs: &Vec<T>, n_party: usize, tape_len: usize) -> bool {
        if self.ch.get_type() == ChannelType::MPCitH {
            self.zk_verify_mith(outputs, n_party, tape_len)
        } else if self.ch.get_type() == ChannelType::VOLEitH {
            self.zk_verify_vith(outputs, n_party, tape_len)
        } else {
            unreachable!()
        }
    }

    pub fn zk_verify_mith(&mut self, outputs: &Vec<T>, n_party: usize, tape_len: usize) -> bool {

        let mut res = true;

        let layer = &self.circ.layers[self.circ.depth - 1];
        let mut bit_len = layer.bit_len;

        let mut tau = vec![T::from(0); bit_len];
        let mut rs: Vec<Vec<Vec<T>>> = vec![];
        let mut scalars: Vec<T> = vec![];
        let mut final_vals: Vec<Vec<T>> = vec![];
        let mut coefs: Vec<Vec<T>> = vec![];

        for i in 0..bit_len {
            tau[i] = T::random();
        }
        self.ch.send(1, Tx::<T>::MultiRandom(tau.clone()));

        for j in 0..(self.circ.depth - 1) {
            let layer = &self.circ.layers[self.circ.depth - 1 - j];
            let back_layer = &self.circ.layers[self.circ.depth - 2 - j];
            bit_len = back_layer.bit_len;

            let mut r_p: Vec<Vec<T>> = Vec::with_capacity(layer.num_phase);
            for _ in 0..layer.num_phase {
                let mut r = vec![T::from(0); bit_len];
                for i in 0..bit_len {
                    let coef: Vec<T> = recv!(self.ch, Sumcheck, 0);
                    mpc_bw_inc(T::SIZE * 2);
                    coefs.push(coef);
                    r[i] = T::random();
                    self.ch.send(1, Tx::<T>::Random(r[i]));
                }
                r_p.push(r);
            }
            rs.push(r_p);


            if layer.num_phase > 1 && j != self.circ.depth - 2 {
                let final_val = recv!(self.ch, Layer, 0);
                mpc_bw_inc(T::SIZE * final_val.len());
                final_vals.push(final_val);

                let scalar = T::random();
                self.ch.send(1, Tx::<T>::Random(scalar));
                // if j == self.circ.depth - 2 {
                    // self.q_scalar = scalar.clone();
                // }
                scalars.push(scalar);
            } else if j != self.circ.depth - 2 {
                let final_val = recv!(self.ch, Layer, 0);
                mpc_bw_inc(T::SIZE * final_val.len());
                final_vals.push(final_val);

                self.ch.send(1, Tx::<T>::Syn);

            // final layer
            } else if layer.num_phase > 1 {
                self.v_ab = recv!(self.ch, FinalLayer, 0);
                mpc_bw_inc(T::SIZE * self.v_ab.len());

                self.q_scalar = T::random();
                self.ch.send(1, Tx::<T>::Random(self.q_scalar));

            } else {
                self.v_ab = recv!(self.ch, FinalLayer, 0);
                mpc_bw_inc(T::SIZE * self.v_ab.len());
                self.ch.send(1, Tx::<T>::Syn);
            }
        }

        let h_final = recv!(self.ch, FinalLayerOpen, 0);
        mpc_bw_inc(SHA256::SIZE);
        mpc_bw_hash_inc(1);
        self.ch.send(1, Tx::<T>::Syn);

        let poly_icr = recv!(self.ch, FLPCPCommit, 0);
        mpc_bw_inc(T::SIZE * poly_icr.len());

        let pr = T::random();
        self.ch.send(1, Tx::<T>::Random(pr));

        let (h, mut a, mut b, mut c) = recv!(self.ch, FLPCPOpen, 0);
        mpc_bw_inc(SHA256::SIZE + T::SIZE * 3);
        mpc_bw_hash_inc(1);
        self.ch.send(1, Tx::<T>::Syn);

        let missing_party: u16 = rand::thread_rng().gen::<u16>() % (n_party as u16);
        let _ = recv!(self.ch, @MPCitH, 0);
        self.ch.send(1, Tx::<T>::Party(vec![missing_party as usize]));

        let (stc, mut sto) = recv!(self.ch, Seed, 0);
        mpc_bw_inc(sto.proof_size);
        self.ch.send(1, Tx::<T>::Syn);


        // begin verify


        let bit_len = self.circ.layers[self.circ.depth - 1].bit_len;
        let log_grp_len = self.circ.layers[self.circ.depth - 1].log_grp_len;
        // not included in verifier time since linear in outputs
        let fr = self.eval_poly_output(&outputs, &tau, bit_len, log_grp_len);

let start = Instant::now();
let mut gkr_time = 0;
let mut gkr_rdl_time = 0; // to compare with virgo's circuit
let mut flpcp_prepare_time = 0;

        assert!(sto.verify(&stc));
        let mut tapes: Vec<Tape> = vec![];
        for p in 0..n_party {
            tapes.push(Tape::new(sto.get_seed(0, p), tape_len * T::SIZE));
        }

        // let mut fr_shares: Vec<T> = vec![T::from(0); n_party];
        let mut fr_shares: Vec<Vec<T>> = vec![vec![]; n_party];
        let mut va_shares: Vec<Vec<T>> = vec![vec![]; n_party];
        let mut vb_shares: Vec<Vec<T>> = vec![vec![]; n_party];
        let mut num_mul = 0;

        for pl in 0..n_party {
            fr_shares[pl].push(T::from(0));
        }
        fr_shares[0][0] = fr;

        let mut idx = 0;
        let mut r_iter = rs.iter();
        let mut scalar_iter = scalars.iter();
        let mut scalar: T = T::from(0);
        let mut rs = &vec![tau];

        for j in 0..(self.circ.depth - 1) {
            let forword_id = self.circ.depth - 1 - j;
            let layer = &self.circ.layers[forword_id];
            let back_layer = &self.circ.layers[forword_id - 1];
            let bit_len = back_layer.bit_len;

            let old_rs = rs;
            rs = r_iter.next().unwrap();
            for p in 0..layer.num_phase {
                for i in 0..bit_len {
                    if coefs[idx].len() == 2 {
                        // let tmp = fr - coefs[idx][0] - coefs[idx][0] - coefs[idx][1];
                        // fr = rs[p][i] * ((rs[p][i] * coefs[idx][1]) + tmp) + coefs[idx][0];
                        for pl in 0..n_party {
                            let mut coef0_share = tapes[pl].get::<T>();
                            let mut coef2_share = tapes[pl].get::<T>();
                            if pl == 0 {
                                coef0_share += coefs[idx][0];
                                coef2_share += coefs[idx][1];
                            }
                            let coef1_share = fr_shares[pl][j] - coef0_share - coef0_share - coef2_share;
                            fr_shares[pl][j] =
                                (coef2_share
                                 * rs[p][i] + coef1_share)
                                * rs[p][i] + coef0_share;
                        }
                    } else {
                        for pl in 0..n_party {
                            let mut coef0_share = tapes[pl].get::<T>();
                            let mut coef2_share = tapes[pl].get::<T>();
                            let mut coef3_share = tapes[pl].get::<T>();
                            if pl == 0 {
                                coef0_share += coefs[idx][0];
                                coef2_share += coefs[idx][1];
                                coef3_share += coefs[idx][2];
                            }
                            let coef1_share = fr_shares[pl][j] - coef0_share - coef0_share - coef2_share - coef3_share;
                            fr_shares[pl][j] =
                                ((coef3_share
                                  * rs[p][i] + coef2_share)
                                 * rs[p][i] + coef1_share)
                                * rs[p][i] + coef0_share;
                        }
                    }
                    idx += 1;
                }
            }

    let start0 = Instant::now();
            let pred =
                if old_rs.len() == 1 {
                    self.predicate_out(forword_id, old_rs, rs)
                } else {
                    self.predicate(forword_id, scalar, old_rs, rs)
                };
    let elapsed = start0.elapsed();
    gkr_rdl_time = elapsed.as_micros() as usize;
    gkr_time += elapsed.as_micros() as usize;

            if layer.num_phase == 1 {
                if j != self.circ.depth - 2 {
                    for pl in 0..n_party {
                        let mut va_share = tapes[pl].get::<T>();
                        if pl == 0 {
                            va_share += final_vals[j][0];
                        }
                        va_shares[pl].push(va_share * pred[0]);
                        vb_shares[pl].push(va_share);
                        fr_shares[pl][j] -= pred[1] * va_share;
                        if pl == 0 {
                            // constant term: should this be pred[3]? XXX
                            fr_shares[pl][j] -= pred[2];
                        }
                        fr_shares[pl].push(va_share);
                    }
                    num_mul += 1;
                } else {
                    let mut h = SHA256::from(0);
                    let mut output_shares: Vec<T> = vec![];
                    let mut missing = T::from(0);
                    for pl in 0..n_party {
                        output_shares.push(fr_shares[pl][num_mul]);
                        if pl != missing_party as usize {
                            missing += output_shares[pl];
                        }
                        fr_shares[pl].truncate(num_mul);
                    }
                    let tmp = pred[0] * self.v_ab[0] * self.v_ab[0] + pred[1] * self.v_ab[0] + pred[2];
                    output_shares[0] -= tmp;
                    if missing_party != 0 {
                        missing -= tmp;
                    }
                    output_shares[missing_party as usize] = T::from(0) - missing;
                    h.commit(&output_shares);
                    if h != h_final {
                        res = false;
                        println!("failed at input layer check");
                    }
                }
            } else {
                if j != self.circ.depth - 2 {
                    scalar = *scalar_iter.next().unwrap();

                    for pl in 0..n_party {
                        let mut va_share = tapes[pl].get::<T>();
                        let mut vb_share = tapes[pl].get::<T>();
                        if pl == 0 {
                            va_share += final_vals[j][0];
                            vb_share += final_vals[j][1];
                        }
                        va_shares[pl].push(va_share * pred[0]);
                        vb_shares[pl].push(vb_share);
                        fr_shares[pl][j] -= pred[1] * va_share + pred[2] * vb_share;
                        if pl == 0 {
                            fr_shares[pl][j] -= pred[3];
                        }
                        fr_shares[pl].push((T::from(1) - scalar) * va_share + scalar * vb_share);
                    }
                    num_mul += 1;
                } else {
                    let mut h = SHA256::from(0);
                    let mut output_shares: Vec<T> = vec![];
                    let mut missing = T::from(0);
                    for pl in 0..n_party {
                        output_shares.push(fr_shares[pl][num_mul]);
                        if pl != missing_party as usize {
                            missing += output_shares[pl];
                        }
                        fr_shares[pl].truncate(num_mul);
                    }
                    let tmp = pred[0] * self.v_ab[0] * self.v_ab[1] + pred[1] * self.v_ab[0] + pred[2] * self.v_ab[1] + pred[3];
                    output_shares[0] -= tmp;
                    if missing_party != 0 {
                        missing -= tmp;
                    }
                    output_shares[missing_party as usize] = T::from(0) - missing;
                    h.commit(&output_shares);
                    if h != h_final {
                        res = false;
                        println!("failed at input layer check");
                    }
                }
            }
        }

        // flpcp

    let start0 = Instant::now();
        let d = num_mul;
        let x_values_to_m = first_n_field_elements::<T>(d + 1);
        let precompute_to_m = precompute_lagrange_polynomials::<T>(&x_values_to_m);
        let x_values_to_2m = first_n_field_elements::<T>(2 * d + 1);
        let precompute_to_2m = precompute_lagrange_polynomials::<T>(&x_values_to_2m);
    let elapsed = start0.elapsed();
    flpcp_prepare_time += elapsed.as_micros() as usize;

        let mut p_shares: Vec<Vec<T>> = vec![vec![T::from(0); d]; n_party];

        for p in 0..n_party {
            va_shares[p].push(tapes[p].get::<T>());
            vb_shares[p].push(tapes[p].get::<T>());
            fr_shares[p].push(tapes[p].get::<T>());

            for _ in d..(2 * d + 1) {
                p_shares[p].push(tapes[p].get::<T>());
            }

            if p == 0 {
                for i in d..(2 * d + 1) {
                    p_shares[p][i] += poly_icr[i - d];
                }
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
                a_shares[p] = dot_product(&lagrange_poly_m_eval_r, &va_shares[p]);
                b_shares[p] = dot_product(&lagrange_poly_m_eval_r, &vb_shares[p]);
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

        if h_ != h {
            res = false;
            println!("fail at flpcp");
        }

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
mpc_vrf_inc(in_ms - gkr_time - flpcp_prepare_time);
if virgo_flag_get() > 0 {
    zkgkr_vrf_inc(gkr_time - gkr_rdl_time);
} else {
    zkgkr_vrf_inc(gkr_time);
}
// assert_eq!(tapes[0].pointer, tapes[0].cache.len());
println!("{} vs {}", tapes[0].pointer, tapes[0].cache.len());
// println!("tape pointer pos {}", tapes[0].pointer);

        res
    }

    pub fn zk_verify_vith(&mut self, outputs: &Vec<T>, n_party: usize, tape_len: usize) -> bool {
        let n = if n_party == 256 { 32 } else if n_party == 16 { 64 } else { unreachable!() };
        let m = n / 2;
        let log_n = ceil_log(n);
        let log_m = ceil_log(m);

        let mut res = true;

        let layer = &self.circ.layers[self.circ.depth - 1];
        let mut bit_len = layer.bit_len;

        let mut tau = vec![T::from(0); bit_len];
        let mut rs: Vec<Vec<Vec<T>>> = vec![];
        let mut scalars: Vec<T> = vec![];
        let mut final_vals: Vec<Vec<T>> = vec![];
        let mut coefs: Vec<Vec<T>> = vec![];

        for i in 0..bit_len {
            tau[i] = T::random();
        }
        self.ch.send(1, Tx::<T>::MultiRandom(tau.clone()));

        for j in 0..(self.circ.depth - 1) {
            let layer = &self.circ.layers[self.circ.depth - 1 - j];
            let back_layer = &self.circ.layers[self.circ.depth - 2 - j];
            bit_len = back_layer.bit_len;

            let mut r_p: Vec<Vec<T>> = Vec::with_capacity(layer.num_phase);
            for _ in 0..layer.num_phase {
                let mut r = vec![T::from(0); bit_len];
                for i in 0..bit_len {
                    let coef: Vec<T> = recv!(self.ch, Sumcheck, 0);
                    mpc_bw_inc(T::SIZE * 2);
                    coefs.push(coef);
                    r[i] = T::random();
                    self.ch.send(1, Tx::<T>::Random(r[i]));
                }
                r_p.push(r);
            }
            rs.push(r_p);


            if layer.num_phase > 1 && j != self.circ.depth - 2 {
                let final_val = recv!(self.ch, Layer, 0);
                mpc_bw_inc(T::SIZE * final_val.len());
                final_vals.push(final_val);

                let scalar = T::random();
                self.ch.send(1, Tx::<T>::Random(scalar));
                // if j == self.circ.depth - 2 {
                    // self.q_scalar = scalar.clone();
                // }
                scalars.push(scalar);
            } else if j != self.circ.depth - 2 {
                let final_val = recv!(self.ch, Layer, 0);
                mpc_bw_inc(T::SIZE * final_val.len());
                final_vals.push(final_val);

                self.ch.send(1, Tx::<T>::Syn);

            // final layer
            } else if layer.num_phase > 1 {
                self.v_ab = recv!(self.ch, FinalLayer, 0);
                mpc_bw_inc(T::SIZE * self.v_ab.len());

                self.q_scalar = T::random();
                self.ch.send(1, Tx::<T>::Random(self.q_scalar));

            } else {
                self.v_ab = recv!(self.ch, FinalLayer, 0);
                mpc_bw_inc(T::SIZE * self.v_ab.len());
                self.ch.send(1, Tx::<T>::Syn);
            }
        }

        let h_final = recv!(self.ch, FinalLayerOpen, 0);
        mpc_bw_inc(SHA256::SIZE);
        mpc_bw_hash_inc(1);
        self.ch.send(1, Tx::<T>::Syn);

        let poly_random = recv!(self.ch, FLPCPRandom, 0);
        mpc_bw_inc(T::SIZE * 3);
        self.ch.send(1, Tx::<T>::Syn);

        let poly_icr = recv!(self.ch, FLPCPCommit, 0);
        mpc_bw_inc(T::SIZE * poly_icr.len());

        let pr = T::random();
        self.ch.send(1, Tx::<T>::Random(pr));

        let (h, a, b, c) = recv!(self.ch, FLPCPOpen, 0);
        mpc_bw_inc(SHA256::SIZE + T::SIZE * 3);
        mpc_bw_hash_inc(1);
        self.ch.send(1, Tx::<T>::Syn);

        let _ = recv!(self.ch, @MPCitH, 0);
        let alpha = T::random();
        self.ch.send(1, Tx::<T>::Random(alpha));

        let (vith_c, vith_s, hv) = recv!(self.ch, VOLEOpen, 0);
        mpc_bw_inc(SHA256::SIZE + T::SIZE * (vith_c.len() + vith_s.len()));
        mpc_bw_hash_inc(1);

        let mut missings: Vec<usize> = vec![0; n];
        for i in 0..n {
            missings[i] = rand::thread_rng().gen::<usize>() % n_party;
        }
        self.ch.send(1, Tx::<T>::Party(missings.clone()));

        let (stc, mut sto) = recv!(self.ch, Seed, 0);
        mpc_bw_inc(sto.proof_size);
        self.ch.send(1, Tx::<T>::Syn);

        // begin verify

        let bit_len = self.circ.layers[self.circ.depth - 1].bit_len;
        let log_grp_len = self.circ.layers[self.circ.depth - 1].log_grp_len;
        // not included in verifier time since linear in outputs
        let mut tmp = self.eval_poly_output(&outputs, &tau, bit_len, log_grp_len);
        let mut tmp_s = tmp;
        let mut ptr_s = 0;

        let coset = Coset::<T>::init(log_n, log_m, 1);

let start = Instant::now();
let mut gkr_time = 0;
let mut gkr_rdl_time = 0; // to compare with virgo's circuit
let mut flpcp_prepare_time = 0;

        assert!(sto.verify(&stc));
        let mut tapes: Vec<Tape> = vec![];
        for i in 0..n {
            for p in 0..n_party {
                tapes.push(Tape::new(sto.get_seed(i, p), tape_len * T::SIZE));
            }
        }

        // verify hv = hash(aQ1+Q2-SGcD)
        let l = vith_s.len() / m;
        assert_eq!(l, tape_len / 2);
        let mut vith_q: Vec<T> = vec![T::from(0); 2 * l * n];

        if T::CHAR == 2 {
            for i in 0..(l * 2) {
                for j in 0..n {
                    for k in 0..n_party {
                        let t = tapes[j * n_party + k].get::<T>();
                        vith_q[i * n + j] += (T::from(missings[j] + 1) - T::from(k + 1)) * t;
                    }
                    if j >= m {
                        vith_q[i * n + j] -= vith_c[i * (n - m) + j - m] * T::from(missings[j] + 1);
                    }
                }
            }
        } else {
            let r = n / m;
            for i in 0..(l * 2) {
                for j in 0..n {
                    for k in 0..n_party {
                        let t = tapes[j * n_party + k].get::<T>();
                        vith_q[i * n + j] += (T::from(missings[j] + 1) - T::from(k + 1)) * t;
                    }
                    if j % r != 0 {
                        vith_q[i * n + j] -= vith_c[i * (n - m) + j - (j / r)] * T::from(missings[j] + 1);
                    }
                }
            }
        }


        let mut vith_v: Vec<T> = vec![T::from(0); l * n];
        for i in 0..l {
            let tmp = (coset.fft)(
                &coset, &(coset.ifft)(
                    &coset,
                    &vith_s[i * m..(i + 1) * m].to_vec(),
                    log_m), log_m, log_n);
            for j in 0..n {
                vith_v[i * n + j] =
                    alpha * vith_q[i * n + j] +
                    vith_q[i * n + j + l * n] -
                    tmp[j] * T::from(missings[j] + 1);
            }
        }

        let mut _hv = SHA256::from(0);
        _hv.commit(&vith_v);
        if _hv != hv {
            res = false;
            println!("failed at vole check");
        }

        let mut fr: Vec<T> = vec![];
        let mut va: Vec<T> = vec![];
        let mut vb: Vec<T> = vec![];
        let mut fr_s: Vec<T> = vec![];
        let mut va_s: Vec<T> = vec![];
        let mut vb_s: Vec<T> = vec![];
        let mut num_mul = 0;

        let mut idx = 0;
        let mut r_iter = rs.iter();
        let mut scalar_iter = scalars.iter();
        let mut scalar: T = T::from(0);
        let mut rs = &vec![tau];

        for j in 0..(self.circ.depth - 1) {
            let forword_id = self.circ.depth - 1 - j;
            let layer = &self.circ.layers[forword_id];
            let back_layer = &self.circ.layers[forword_id - 1];
            let bit_len = back_layer.bit_len;

            let old_rs = rs;
            rs = r_iter.next().unwrap();
            for p in 0..layer.num_phase {
                for i in 0..bit_len {
                    if coefs[idx].len() == 2 {
                        tmp = tmp - coefs[idx][0] - coefs[idx][0] - coefs[idx][1];
                        tmp = (coefs[idx][1] * rs[p][i] + tmp) * rs[p][i] + coefs[idx][0];
                        tmp_s = tmp_s - vith_s[ptr_s] - vith_s[ptr_s] - vith_s[ptr_s + 1];
                        tmp_s = (vith_s[ptr_s + 1] * rs[p][i] + tmp_s) * rs[p][i] + vith_s[ptr_s];
                        ptr_s += 2;
                    } else {
                        tmp = tmp - coefs[idx][0] - coefs[idx][0] - coefs[idx][1] - coefs[idx][2];
                        tmp = ((coefs[idx][2] * rs[p][i] + coefs[idx][1]) * rs[p][i] + tmp) * rs[p][i] + coefs[idx][0];
                        tmp_s = tmp_s - vith_s[ptr_s] - vith_s[ptr_s] - vith_s[ptr_s + 1] - vith_s[ptr_s + 2];
                        tmp_s = ((vith_s[ptr_s + 2] * rs[p][i] + vith_s[ptr_s + 1]) * rs[p][i] + tmp_s) * rs[p][i] + vith_s[ptr_s];
                        ptr_s += 3;
                    }
                    idx += 1;
                }
            }

    let start0 = Instant::now();
            let pred =
                if old_rs.len() == 1 {
                    self.predicate_out(forword_id, old_rs, rs)
                } else {
                    self.predicate(forword_id, scalar, old_rs, rs)
                };
    let elapsed = start0.elapsed();
    gkr_rdl_time = elapsed.as_micros() as usize;
    gkr_time += elapsed.as_micros() as usize;

            if layer.num_phase == 1 {
                if j != self.circ.depth - 2 {
                    va.push(final_vals[j][0] * pred[0]);
                    vb.push(final_vals[j][0]);
                    va_s.push(vith_s[ptr_s] * pred[0]);
                    vb_s.push(vith_s[ptr_s]);
                    fr.push(tmp - final_vals[j][0] * pred[1] + pred[2]);
                    fr_s.push(tmp_s - vith_s[ptr_s] * pred[1] + pred[2]);
                    tmp = final_vals[j][0];
                    tmp_s = vith_s[ptr_s];
                    ptr_s += 1;
                    num_mul += 1;
                } else {
                    let mut h = SHA256::from(0);

                    let real = pred[0] * self.v_ab[0] * self.v_ab[0] + pred[1] * self.v_ab[0] + pred[2];

                    tmp -= real;
                    tmp_s -= real;
                    tmp = tmp_s + alpha * tmp;

                    h.commit(&[tmp]);
                    if h != h_final {
                        res = false;
                        println!("failed at input layer check");
                    }
                }
            } else {
                if j != self.circ.depth - 2 {
                    scalar = *scalar_iter.next().unwrap();

                    va.push(final_vals[j][0] * pred[0]);
                    vb.push(final_vals[j][1]);
                    fr.push(tmp - pred[1] * final_vals[j][0] - pred[2] * final_vals[j][1] - pred[3]);
                    va_s.push(vith_s[ptr_s] * pred[0]);
                    vb_s.push(vith_s[ptr_s + 1]);
                    fr_s.push(tmp_s - pred[1] * vith_s[ptr_s] - pred[2] * vith_s[ptr_s + 1] - pred[3]);
                    tmp = (T::from(1) - scalar) * final_vals[j][0] + scalar * final_vals[j][1];
                    tmp_s = (T::from(1) - scalar) * vith_s[ptr_s] + scalar * vith_s[ptr_s + 1];
                    ptr_s += 2;
                    num_mul += 1;
                } else {
                    let mut h = SHA256::from(0);
                    let real = pred[0] * self.v_ab[0] * self.v_ab[1] + pred[1] * self.v_ab[0] + pred[2] * self.v_ab[1] + pred[3];
                    tmp -= real;
                    tmp_s -= real;
                    tmp = tmp_s + alpha * tmp;
                    h.commit(&[tmp]);
                    if h != h_final {
                        res = false;
                        println!("failed at input layer check");
                    }
                }
            }
        }

        // flpcp
        // let poly_random = recv!(self.ch, FLPCPRandom, 0);
        // let poly_icr = recv!(self.ch, FLPCPCommit, 0);
        // let (h, mut a, mut b, mut c) = recv!(self.ch, FLPCPOpen, 0);

    let start0 = Instant::now();
        let d = num_mul;
        let x_values_to_m = first_n_field_elements::<T>(d + 1);
        let precompute_to_m = precompute_lagrange_polynomials::<T>(&x_values_to_m);
        let x_values_to_2m = first_n_field_elements::<T>(2 * d + 1);
        let precompute_to_2m = precompute_lagrange_polynomials::<T>(&x_values_to_2m);
    let elapsed = start0.elapsed();
    flpcp_prepare_time += elapsed.as_micros() as usize;

        let mut poly: Vec<T> = vec![T::from(0); 2 * d + 1];
        let mut poly_s: Vec<T> = vec![T::from(0); 2 * d + 1];

        fr.push(poly_random[0]);
        va.push(poly_random[1]);
        vb.push(poly_random[2]);
        fr_s.push(vith_s[ptr_s]);
        va_s.push(vith_s[ptr_s + 1]);
        vb_s.push(vith_s[ptr_s + 2]);
        ptr_s += 3;

        for i in d..(2 * d + 1) {
            poly[i] = poly_icr[i - d];
            poly_s[i] = vith_s[ptr_s + i - d];
            // println!("v poly {} {:?}", i, poly[i] + poly_s[i]);
            // println!("v poly {} {:?}", i, poly_s[i]);
        }
        // ptr_s += d + 1;

        let mut lagrange_poly_m_eval_r = vec![T::from(0); d + 1];
        let mut lagrange_poly_2m_eval_r = vec![T::from(0); 2 * d + 1];
        for k in 0..(d + 1) {
            lagrange_poly_m_eval_r[k] = poly_eval(&precompute_to_m[k], pr);
        }
        for k in 0..(2 * d + 1) {
            lagrange_poly_2m_eval_r[k] = poly_eval(&precompute_to_2m[k], pr);
        }

        let mut output = vec![T::from(0); 4];
        let mut output_s = vec![T::from(0); 4];
        output[0] = dot_product(&lagrange_poly_m_eval_r, &va) - a;
        output_s[0] = dot_product(&lagrange_poly_m_eval_r, &va_s) - a;
        output[1] = dot_product(&lagrange_poly_m_eval_r, &vb) - b;
        output_s[1] = dot_product(&lagrange_poly_m_eval_r, &vb_s) - b;
        output[2] = dot_product(&lagrange_poly_m_eval_r, &fr) - c;
        output_s[2] = dot_product(&lagrange_poly_m_eval_r, &fr_s) - c;
        output[3] = dot_product(&lagrange_poly_2m_eval_r, &poly) - a * b + c;
        output_s[3] = dot_product(&lagrange_poly_2m_eval_r, &poly_s) - a * b + c;

        for i in 0..4 {
            output[i] = alpha * output[i] + output_s[i];
        }

        let mut h_ = SHA256::from(0);
        h_.commit(&output);

        if h_ != h {
            res = false;
            println!("fail at flpcp");
        }

let elapsed = start.elapsed();
let in_ms = elapsed.as_micros() as usize;
mpc_vrf_inc(in_ms - gkr_time - flpcp_prepare_time);
if virgo_flag_get() > 0 {
    zkgkr_vrf_inc(gkr_time - gkr_rdl_time);
} else {
    zkgkr_vrf_inc(gkr_time);
}
// assert_eq!(tapes[0].pointer, tapes[0].cache.len());
// println!("{} vs {}", tapes[0].pointer, tapes[0].cache.len());
// println!("tape pointer pos {}", tapes[0].pointer);

        res
    }
}

