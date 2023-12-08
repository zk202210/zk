use rand::Rng;
use std::sync::Arc;

use crate::field::Field;
use crate::field::Coset;
use crate::gkr::BKT;
use crate::util::*;

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum GateType {
    ADD,        // x₀ + x₁
    MUL,        // x₀ × x₁
    SUB,
    XOR,        // x₀ + x₁ - 2x₀x₁
    // OR,
    NOT,        // 1 - x₀
    NAAB,       // (1 - x₀) × x₁
    EXPSUM,     // x₀ + x₁ + x₂ + ... + xₙ₋₁
    SUM,        // x₀ + x₁ + x₂ + ... + xₙ₋₁
    SUMMUL,     // (x₀ + x₁ + x₂ + ... + xₙ₋₂) × xₙ₋₁
    SUMLIN,     // (x₀ + x₁ + x₂ + ... + xₙ₋₂ + s₂) × xₙ₋₁ × s₀ + s₁
    SUMCONS,    // x₀ + x₁ + x₂ + ... + xₙ₋₂ + xₙ₋₁ + s₀
    ADDCONS,    // x₀ + s₀
    ADDMUL,     // (x₀ + s₀) × (x₁ + s₁)
    P1D2C0,     // x(x + s), phase 1, degree 2, const term 0
    SELECTOR,   // (1 - x₀)x₁ + x₀x₂
    KECCAKX,   // (1 - x₀)x₁ + x₂
    NKECCAKX, // 1 - ((1 - x₀)x₁ + x₂)
    // EXPSUM,
    // GEXPSUM,
    RELAY,      // x₀
    LINEAR,     // s₀x₀ + s₁x₁
    FLINEAR,    // s₀x₀ + s₁x₁ + s₂
    SCALAR,     // s₀x₀
    SCAMUL,     // s₀x₀x₁
    IN,
    DUMMY,
}

impl GateType {
    pub const AND: GateType = GateType::MUL;
    pub const LEN: usize = GateType::DUMMY as usize + 1;
    // pub const BITTEST: GateType = GateType::NAAB;
    // pub const RELAY: GateType = GateType::SUM;
    // pub const ONESUB: GateType = GateType::NOT;
}

#[derive(Clone)]
pub struct Gate {
    pub gtype: GateType,
    pub inputs: Vec<u32>,
    pub utils: Vec<u32>,
}

impl Gate {
    pub fn default() -> Self {
        Self {
            gtype: GateType::DUMMY,
            inputs:  vec![],
            utils: vec![],
        }
    }
}

type Predicate<'a, T> = Arc<dyn Fn(&Vec<T>, &Vec<T>, &Vec<T>) -> [T; 4] + 'a + Sync + Send>;
type SumcheckInit<'a, T> = Arc<dyn Fn(&mut BKT<T>, &Layer<'a, T>, usize) + 'a + Sync + Send>;
type SumcheckEval<'a, T> = Arc<dyn Fn(&mut BKT<T>, usize, usize) -> Vec<T> + 'a + Sync + Send>;
type SumcheckUpdate<'a, T> = Arc<dyn Fn(&mut BKT<T>, usize, usize, T) + 'a + Sync + Send>;
type Evaluate<'a, T> = Arc<dyn Fn(&Layer<'a, T>, &Vec<T>) -> Vec<T> + 'a + Sync + Send>;

#[derive(Clone)]
pub struct Layer<'a, T: Field> {
    pub gates: Vec<Gate>,
    pub values: Vec<T>,
    pub bit_len: usize,
    pub log_grp_len: usize,
    pub utils: Vec<T>,
    pub num_table: Vec<usize>,
    pub val_degree: Vec<usize>,
    pub num_phase: usize,
    pub sumcheck_init: Vec<SumcheckInit<'a, T>>,
    pub sumcheck_eval: Vec<SumcheckEval<'a, T>>,
    pub sumcheck_update: Vec<SumcheckUpdate<'a, T>>,
    // pub generate_point_table: Vec<GeneratePointTable<T>>,
    pub predicate: Vec<Predicate<'a, T>>,
    pub evaluate: Vec<Evaluate<'a, T>>,
}

impl<'a, T: Field> Layer<'a, T> {
    pub fn default() -> Self {
        Self {
            gates: vec![],
            values: vec![],
            bit_len: 0,
            log_grp_len: 0,
            utils: vec![],
            num_table: vec![2, 4],
            val_degree: vec![1, 1],
            num_phase: 2,
            sumcheck_init: vec![
                Arc::new(BKT::<T>::sumcheck_phase1_init),
                Arc::new(BKT::<T>::sumcheck_phase2_init)],
            sumcheck_eval: vec![
                Arc::new(BKT::<T>::sumcheck_phase1_eval),
                Arc::new(BKT::<T>::sumcheck_phase2_eval)],
            sumcheck_update: vec![
                Arc::new(BKT::<T>::sumcheck_phase1_update),
                Arc::new(BKT::<T>::sumcheck_phase2_update)],
            // generate_point_table: vec![
            //     Arc::new(BKT::<T>::generate_point_table),
            //     Arc::new(BKT::<T>::generate_point_table)],
            predicate: vec![],
            evaluate: vec![Arc::new(Layer::<'a, T>::default_evaluate)],
        }
    }

    pub fn new(bit_len: usize) -> Self {
        Self {
            gates: vec![Gate::default(); 1 << bit_len],
            values: vec![T::from(0); 1 << bit_len],
            bit_len,
            log_grp_len: 0,
            utils: vec![],
            num_table: vec![2, 4],
            val_degree: vec![1, 1],
            num_phase: 2,
            sumcheck_init: vec![
                Arc::new(BKT::<T>::sumcheck_phase1_init),
                Arc::new(BKT::<T>::sumcheck_phase2_init)],
            sumcheck_eval: vec![
                Arc::new(BKT::<T>::sumcheck_phase1_eval),
                Arc::new(BKT::<T>::sumcheck_phase2_eval)],
            sumcheck_update: vec![
                Arc::new(BKT::<T>::sumcheck_phase1_update),
                Arc::new(BKT::<T>::sumcheck_phase2_update)],
            // generate_point_table: vec![
            //     Arc::new(BKT::<T>::generate_point_table),
            //     Arc::new(BKT::<T>::generate_point_table)],
            predicate: vec![],
            evaluate: vec![Arc::new(Layer::<T>::default_evaluate)],
        }
    }

    pub fn default_evaluate(&self, prev_vals: &Vec<T>) -> Vec<T> {
        let mut res = Vec::with_capacity(1 << self.bit_len);
        for j in 0..(1 << self.bit_len) {
            let gate = &self.gates[j];
            // let left = gate.inputs[0] as usize;
            // let right = gate.inputs[1] as usize;
            let inp = &gate.inputs;
            match gate.gtype {
                GateType::ADD =>
                    res.push(
                        prev_vals[inp[0] as usize] +
                        prev_vals[inp[1] as usize]),
                GateType::SUB =>
                    res.push(
                        prev_vals[inp[0] as usize] -
                        prev_vals[inp[1] as usize]),
                GateType::XOR =>
                    res.push(
                        prev_vals[inp[0] as usize] +
                        prev_vals[inp[1] as usize] -
                        T::from(2) *
                        prev_vals[inp[0] as usize] *
                        prev_vals[inp[1] as usize]),
                GateType::LINEAR =>
                    res.push(
                        prev_vals[inp[0] as usize] * self.utils[gate.utils[0] as usize] +
                        prev_vals[inp[1] as usize] * self.utils[gate.utils[1] as usize]),
                GateType::ADDMUL =>
                    res.push(
                        (prev_vals[inp[0] as usize] + self.utils[gate.utils[0] as usize]) *
                        (prev_vals[inp[1] as usize] + self.utils[gate.utils[1] as usize])),
                GateType::FLINEAR =>
                    res.push(
                        prev_vals[inp[0] as usize] * self.utils[gate.utils[0] as usize] +
                        prev_vals[inp[1] as usize] * self.utils[gate.utils[1] as usize] +
                        self.utils[gate.utils[2] as usize]),
                GateType::MUL =>
                    res.push(
                        prev_vals[inp[0] as usize] *
                        prev_vals[inp[1] as usize]),
                GateType::NAAB =>
                    res.push(
                        (T::from(1) - prev_vals[inp[0] as usize]) *
                        prev_vals[inp[1] as usize]),
                GateType::NOT =>
                    res.push(
                        T::from(1) - prev_vals[inp[0] as usize]),
                GateType::ADDCONS =>
                    res.push(
                        self.utils[gate.utils[0] as usize] + prev_vals[inp[0] as usize]),
                GateType::SUMCONS => {
                    let mut sum = T::from(0);
                    for i in 0..gate.inputs.len() {
                        sum += prev_vals[gate.inputs[i] as usize];
                    }
                    res.push(self.utils[gate.utils[0] as usize] + sum);
                }
                GateType::EXPSUM => {
                    let mut sum = T::from(0);
                    let u = gate.inputs[0] as usize;
                    let v = gate.inputs[1] as usize;
                    let mut pow = T::from(1);
                    for i in u..=v {
                        sum += pow * prev_vals[i as usize];
                        pow = pow + pow;
                    }
                    res.push(sum);
                }
                GateType::SUM => {
                    let mut sum = T::from(0);
                    for i in 0..gate.inputs.len() {
                        sum += prev_vals[gate.inputs[i] as usize];
                    }
                    res.push(sum);
                }
                GateType::SUMLIN => {
                    let mut sum = T::from(0);
                    for i in 0..(gate.inputs.len() - 1) {
                        sum += prev_vals[gate.inputs[i] as usize];
                    }
                    let a = self.utils[gate.utils[0] as usize];
                    let b = self.utils[gate.utils[1] as usize];
                    let c = self.utils[gate.utils[2] as usize];
                    res.push((sum + c) * prev_vals[*gate.inputs.last().unwrap() as usize] * a + b);
                }
                GateType::SUMMUL => {
                    let mut sum = T::from(0);
                    for i in 0..(gate.inputs.len() - 1) {
                        sum += prev_vals[gate.inputs[i] as usize];
                    }
                    res.push(sum * prev_vals[*gate.inputs.last().unwrap() as usize]);
                }
                GateType::RELAY =>
                    res.push(
                        prev_vals[inp[0] as usize]),
                GateType::SCALAR =>
                    res.push(
                        prev_vals[inp[0] as usize]
                        * self.utils[
                        gate.utils[0] as usize]),
                GateType::SCAMUL =>
                    res.push(
                        prev_vals[inp[0] as usize] *
                        prev_vals[inp[1] as usize] * self.utils[gate.utils[0] as usize]),
                GateType::P1D2C0 =>
                    res.push(prev_vals[inp[0] as usize] *
                             (prev_vals[inp[0] as usize] + self.utils[gate.utils[0] as usize])),
                GateType::SELECTOR =>
                    res.push(
                        (T::from(1) - prev_vals[inp[0] as usize]) *
                        prev_vals[inp[1] as usize] +
                        prev_vals[inp[0] as usize] *
                        prev_vals[inp[2] as usize]),
                GateType::KECCAKX =>
                    res.push(
                        (T::from(1) - prev_vals[inp[0] as usize]) *
                        prev_vals[inp[1] as usize] +
                        prev_vals[inp[2] as usize]),
                GateType::NKECCAKX =>
                    res.push(T::from(1) - (
                        (T::from(1) - prev_vals[inp[0] as usize]) *
                        prev_vals[inp[1] as usize] +
                        prev_vals[inp[2] as usize])),
                GateType::DUMMY =>
                    res.push(T::from(0)),
                GateType::IN =>
                    res.push(prev_vals[inp[0] as usize]),
            }
        }
        res
    }
}


#[derive(Clone)]
pub struct Circuit<'a, T: Field> {
    pub layers: Vec<Layer<'a, T>>,
    pub depth: usize,
    pub inputs: Vec<T>,
    // pub predicate: Vec<Predicate<'a, T>>,
}

impl<'a, T: Field> Circuit<'a, T> {
    pub fn new(depth: usize) -> Self {
        Self {
            layers: vec![],
            inputs: vec![],
            // predicate: vec![],
            depth,
        }
    }

    pub fn read(file_path: &str) -> Self {
        use std::fs::File;
        use std::io::{self, BufRead};

        let file = File::open(file_path).unwrap();

        // let contents = file.read_to_string(& mut s)?;

        let contents = io::BufReader::new(file).lines();

        let mut i = 0;
        let mut c = Self::new(0);
        for line in contents {
            let line = line.unwrap();
            if i == 0 {
                let depth = line.parse::<usize>().unwrap();
                c = Self::new(depth);
            } else {
                let mut j = 0;
                let mut gtype_num = 0;
                let mut idx = 0;
                let mut l = 0;
                for word in line.split_whitespace() {
                    let num = word.parse::<usize>().unwrap();
                    if j == 0 {
                        c.layers.push(Layer::<T>::default());
                        assert_eq!(c.layers.len(), i);
                        let layer = &mut c.layers[i - 1];
                        layer.bit_len = ceil_log(num);
                        layer.gates = vec![Gate::default(); num];
                        layer.values = vec![T::from(0); num];
                    } else if j % 4 == 0 {
                        let r = num;
                        let layer = &mut c.layers[i - 1];
                        match gtype_num {
                            0 => {
                                layer.gates[idx].gtype = GateType::ADD;
                                layer.gates[idx].inputs = vec![l as u32, r as u32];
                            }
                            1 => {
                                layer.gates[idx].gtype = GateType::MUL;
                                layer.gates[idx].inputs = vec![l as u32, r as u32];
                            }
                            2 => {
                                layer.gates[idx].gtype = GateType::DUMMY;
                            }
                            3 => {
                                layer.gates[idx].gtype = GateType::IN;
                                layer.gates[idx].inputs.push(idx as u32);
                                layer.values[idx] = T::from(l);
                            }
                            4 => {
                                layer.gates[idx].gtype = GateType::RELAY;
                                layer.gates[idx].inputs = vec![l as u32];
                            }
                            5 => {
                                layer.gates[idx].gtype = GateType::SUM;
                                for i in l..r {
                                    layer.gates[idx].inputs.push(i as u32);
                                }
                            }
                            6 => {
                                layer.gates[idx].gtype = GateType::NOT;
                                layer.gates[idx].inputs = vec![l as u32];
                            }
                            7 => {
                                layer.gates[idx].gtype = GateType::SUB;
                                layer.gates[idx].inputs = vec![l as u32, r as u32];
                            }
                            8 => {
                                layer.gates[idx].gtype = GateType::XOR;
                                // layer.gates[idx].gtype = GateType::MUL;
                                layer.gates[idx].inputs = vec![l as u32, r as u32];
                            }
                            9 => {
                                layer.gates[idx].gtype = GateType::NAAB;
                                layer.gates[idx].inputs = vec![l as u32, r as u32];
                            }
                            10 => {
                                layer.gates[idx].gtype = GateType::RELAY;
                                layer.gates[idx].inputs = vec![l as u32];
                            }
                            // 12 => {
                            //     layer.gates[idx].gtype = GateType::SUM;
                            //     for i in l..r {
                            //         layer.gates[idx].inputs.push(i as u32);
                            //     }
                            // }
                            12 => {
                                layer.gates[idx].gtype = GateType::EXPSUM;
                                layer.gates[idx].inputs = vec![l as u32, r as u32];
                            }
                            13 => {
                                layer.gates[idx].gtype = GateType::NAAB;
                                layer.gates[idx].inputs = vec![l as u32, l as u32];
                            }
                            _ => {
                                println!("unknown type {}", gtype_num);
                                unreachable!();
                            }
                        }
                    } else if j % 4 == 1 {
                        gtype_num = num;
                    } else if j % 4 == 2 {
                        idx = num;
                    } else if j % 4 == 3 {
                        l = num;
                    }
                    j += 1;
                }
                assert_eq!(j % 4, 1);
            }
            i += 1;
        }
        c
    }

    pub fn read_with_meta(cfile_path: &str, mfile_path: &str) -> Self {
        let mut c = Self::read(cfile_path);

        use std::fs::File;
        use std::io::{self, BufRead};

        let file = File::open(mfile_path).unwrap();

        let contents = io::BufReader::new(file).lines();

        let mut i = 0;
        for line in contents {
            let line = line.unwrap();
            if i != 0 {
                let mut j = 0;
                for word in line.split_whitespace() {
                    if j == 4 {
                        let num = word.parse::<usize>().unwrap();
                        c.layers[i].log_grp_len = num;
                    }
                    j += 1;
                }
            }
            i += 1;
        }
        c
    }

    // TODO specify phase number
    pub fn random(depth: usize, minlen: usize, maxlen: usize, log_input: usize) -> Self {
        let mut rng = rand::thread_rng();
        let mut c = Self::new(depth);
        c.layers.push(Layer::<T>::default());
        c.layers[0].bit_len = log_input;
        let size = 1 << log_input;
        // let len: usize = rng.gen_range(minlen..=maxlen);
        // let size: usize = rng.gen_range((1 << (len - 1))..=(1 << len));
        for i in 0..size {
            c.layers[0].gates.push(Gate::default());
            c.layers[0].gates[i].gtype = GateType::IN;
            c.layers[0].gates[i].inputs.push(i as u32);
            c.layers[0].values.push(T::from(0));
        }
        // for _ in size..(1 << len) {
        //     c.layers[0].gates.push(Gate::default());
        //     c.layers[0].values.push(T::from(0));
        // }
        let mut size0 = size;
        for d in 1..(c.depth) {
            c.layers.push(Layer::<T>::default());
            let len: usize = rng.gen_range(minlen..=maxlen);
            let size: usize = rng.gen_range((1 << (len - 1))..=(1 << len));
            c.layers[d].bit_len = len;
            for i in 0..size {
                c.layers[d].gates.push(Gate::default());
                match rng.gen_range(0..=19) as usize {
                    0 => {
                        c.layers[d].gates[i].gtype = GateType::ADD;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    1 => {
                        c.layers[d].gates[i].gtype = GateType::SUB;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    2 => {
                        c.layers[d].gates[i].gtype = GateType::XOR;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    3 => {
                        c.layers[d].gates[i].gtype = GateType::MUL;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    4 => {
                        c.layers[d].gates[i].gtype = GateType::NOT;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    5 => {
                        c.layers[d].gates[i].gtype = GateType::EXPSUM;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    }
                    6 => {
                        c.layers[d].gates[i].gtype = GateType::RELAY;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    7 => {
                        c.layers[d].gates[i].gtype = GateType::LINEAR;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[1] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                    },
                    8 => {
                        c.layers[d].gates[i].gtype = GateType::NAAB;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    9 => {
                        c.layers[d].gates[i].gtype = GateType::SCALAR;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                    },
                    10 => {
                        c.layers[d].gates[i].gtype = GateType::SCAMUL;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                    }
                    11 => {
                        c.layers[d].gates[i].gtype = GateType::SUM;
                        let n_sum: usize = rng.gen_range(4..8);
                        for _ in 0..n_sum {
                            c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        }
                    }
                    12 => {
                        c.layers[d].gates[i].gtype = GateType::FLINEAR;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[1] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[2] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                    },
                    13 => {
                        c.layers[d].gates[i].gtype = GateType::SUMMUL;
                        let n_sum: usize = rng.gen_range(4..8);
                        for _ in 0..n_sum {
                            c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        }
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    }
                    14 => {
                        c.layers[d].gates[i].gtype = GateType::SUMLIN;
                        let n_sum: usize = rng.gen_range(4..8);
                        for _ in 0..n_sum {
                            c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        }
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);

                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[1] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[2] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                    }
                    15 => {
                        c.layers[d].gates[i].gtype = GateType::ADDCONS;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                    },
                    16 => {
                        c.layers[d].gates[i].gtype = GateType::ADDMUL;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[1] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                    },
                    17 => {
                        c.layers[d].gates[i].gtype = GateType::SUMCONS;
                        let n_sum: usize = rng.gen_range(4..8);
                        for _ in 0..n_sum {
                            c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        }
                        let util_len = c.layers[d].utils.len();
                        c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                        if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                            c.layers[d].utils.push(T::random());
                        }
                    },
                    18 => {
                        c.layers[d].gates[i].gtype = GateType::KECCAKX;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    19 => {
                        c.layers[d].gates[i].gtype = GateType::NKECCAKX;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    20 => {
                        c.layers[d].gates[i].gtype = GateType::SELECTOR;
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                        c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    },
                    // 15 => {
                    //     c.layers[d].gates[i].gtype = GateType::P1D2C0;
                    //     c.layers[d].gates[i].inputs.push(rng.gen_range(0..size0) as u32);
                    //     let util_len = c.layers[d].utils.len();
                    //     c.layers[d].gates[i].utils.push(rng.gen_range(0..=util_len) as u32);
                    //     if c.layers[d].gates[i].utils[0] >= c.layers[d].utils.len() as u32 {
                    //         c.layers[d].utils.push(T::random());
                    //     }
                    // },
                    _ => unreachable!(),
                }
                c.layers[d].values.push(T::from(0));
            }
            for _ in size..(1 << len) {
                c.layers[d].gates.push(Gate::default());
                c.layers[d].values.push(T::from(0));
            }
            size0 = size;
        }
        c
    }

    pub fn continue_eval(&mut self, start_depth: usize) {
        for i in start_depth..self.depth {
            let (prev_layer, layers) = self.layers.split_at_mut(i);
            let prev_values = &prev_layer.last().unwrap().values;
            let mut layer = &mut layers[0];
            layer.values = layer.evaluate[0](&layer, prev_values);
        }
    }

    pub fn continue_eval_only(&mut self, start_depth: usize) {
        let mut values = self.layers[start_depth - 1].values.to_vec();
        for i in start_depth..self.depth {
            let layer = &self.layers[i];
            values = layer.evaluate[0](&layer, &values);
        }
        let depth = self.depth;
        self.layers[depth - 1].values = values;
    }

    pub fn eval(&mut self, inputs: &Vec<T>) {
        self.inputs = inputs.to_vec();
        self.layers[0].values = self.inputs.to_vec();
        self.continue_eval(1);
    }

    pub fn clean(&mut self) {
        self.inputs = vec![];
        self.layers = vec![];
    }
}

impl std::fmt::Display for Gate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.gtype {
            GateType::ADD => write!(f, "add({}, {})", self.inputs[0], self.inputs[1]),
            GateType::SUB => write!(f, "sub({}, {})", self.inputs[0], self.inputs[1]),
            GateType::MUL => write!(f, "mul({}, {})", self.inputs[0], self.inputs[1]),
            GateType::IN => write!(f, "in({})", self.inputs[0]),
            GateType::SCALAR => write!(f, "scalar({}*{})", self.utils[0], self.inputs[0]),
            GateType::RELAY => write!(f, "relay({})", self.inputs[0]),
            GateType::NOT => write!(f, "not({})", self.inputs[0]),
            GateType::LINEAR => write!(f, "linear({}*{}, {}*{})", self.utils[0], self.inputs[0], self.utils[1], self.inputs[1]),
            GateType::DUMMY => write!(f, "dummy"),
            GateType::ADDMUL => write!(f, "addmul({}+{}, {}+{})", self.utils[0], self.inputs[0], self.utils[1], self.inputs[1]),
            GateType::SUM => write!(f, "sum({:?})", self.inputs),
            GateType::SUMCONS => write!(f, "sumcons(({:?}) + u{})", self.inputs, self.utils[0]),
            GateType::SUMLIN => {
                let s = self.inputs.len();
                write!(f, "sumlin(({:?} + u{}) * {} * u{} + u{}", self.inputs[0..(s-1)].to_vec(), self.inputs[s - 1], self.utils[2], self.utils[0], self.utils[1])
            },
            _ => write!(f, ""),
        }
    }
}

impl<'a, T: Field> std::fmt::Display for Layer<'a, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..(1 << self.bit_len) {
            write!(f, "{}:{}; ", self.gates[i], self.values[i]).unwrap();
        }
        write!(f, "")
    }
}

impl<'a, T: Field> std::fmt::Display for Circuit<'a, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for d in 0..(self.depth) {
            write!(f, "{}\n", self.layers[d]).unwrap();
        }
        write!(f, "")
    }
}

fn pred_ifft<T: Field>(len: usize, factors: Vec<T>, tau: &Vec<T>, rho: &Vec<T>, _sgm: &Vec<T>) -> [T; 4] {
    // rho[hi | 0 | lo: len], sgm[hi | 1 | lo: len]
    // -> tau[hi | 0 | lo: len] : LINEAR gate w/ vanishing(hi)
    // -> tau[hi | 1 | lo: len] : ADD gate

    let mut lo = T::from(1);
    for i in 0..len {
        lo *= T::from(1) + rho[i] + tau[i];
    }
    // rho=0, sgm=1
    // lo *= (T::from(1) - rho[len]) * sgm[len];

    // multiplication tree
    let log_len = ceil_log(rho.len() - len - 1);
    let mut cache = vec![T::from(1); 2 << log_len];
    // count number of 0s
    let mut cnt: isize = -1;
    for i in (len + 1)..rho.len() {
        let idx = (i - len - 1) + (1 << log_len);
        cache[idx] = T::from(1) + rho[i] + tau[i];

        if cache[idx] == T::from(0) && cnt == -1 {
            cnt = i as isize;
        } else if cache[idx] == T::from(0) && cnt > 0 {
            return [T::from(0); 4];
        }
    }

    if cnt > 0 {
        // XXX
        let cnt = cnt as usize;
        let mut hi_lin1 = factors[cnt - len - 1] * rho[cnt] * tau[cnt];
        for i in (len + 1)..rho.len() {
            let idx = (i - len - 1) + (1 << log_len);
            if i != cnt {
                hi_lin1 *= cache[idx];
            }
        }
        lo *= T::from(1) - tau[len];
        return [T::from(0), hi_lin1 * lo, hi_lin1 * lo, T::from(0)];
    }

    // compute mul tree
    for i in (1..(1 << log_len)).rev() {
        cache[i] = cache[2 * i] * cache[2 * i + 1];
    }

    let add = lo * tau[len] * cache[1];
    lo *= T::from(1) - tau[len];
    let mut hi_lin0 = cache[1];
    let mut hi_lin1 = T::from(0);

    // compute products with missing element
    let mut last: usize = (1 << log_len) + (rho.len() - len - 1) - 1;
    last += 1 - (last & 1);
    for i in 4..(last + 1) {
        cache[i] = cache[i] * cache[(i / 2) ^ 1];
    }

    for i in (len + 1)..rho.len() {
        let idx = (i - len - 1) + (1 << log_len);
        hi_lin1 += factors[i - len - 1] * rho[i] * tau[i] *
            cache[idx ^ 1];
    }

    hi_lin0 += hi_lin1;

    [T::from(0), (hi_lin0 * lo + add) * (T::from(1) - rho[len]) +  (hi_lin1 * lo + add) * rho[len], T::from(0), T::from(0)]
}

// fn _pred_ifft<T: Field>(len: usize, factors: Vec<T>, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
//     // rho[hi | 0 | lo: len], sgm[hi | 1 | lo: len]
//     // -> tau[hi | 0 | lo: len] : LINEAR gate w/ vanishing(hi)
//     // -> tau[hi | 1 | lo: len] : ADD gate
// 
//     let mut lo = T::from(1);
//     for i in 0..len {
//         lo *=
//             (T::from(1) - rho[i]) *
//             (T::from(1) - sgm[i]) *
//             (T::from(1) - tau[i]) +
//             rho[i] * sgm[i] * tau[i];
//     }
//     // rho=0, sgm=1
//     lo *= (T::from(1) - rho[len]) * sgm[len];
// 
//     // multiplication tree
//     let log_len = ceil_log(rho.len() - len - 1);
//     let mut cache = vec![T::from(1); 2 << log_len];
//     // count number of 0s
//     let mut cnt: isize = -1;
//     for i in (len + 1)..rho.len() {
//         let idx = (i - len - 1) + (1 << log_len);
//         cache[idx] =
//             (T::from(1) - rho[i]) *
//             (T::from(1) - sgm[i]) *
//             (T::from(1) - tau[i]) +
//             rho[i] * sgm[i] * tau[i];
// 
//         if cache[idx] == T::from(0) && cnt == -1 {
//             cnt = i as isize;
//         } else if cache[idx] == T::from(0) && cnt > 0 {
//             return [T::from(0); 4];
//         }
//     }
// 
//     if cnt > 0 {
//         let cnt = cnt as usize;
//         let mut hi_lin1 = factors[cnt - len - 1] *
//             rho[cnt] * sgm[cnt] * tau[cnt];
//         for i in (len + 1)..rho.len() {
//             let idx = (i - len - 1) + (1 << log_len);
//             if i != cnt {
//                 hi_lin1 *= cache[idx];
//             }
//         }
//         lo *= T::from(1) - tau[len];
//         return [T::from(0), hi_lin1 * lo, hi_lin1 * lo, T::from(0)];
//     }
// 
//     // compute mul tree
//     for i in (1..(1 << log_len)).rev() {
//         cache[i] = cache[2 * i] * cache[2 * i + 1];
//     }
// 
//     let add = lo * tau[len] * cache[1];
//     lo *= T::from(1) - tau[len];
//     let mut hi_lin0 = cache[1];
//     let mut hi_lin1 = T::from(0);
// 
//     // compute products with missing element
//     let mut last: usize = (1 << log_len) + (rho.len() - len - 1) - 1;
//     last += 1 - (last & 1);
//     for i in 4..(last + 1) {
//         cache[i] = cache[i] * cache[(i / 2) ^ 1];
//     }
// 
//     for i in (len + 1)..rho.len() {
//         let idx = (i - len - 1) + (1 << log_len);
//         hi_lin1 += factors[i - len - 1] * rho[i] * sgm[i] * tau[i] *
//             cache[idx ^ 1];
//     }
// 
//     hi_lin0 += hi_lin1;
// 
//     [T::from(0), hi_lin0 * lo + add, hi_lin1 * lo + add, T::from(0)]
// }

fn pred_fft_first<T: Field>(len: usize, factors: Vec<T>, tau: &Vec<T>, rho: &Vec<T>, _sgm: &Vec<T>) -> [T; 4] {
    // rest fft layer: rho[hi | 0 | lo: len],
    // first fft layer: rho[0 | lo:len] sigme[1 | lo:len]
    // -> tau[hi | 0 | lo: len] : LINEAR gate w/ 1, factor(hi)
    // -> tau[hi | 1 | lo: len] : LINEAR gate w/ 1, 1 + factor(hi)

    let mut lo = T::from(1);
    for i in 0..len {
        lo *= T::from(1) + rho[i] + tau[i];
    }

    let mut coeff = T::from(0);
    for i in (len + 1)..tau.len() {
        coeff += factors[i - len - 1] * tau[i];
    }

    // tau = 0 && rho = 0 => 1
    // tau = 1 && rho = 0 => 1
    // tau = 0 && rho = 1 => coeff
    // tau = 1 && rho = 1 => coeff + 1

    let res = zer(rho[len]) + one(rho[len]) * (coeff + one(tau[len]));

    [T::from(0), res * lo, T::from(0), T::from(0)]
}

fn pred_fft<T: Field>(len: usize, factors: Vec<T>, tau: &Vec<T>, rho: &Vec<T>, _sgm: &Vec<T>) -> [T; 4] {
    // rho[hi | 0 | lo: len], sgm[hi | 1 | lo: len]
    // -> tau[hi | 0 | lo: len] : LINEAR gate w/ 1, factor(hi)
    // -> tau[hi | 1 | lo: len] : LINEAR gate w/ 1, 1 + factor(hi)

    let mut lo = T::from(1);
    for i in 0..len {
        lo *= T::from(1) + rho[i] + tau[i];
    }
    // rho=0, sgm=1
    // lo *= (T::from(1) - rho[len]) * sgm[len];

    // multiplication tree
    let log_len = ceil_log(tau.len() - len - 1);
    let mut cache = vec![T::from(1); 2 << log_len];
    // count number of 0s
    let mut cnt: isize = -1;
    for i in (len + 1)..tau.len() {
        let idx = (i - len - 1) + (1 << log_len);
        cache[idx] = T::from(1) + rho[i] + tau[i];

        if cache[idx] == T::from(0) && cnt == -1 {
            cnt = i as isize;
        } else if cache[idx] == T::from(0) && cnt > 0 {
            return [T::from(0); 4];
        }
    }

    let mut evn_hi1 = T::from(0);

    if cnt > 0 {
        // XXX not tested
        let cnt = cnt as usize;
        evn_hi1 = factors[cnt - len - 1] * rho[cnt] * tau[cnt];
        for i in (len + 1)..tau.len() {
            let idx = (i - len - 1) + (1 << log_len);
            if i != cnt {
                evn_hi1 *= cache[idx];
            }
        }
        return
            [T::from(0), T::from(0), evn_hi1 * lo, T::from(0)];
    }

    // compute mul tree
    for i in (1..(1 << log_len)).rev() {
        cache[i] = cache[2 * i] * cache[2 * i + 1];
    }

    // compute products with missing element
    let mut last: usize = (1 << log_len) + (tau.len() - len - 1) - 1;
    last += 1 - (last & 1);
    for i in 4..(last + 1) {
        cache[i] = cache[i] * cache[(i / 2) ^ 1];
    }

    for i in (len + 1)..tau.len() {
        let idx = (i - len - 1) + (1 << log_len);
        evn_hi1 += factors[i - len - 1] * rho[i] * tau[i] *
            cache[idx ^ 1];
    }

    [T::from(0), cache[1] * lo * (T::from(1) - rho[len]) + (evn_hi1 + cache[1] * tau[len]) * lo * rho[len], T::from(0), T::from(0)]
}


// fn _pred_fft<T: Field>(len: usize, factors: Vec<T>, tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>) -> [T; 4] {
//     // rho[hi | 0 | lo: len], sgm[hi | 1 | lo: len]
//     // -> tau[hi | 0 | lo: len] : LINEAR gate w/ 1, factor(hi)
//     // -> tau[hi | 1 | lo: len] : LINEAR gate w/ 1, 1 + factor(hi)
// 
//     let mut lo = T::from(1);
//     for i in 0..len {
//         lo *=
//             (T::from(1) - rho[i]) *
//             (T::from(1) - sgm[i]) *
//             (T::from(1) - tau[i]) +
//             rho[i] * sgm[i] * tau[i];
//     }
//     // rho=0, sgm=1
//     lo *= (T::from(1) - rho[len]) * sgm[len];
// 
//     // multiplication tree
//     let log_len = ceil_log(tau.len() - len - 1);
//     let mut cache = vec![T::from(1); 2 << log_len];
//     // count number of 0s
//     let mut cnt: isize = -1;
//     for i in (len + 1)..tau.len() {
//         let idx = (i - len - 1) + (1 << log_len);
//         cache[idx] =
//             (T::from(1) - rho[i]) *
//             (T::from(1) - sgm[i]) *
//             (T::from(1) - tau[i]) +
//             rho[i] * sgm[i] * tau[i];
// 
//         if cache[idx] == T::from(0) && cnt == -1 {
//             cnt = i as isize;
//         } else if cache[idx] == T::from(0) && cnt > 0 {
//             return [T::from(0); 4];
//         }
//     }
// 
//     let mut evn_hi1 = T::from(0);
// 
//     if cnt > 0 {
//         let cnt = cnt as usize;
//         evn_hi1 = factors[cnt - len - 1] *
//             rho[cnt] * sgm[cnt] * tau[cnt];
//         for i in (len + 1)..tau.len() {
//             let idx = (i - len - 1) + (1 << log_len);
//             if i != cnt {
//                 evn_hi1 *= cache[idx];
//             }
//         }
//         return
//             [T::from(0), T::from(0), evn_hi1 * lo, T::from(0)];
//     }
// 
//     // compute mul tree
//     for i in (1..(1 << log_len)).rev() {
//         cache[i] = cache[2 * i] * cache[2 * i + 1];
//     }
// 
//     // compute products with missing element
//     let mut last: usize = (1 << log_len) + (tau.len() - len - 1) - 1;
//     last += 1 - (last & 1);
//     for i in 4..(last + 1) {
//         cache[i] = cache[i] * cache[(i / 2) ^ 1];
//     }
// 
//     for i in (len + 1)..tau.len() {
//         let idx = (i - len - 1) + (1 << log_len);
//         evn_hi1 += factors[i - len - 1] * rho[i] * sgm[i] * tau[i] *
//             cache[idx ^ 1];
//     }
// 
//     [T::from(0), cache[1] * lo,
//     (evn_hi1 + cache[1] * tau[len]) * lo,
//     T::from(0)]
// }

fn pred_selector<T: Field>(pts: Vec<usize>, tau: &Vec<T>, rho: &Vec<T>, _sgm: &Vec<T>) -> [T; 4] {
    // let mut s = T::from(1);
    // for i in 0..sgm.len() {
    //     s *= T::from(1) - sgm[i];
    // }

    let mut tr = T::from(0);
    for t in 0..pts.len() {
        let mut t = t;
        let mut r = pts[t];
        let mut tmp = T::from(1);
        for i in 0..tau.len() {
            if (t & 1) == 1 {
                tmp *= tau[i];
            } else {
                tmp *= T::from(1) - tau[i];
            }
            t >>= 1;
        }
        for i in 0..rho.len() {
            if (r & 1) == 1 {
                tmp *= rho[i];
            } else {
                tmp *= T::from(1) - rho[i];
            }
            r >>= 1;
        }
        tr += tmp;
    }
    [T::from(0), tr, T::from(0), T::from(0)]
}


// pub fn generate_point_table_fft<T: Field>(bkt: &mut BKT<T>, r: &Vec<T>, j: usize) {
//     // gen pts skipping the j-th in r
//     let mut tbl: Vec<T> = Vec::with_capacity((1 << (r.len() - 1)) + 1);
//     tbl.push(T::from(1));
//     let mut mask = 1;
//     for k in 0..r.len() {
//         if k == j {
//             continue;
//         }
//         for i in 0..mask {
//             let tbl_i_mask = tbl[i] * r[k];
//             // tbl[i + mask] = tbl[i] * rho[k];
//             tbl.push(tbl_i_mask);
//             tbl[i] = tbl[i] - tbl_i_mask;
//         }
//         mask <<= 1;
//     }
//     tbl.push(r[j]);
//     bkt.tbl_point.push(tbl);
// }

fn sumcheck_init_ifft<T: Field>(i: usize, log_len: usize, coeff: &Vec<T>, tbl: &mut BKT<T>, back_bit_len: usize) {
    let back_size = 1 << back_bit_len;
    tbl.tbl_coeff = vec![];
    tbl.tbl_coeff.push(vec![T::from(0); back_size]);

    let tbl_prev_p = &tbl.tbl_prev_point[0];
    // let size = 1 << layer.bit_len;
    for j in 0..(1 << log_len) {
        let m = j % (1usize << i);
        let p = (j - m) * 2;
        let idx0 = p + m;
        let idx1 = idx0 + (1usize << i);
        let coeff_idx = p >> (i + 1);

        tbl.tbl_coeff[0][idx0] = tbl_prev_p[idx1] +
            tbl_prev_p[idx0] * (T::from(1) + coeff[coeff_idx]);

        tbl.tbl_coeff[0][idx1] = tbl_prev_p[idx1] +
            tbl_prev_p[idx0] * coeff[coeff_idx];
    }
}

fn sumcheck_init_fft<T: Field>(i: usize, log_len: usize, log_rho: usize, coeff: &Vec<T>, tbl: &mut BKT<T>, back_bit_len: usize) {
    let back_size = 1 << back_bit_len;
    tbl.tbl_coeff = vec![];
    tbl.tbl_coeff.push(vec![T::from(0); back_size]);

    let tbl_prev_p = &tbl.tbl_prev_point[0];
    // let size = 1 << layer.bit_len;
    for k in 0..(1 << log_rho) {
        let offset = k << (log_len + 1);
        for j in 0..(1 << log_len) {
            let m = j % (1usize << i);
            let p = (j - m) * 2;
            let idx0 = p + m + offset;
            let idx1 = idx0 + (1usize << i);
            // let coeff_idx = (k << (log_l - i)) + 2 * (p >> (i + 1));
            // let coeff_idx = (p >> (i + 1)) + (1 << log_rho);
            let coeff_idx = (p >> (i + 1)) ^ (k << (log_len - i));

            tbl.tbl_coeff[0][idx0] = tbl_prev_p[idx1] + tbl_prev_p[idx0];
            tbl.tbl_coeff[0][idx1] = tbl.tbl_coeff[0][idx0] * coeff[coeff_idx] + tbl_prev_p[idx1];

            // tbl.tbl_coeff[0][idx1] =
            //     // tbl_prev_p[idx0] * (layer.utils[coeff_idx] + layer.utils[k]) +
            //     // tbl_prev_p[idx1] * (layer.utils[coeff_idx] + layer.utils[k] + T::from(1));
            //     tbl_prev_p[idx0] * coeff[coeff_idx] +
            //     tbl_prev_p[idx1] * (coeff[coeff_idx] + T::from(1));
        }
    }
}

fn forward_ifft<T: Field>(i: usize, log_len: usize, coeff: &Vec<T>, inputs: &Vec<T>) -> Vec<T> {
    let mut res = vec![T::from(0); 2 << log_len];
    for j in 0..(1 << log_len) {
        let m = j % (1usize << i);
        let p = (j - m) * 2;
        let idx0 = p + m;
        let idx1 = idx0 + (1usize << i);
        let coeff_idx = p >> (i + 1);

        // assert_eq!(coeff[coeff_idx], layer.utils[coeff_idx]);
        res[idx1] = inputs[idx0] + inputs[idx1];
        res[idx0] = inputs[idx0] - res[idx1] * coeff[coeff_idx];
        // res[idx0] = inputs[idx0] * (T::from(1) - coeff[coeff_idx]) + inputs[idx1] * coeff[coeff_idx];
    }
    res
}

fn backward_ifft<T: Field>(i: usize, log_len: usize, coeff: &Vec<T>, inputs: &Vec<T>) -> Vec<T> {
    let mut res = vec![T::from(0); 2 << log_len];
    for j in 0..(1 << log_len) {
        let m = j % (1usize << i);
        let p = (j - m) * 2;
        let idx0 = p + m;
        let idx1 = idx0 + (1usize << i);
        let coeff_idx = p >> (i + 1);

        let tmp = inputs[idx1] * coeff[coeff_idx];
        res[idx0] = inputs[idx0] + tmp;
        res[idx1] = inputs[idx1] + res[idx0];
        // res[idx0] = inputs[idx0] + inputs[idx1] * coeff[coeff_idx];
        // res[idx1] = inputs[idx0] + inputs[idx1] * (coeff[coeff_idx] + T::from(1));
    }
    res
}

fn forward_fft<T: Field>(i: usize, log_len: usize, log_rho: usize, log_l: usize, coeff: &Vec<T>, inputs: &Vec<T>) -> Vec<T> {
    let mut res = vec![T::from(0); 1 << log_l];

    for k in 0..(1 << log_rho) {
        let offset = k << (log_len + 1);
        for j in 0..(1 << log_len) {
            let m = j % (1usize << i);
            let p = (j - m) * 2;
            let idx0 = p + m + offset;
            let idx1 = idx0 + (1usize << i);
            // let coeff_idx = (k << (log_l - i)) + 2 * (p >> (i + 1));
            // let coeff_idx = (p >> (i + 1)) + (1 << log_rho);
            let coeff_idx = (p >> (i + 1)) ^ (k << (log_len - i));

            let tmp = inputs[idx1] * coeff[coeff_idx];
            res[idx0] = inputs[idx0] + tmp;
            res[idx1] = res[idx0] + inputs[idx1];
            // res[idx0] = inputs[idx0] + inputs[idx1] * coeff[coeff_idx];
            // res[idx1] = inputs[idx0] + inputs[idx1] * (coeff[coeff_idx] + T::from(1));
            //
            // res[idx0] = inputs[idx0] + inputs[idx1] * (layer.utils[coeff_idx] + layer.utils[k]);
            // res[idx1] = inputs[idx0] + inputs[idx1] * (layer.utils[coeff_idx] + layer.utils[k] + T::from(1));
            // res[idx0] = inputs[idx0] + inputs[idx1] * layer.utils[coeff_idx];
            // res[idx1] = inputs[idx0] + inputs[idx1] * layer.utils[coeff_idx + 1];
        }
    }
    res
}

fn backward_fft<T: Field>(i: usize, log_len: usize, log_rho: usize, log_l: usize, coeff: &Vec<T>, inputs: &Vec<T>) -> Vec<T> {
    let mut res = vec![T::from(0); 1 << log_l];

    for k in 0..(1 << log_rho) {
        let offset = k << (log_len + 1);
        for j in 0..(1 << log_len) {
            let m = j % (1usize << i);
            let p = (j - m) * 2;
            let idx0 = p + m + offset;
            let idx1 = idx0 + (1usize << i);
            let coeff_idx = (p >> (i + 1)) ^ (k << (log_len - i));

            res[idx1] = inputs[idx0] + inputs[idx1];
            res[idx0] = inputs[idx0] + res[idx1] * coeff[coeff_idx];
            // res[idx0] = inputs[idx0] * (coeff[coeff_idx] + T::from(1)) + inputs[idx1] * coeff[coeff_idx];

            // res[idx0] = inputs[idx0] * (coeff[coeff_idx] + T::from(1)) + inputs[idx1] * coeff[coeff_idx];
            // res[idx1] = inputs[idx0] + inputs[idx1];
        }
    }
    res
}

pub fn append_ifft_fprime<'a, T: Field>(c: &mut Circuit<'a, T>, log_l: usize, coset: &'a Coset<T>) {
    let mut depth = c.depth;
    let log_len = c.layers[depth - 1].bit_len;

    // first ifft layer
    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_len;
        // let size_inv = T::from(1 << log_len).inv();
        let size_inv = coset.ksi;

        layer.num_phase = 1;
        layer.num_table = vec![1];
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _back_bit_len: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_len]];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in (0..(1 << log_len)).step_by(2) {
                    let idx0 = k;
                    let idx1 = k + 1;
                    let u = tbl_prev_p[idx0] * size_inv;
                    let v = tbl_prev_p[idx1] * size_inv;
                    tbl.tbl_coeff[0][swap_bits(idx0, log_len)] = u + v;
                    tbl.tbl_coeff[0][swap_bits(idx1, log_len)] = u - v;
                }
            });
        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_len];
                for k in (0..(1 << log_len)).step_by(2) {
                    let idx0 = k;
                    let idx1 = k + 1;
                    let u = inputs[swap_bits(idx0, log_len)] * size_inv;
                    let v = inputs[swap_bits(idx1, log_len)] * size_inv;
                    res[idx0] = u + v;
                    res[idx1] = u - v;
                }
                res
            });
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                // lo log_len - 1 bits are equal
                let mut eq_lo = T::from(1);
                for j in 1..log_len {
                    eq_lo *= T::xor(rho[log_len - j - 1], tau[j]);
                }
                let u = zer(rho[log_len - 1]);
                let v = one(rho[log_len - 1]);
                let res = zer(tau[0]) * (u + v) + one(tau[0]) * (u - v);
                [T::from(0), res * eq_lo * size_inv, T::from(0), T::from(0)]
            }));
        // layer.gates = vec![Gate::default(); 1 << log_len];
        // layer.values = vec![T::from(0); 1 << log_len];
        // layer.utils = vec![size_inv, - size_inv];
        // for k in (0..(1 << log_len)).step_by(2) {
        //     let idx0 = k;
        //     let idx1 = k + 1;
        //     layer.gates[idx0].gtype = GateType::LINEAR;
        //     layer.gates[idx0].inputs.push(swap_bits(idx0, log_len) as u32);
        //     layer.gates[idx0].inputs.push(swap_bits(idx1, log_len) as u32);
        //     layer.gates[idx0].utils.push(0);
        //     layer.gates[idx0].utils.push(0);

        //     layer.gates[idx1].gtype = GateType::LINEAR;
        //     layer.gates[idx1].inputs.push(swap_bits(idx0, log_len) as u32);
        //     layer.gates[idx1].inputs.push(swap_bits(idx1, log_len) as u32);
        //     layer.gates[idx1].utils.push(0);
        //     layer.gates[idx1].utils.push(1);
        // }
        depth += 1;
    }

    // ifft layers
    for i in 1..log_len {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_len;

        let m = 2 << i;
        let lvl = log_l - i - 1;

        layer.num_phase = 1;
        layer.num_table = vec![1];
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _back_bit_len: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_len]];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in (0..(1 << log_len)).step_by(m) {
                    for j in 0..(m / 2) {
                        let idx0 = k + j;
                        let idx1 = k + j + m / 2;
                        let w = coset.element[lvl][(m - j) % m];
                        let u = tbl_prev_p[idx0];
                        let v = tbl_prev_p[idx1];
                        tbl.tbl_coeff[0][idx0] = u + v;
                        tbl.tbl_coeff[0][idx1] = (u - v) * w;
                    }
                }
            });
        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_len];
                for k in (0..(1 << log_len)).step_by(m) {
                    for j in 0..(m / 2) {
                        let idx0 = k + j;
                        let idx1 = k + j + m / 2;
                        let w = coset.element[lvl][(m - j) % m];
                        let u = inputs[idx0];
                        let v = inputs[idx1] * w;
                        res[idx0] = u + v;
                        res[idx1] = u - v;
                    }
                }
                res
            });
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                // high log_len-(i+1) bits and low i bits are equal
                // low i bits indicate fft coeff
                let mut eq_hi = T::from(1);
                for j in (i + 1)..log_len {
                    eq_hi *= T::xor(rho[j], tau[j]);
                }
                let mut eq_lo = T::from(1);
                for j in 0..i {
                    eq_lo *= T::xor(rho[j], tau[j]);
                }
                let mut coeff = T::from(1);
                for j in 0..i {
                    // 0 if rho[j] != tau[j], 1 if rho[j] == tau[j] == 0
                    // w^{2^j} if rho[j] == tau[j] == 1
                    coeff *= zer(rho[j]) * zer(tau[j]) +
                        one(rho[j]) * one(tau[j]) * coset.element[lvl][(m - (1 << j)) % m];
                }
                let u = zer(rho[i]) * eq_lo;
                let v = one(rho[i]) * coeff;
                let res = zer(tau[i]) * (u + v) + one(tau[i]) * (u - v);
                [T::from(0), res * eq_hi, T::from(0), T::from(0)]
            }));

        // layer.gates = vec![Gate::default(); 1 << log_len];
        // layer.values = vec![T::from(0); 1 << log_len];
        // layer.utils = coset.element[lvl].to_vec();
        // for k in (0..(1 << log_len)).step_by(m) {
        //     for j in 0..(m / 2) {
        //         let idx0 = k + j;
        //         let idx1 = k + j + m / 2;
        //         // let w = coset.element[lvl][(m - swap_bits(j, i)) % m];

        //         layer.gates[idx0].gtype = GateType::LINEAR;
        //         layer.gates[idx0].inputs.push(idx0 as u32);
        //         layer.gates[idx0].inputs.push(idx1 as u32);
        //         layer.gates[idx0].utils.push(0);
        //         layer.gates[idx0].utils.push(((m - j) % m) as u32);

        //         layer.gates[idx1].gtype = GateType::LINEAR;
        //         layer.gates[idx1].inputs.push(idx0 as u32);
        //         layer.gates[idx1].inputs.push(idx1 as u32);
        //         layer.gates[idx1].utils.push(0);
        //         layer.gates[idx1].utils.push(((m / 2 - j) % m) as u32);
        //     }
        // }

        depth += 1;
    }
    c.depth = depth;
}

pub fn append_fft_fprime<'a, T: Field>(c: &mut Circuit<'a, T>, log_l: usize, coset: &'a Coset<T>) {
    let mut depth = c.depth;
    let log_len = c.layers[depth - 1].bit_len;
    let log_rho = log_l - log_len;

    // first fft layer
    {
        let i = log_rho;
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_l;
        layer.num_phase = 1;
        layer.num_table = vec![1];
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);

        let lvl = log_l - i - 1;

        // layer.sumcheck_init[0] = Arc::new(
        //     move |tbl: &mut BKT<T>, _: &Layer<T>, _back_bit_len: usize| {
        //         tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_len]];
        //         let tbl_prev_p = &tbl.tbl_prev_point[0];
        //         for k in (0..(1 << log_l)).step_by(2 << log_rho) {
        //             for j in 0..(1 << log_rho) {
        //                 let idx0 = k + j;
        //                 let idx1 = k + j + (1 << log_rho);
        //                 let w = coset.element[lvl][j];
        //                 let u = tbl_prev_p[idx0];
        //                 let v = tbl_prev_p[idx1];
        //                 tbl.tbl_coeff[0][swap_bits(idx0 >> log_rho, log_len)] += u + v;
        //                 tbl.tbl_coeff[0][swap_bits(idx1 >> log_rho, log_len)] += (u - v) * w;
        //             }
        //         }
        //     });
        // layer.evaluate[0] = Arc::new(
        //     move |_: &Layer<T>, inputs: &Vec<T>| {
        //         let mut res = vec![T::from(0); 1 << log_l];
        //         for k in (0..(1 << log_l)).step_by(2 << log_rho) {
        //             for j in 0..(1 << log_rho) {
        //                 let idx0 = k + j;
        //                 let idx1 = k + j + (1 << log_rho);
        //                 let w = coset.element[lvl][j];
        //                 let u = inputs[swap_bits(idx0 >> log_rho, log_len)];
        //                 let v = inputs[swap_bits(idx1 >> log_rho, log_len)] * w;
        //                 res[idx0] = u + v;
        //                 res[idx1] = u - v;
        //             }
        //         }
        //         res
        //     });
        // layer.predicate.push(Arc::new(
        //     move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
        //         // high log_len-1 bits of tau and rho are equal
        //         // low log_rho bits indicate fft coeff
        //         let mut eq_hi = T::from(1);
        //         for j in (log_rho + 1)..log_l {
        //             eq_hi *= T::xor(rho[log_len - 1 - (j - log_rho)], tau[j]);
        //         }
        //         let mut coeff = T::from(1);
        //         for j in 0..log_rho {
        //             // 1 if tau[j] == 0
        //             // w^{2^j} if tau[j] == 1
        //             coeff *= zer(tau[j]) +
        //                 one(tau[j]) * coset.element[lvl][1 << j];
        //         }
        //         let u = zer(rho[log_len - 1]);
        //         let v = one(rho[log_len - 1]) * coeff;
        //         let res = zer(tau[log_rho]) * (u + v) + one(tau[log_rho]) * (u - v);
        //         [T::from(0), res * eq_hi, T::from(0), T::from(0)]
        //     }));

        layer.gates = vec![Gate::default(); 1 << log_l];
        layer.values = vec![T::from(0); 1 << log_l];
        layer.utils = coset.element[lvl].to_vec();
        for k in (0..(1 << log_l)).step_by(2 << log_rho) {
            for j in 0..(1 << log_rho) {
                let idx0 = k + j;
                let idx1 = k + j + (1 << log_rho);
                layer.gates[idx0].gtype = GateType::LINEAR;
                layer.gates[idx0].inputs.push(swap_bits(idx0 >> log_rho, log_len) as u32);
                layer.gates[idx0].inputs.push(swap_bits(idx1 >> log_rho, log_len) as u32);
                layer.gates[idx0].utils.push(0);
                layer.gates[idx0].utils.push(j as u32);

                layer.gates[idx1].gtype = GateType::LINEAR;
                layer.gates[idx1].inputs.push(swap_bits(idx0 >> log_rho, log_len) as u32);
                layer.gates[idx1].inputs.push(swap_bits(idx1 >> log_rho, log_len) as u32);
                layer.gates[idx1].utils.push(0);
                layer.gates[idx1].utils.push((j + (1 << log_rho)) as u32);
            }
        }

        depth += 1;
    }

    // fft layers
    for i in (log_rho + 1)..log_l {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_l;
        layer.num_phase = 1;
        layer.num_table = vec![1];
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);

        let m = 2 << i;
        let lvl = log_l - i - 1;

        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _back_bit_len: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_l]];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in (0..(1 << log_l)).step_by(m) {
                    for j in 0..(m / 2) {
                        let idx0 = k + j;
                        let idx1 = k + j + m / 2;
                        let w = coset.element[lvl][j];
                        let u = tbl_prev_p[idx0];
                        let v = tbl_prev_p[idx1];
                        tbl.tbl_coeff[0][idx0] = u + v;
                        tbl.tbl_coeff[0][idx1] = (u - v) * w;
                    }
                }
            });
        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_l];
                for k in (0..(1 << log_l)).step_by(m) {
                    for j in 0..(m / 2) {
                        let idx0 = k + j;
                        let idx1 = k + j + m / 2;
                        let w = coset.element[lvl][j];
                        let u = inputs[idx0];
                        let v = inputs[idx1] * w;
                        res[idx0] = u + v;
                        res[idx1] = u - v;
                    }
                }
                res
            });
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                // high log_l-(i+1) bits and low i bits are equal
                // low i bits indicate fft coeff
                let mut eq_hi = T::from(1);
                for j in (i + 1)..log_l {
                    eq_hi *= T::xor(rho[j], tau[j]);
                }
                let mut eq_lo = T::from(1);
                for j in 0..i {
                    eq_lo *= T::xor(rho[j], tau[j]);
                }
                let mut coeff = T::from(1);
                for j in 0..i {
                    // 0 if rho[j] != tau[j], 1 if rho[j] == tau[j] == 0
                    // w^{2^j} if rho[j] == tau[j] == 1
                    coeff *= zer(rho[j]) * zer(tau[j]) +
                        one(rho[j]) * one(tau[j]) * coset.element[lvl][1 << j];
                }
                let u = zer(rho[i]) * eq_lo;
                let v = one(rho[i]) * coeff;
                let res = zer(tau[i]) * (u + v) + one(tau[i]) * (u - v);
                [T::from(0), res * eq_hi, T::from(0), T::from(0)]
            }));

        // layer.gates = vec![Gate::default(); 1 << log_l];
        // layer.values = vec![T::from(0); 1 << log_l];
        // layer.utils = coset.element[lvl].to_vec();
        // for k in (0..(1 << log_l)).step_by(m) {
        //     for j in 0..(m / 2) {
        //         let idx0 = k + j;
        //         let idx1 = k + j + m / 2;
        //         layer.gates[idx0].gtype = GateType::LINEAR;
        //         layer.gates[idx0].inputs.push(idx0 as u32);
        //         layer.gates[idx0].inputs.push(idx1 as u32);
        //         layer.gates[idx0].utils.push(0);
        //         layer.gates[idx0].utils.push(j as u32);

        //         layer.gates[idx1].gtype = GateType::LINEAR;
        //         layer.gates[idx1].inputs.push(idx0 as u32);
        //         layer.gates[idx1].inputs.push(idx1 as u32);
        //         layer.gates[idx1].utils.push(0);
        //         layer.gates[idx1].utils.push((j + m / 2) as u32);
        //     }
        // }

        depth += 1;
    }
    c.depth = depth;
}


pub fn append_ifft<'a, T: Field>(c: &mut Circuit<'a, T>, _log_l: usize, fft_coeff: &'a Vec<Vec<T>>) {
    let mut depth = c.depth;
    let log_len = c.layers[depth - 1].bit_len - 1;

    // ifft
    for i in 0..=log_len {
        c.layers.push(Layer::<T>::default());
        c.layers[depth].bit_len = log_len + 1;
        c.layers[depth].num_phase = 1;
        c.layers[depth].num_table = vec![1];

        // constants
        // c.layers[depth].utils = fft_coeff[i].to_vec();

        // for j in 0..(1 << (log_l - i - 1)) {
        //     c.layers[depth].utils.push(factors0[i][j]);
        //     c.layers[depth].utils.push(T::from(1) - factors0[i][j]);
        // }

        // c.layers[depth].gates = vec![Gate::default(); 2 << log_len];
        // c.layers[depth].values = vec![T::from(0); 2 << log_len];
        // for j in 0..(1 << log_len) {
        //     let m = j % (1usize << i);
        //     let p = (j - m) * 2;
        //     let idx0 = p + m;
        //     let idx1 = idx0 + (1usize << i);
        //     c.layers[depth].gates[idx1].gtype = GateType::ADD;
        //     c.layers[depth].gates[idx1].inputs[0] = idx0 as u32;
        //     c.layers[depth].gates[idx1].inputs[1] = idx1 as u32;

        //     c.layers[depth].gates[idx0].gtype = GateType::LINEAR;
        //     c.layers[depth].gates[idx0].inputs[0] = idx0 as u32;
        //     c.layers[depth].gates[idx0].inputs[1] = idx1 as u32;
        //     c.layers[depth].gates[idx0].utils[0] = (2 * (p >> (i + 1)) + 1) as u32;
        //     c.layers[depth].gates[idx0].utils[1] = (2 * (p >> (i + 1))) as u32;
        //     // c.layers[depth].gates[idx0].utils[0] = T::from(1) + fact;
        //     // c.layers[depth].gates[idx0].utils[1] = fact;
        // }

        c.layers[depth].sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, back_bit_len: usize|
            sumcheck_init_ifft(i, log_len, &fft_coeff[i], tbl, back_bit_len));
        c.layers[depth].sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
        c.layers[depth].evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>|
            forward_ifft(i, log_len, &fft_coeff[i], inputs));
        c.layers[depth].evaluate.push(Arc::new(
                move |_: &Layer<T>, inputs: &Vec<T>|
                backward_ifft(i, log_len, &fft_coeff[i], inputs)));

        // Remark: vanishing is a linear function!
        let mut fft_coeff_base = Vec::with_capacity(log_len - i);
        for j in 0..(log_len - i) {
            let base = fft_coeff[i][1 << j];
            fft_coeff_base.push(base);
        }
        c.layers[depth].predicate.push(Arc::new(
                move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                pred_ifft(i, fft_coeff_base.to_vec(), tau, rho, sgm)));
        depth += 1;
    }
    c.depth = depth;
}

pub fn append_fft<'a, T: Field>(c: &mut Circuit<'a, T>, log_l: usize, fft_coeff: &'a Vec<Vec<T>>) {
    let mut depth = c.depth;
    let log_len = c.layers[depth - 1].bit_len - 1;
    let log_rho = log_l - log_len - 1;
    // first fft layer
    {
        let i = log_len;
        c.layers.push(Layer::<T>::default());
        c.layers[depth].bit_len = log_l;
        c.layers[depth].num_phase = 1;
        c.layers[depth].num_table = vec![1];

        c.layers[depth].sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << (log_len + 1)]];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in 0..(1 << log_rho) {
                    let offset = k << (log_len + 1);
                    for j in 0..(1 << log_len) {
                        let m = j % (1usize << i);
                        let p = (j - m) * 2;
                        let idx0 = p + m;
                        let idx1 = idx0 + (1usize << i);
                        let coeff_idx = (p >> (i + 1)) ^ (k << (log_len - i));

                        let tmp = tbl_prev_p[idx1 + offset] + tbl_prev_p[idx0 + offset];
                        tbl.tbl_coeff[0][idx0] += tmp;
                        tbl.tbl_coeff[0][idx1] += tmp * fft_coeff[i][coeff_idx] + tbl_prev_p[idx1 + offset];
                    }
                }
            });
        c.layers[depth].sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
        c.layers[depth].evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_l];
                for k in 0..(1 << log_rho) {
                    let offset = k << (log_len + 1);
                    for j in 0..(1 << log_len) {
                        let m = j % (1usize << i);
                        let p = (j - m) * 2;
                        let idx0 = p + m;
                        let idx1 = idx0 + (1usize << i);
                        // let coeff_idx = (k << (log_l - i)) + 2 * (p >> (i + 1));
                        // let coeff_idx = (p >> (i + 1)) + (1 << log_rho);
                        let coeff_idx = (p >> (i + 1)) ^ (k << (log_len - i));

                        let tmp = inputs[idx1] * fft_coeff[i][coeff_idx];
                        res[idx0 + offset] = inputs[idx0] + tmp;
                        res[idx1 + offset] = res[idx0 + offset] + inputs[idx1];
                    }
                }
                res
            });
        c.layers[depth].evaluate.push(Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 2 << log_len];
                for j in 0..(1 << log_len) {
                    let m = j % (1usize << i);
                    let p = (j - m) * 2;
                    let idx0 = p + m;
                    let idx1 = idx0 + (1usize << i);
                    let coeff_idx = p >> (i + 1);
                    res[idx1] = inputs[idx0] + inputs[idx1];
                    res[idx0] = inputs[idx0] + res[idx1] * fft_coeff[i][coeff_idx];
                }
                res
            }));

        // c.layers[depth].utils.push(T::from(1));
        // for k in 0..fft_coeff[i].len() {
        //     c.layers[depth].utils.push(fft_coeff[i][k]);
        //     c.layers[depth].utils.push(T::from(1) + fft_coeff[i][k]);
        // }

        // c.layers[depth].gates.resize(1 << log_l, Gate::default());
        // c.layers[depth].values.resize(1 << log_l, T::from(0));
        // for k in 0..(1 << log_rho) {
        //     let offset = k << (log_len + 1);
        //     for j in 0..(1 << log_len) {
        //         let m = j % (1usize << i);
        //         let p = (j - m) * 2;
        //         let idx0 = p + m;
        //         let idx1 = idx0 + (1usize << i);
        //         let coeff_idx = ((p >> (i + 1)) ^ (k << (log_len - i))) as u32;

        //         c.layers[depth].gates[idx0 + offset].gtype = GateType::LINEAR;
        //         c.layers[depth].gates[idx0 + offset].inputs.push(idx0 as u32);
        //         c.layers[depth].gates[idx0 + offset].inputs.push(idx1 as u32);
        //         c.layers[depth].gates[idx0 + offset].utils.push(0);
        //         c.layers[depth].gates[idx0 + offset].utils.push(2 * coeff_idx + 1);
        //         // c.layers[depth].gates[idx0].utils[0] = T::from(1);
        //         // c.layers[depth].gates[idx0].utils[1] = fact;

        //         c.layers[depth].gates[idx1 + offset].gtype = GateType::LINEAR;
        //         c.layers[depth].gates[idx1 + offset].inputs.push(idx0 as u32);
        //         c.layers[depth].gates[idx1 + offset].inputs.push(idx1 as u32);
        //         c.layers[depth].gates[idx1 + offset].utils.push(0);
        //         c.layers[depth].gates[idx1 + offset].utils.push(2 * coeff_idx + 2);
        //         // c.layers[depth].gates[idx1].utils[0] = T::from(1);
        //         // c.layers[depth].gates[idx1].utils[1] = T::from(1) + fact;
        //     }
        // }

        let mut fft_coeff_base = Vec::with_capacity(log_l - i - 1);
        for j in 0..(log_len - i) {
            let base = fft_coeff[i][1 << j];
            fft_coeff_base.push(base);
        }
        for j in 0..log_rho {
            let base = fft_coeff[i][1 << (log_len - i + j)];
            fft_coeff_base.push(base);
        }
        c.layers[depth].predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>| {
                pred_fft_first(i, fft_coeff_base.to_vec(), tau, rho, sgm)
            }));
        depth += 1;
    }

    // fft
    for i in (0..log_len).rev() {
        c.layers.push(Layer::<T>::default());
        c.layers[depth].bit_len = log_l;
        c.layers[depth].num_phase = 1;
        c.layers[depth].num_table = vec![1];

        // constants
        // c.layers[depth].utils = fft_coeff[i].to_vec();
        // c.layers[depth].utils = Vec::with_capacity((1 << log_rho) + (1 << (log_l - i - 1)));
        // for k in 0..(1 << log_rho) {
        //     c.layers[depth].utils.push(factors1[i][k]);
        // }
        // for j in 0..(1 << (log_l - i - 1)) {
        //     c.layers[depth].utils.push(factors0[i][j]);
        // }
        // each gate need 1 as the first coeff,
        // omit for non-gate, function only approach
        // c.layers[depth].utils.push(T::from(1));
        // for k in 0..(1 << log_rho) {
        //     for j in 0..(1 << (log_l - i - 1)) {
        //         let factor = factors0[i][j] + factors1[i][k];
        //         c.layers[depth].utils.push(factor);
        //         // c.layers[depth].utils.push(T::from(1) + factor);
        //     }
        // }

        // let i = log_len - d;
        // c.layers[depth].gates.resize(1 << log_l, Gate::default());
        // c.layers[depth].values.resize(1 << log_l, T::from(0));
        // for k in 0..(1 << log_rho) {
        //     let offset = k << (log_len + 1);
        //     for j in 0..(1 << log_len) {
        //         let m = j % (1usize << i);
        //         let p = (j - m) * 2;
        //         let idx0 = p + m + offset;
        //         let idx1 = idx0 + (1usize << i);

        //         c.layers[depth].gates[idx0].gtype = GateType::LINEAR;
        //         c.layers[depth].gates[idx0].inputs[0] = idx0 as u32;
        //         c.layers[depth].gates[idx0].inputs[1] = idx1 as u32;
        //         c.layers[depth].gates[idx0].utils[0] = 0;
        //         c.layers[depth].gates[idx0].utils[1] = ((k << (log_l - i)) + 2 * (p >> (i + 1)) + 1) as u32;
        //         // c.layers[depth].gates[idx0].utils[0] = T::from(1);
        //         // c.layers[depth].gates[idx0].utils[1] = fact;

        //         c.layers[depth].gates[idx1].gtype = GateType::LINEAR;
        //         c.layers[depth].gates[idx1].inputs[0] = idx0 as u32;
        //         c.layers[depth].gates[idx1].inputs[1] = idx1 as u32;
        //         c.layers[depth].gates[idx1].utils[0] = 0;
        //         c.layers[depth].gates[idx1].utils[1] = ((k << (log_l - i)) + 2 * (p >> (i + 1)) + 2) as u32;
        //         // c.layers[depth].gates[idx1].utils[0] = T::from(1);
        //         // c.layers[depth].gates[idx1].utils[1] = T::from(1) + fact;
        //     }
        // }

        c.layers[depth].sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, back_bit_len: usize|
            sumcheck_init_fft(i, log_len, log_rho, &fft_coeff[i], tbl, back_bit_len));
        c.layers[depth].sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);
        c.layers[depth].evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>|
            forward_fft(i, log_len, log_rho, log_l, &fft_coeff[i], inputs));
        c.layers[depth].evaluate.push(Arc::new(
                move |_: &Layer<T>, inputs: &Vec<T>|
                backward_fft(i, log_len, log_rho, log_l, &fft_coeff[i], inputs)));

        // let mut factors = Vec::with_capacity(log_l - i - 1);
        // for j in 0..(log_len - i) {
        //     let fact = factors0[i][1 << j];
        //     factors.push(fact);
        // }
        // for j in 0..log_rho {
        //     let fact = factors1[i][1 << j];
        //     factors.push(fact);
        // }
        // let mut fft_coeff_base = Vec::with_capacity(log_len - i);
        let mut fft_coeff_base = Vec::with_capacity(log_l - i - 1);
        for j in 0..(log_len - i) {
            let base = fft_coeff[i][1 << j];
            fft_coeff_base.push(base);
        }
        for j in 0..log_rho {
            let base = fft_coeff[i][1 << (log_len - i + j)];
            fft_coeff_base.push(base);
        }
        c.layers[depth].predicate.push(Arc::new(
                move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
                pred_fft(i, fft_coeff_base.to_vec(), tau, rho, sgm)));
        depth += 1;
    }
    c.depth = depth;
}

pub fn append_selector<'a, T: Field>(c: &mut Circuit<'a, T>, query_pts: &Vec<usize>) {
    let kappa = query_pts.len();
    let log_kappa = ceil_log(kappa);
    let depth = c.depth;
    c.layers.push(Layer::<T>::default());
    c.layers[depth].bit_len = log_kappa;
    c.layers[depth].num_phase = 1;
    c.layers[depth].num_table = vec![1];
    c.layers[depth].gates.resize(1 << log_kappa, Gate::default());
    c.layers[depth].values.resize(1 << log_kappa, T::from(0));

    let tmp = query_pts.to_vec();
    c.layers[depth].sumcheck_init[0] = Arc::new(
        move |tbl: &mut BKT<T>, _layer: &Layer<T>, back_bit_len: usize| {
            let back_size = 1 << back_bit_len;
            tbl.tbl_coeff = vec![];
            tbl.tbl_coeff.push(vec![T::from(0); back_size]);
            let tbl_prev_p = &tbl.tbl_prev_point[0];
            for i in 0..kappa {
                tbl.tbl_coeff[0][tmp[i]] += tbl_prev_p[i];
            }
        });
    c.layers[depth].sumcheck_eval[0] = Arc::new(BKT::<T>::sumcheck_phase1_eval_no_const);

    for i in 0..kappa {
        c.layers[depth].gates[i].gtype = GateType::RELAY;
        c.layers[depth].gates[i].inputs.push(query_pts[i] as u32);
    }

    let tmp = query_pts.to_vec();
    c.layers[depth].predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, sgm: &Vec<T>|
            pred_selector(tmp.to_vec(), tau, rho, sgm)));
    c.depth += 1;
}

pub fn append_poly_eval_fprime<'a, T: Field>(c: &mut Circuit<'a, T>, power_series: &Vec<Vec<T>>) {
    let kappa = power_series.len();
    let log_kappa = ceil_log(kappa);
    let log_len = power_series[0].len();
    let mut depth = c.depth;

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        layer.bit_len = log_kappa;
        layer.num_phase = 1;
        layer.num_table = vec![1];
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);

        let p = power_series.clone();
        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_len]];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in 0..kappa {
                    let mut power = T::from(1);
                    for i in 0..(1 << log_len) {
                        tbl.tbl_coeff[0][i] += power * tbl_prev_p[k];
                        power *= p[k][0];
                    }
                }
            });
        let p = power_series.clone();
        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_kappa];
                for k in 0..kappa {
                    let mut val = inputs[(1 << log_len) - 1];
                    for i in (0..((1 << log_len) - 1)).rev() {
                        val = val * p[k][0] + inputs[i];
                    }
                    res[k] = val;
                }
                res
            });
        let p = power_series.clone();
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                let mut res = T::from(0);
                for k in 0..kappa {
                    let mut out = expand_binary_eq(k, tau);
                    for i in 0..log_len {
                        out *= p[k][i] * one(rho[i]) + zer(rho[i]);
                    }
                    res += out;
                }
                let zero = T::from(0);
                [zero, res, zero, zero]
            }));

        // layer.gates.resize(1 << layer.bit_len, Gate::default());
        // layer.values.resize(1 << layer.bit_len, T::from(0));
        // for k in 0..kappa {
        //     let mut power = T::from(1);
        //     for i in 0..(1 << log_len) {
        //         layer.utils.push(power);
        //         power *= power_series[k][0];

        //         let idx = (k << log_len) + i;
        //         layer.gates[idx].gtype = GateType::SCALAR;
        //         layer.gates[idx].inputs.push(i as u32);
        //         layer.gates[idx].utils.push(idx as u32);
        //     }
        // }
        depth += 1;
    }

    // {
    //     c.layers.push(Layer::<T>::default());
    //     let layer = &mut c.layers[depth];
    //     layer.bit_len = log_kappa;
    //     layer.num_phase = 1;
    //     layer.num_table = vec![1];
    //     layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);

    //     layer.gates.resize(1 << layer.bit_len, Gate::default());
    //     layer.values.resize(1 << layer.bit_len, T::from(0));
    //     for k in 0..kappa {
    //         layer.gates[k].gtype = GateType::SUM;
    //         for i in 0..(1 << log_len) {
    //             layer.gates[k].inputs.push(((k << log_len) + i) as u32);
    //         }
    //     }
    //     depth += 1;
    // }
    c.depth = depth;
}

pub fn append_poly_eval<'a, T: Field>(c: &mut Circuit<'a, T>, power_series: &Vec<Vec<T>>) {
    let kappa = power_series.len();
    let log_kappa = ceil_log(kappa);
    // let log_len = power_series[0].len();
    // assert_eq!(c.layers[c.depth - 1].bit_len, log_len);
    let log_len = c.layers[c.depth - 1].bit_len;
    let mut depth = c.depth;

    {
        c.layers.push(Layer::<T>::default());
        let layer = &mut c.layers[depth];
        // layer.bit_len = log_kappa;
        layer.bit_len = log_kappa;
        layer.num_phase = 1;
        layer.num_table = vec![1];
        layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);

        let p = power_series.clone();
        layer.sumcheck_init[0] = Arc::new(
            move |tbl: &mut BKT<T>, _: &Layer<T>, _: usize| {
                tbl.tbl_coeff = vec![vec![T::from(0); 1 << log_len]];
                let tbl_prev_p = &tbl.tbl_prev_point[0];
                for k in 0..kappa {
                    for i in 0..(1 << log_len) {
                        let mut power = T::from(1);
                        let mut x = i;
                        for j in 0..log_len {
                            if (x & 1) == 1 {
                                power *= p[k][j];
                            }
                            x >>= 1;
                        }
                        tbl.tbl_coeff[0][i] += power * tbl_prev_p[k];
                    }
                }
            });
        let p = power_series.clone();
        layer.evaluate[0] = Arc::new(
            move |_: &Layer<T>, inputs: &Vec<T>| {
                let mut res = vec![T::from(0); 1 << log_kappa];
                for k in 0..kappa {
                    let mut vals = inputs.to_vec();
                    for i in 0..(1 << log_len) {
                        let mut x = i;
                        for j in 0..log_len {
                            if (x & 1) == 1 {
                                vals[i] *= p[k][j];
                            }
                            x >>= 1;
                        }
                    }
                    for v in vals.iter() {
                        res[k] += *v;
                    }
                }
                res
            });
        let p = power_series.clone();
        layer.predicate.push(Arc::new(
            move |tau: &Vec<T>, rho: &Vec<T>, _: &Vec<T>| {
                let mut res = T::from(0);
                for k in 0..kappa {
                    let mut out = expand_binary_eq(k, tau);
                    for i in 0..log_len {
                        out *= p[k][i] * one(rho[i]) + zer(rho[i]);
                    }
                    res += out;
                }
                let zero = T::from(0);
                [zero, res, zero, zero]
            }));

        // layer.gates.resize(1 << layer.bit_len, Gate::default());
        // layer.values.resize(1 << layer.bit_len, T::from(0));
        // for k in 0..kappa {
        //     for i in 0..(1 << log_len) {
        //         let mut power = T::from(1);
        //         let mut x = i;
        //         for j in 0..log_len {
        //             if (x & 1) == 1 {
        //                 power *= power_series[k][j];
        //             }
        //             x >>= 1;
        //         }
        //         layer.utils.push(power);

        //         let idx = (k << log_len) + i;
        //         layer.gates[idx].gtype = GateType::SCALAR;
        //         layer.gates[idx].inputs.push(i as u32);
        //         layer.gates[idx].utils.push(idx as u32);
        //     }
        // }
        depth += 1;
    }

    // {
    //     c.layers.push(Layer::<T>::default());
    //     let layer = &mut c.layers[depth];
    //     layer.bit_len = log_kappa;
    //     layer.num_phase = 1;
    //     layer.num_table = vec![1];
    //     layer.sumcheck_eval[0] = Arc::new(BKT::sumcheck_phase1_eval_no_const);

    //     layer.gates.resize(1 << layer.bit_len, Gate::default());
    //     layer.values.resize(1 << layer.bit_len, T::from(0));
    //     for k in 0..kappa {
    //         layer.gates[k].gtype = GateType::SUM;
    //         for i in 0..(1 << log_len) {
    //             layer.gates[k].inputs.push(((k << log_len) + i) as u32);
    //         }
    //     }
    //     depth += 1;
    // }
    c.depth = depth;
}
