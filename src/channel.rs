extern crate mio;

use std::ops::{Deref, DerefMut};
use std::time::Instant;

use crate::field::*;
use crate::tape::*;
use crate::merkle::*;
use crate::hash::*;
use crate::poly::*;
use crate::eventfd::*;
use crate::statistic::*;
use crate::util::*;

use mio::{Events, Poll, Interest, Token};
use mio::unix::SourceFd;
use queues::*;


#[derive(Clone)]
pub enum Tx<T: Field> {
    Syn,
    Random(T),
    OLayer(T),
    Input(Vec<T>),
    InputZK(Vec<Vec<T>>),
    Output(Vec<T>),
    MultiRandom(Vec<T>),
    Sumcheck(Vec<T>),
    Layer(Vec<T>),
    FinalLayer(Vec<T>),
    FinalLayerPred([T; 4]),
    FinalLayerOpen(SHA256),
    LayerPred([T; 4]),
    LayerScalar(T),
    MPCitH0,
    FLPCP,
    FLPCPRandom(Vec<T>),
    FLPCPCommit(Vec<T>),
    FLPCPOpen((SHA256, T, T, T)),
    MPCitH,
    VOLEOpen((Vec<T>, Vec<T>, SHA256)),
    Party(Vec<usize>),
    SeedTree(SHA256),
    Seed((SHA256, SeedTreeOpen)),
    VPDOpen0((T, T, SHA256)),
    VPDOpen1(SHA256),
    Hash(SHA256),
    LDTCommitHash(SHA256),
    LDTCommitValue(Vec<T>),
    LDTQuery(Vec<MTO<T>>),
    LDTQueryPoint(Vec<usize>),
}


#[derive(Clone, PartialEq)]
pub enum ChannelType {
    Plain,
    MPCitH,
    VOLEitH
}


impl<T: Field> std::fmt::Display for Tx<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Output(_) => write!(f, "output"),
            Self::FinalLayerPred(_) => write!(f, "final layer pred"),
            Self::LayerPred(_) => write!(f, "layer pred"),
            Self::LayerScalar(_) => write!(f, "layer scalar"),
            Self::VPDOpen0(_) => write!(f, "vpd open0"),
            Self::VPDOpen1(_) => write!(f, "vpd open1"),
            Self::LDTCommitHash(_) => write!(f, "ldt commit hash"),
            Self::LDTCommitValue(_) => write!(f, "ldt commit value"),
            Self::LDTQuery(_) => write!(f, "ldt query"),
            Self::LDTQueryPoint(_) => write!(f, "ldt query point"),
            Self::Input(_) => write!(f, "input"),
            Self::InputZK(_) => write!(f, "input_zk"),
            Self::Random(_) => write!(f, "random"),
            Self::MultiRandom(_) => write!(f, "multirandom"),
            Self::Sumcheck(_) => write!(f, "sumcheck"),
            Self::Layer(_) => write!(f, "layer"),
            Self::OLayer(_) => write!(f, "olayer"),
            Self::FinalLayer(_) => write!(f, "finallayer"),
            Self::FinalLayerOpen(_) => write!(f, "finallayeropen"),
            Self::Syn => write!(f, "syn"),
            Self::FLPCPRandom(_) => write!(f, "flpcprandom"),
            Self::FLPCPOpen(_) => write!(f, "flpcpopen"),
            Self::FLPCPCommit(_) => write!(f, "flpcpcommit"),
            Self::Seed(_) => write!(f, "seed"),
            Self::SeedTree(_) => write!(f, "seedtree"),
            Self::Party(_) => write!(f, "party"),
            Self::MPCitH0 => write!(f, "mpcith0"),
            Self::MPCitH => write!(f, "mpcith"),
            Self::VOLEOpen(_) => write!(f, "voleopen"),
            Self::FLPCP => write!(f, "flpcp"),
            Self::Hash(_) => write!(f, "hash"),
            // _ => write!(f, "unknown"),
        }
    }
}

pub trait Channel {
    type Output: Field;
    const CHANNEL_TYPE: ChannelType;
    fn get_type(&mut self) -> ChannelType;
    fn force_plain(&mut self);
    fn send(&mut self, id: usize, tx: Tx::<Self::Output>);
    fn recv(&mut self, id: usize) -> Tx::<Self::Output>;
}


#[macro_export]
macro_rules! recv {
    ($ch:expr, @$msg_t:ident, $id:expr) => {
        match $ch.recv($id) {
            Tx::$msg_t => (),
            tx => {
                println!("{}", tx);
                unreachable!()
            },
        }
    };
    ($ch:expr, $msg_t:ident, $T:ty, $id:expr) => {
        match $ch.recv($id) {
            Tx::<$T>::$msg_t(payload) => payload,
            tx => {
                println!("{}", tx);
                unreachable!()
            },
        }
    };
    ($ch:expr, $msg_t:ident, $id:expr) => {
        match $ch.recv($id) {
            Tx::<T>::$msg_t(payload) => payload,
            tx => {
                println!("{}", tx);
                unreachable!()
            },
        }
    };
}

pub use recv;


pub struct LocalChannel<T: Field> {
    buffer: [Queue<Tx<T>>; 2],
    event: [i32; 2],
}

impl<T: Field> LocalChannel<T> {
    pub fn new() -> Self {
        Self {
            buffer: [queue![], queue![]],
            event: [unsafe { eventfd(0, 0) }, unsafe { eventfd(0, 0) }],
        }
    }
    pub fn close(&self) {
        unsafe {
            close(self.event[0]);
            close(self.event[1]);
        }
    }
}

impl<T: Field> Channel for LocalChannel<T> {
    type Output = T;

    const CHANNEL_TYPE: ChannelType = ChannelType::Plain;

    fn get_type(&mut self) -> ChannelType {
        ChannelType::Plain
    }

    fn force_plain(&mut self) {}

    fn send(&mut self, id: usize, tx: Tx::<Self::Output>) {
        let start = Instant::now();
        // println!("{} send {}", id, tx);
        self.buffer[id].add(tx).unwrap();
        unsafe { eventfd_write(self.event[id], 1) };

        let elapsed = start.elapsed();
        let in_ms = elapsed.as_micros() as usize;
        if id == 0 {
            channel0_inc(in_ms);
        } else {
            channel1_inc(in_ms);
        }
    }

    fn recv(&mut self, id: usize) -> Tx::<Self::Output> {
        // println!("{} recving...", 1 - id);
        const TOKEN: Token = Token(0);
        let mut poll = Poll::new().unwrap();
        poll.registry().register(&mut SourceFd(&self.event[id]), TOKEN, Interest::READABLE).unwrap();
        let mut events = Events::with_capacity(1);

        poll.poll(&mut events, None).unwrap();
        for event in &events {
            match event.token() {
                TOKEN => {
                    let mut tmp: u64 = 0;
                    let tmp_ptr: *mut u64 = &mut tmp;
                    unsafe { eventfd_read(self.event[id], tmp_ptr) };
                    return self.buffer[id].remove().unwrap();
                }
                _ => unreachable!(),
            }
        }
        unreachable!()
    }
}

pub struct PLocalChannel<T: Field>(*mut LocalChannel<T>);
impl<T: Field> PLocalChannel<T> {
    pub fn new(t: *mut LocalChannel<T>) -> Self {
        PLocalChannel::<T>(t)
    }
}
unsafe impl<T: Field> Send for PLocalChannel<T> {}
unsafe impl<T: Field> Sync for PLocalChannel<T> {}
impl<T: Field> Deref for PLocalChannel<T> {
    type Target = LocalChannel<T>;

    fn deref(&self) -> &Self::Target {
        unsafe { &*self.0 }
    }
}
impl<T: Field> DerefMut for PLocalChannel<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        unsafe { &mut *self.0 }
    }
}


pub struct LocalZKChannel<T: Field> {
    plain: LocalChannel<T>,
    seedtree: SeedTree,
    tapes: Vec<Tape>,
    n_party: usize,
    _fr: Vec<T>,
    _Vr: Vec<T>,
    _Vs: Vec<T>,
    _fr_shares: Vec<Vec<T>>,
    _Vr_shares: Vec<Vec<T>>,
    _Vs_shares: Vec<Vec<T>>,
    _p_shares: Vec<Vec<T>>,
    _fr_poly_shares: Vec<Vec<T>>,
    _Vr_poly_shares: Vec<Vec<T>>,
    _Vs_poly_shares: Vec<Vec<T>>,
    _scalar: T,
    _pred: [T; 4],
    _r: T,
    _layer_id: usize,
    flag: bool,
    local_msg: Tx::<T>,
    plain_flag: bool,
}

impl<T: Field> LocalZKChannel<T> {
    pub fn new(mseed: &Seed, tape_len: usize, n_party: usize) -> Self {
        let start = Instant::now();

        let st = SeedTree::new(mseed, n_party, 1);
        let mut tapes: Vec<Tape> = vec![];
        for p in 0..n_party {
            let seed = st.get_seed(0, p);
            tapes.push(Tape::new(seed, tape_len * T::SIZE));
        }

        let elapsed = start.elapsed();
        let in_ms = elapsed.as_micros() as usize;
        mpc_prv_inc(in_ms);
        Self {
            plain: LocalChannel::<T>::new(),
            seedtree: st,
            tapes: tapes,
            n_party: n_party,
            _fr: vec![],
            _Vr: vec![],
            _Vs: vec![],
            _fr_shares: vec![vec![]; n_party],
            _Vr_shares: vec![vec![]; n_party],
            _Vs_shares: vec![vec![]; n_party],
            _p_shares: vec![vec![]; n_party],
            _fr_poly_shares: vec![],
            _Vr_poly_shares: vec![],
            _Vs_poly_shares: vec![],
            _r: T::from(0),
            _scalar: T::from(0),
            _pred: [T::from(0); 4],
            _layer_id: 0,
            flag: false,
            local_msg: Tx::<T>::Syn,
            plain_flag: false,
        }
    }
    pub fn close(&self) {
        self.plain.close()
    }
}

impl<T: Field> Channel for LocalZKChannel<T> {
    type Output = T;

    const CHANNEL_TYPE: ChannelType = ChannelType::MPCitH;

    fn get_type(&mut self) -> ChannelType {
        ChannelType::MPCitH
    }

    fn force_plain(&mut self) {
        self.plain_flag = true;
    }

    fn send(&mut self, id: usize, tx: Tx::<Self::Output>) {
        if self.plain_flag {
            return self.plain.send(id, tx);
        }

        if id == 0 {
            channel0_rst();
        }
        let start0 = Instant::now();
        let flp0_ms = flp_prv_get();
        let mut recv_wait_ms = 0;
        match tx {
            Tx::<T>::OLayer(r) => {
                let mut fr_share = vec![T::from(0); self.n_party];
                fr_share[0] = r;
                for p in 0..self.n_party {
                    self._fr_shares[p].push(fr_share[p]);
                }
                self._fr.push(r);
                self.flag = true;
                self.local_msg = Tx::Syn;
            },
            Tx::<T>::LayerPred(r) => {
                self._pred = r;
                if self._layer_id > 0 {
                    let lid = self._layer_id - 1;
                    for p in 0..self.n_party {
                        self._fr_shares[p][lid] -= self._pred[1] * self._Vr_shares[p][lid] + self._pred[2] * self._Vs_shares[p][lid];
                    }
                    self._fr_shares[0][lid] -= self._pred[3];
                    self._fr[lid] -= self._pred[1] * self._Vr[lid] + self._pred[2] * self._Vs[lid] + self._pred[3];
                }
                self.flag = true;
                self.local_msg = Tx::Syn;
            },
            Tx::<T>::LayerScalar(r) => {
                self._scalar = r;
                let scalar0 = T::from(1) - r;
                let scalar1 = r;
                let mut fr_share = vec![T::from(0); self.n_party];
                let lid = self._layer_id - 1;
                for p in 0..self.n_party {
                    fr_share[p] =
                        scalar0 * self._Vr_shares[p][lid] +
                        scalar1 * self._Vs_shares[p][lid];
                    self._Vr_shares[p][lid] *= self._pred[0];
                }
                for p in 0..self.n_party {
                    self._fr_shares[p].push(fr_share[p]);
                }
                self._fr.push(
                    scalar0 * self._Vr[lid] +
                    scalar1 * self._Vs[lid]);
                self._Vr[lid] *= self._pred[0];
                self.flag = true;
                self.local_msg = Tx::Syn;
            },
            Tx::<T>::Input(r) => {
                let mut msg: Vec<T> = vec![];
                for w in r.iter() {
                    let mut icr = T::from(0);
                    for t in self.tapes.iter_mut() {
                        icr += t.get::<T>();
                    }
                    icr = *w - icr;
                    msg.push(icr);
                }
                self.plain.send(id, Tx::<T>::Input(msg));
            },
            Tx::<T>::Sumcheck(r) => {
                let mut coef_share = vec![vec![]; r.len()];

                let mut msg: Vec<T> = vec![];
                for i in 0..r.len() {
                    let w = r[i];
                    let mut icr = T::from(0);
                    for t in self.tapes.iter_mut() {
                        let share = t.get::<T>();
                        coef_share[i].push(share);
                        icr += share;
                    }
                    icr = w - icr;
                    msg.push(icr);
                }
                for i in 0..r.len() {
                    coef_share[i][0] += msg[i];
                }

                self.plain.send(id, Tx::<T>::Sumcheck(msg));
        let start1 = Instant::now();
                let _r = recv!(self.plain, Random, 1 - id);
        let elapsed1 = start1.elapsed();
        recv_wait_ms += elapsed1.as_micros() as usize;

                if r.len() == 2 {
                    for p in 0..self.n_party {
                        let coef2 = self._fr_shares[p][self._layer_id]
                            - coef_share[0][p] - coef_share[0][p]
                            - coef_share[1][p];
                        self._fr_shares[p][self._layer_id] =
                            (coef_share[1][p]
                             * _r + coef2)
                            * _r + coef_share[0][p];
                    }
                    let d = self._layer_id;
                    let tmp = self._fr[d] - r[0] - r[0] - r[1];
                    self._fr[d] = _r * ((_r * r[1]) + tmp) + r[0];
                } else {
                    for p in 0..self.n_party {
                        let coef2 = self._fr_shares[p][self._layer_id]
                            - coef_share[0][p] - coef_share[0][p]
                            - coef_share[1][p] - coef_share[2][p];
                        self._fr_shares[p][self._layer_id] =
                            ((coef_share[2][p]
                              * _r + coef_share[1][p])
                             * _r + coef2)
                            * _r + coef_share[0][p];
                    }
                    let d = self._layer_id;
                    let tmp = self._fr[d] - r[0] - r[0] - r[1] - r[2];
                    self._fr[d] = ((r[2] * _r + r[1]) * _r + tmp) * _r + r[0];
                }
                self.flag = true;
                self.local_msg = Tx::<T>::Random(_r);
            },
            Tx::<T>::Layer(r) => {
                let mut Vr_share = vec![T::from(0); self.n_party];
                let mut Vs_share = vec![T::from(0); self.n_party];
                let mut msg: Vec<T> = vec![];

                {
                    let mut icr = T::from(0);
                    for p in 0..self.n_party {
                        let share = self.tapes[p].get::<T>();
                        Vr_share[p] = share;
                        icr += share
                    }
                    icr = r[0] - icr;
                    msg.push(icr);
                    Vr_share[0] += icr;
                }
                for p in 0..self.n_party {
                    self._Vr_shares[p].push(Vr_share[p]);
                }
                self._Vr.push(r[0]);

                if r.len() > 1 {
                    {
                        let mut icr = T::from(0);
                        for p in 0..self.n_party {
                            let share = self.tapes[p].get::<T>();
                            Vs_share[p] = share;
                            icr += share
                        }
                        icr = r[1] - icr;
                        msg.push(icr);
                        Vs_share[0] += icr;
                    }
                    for p in 0..self.n_party {
                        self._Vs_shares[p].push(Vs_share[p]);
                    }
                    self._Vs.push(r[1]);
                } else {
                    for p in 0..self.n_party {
                        self._Vs_shares[p].push(*self._Vr_shares[p].last().unwrap());
                    }
                    self._Vs.push(r[0]);
                }
                self._layer_id += 1;
                self.plain.send(id, Tx::<T>::Layer(msg));
            },
            Tx::<T>::FinalLayer(r) => {
                self._scalar = r[0];
                if r.len() > 1 {
                    self._r = r[1];
                } else {
                    self._r = r[0];
                }
                // let prod = r[0] * r[1];
                // self._r = prod;
                self.plain.send(id, Tx::<T>::FinalLayer(r));
            },

            Tx::<T>::FinalLayerPred(r) => {
                // self._scalar stores last layer v_a, v_b
                let tmp = r[0] * self._scalar * self._r +
                    r[1] * self._scalar + r[2] * self._r + r[3];
                self._r = tmp;
                self.flag = true;
                self.local_msg = Tx::Syn;
            }

            Tx::<T>::MPCitH0 => {
                // TODO XXX merge hash with flpcp hash
                // TODO XXX merge hash of non secmul layer opens
                let mut h = SHA256::from(0);
                let mut output_shares: Vec<T> = vec![];
                // assert_eq!(self._r, *self._fr.last().unwrap());
                // let layer_id = self._fr_shares[0].len() - 1;
                self._fr.truncate(self._layer_id);
                for p in 0..self.n_party {
                    output_shares.push(self._fr_shares[p][self._layer_id]);
                    self._fr_shares[p].truncate(self._layer_id);
                }
                output_shares[0] -= self._r;
                h.commit(&output_shares);

                self.plain.send(id, Tx::<T>::FinalLayerOpen(h));
                // self.flag = true;
                // self.local_msg = Tx::Syn;
            },
            Tx::<T>::FLPCP => {
            let start = Instant::now();

                let d = self._layer_id;
                let x_values_to_m = first_n_field_elements::<T>(d + 1);
                let precompute_to_m = precompute_lagrange_polynomials::<T>(&x_values_to_m);
                let x_values_to_2m = first_n_field_elements::<T>(2 * d + 1);
                let precompute_to_2m = precompute_lagrange_polynomials::<T>(&x_values_to_2m);

            let elapsed = start.elapsed();
            let in_ms = elapsed.as_micros() as usize;
            flp_prv_inc(in_ms);

                self._fr.push(T::from(0));
                self._Vr.push(T::from(0));
                self._Vs.push(T::from(0));

                let mut poly_icr = vec![T::from(0); d + 1];

                for p in 0..self.n_party {
                    // let mut a = vec![T::from(0); d + 1];
                    // let mut b = vec![T::from(0); d + 1];
                    // let mut c = vec![T::from(0); d + 1];
                    // for i in 0..d {
                    //     a[i] = self._Vr_shares[p][i];
                    //     b[i] = self._Vs_shares[p][i];
                    //     c[i] = self._fr_shares[p][i];
                    // }
                    let a = self.tapes[p].get::<T>();
                    let b = self.tapes[p].get::<T>();
                    let c = self.tapes[p].get::<T>();
                    self._Vr_shares[p].push(a);
                    self._Vs_shares[p].push(b);
                    self._fr_shares[p].push(c);

                    self._Vr[d] += a;
                    self._Vs[d] += b;
                    self._fr[d] += c;
                    self._Vr_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &self._Vr_shares[p]));
                    self._Vs_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &self._Vs_shares[p]));
                    self._fr_poly_shares.push(interpolate_with_precomputation(&precompute_to_m, &self._fr_shares[p]));
                    for _ in 0..d {
                        self._p_shares[p].push(T::from(0));
                    }
                    for i in 0..(d + 1) {
                        let share = self.tapes[p].get::<T>();
                        self._p_shares[p].push(share);
                        poly_icr[i] += share;
                    }
                }
                let start = Instant::now();

                let Vr_poly = interpolate_with_precomputation(&precompute_to_m, &self._Vr);
                let Vs_poly = interpolate_with_precomputation(&precompute_to_m, &self._Vs);
                let fr_poly = interpolate_with_precomputation(&precompute_to_m, &self._fr);
                let p_poly = poly_mul_sub(&Vr_poly, &Vs_poly, &fr_poly);
                for k in 0..(d + 1) {
                    let x = x_values_to_2m[k + d];
                    poly_icr[k] = poly_eval(&p_poly, x) - poly_icr[k];
                    self._p_shares[0][d + k] += poly_icr[k];
                }

                let elapsed = start.elapsed();
                let in_ms = elapsed.as_micros() as usize;
                flp_prv_inc(in_ms);

                self.plain.send(id, Tx::<T>::FLPCPCommit(poly_icr));

        let start1 = Instant::now();
                let _r = recv!(self.plain, Random, 1 - id);
        let elapsed1 = start1.elapsed();
        recv_wait_ms += elapsed1.as_micros() as usize;

            // },
            // Tx::<T>::FLPCP1 => {
            //     let d = self._layer_id;
            //     let x_values_to_m = first_n_field_elements::<T>(d + 1);
            //     let precompute_to_m = precompute_lagrange_polynomials::<T>(&x_values_to_m);
            //     let x_values_to_2m = first_n_field_elements::<T>(2 * d + 1);
            //     let precompute_to_2m = precompute_lagrange_polynomials::<T>(&x_values_to_2m);

                // let start = Instant::now();

                let mut lagrange_poly_m_eval_r = vec![T::from(0); d + 1];
                let mut lagrange_poly_2m_eval_r = vec![T::from(0); 2 * d + 1];

                let mut a_shares = vec![];
                let mut b_shares = vec![];
                let mut c_shares = vec![];
                let mut d_shares = vec![];
                let mut a = T::from(0);
                let mut b = T::from(0);
                let mut c = T::from(0);

                for k in 0..(d + 1) {
                    lagrange_poly_m_eval_r[k] = poly_eval(&precompute_to_m[k], _r);
                }
                for k in 0..(2 * d + 1) {
                    lagrange_poly_2m_eval_r[k] = poly_eval(&precompute_to_2m[k], _r);
                }

                // let elapsed = start.elapsed();
                // let in_ms = elapsed.as_micros() as usize;
                // flp_prv_inc(in_ms);

                for p in 0..self.n_party {
                    a_shares.push(dot_product(&lagrange_poly_m_eval_r, &self._Vr_shares[p]));
                    b_shares.push(dot_product(&lagrange_poly_m_eval_r, &self._Vs_shares[p]));
                    c_shares.push(dot_product(&lagrange_poly_m_eval_r, &self._fr_shares[p]));
                    d_shares.push(dot_product(&lagrange_poly_2m_eval_r, &self._p_shares[p]));
                    a += a_shares[p];
                    b += b_shares[p];
                    c += c_shares[p];
                }
                // assert_eq!(a*b-c, d_shares[0]+d_shares[1]+d_shares[2]+d_shares[3]);

                let mut h = SHA256::from(0);
                a_shares.append(&mut b_shares);
                a_shares.append(&mut c_shares);
                a_shares.append(&mut d_shares);
                h.commit(&a_shares);

                self.plain.send(id, Tx::<T>::FLPCPOpen((h, a, b, c)));
                // self.flag = true;
                // self.local_msg = Tx::Syn;
            },
            Tx::<T>::MPCitH => {
                self.plain.send(id, Tx::<T>::MPCitH);

                // TODO seedtree commit root should be sent in the very first round
                let root = self.seedtree.commit();
                let mut sto = SeedTreeOpen::new(&self.seedtree);
        let start1 = Instant::now();
                sto.missings = recv!(self.plain, Party, 1 - id);
        let elapsed1 = start1.elapsed();
        recv_wait_ms += elapsed1.as_micros() as usize;

                self.seedtree.open(&mut sto);

                self.plain.send(id, Tx::<T>::Seed((root, sto.clone())));
            },
            // Tx::<T>::Random(_) => self.plain.send(id, tx),
            // Tx::<T>::MultiRandom(_) => self.plain.send(id, tx),
            // Tx::<T>::Party(_) => self.plain.send(id, tx),
            // Tx::<T>::VPDOpen0(_) => self.plain.send(id, tx),
            // Tx::<T>::VPDOpen1(_) => self.plain.send(id, tx),
            // Tx::<T>::Syn => self.plain.send(id, tx),
            _ => {
                self.plain.send(id, tx)
                // println!("{}", tx);
                // unreachable!()
            },
        }
        if id == 0 {
            let elapsed = start0.elapsed();
            let in_ms = elapsed.as_micros() as usize;
            let wait_ms = channel0_get();
            let flp_ms = flp_prv_get() - flp0_ms;
            channel0_rst();
            mpc_prv_inc(in_ms - wait_ms - flp_ms - recv_wait_ms);
        }
    }

    fn recv(&mut self, id: usize) -> Tx::<Self::Output> {
        let start = Instant::now();

        if self.flag {
            self.flag = false;
            return self.local_msg.clone();
        }
        let res = self.plain.recv(id);

        let elapsed = start.elapsed();
        let in_ms = elapsed.as_micros() as usize;
        if id == 1 {
            channel0_inc(in_ms);
        } else {
            channel1_inc(in_ms);
        }
        res
        // match res {
        //     Tx::<T>::Random(r) => {
        //         self._r = r;
        //         res
        //     },
        //     _ => res,
        // }
    }
}

pub struct PLocalZKChannel<T: Field>(*mut LocalZKChannel<T>);
impl<T: Field> PLocalZKChannel<T> {
    pub fn new(t: *mut LocalZKChannel<T>) -> Self {
        PLocalZKChannel::<T>(t)
    }
}
unsafe impl<T: Field> Send for PLocalZKChannel<T> {}
unsafe impl<T: Field> Sync for PLocalZKChannel<T> {}
impl<T: Field> Deref for PLocalZKChannel<T> {
    type Target = LocalZKChannel<T>;

    fn deref(&self) -> &Self::Target {
        unsafe { &*self.0 }
    }
}
impl<T: Field> DerefMut for PLocalZKChannel<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        unsafe { &mut *self.0 }
    }
}


pub struct LocalVithChannel<T: Field> {
    plain: LocalChannel<T>,
    seedtree: SeedTree,
    // tapes: Vec<Tape>,
    // n_party: usize,
    n: usize,
    l: usize,
    m: usize,
    // log_m: usize,
    // log_n: usize,
    vith_c: Vec<T>,
    vith_u0: Vec<T>,
    vith_v0: Vec<T>,
    vith_u1_ptr: usize,
    vith_u2_ptr: usize,
    _fr: Vec<T>,
    _Vr: Vec<T>,
    _Vs: Vec<T>,
    _fr_mask: Vec<T>,
    _Vr_mask: Vec<T>,
    _Vs_mask: Vec<T>,
    // _fr_shares: Vec<Vec<T>>,
    // _Vr_shares: Vec<Vec<T>>,
    // _Vs_shares: Vec<Vec<T>>,
    // _p_shares: Vec<Vec<T>>,
    // _fr_poly_shares: Vec<Vec<T>>,
    // _Vr_poly_shares: Vec<Vec<T>>,
    // _Vs_poly_shares: Vec<Vec<T>>,
    _scalar: T,
    _pred: [T; 4],
    _r: T,
    _layer_id: usize,
    flag: bool,
    local_msg: Tx::<T>,
    plain_flag: bool,
}

impl<T: Field> LocalVithChannel<T> {
    pub fn new(mseed: &Seed, tape_len: usize, n_iter: usize, n_party: usize) -> Self {
        // tape_len is 2\ell in Vith-linear protocol
        // n_iter   is n     in Vith-linear protocol
        // n_party  is N     in Vith-linear protocol

        let log_n = ceil_log(n_iter);
        // TODO: configurable RS code rate
        let log_m = log_n - 1;
        let m = 1 << log_m;
        let n = n_iter;
        let r = n / m;

        let coset = Coset::<T>::init(log_n, log_m, 1);

        let start = Instant::now();

        let st = SeedTree::new(mseed, n_party, n_iter);
        let mut tapes: Vec<Tape> = vec![];
        for i in 0..n_iter {
            for p in 0..n_party {
                let seed = st.get_seed(i, p);
                // println!("{:?}", seed);
                tapes.push(Tape::new(seed, tape_len * T::SIZE));
            }
        }

        let mut vith_uu: Vec<T> = vec![T::from(0); tape_len * n];

        let mut vith_u0: Vec<T> = vec![T::from(0); tape_len * m];
        let mut vith_c: Vec<T> = vec![T::from(0); tape_len * (n - m)];
        let mut vith_v0: Vec<T> = vec![T::from(0); tape_len * n];

        if T::CHAR == 2 {
            // additive coset
            for i in 0..tape_len {
                // compute u0, v0
                for j in 0..n {
                    for k in 0..n_party {
                        let t = tapes[j * n_party + k].get::<T>();
                        // if i == tape_len / 2 && j < 2 && k != 0 {
                        //     println!("P tape{} party{}: {:?}", j, k, t);
                        // }
                        vith_uu[i * n + j] += t;
                        if j < m {
                            vith_u0[i * m + j] += t;
                        } else {
                            vith_c[i * (n - m) + (j - m)] += t;
                        }
                        // TODO configurable set S,
                        // here S = {1, 2, 3, ..., N}
                        vith_v0[i * n + j] += T::from(k + 1) * t;
                    }
                }
                // compute correction
                let tmp = (coset.fft)(
                    &coset, &(coset.ifft)(
                        &coset,
                        &vith_u0[i * m..(i + 1) * m].to_vec(),
                        log_m), log_m, log_n);
                for j in 0..m {
                    assert_eq!(tmp[j], vith_u0[i * m + j]);
                }
                for j in m..n {
                    vith_c[i * (n - m) + j - m] -= tmp[j];
                }
            }
        } else {
            // multiplicative coset
            for i in 0..tape_len {
                // compute u0, v0
                for j in 0..n {
                    for k in 0..n_party {
                        let t =  tapes[j * n_party + k].get::<T>();
                        if j % r == 0 {
                            vith_u0[i * m + (j / r)] += t;
                        } else {
                            vith_c[i * (n - m) + j - (j / r)] += t;
                        }
                        // TODO configurable set S,
                        // here S = {1, 2, 3, ..., N}
                        vith_v0[i * n + j] += T::from(k + 1) * t;
                    }
                }
                // compute correction
                let tmp = (coset.fft)(
                    &coset, &(coset.ifft)(
                        &coset,
                        &vith_u0[i * m..(i + 1) * m].to_vec(),
                        log_m), log_m, log_n);
                for j in 0..n {
                    if j % r != 0 {
                        vith_c[i * (n - m) + j - (j / r)] -= tmp[j];
                    } else {
                        assert_eq!(tmp[j], vith_u0[i * m + (j / r)]);
                    }
                }
            }
        }

        let elapsed = start.elapsed();
        let in_ms = elapsed.as_micros() as usize;
        mpc_prv_inc(in_ms);

        Self {
            plain: LocalChannel::<T>::new(),
            seedtree: st,
            // tapes: tapes,
            // n_party: n_party,
            n: n,
            l: tape_len,
            m: m,
            // log_m: log_m,
            // log_n: log_n,
            vith_c: vith_c,
            vith_u0: vith_u0,
            vith_v0: vith_v0,
            vith_u1_ptr: 0,
            vith_u2_ptr: tape_len / 2 * m,
            _fr: vec![],
            _fr_mask: vec![],
            _Vr: vec![],
            _Vr_mask: vec![],
            _Vs: vec![],
            _Vs_mask: vec![],
            // _fr_shares: vec![vec![]; n_party],
            // _Vr_shares: vec![vec![]; n_party],
            // _Vs_shares: vec![vec![]; n_party],
            // _p_shares: vec![vec![]; n_party],
            // _fr_poly_shares: vec![],
            // _Vr_poly_shares: vec![],
            // _Vs_poly_shares: vec![],
            _r: T::from(0),
            _scalar: T::from(0),
            _pred: [T::from(0); 4],
            _layer_id: 0,
            flag: false,
            local_msg: Tx::<T>::Syn,
            plain_flag: false,
        }
    }
    pub fn close(&self) {
        self.plain.close()
    }
}

impl<T: Field> Channel for LocalVithChannel<T> {
    type Output = T;

    const CHANNEL_TYPE: ChannelType = ChannelType::VOLEitH;

    fn get_type(&mut self) -> ChannelType {
        ChannelType::VOLEitH
    }

    fn force_plain(&mut self) {
        self.plain_flag = true;
    }

    fn send(&mut self, id: usize, tx: Tx::<Self::Output>) {
        if self.plain_flag {
            return self.plain.send(id, tx);
        }

        if id == 0 {
            channel0_rst();
            // println!("{}", tx);
        }
        let start0 = Instant::now();
        let flp0_ms = flp_prv_get();
        let mut recv_wait_ms = 0;
        match tx {
            Tx::<T>::OLayer(r) => {
                self._fr.push(r);
                // TODO???: remove constant terms to mask
                self._fr_mask.push(r);
                self.flag = true;
                self.local_msg = Tx::Syn;
            },
            Tx::<T>::Sumcheck(r) => {
                let mut msg: Vec<T> = vec![];
                for i in 0..r.len() {
                    msg.push(r[i] - self.vith_u0[self.vith_u1_ptr + i]);
                }
                self.plain.send(id, Tx::<T>::Sumcheck(msg));
        let start1 = Instant::now();
                let _r = recv!(self.plain, Random, 1 - id);
        let elapsed1 = start1.elapsed();
        recv_wait_ms += elapsed1.as_micros() as usize;

                let d = self._layer_id;

                if r.len() == 2 {
                    let tmp = self._fr[d] - r[0] - r[0] - r[1];
                    self._fr[d] = _r * ((_r * r[1]) + tmp) + r[0];

                    let r0 = self.vith_u0[self.vith_u2_ptr];
                    let r1 = self.vith_u0[self.vith_u2_ptr + 1];
                    let tmp = self._fr_mask[d] - r0 - r0 - r1;
                    self._fr_mask[d] = _r * ((_r * r1) + tmp) + r0;
                } else {
                    let tmp = self._fr[d] - r[0] - r[0] - r[1] - r[2];
                    self._fr[d] = ((r[2] * _r + r[1]) * _r + tmp) * _r + r[0];

                    let r0 = self.vith_u0[self.vith_u2_ptr];
                    let r1 = self.vith_u0[self.vith_u2_ptr + 1];
                    let r2 = self.vith_u0[self.vith_u2_ptr + 2];
                    let tmp = self._fr_mask[d] - r0 - r0 - r1 - r2;
                    self._fr_mask[d] = ((r2 * _r + r1) * _r + tmp) * _r + r0;
                }
                self.vith_u1_ptr += r.len();
                self.vith_u2_ptr += r.len();

                self.flag = true;
                self.local_msg = Tx::<T>::Random(_r);
            },
            Tx::<T>::Layer(r) => {
                let mut msg: Vec<T> = vec![];
                for i in 0..r.len() {
                    msg.push(r[i] - self.vith_u0[self.vith_u1_ptr + i]);
                }

                self._Vr.push(r[0]);
                self._Vr_mask.push(self.vith_u0[self.vith_u2_ptr]);

                if r.len() > 1 {
                    self._Vs.push(r[1]);
                    self._Vs_mask.push(self.vith_u0[self.vith_u2_ptr + 1])
                } else {
                    // XXX???: maybe some problems here when r.len() == 1
                    self._Vs.push(r[0]);
                    self._Vs_mask.push(self.vith_u0[self.vith_u2_ptr])
                }
                self.vith_u1_ptr += r.len();
                self.vith_u2_ptr += r.len();
                self._layer_id += 1;
                self.plain.send(id, Tx::<T>::Layer(msg));
            },
            Tx::<T>::LayerPred(r) => {
                self._pred = r;
                if self._layer_id > 0 {
                    let lid = self._layer_id - 1;
                    self._fr[lid] -= self._pred[1] * self._Vr[lid] + self._pred[2] * self._Vs[lid] + self._pred[3];
                    // XXX??? remove constant terms to mask
                    self._fr_mask[lid] -= self._pred[1] * self._Vr_mask[lid] + self._pred[2] * self._Vs_mask[lid] + self._pred[3];
                }
                self.flag = true;
                self.local_msg = Tx::Syn;
            },
            Tx::<T>::LayerScalar(r) => {
                self._scalar = r;
                let scalar0 = T::from(1) - r;
                let scalar1 = r;

                let lid = self._layer_id - 1;

                self._fr.push(
                    scalar0 * self._Vr[lid] +
                    scalar1 * self._Vs[lid]);
                self._fr_mask.push(
                    scalar0 * self._Vr_mask[lid] +
                    scalar1 * self._Vs_mask[lid]);
                self._Vr[lid] *= self._pred[0];
                self._Vr_mask[lid] *= self._pred[0];
                // println!("Vr_mask: {:?}", self._Vr_mask[lid]);
                // println!("Vs_mask: {:?}", self._Vs_mask[lid]);
                // println!("fr_mask: {:?}", self._fr_mask[lid]);
                self.flag = true;
                self.local_msg = Tx::Syn;
            },
            Tx::<T>::Input(r) => {
                // XXX
                // let mut msg: Vec<T> = vec![];
                // for w in r.iter() {
                //     let mut icr = T::from(0);
                //     for t in self.tapes.iter_mut() {
                //         icr += t.get::<T>();
                //     }
                //     icr = *w - icr;
                //     msg.push(icr);
                // }
                // self.plain.send(id, Tx::<T>::Input(msg));
                self.plain.send(id, Tx::<T>::Input(r));
            },

            Tx::<T>::FinalLayer(r) => {
                self._scalar = r[0];
                if r.len() > 1 {
                    self._r = r[1];
                } else {
                    self._r = r[0];
                }
                // let prod = r[0] * r[1];
                // self._r = prod;
                self.plain.send(id, Tx::<T>::FinalLayer(r));
            },

            Tx::<T>::FinalLayerPred(r) => {
                // self._scalar stores last layer v_a, v_b
                let tmp = r[0] * self._scalar * self._r +
                    r[1] * self._scalar + r[2] * self._r + r[3];
                self._r = tmp;
                self.flag = true;
                self.local_msg = Tx::Syn;
            }

            // XXX: weird
            Tx::<T>::MPCitH0 => {
                // TODO XXX merge hash with flpcp hash
                // TODO XXX merge hash of non secmul layer opens
                let mut h = SHA256::from(0);
                // assert_eq!(self._r, *self._fr.last().unwrap());
                // let layer_id = self._fr_shares[0].len() - 1;
                self._fr[self._layer_id] -= self._r;
                self._fr_mask[self._layer_id] -= self._r;

                h.commit(&[self._fr_mask[self._layer_id]]);

                self._fr.truncate(self._layer_id);
                self._fr_mask.truncate(self._layer_id);

                self.plain.send(id, Tx::<T>::FinalLayerOpen(h));
                // self.flag = true;
                // self.local_msg = Tx::Syn;
            },

            Tx::<T>::FLPCP => {
            flp_prv_rst();
            let start = Instant::now();

                let d = self._layer_id;
                let x_values_to_m = first_n_field_elements::<T>(d + 1);
                let precompute_to_m = precompute_lagrange_polynomials::<T>(&x_values_to_m);
                let x_values_to_2m = first_n_field_elements::<T>(2 * d + 1);
                let precompute_to_2m = precompute_lagrange_polynomials::<T>(&x_values_to_2m);

            let elapsed = start.elapsed();
            let in_ms = elapsed.as_micros() as usize;
            flp_prv_inc(in_ms);

                self._fr.push(T::random());
                self._Vr.push(T::random());
                self._Vs.push(T::random());
                self._fr_mask.push(self.vith_u0[self.vith_u2_ptr + 0]);
                self._Vr_mask.push(self.vith_u0[self.vith_u2_ptr + 1]);
                self._Vs_mask.push(self.vith_u0[self.vith_u2_ptr + 2]);
                assert_eq!(self._fr.len(), d + 1);
                let random_commit = vec![
                    self._fr[d] - self.vith_u0[self.vith_u1_ptr + 0],
                    self._Vr[d] - self.vith_u0[self.vith_u1_ptr + 1],
                    self._Vs[d] - self.vith_u0[self.vith_u1_ptr + 2]];
                self.vith_u1_ptr += 3;
                self.vith_u2_ptr += 3;

                self.plain.send(id, Tx::<T>::FLPCPRandom(random_commit.to_vec()));
        let start1 = Instant::now();
                let _ = recv!(self.plain, @Syn, 1 - id);
        let elapsed1 = start1.elapsed();
        recv_wait_ms += elapsed1.as_micros() as usize;

                let start = Instant::now();
                let mut poly = vec![T::from(0); 2 * d + 1];
                let mut poly_mask = vec![T::from(0); 2 * d + 1];
                let mut poly_commit = vec![T::from(0); d + 1];

                let Vr_poly = interpolate_with_precomputation(&precompute_to_m, &self._Vr);
                let Vs_poly = interpolate_with_precomputation(&precompute_to_m, &self._Vs);
                let fr_poly = interpolate_with_precomputation(&precompute_to_m, &self._fr);
                // let Vr_mask_poly = interpolate_with_precomputation(&precompute_to_m, &self._Vr_mask);
                // let Vs_mask_poly = interpolate_with_precomputation(&precompute_to_m, &self._Vs_mask);
                // let fr_mask_poly = interpolate_with_precomputation(&precompute_to_m, &self._fr_mask);
                let p_poly = poly_mul_sub(&Vr_poly, &Vs_poly, &fr_poly);
                for k in 0..(d + 1) {
                    let x = x_values_to_2m[k + d];
                    poly[k + d] = poly_eval(&p_poly, x);
                    poly_mask[k + d] = self.vith_u0[self.vith_u2_ptr + k];
                    poly_commit[k] = poly[k + d] - self.vith_u0[self.vith_u1_ptr + k];
                    // println!("P poly {} {:?}", k, poly_mask[k+d]+poly[k+d]);
                    // println!("a=1, poly_mask {} {:?}", k, self.vith_u0[self.vith_u2_ptr + k] + self.vith_u0[self.vith_u1_ptr + k])
                }
                self.vith_u1_ptr += d + 1;
                self.vith_u2_ptr += d + 1;

                let elapsed = start.elapsed();
                let in_ms = elapsed.as_micros() as usize;
                flp_prv_inc(in_ms);

                self.plain.send(id, Tx::<T>::FLPCPCommit(poly_commit.to_vec()));

        let start1 = Instant::now();
                let _r = recv!(self.plain, Random, 1 - id);
        let elapsed1 = start1.elapsed();
        recv_wait_ms += elapsed1.as_micros() as usize;

                // let start = Instant::now();

                let mut lagrange_poly_m_eval_r = vec![T::from(0); d + 1];
                let mut lagrange_poly_2m_eval_r = vec![T::from(0); 2 * d + 1];

                for k in 0..(d + 1) {
                    lagrange_poly_m_eval_r[k] = poly_eval(&precompute_to_m[k], _r);
                }
                for k in 0..(2 * d + 1) {
                    lagrange_poly_2m_eval_r[k] = poly_eval(&precompute_to_2m[k], _r);
                }

                // let elapsed = start.elapsed();
                // let in_ms = elapsed.as_micros() as usize;
                // flp_prv_inc(in_ms);

                let a = dot_product(&lagrange_poly_m_eval_r, &self._Vr);
                let b = dot_product(&lagrange_poly_m_eval_r, &self._Vs);
                let c = dot_product(&lagrange_poly_m_eval_r, &self._fr);
                // let d = dot_product(&lagrange_poly_2m_eval_r, &p);
                let mut output = vec![T::from(0); 4];
                output[0] = dot_product(&lagrange_poly_m_eval_r, &self._Vr_mask) - a;
                output[1] = dot_product(&lagrange_poly_m_eval_r, &self._Vs_mask) - b;
                output[2] = dot_product(&lagrange_poly_m_eval_r, &self._fr_mask) - c;
                output[3] = dot_product(&lagrange_poly_2m_eval_r, &poly_mask) - a * b + c;
                // println!("p output: {:?}", output);

                // assert_eq!(a*b-c, d);

                let mut h = SHA256::from(0);
                h.commit(&output);

                self.plain.send(id, Tx::<T>::FLPCPOpen((h, a, b, c)));
                // self.flag = true;
                // self.local_msg = Tx::Syn;
            },
            Tx::<T>::MPCitH => {
                self.plain.send(id, Tx::<T>::MPCitH);
        let start1 = Instant::now();
                let alpha = recv!(self.plain, Random, 1 - id);
        let elapsed1 = start1.elapsed();
        recv_wait_ms += elapsed1.as_micros() as usize;

                let lm = self.l / 2 * self.m;
                let ln = self.l / 2 * self.n;
                let mut vith_s: Vec<T> = vec![T::from(0); lm];
                let mut vith_v: Vec<T> = vec![T::from(0); ln];
                for i in 0..lm {
                    vith_s[i] = alpha * self.vith_u0[i] + self.vith_u0[i + lm];
                }
                for i in 0..ln {
                    vith_v[i] = alpha * self.vith_v0[i] + self.vith_v0[i + ln];
                }
                let mut h = SHA256::from(0);
                h.commit(&vith_v);
                // TODO: vith_c at begining
                self.plain.send(id, Tx::<T>::VOLEOpen((self.vith_c.clone(), vith_s, h)));

                // TODO seedtree commit root should be sent in the very first round
                let root = self.seedtree.commit();
                let mut sto = SeedTreeOpen::new(&self.seedtree);
        let start1 = Instant::now();
                sto.missings = recv!(self.plain, Party, 1 - id);
        let elapsed1 = start1.elapsed();
        recv_wait_ms += elapsed1.as_micros() as usize;

                self.seedtree.open(&mut sto);

                self.plain.send(id, Tx::<T>::Seed((root, sto.clone())));
            },
            // Tx::<T>::Random(_) => self.plain.send(id, tx),
            // Tx::<T>::MultiRandom(_) => self.plain.send(id, tx),
            // Tx::<T>::Party(_) => self.plain.send(id, tx),
            // Tx::<T>::VPDOpen0(_) => self.plain.send(id, tx),
            // Tx::<T>::VPDOpen1(_) => self.plain.send(id, tx),
            // Tx::<T>::Syn => self.plain.send(id, tx),
            _ => {
                self.plain.send(id, tx)
                // println!("{}", tx);
                // unreachable!()
            },
        }
        if id == 0 {
            let elapsed = start0.elapsed();
            let in_ms = elapsed.as_micros() as usize;
            let wait_ms = channel0_get();
            let flp_ms = flp_prv_get() - flp0_ms;
            // println!("{} - {} - {} - {} = {}", in_ms, wait_ms, flp_ms, recv_wait_ms, in_ms - wait_ms - flp_ms - recv_wait_ms);
            channel0_rst();
            mpc_prv_inc(in_ms - wait_ms - flp_ms - recv_wait_ms);
        }
    }

    fn recv(&mut self, id: usize) -> Tx::<Self::Output> {
        let start = Instant::now();

        if self.flag {
            self.flag = false;
            return self.local_msg.clone();
        }
        let res = self.plain.recv(id);

        let elapsed = start.elapsed();
        let in_ms = elapsed.as_micros() as usize;
        if id == 1 {
            channel0_inc(in_ms);
        } else {
            channel1_inc(in_ms);
        }
        res
        // match res {
        //     Tx::<T>::Random(r) => {
        //         self._r = r;
        //         res
        //     },
        //     _ => res,
        // }
    }
}

pub struct PLocalVithChannel<T: Field>(*mut LocalVithChannel<T>);
impl<T: Field> PLocalVithChannel<T> {
    pub fn new(t: *mut LocalVithChannel<T>) -> Self {
        PLocalVithChannel::<T>(t)
    }
}
unsafe impl<T: Field> Send for PLocalVithChannel<T> {}
unsafe impl<T: Field> Sync for PLocalVithChannel<T> {}
impl<T: Field> Deref for PLocalVithChannel<T> {
    type Target = LocalVithChannel<T>;

    fn deref(&self) -> &Self::Target {
        unsafe { &*self.0 }
    }
}
impl<T: Field> DerefMut for PLocalVithChannel<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        unsafe { &mut *self.0 }
    }
}
