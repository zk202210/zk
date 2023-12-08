use std::time::Instant;
use std::thread;
use std::sync::Arc;
use std::fs::OpenOptions;
use std::io::Write;

use rustzk::field::*;
use rustzk::channel::*;
use rustzk::circuit::*;
use rustzk::statistic::*;
use rustzk::zkp::*;
use rustzk::tape::*;
use rustzk::util::*;


const N_PARTY: usize = 16;
const LOG_RHO: usize = 4;
const LOG_LAYER_RATIO: usize = 1;
// 0 for plain_VDP, 1 for mpc_VDP, 2 for Virgo's VPD
const VPD_MODE: usize = 2;
const N_MPC_ITER: usize = 33;


fn type_name<'a, T>() -> &'a str {
    let s = std::any::type_name::<T>();
    let parts: Vec<&str> = s.rsplit(":").collect();
    parts[0]
}


fn zk_vpd<T: Field>(log_len: usize) {

    stats_reset_all();
    let n_ldt_iter = get_n_ldt_iter(
        log_len, T::SIZE * 8,
        (1 << LOG_RHO) / (if VPD_MODE == 3 { 2 } else { 1 })
        - if VPD_MODE == 1 { 0 } else { 1 });
    // let n_ldt_iter = 33;
    println!("LDT num of queries: {}", n_ldt_iter);

    // let tape_len: usize = 29 * (log_len + 11) + 13;
    let tape_len: usize = 4 * log_len + 11;

    let mseed = Seed::random();

    let mut lc = LocalZKChannel::<T>::new(&mseed, tape_len, N_PARTY);
    let mut lc0 = PLocalZKChannel::<T>::new(&mut lc);
    let mut lc1 = PLocalZKChannel::<T>::new(&mut lc);

    let start = Instant::now();

    let mut c = Circuit::<T>::random(3, 1, 1, log_len);
    let mut inputs: Vec<T> = vec![];
    for i in 0..(1 << c.layers[0].bit_len) {
        match c.layers[0].gates[i].gtype {
            GateType::IN => {
                let r = T::random();
                inputs.push(r);
            },
            GateType::DUMMY => {
                inputs.push(T::from(0));
            },
            _ => unreachable!(),
        }
    }
    c.eval(&inputs);

    let mut c_v = c.clone();
    let outputs = c.layers.last().unwrap().values.clone();

    let elapsed = start.elapsed();
    let in_ms = elapsed.as_micros();
    println!("random circuit generated in {} ms", in_ms as f64 / 1000.0);

    let start = Instant::now();
    let log_l = log_len + LOG_RHO;
    let mut coset = Coset::<T>::init(log_l, LOG_RHO, LOG_LAYER_RATIO);
    coset.precompute(n_ldt_iter, VPD_MODE < 2);

    let coset = Arc::new(coset);
    // let coset_v = coset.clone();

    let elapsed = start.elapsed();
    let in_ms = elapsed.as_micros();
    println!("coset initialized in {} ms", in_ms as f64 / 1000.0);

    let coset_p = coset.clone();
    let p_thread = thread::spawn(move || {
        let lc = &mut *lc0;
        let mut prover = Prover::<T, LocalZKChannel<T>>::new(
            log_len, &mut c, &inputs, lc);
        if VPD_MODE == 0 {
            prover.prove_plain(n_ldt_iter, &coset_p);
        } else if VPD_MODE == 1 {
            prover.prove_mpcvpd(&coset_p, n_ldt_iter, N_PARTY, N_MPC_ITER);
        } else if VPD_MODE == 2 {
            prover.prove_v2(n_ldt_iter, &coset_p);
        } else {
            prover.prove(n_ldt_iter, &coset_p);
        }
    });

    let coset_v = coset.clone();
    let v_thread = thread::spawn(move || {
        let lc = &mut *lc1;
        let mut verifier = Verifier::<T, LocalZKChannel<T>>::new(
            log_len, &mut c_v, lc);
        let res =
            if VPD_MODE == 0 {
                verifier.verify_plain(&outputs, N_PARTY, tape_len, n_ldt_iter, &coset_v)
            } else if VPD_MODE == 1 {
                verifier.verify_mpcvpd(&outputs, &coset_v, N_PARTY, tape_len, n_ldt_iter, N_MPC_ITER)
            } else if VPD_MODE == 2 {
                verifier.verify_v2(&outputs, N_PARTY, tape_len, n_ldt_iter, &coset_v)
            } else {
                verifier.verify(&outputs, N_PARTY, tape_len, n_ldt_iter, &coset_v)
            };
        assert!(res);
    });

    v_thread.join().unwrap();
    p_thread.join().unwrap();
    lc.close();

    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(std::format!("log/zk_vpd_{}_{}.log", VPD_MODE, type_name::<T>()))
        .unwrap();
    statistic(&mut file, log_len, n_ldt_iter, N_MPC_ITER);
}

fn test_vpd<T: Field>() {
    println!("Field: {}", type_name::<T>());
    println!("RS code rate: 1/{}", 1 << LOG_RHO);
    println!("Coset ratio: {}", 1 << LOG_LAYER_RATIO);
    println!("MPC num of parties: {}", N_PARTY);
    let vpd_mode = if VPD_MODE == 0 {
        "plain VPD"
    } else if VPD_MODE == 1 {
        "mpc VPD"
    } else if VPD_MODE == 2 {
        "simplified Virgo's VPD"
    } else {
        "Virgo's VPD"
    };
    println!("VPD mode: {}", vpd_mode);
    println!("--------------------------------------------------------------------------");

    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(false)
        .truncate(true)
        .open(std::format!("log/zk_vpd_{}_{}.log", VPD_MODE, type_name::<T>()))
        .unwrap();
    writeln!(file, "Field: {}", type_name::<T>()).unwrap();
    writeln!(file, "RS code rate: 1/{}", 1 << LOG_RHO).unwrap();
    writeln!(file, "Coset ratio: {}", 1 << LOG_LAYER_RATIO).unwrap();
    writeln!(file, "MPC num of parties: {}", N_PARTY).unwrap();
    writeln!(file, "VPD_MODE: {}", vpd_mode).unwrap();

    for log_len in 14..=21 {
    // for log_len in 10..11 {
        println!("log inputs: {0}", log_len);
        zk_vpd::<T>(log_len);
    }
}

#[test]
fn test_vpd_all() {
    println!("");
    // test_zk_aes::<GF2p16>();
    // test_zk_aes::<GF2p32>();
    // test_zk_aes::<GF2p64>();
    // test_zk_aes::<GF2p128>();
    // test_vpd::<GF2p192>();
    test_vpd::<GF2p256>();
    // test_vpd::<M127p2>();
    // test_vpd::<M61p2>();
}
