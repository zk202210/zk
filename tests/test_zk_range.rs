use std::time::Instant;
use std::thread;
use std::sync::Arc;
use std::fs::OpenOptions;
use std::io::Write;

use rand::Rng;

use rustzk::field::*;
use rustzk::channel::*;
use rustzk::statistic::*;
use rustzk::zkp::*;
use rustzk::tape::*;

const N_PARTY: usize = 16;
const LOG_RHO: usize = 4;
const LOG_LAYER_RATIO: usize = 1;
// 0 for plain_VDP, 1 for mpc_VDP, 2 for Virgo's VPD
const VPD_MODE: usize = 0;
const N_MPC_ITER: usize = 33;
const FAST_MODE: bool = true;

fn type_name<'a, T>() -> &'a str {
    let s = std::any::type_name::<T>();
    let parts: Vec<&str> = s.rsplit(":").collect();
    parts[0]
}


fn zk_range_w_vpd<T: Field>(log_range: usize) {
    stats_reset_all();
    let n_ldt_iter = 1;

    let tape_len: usize =
        if FAST_MODE {
            2 * log_range + 2
        } else {
            (2 * log_range + 1) * log_range + 1
        };

    let mseed = Seed::random();

    let mut lc = LocalZKChannel::<T>::new(&mseed, tape_len, N_PARTY);
    let mut lc0 = PLocalZKChannel::<T>::new(&mut lc);
    let mut lc1 = PLocalZKChannel::<T>::new(&mut lc);

    let start = Instant::now();

    let r: u64 = rand::thread_rng().gen::<u64>() % (1 << log_range);
    let inputs = vec![T::from(r as usize), T::from(0), T::from(0), T::from(0)];

    let mut c = rustzk::app::range::build_circuit(log_range, FAST_MODE);
    c.eval(&inputs);

    let log_len = 2;

    let mut c_v = c.clone();
    let outputs = c.layers.last().unwrap().values.clone();
    for i in 0..outputs.len() {
        assert_eq!(outputs[i], T::from(0));
    }

    let elapsed = start.elapsed();
    let in_ms = elapsed.as_micros();
    println!("range circuit evaluated in {} ms", in_ms as f64 / 1000.0);

    let start = Instant::now();
    let log_l = log_len + LOG_RHO;
    let mut coset = Coset::<T>::init(log_l, LOG_RHO, LOG_LAYER_RATIO);
    coset.precompute(n_ldt_iter, VPD_MODE != 2);

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
        .open(std::format!("log/zk_range_{}_{}.log", if FAST_MODE { "fast" } else { "slow" }, type_name::<T>()))
        .unwrap();
    statistic(&mut file, log_range, n_ldt_iter, N_MPC_ITER);
}


fn test_zk_range<T: Field>() {
    println!("Field: {}", type_name::<T>());
    println!("RS code rate: 1/{}", 1 << LOG_RHO);
    println!("Coset ratio: {}", 1 << LOG_LAYER_RATIO);
    println!("MPC num of parties: {}", N_PARTY);
    let vpd_mode = if VPD_MODE == 0 {
        "plain VPD"
    } else if VPD_MODE == 1 {
        "mpc VPD"
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
        .open(std::format!("log/zk_range_{}_{}.log", if FAST_MODE { "fast" } else { "slow" }, type_name::<T>()))
        .unwrap();
    writeln!(file, "Field: {}", type_name::<T>()).unwrap();
    writeln!(file, "RS code rate: 1/{}", 1 << LOG_RHO).unwrap();
    writeln!(file, "Coset ratio: {}", 1 << LOG_LAYER_RATIO).unwrap();
    writeln!(file, "MPC num of parties: {}", N_PARTY).unwrap();
    writeln!(file, "VPD_MODE: {}", vpd_mode).unwrap();

    for log_range in (4..=24).step_by(4) {
        println!("range size: 2^{0}", log_range);
        // println!("log inputs: {0}", log_range + 10);
        zk_range_w_vpd::<T>(log_range);

        use std::time;
        let one_sec = time::Duration::from_secs(1);
        thread::sleep(one_sec);
    }
}

#[test]
fn test_zk_range_all() {
    // test_range_circuit::<GF2p16>();
    println!("");
    // test_zk_range::<GF2p16>();
    // test_zk_range::<GF2p32>();
    // test_zk_range::<GF2p64>();
    // test_zk_range::<GF2p128>();
    test_zk_range::<GF2p192>();
    // test_zk_range::<GF2p256>();
}
