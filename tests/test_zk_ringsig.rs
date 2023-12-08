use std::time::Instant;
use std::thread;
use std::sync::Arc;
use std::fs::OpenOptions;
use std::io::Write;

use rand::Rng;

use rustzk::app::ringsig::*;
use rustzk::field::*;
use rustzk::channel::*;
use rustzk::statistic::*;
use rustzk::zkp::*;
use rustzk::tape::*;
use rustzk::util::*;

const N_PARTY: usize = 16;
const LOG_RHO: usize = 4;
const LOG_LAYER_RATIO: usize = 1;
// 0 for plain_VDP, 1 for mpc_VDP, 2 for Virgo's VPD
const VPD_MODE: usize = 0;
const N_MPC_ITER: usize = 33;

// input + (round * 2 + 2) * depth + 3 * depth + (depth + 1)
// const tape_len: usize = (1 << (2 * log_len + 1)) + log_len * 2 + 2 + 3 + 2;
//


fn type_name<'a, T>() -> &'a str {
    let s = std::any::type_name::<T>();
    let parts: Vec<&str> = s.rsplit(":").collect();
    parts[0]
}

fn test_ringsig_circuit<T: Field>() {
    use rustzk::gkr;

    for log_ring_size in 1..12 {
        let path_size = log_ring_size + 1;
        let log_num_hash = ceil_log(path_size);

        let mut lc = LocalChannel::<T>::new();
        let mut lc0 = PLocalChannel::<T>::new(&mut lc);
        let mut lc1 = PLocalChannel::<T>::new(&mut lc);

        let mut inputs = vec![];
        for _ in 0..path_size {
            let mut tmp = [0u64; 4];
            for i in 0..4 {
                tmp[i] = rand::thread_rng().gen();
            }
            inputs.push(tmp);
        }
        let inputs = gen_ringsig_witness::<T>(&inputs, log_num_hash);


        let mut c = rustzk::app::ringsig::build_circuit(path_size);
        println!("depth {}", c.depth);
        c.eval(&inputs);

        // let log_len = log_num_hash + 9;

        let mut c_v = c.clone();
        let outputs = c.layers.last().unwrap().values.clone();
        for i in 0..outputs.len() {
            assert_eq!(outputs[i], T::from(0));
        }

        let p_thread = thread::spawn(move || {
            let lc = &mut *lc0;
            let mut prover = gkr::Prover::<T, LocalChannel<T>>::new(&mut c, lc);
            prover.prove();
        });

        let v_thread = thread::spawn(move || {
            let lc = &mut *lc1;
            let mut verifier = gkr::Verifier::<T, LocalChannel<T>>::new(&mut c_v, lc);
            let res = verifier.verify(&outputs);
            assert!(res);
        });
        v_thread.join().unwrap();
        p_thread.join().unwrap();
        lc.close();
    }
}


fn zk_ringsig_w_vpd<T: Field>(log_ring_size: usize) {
    let path_size = log_ring_size + 1;
    let log_num_hash = ceil_log(path_size);
    stats_reset_all();
    let n_ldt_iter = get_n_ldt_iter(
        log_num_hash + 9, T::SIZE * 8,
        (1 << LOG_RHO) * (if VPD_MODE == 2 { 2 } else { 1 })
        - if VPD_MODE == 1 { 0 } else { 1 });
    println!("LDT num of queries: {}", n_ldt_iter);

    println!("log_num_hash {}", log_num_hash);
    let tape_len: usize = 3672 / 24 * log_num_hash + 46488 / 24;

    let mseed = Seed::random();

    let mut lc = LocalZKChannel::<T>::new(&mseed, tape_len, N_PARTY);
    let mut lc0 = PLocalZKChannel::<T>::new(&mut lc);
    let mut lc1 = PLocalZKChannel::<T>::new(&mut lc);

    let start = Instant::now();

    let mut inputs = vec![];
    for _ in 0..path_size {
        let mut tmp = [0u64; 4];
        for i in 0..4 {
            tmp[i] = rand::thread_rng().gen();
        }
        inputs.push(tmp);
    }
    let inputs = gen_ringsig_witness::<T>(&inputs, log_num_hash);

    let mut c = rustzk::app::ringsig::build_circuit(path_size);
    c.eval(&inputs);

    let log_len = log_num_hash + 9;

    let mut c_v = c.clone();
    let outputs = c.layers.last().unwrap().values.clone();
    // for i in 0..outputs.len() {
    //     assert_eq!(outputs[i], T::from(0));
    // }

    let elapsed = start.elapsed();
    let in_ms = elapsed.as_micros();
    println!("ring signature mt & aes witnesses generated in {} ms", in_ms as f64 / 1000.0);

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
        .open(std::format!("log/zk_ringsig_{}.log", type_name::<T>()))
        .unwrap();
    statistic(&mut file, log_ring_size, n_ldt_iter, N_MPC_ITER);
}


fn test_zk_ringsig<T: Field>() {
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
        .open(std::format!("log/zk_ringsig_{}.log", type_name::<T>()))
        .unwrap();
    writeln!(file, "Field: {}", type_name::<T>()).unwrap();
    writeln!(file, "RS code rate: 1/{}", 1 << LOG_RHO).unwrap();
    writeln!(file, "Coset ratio: {}", 1 << LOG_LAYER_RATIO).unwrap();
    writeln!(file, "MPC num of parties: {}", N_PARTY).unwrap();
    writeln!(file, "VPD_MODE: {}", vpd_mode).unwrap();

    for log_ring_size in 1..41 {
        println!("ring size: 2^{0}", log_ring_size);
        // println!("log inputs: {0}", log_ring_size + 10);
        zk_ringsig_w_vpd::<T>(log_ring_size);

        use std::time;
        let one_sec = time::Duration::from_secs(2);
        thread::sleep(one_sec);
    }
}

#[test]
fn test_zk_ringsig_all() {
    println!("");
    // test_ringsig_circuit::<GF2p16>();
    // test_zk_ringsig::<GF2p16>();
    // test_zk_ringsig::<GF2p32>();
    // test_zk_ringsig::<GF2p64>();
    // test_zk_ringsig::<GF2p128>();
    test_zk_ringsig::<GF2p192>();
    // test_zk_ringsig::<GF2p256>();
}
