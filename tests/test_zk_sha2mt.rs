use std::time::Instant;
use std::thread;
use std::sync::Arc;
use std::fs::OpenOptions;
use std::io::Write;

use rustzk::field::*;
use rustzk::circuit::*;
use rustzk::channel::*;
use rustzk::statistic::*;
use rustzk::zkp::*;
use rustzk::tape::*;
use rustzk::util::*;

const N_PARTY: usize = 16;
const LOG_RHO: usize = 4;
// const N_LDT_ITER: usize = 33;
const LOG_LAYER_RATIO: usize = 3;
const MPC_VPD: bool = true;
const N_MPC_ITER: usize = 33;

// input + (round * 2 + 2) * depth + 3 * depth + (depth + 1)
// const tape_len: usize = (1 << (2 * log_len + 1)) + log_len * 2 + 2 + 3 + 2;
//


fn type_name<'a, T>() -> &'a str {
    let s = std::any::type_name::<T>();
    let parts: Vec<&str> = s.rsplit(":").collect();
    parts[0]
}

fn test_sha2mt<T: Field>() {
    use rustzk::gkr;

    let mut lc = LocalChannel::<T>::new();
    let mut lc0 = PLocalChannel::<T>::new(&mut lc);
    let mut lc1 = PLocalChannel::<T>::new(&mut lc);

    let mut c = Circuit::<T>::read("resources/sha256mt/SHA256_64_merkle_1_circuit.txt");
    let inputs = c.layers[0].values.clone();
    c.eval(&inputs);

    let mut cv = c.clone();

    c.eval(&inputs);

    let mut outputs: Vec<T> = vec![];
    let output_layer = c.layers.last().unwrap();
    for i in 0..(1 << (output_layer.bit_len)) {
        outputs.push(output_layer.values[i]);
    }

    let p_thread = thread::spawn(move || {
        let lc = &mut *lc0;
        let mut prover = gkr::Prover::<T, LocalChannel<T>>::new(&mut c, lc);
        prover.prove();
    });

    let v_thread = thread::spawn(move || {
        let lc = &mut *lc1;
        let mut verifier = gkr::Verifier::<T, LocalChannel<T>>::new(&mut cv, lc);
        let res = verifier.verify(&outputs);
        assert!(res);
    });

    v_thread.join().unwrap();
    p_thread.join().unwrap();
    lc.close();
}

fn zk_sha2mt_w_vpd<T: Field>(log_leaf_size: usize) {
    stats_reset_all();
    virgo_flag_inc(1);

    let tape_len: usize = 44 * log_leaf_size + 578;

    let mseed = Seed::random();

    let mut lc = LocalZKChannel::<T>::new(&mseed, tape_len, N_PARTY);
    let mut lc0 = PLocalZKChannel::<T>::new(&mut lc);
    let mut lc1 = PLocalZKChannel::<T>::new(&mut lc);

    let start = Instant::now();

    // let mut inputs = vec![];
    // for _ in 0..(1 << log_leaf_size) {
    //     let mut tmp = [0u64; 8];
    //     for i in 0..8 {
    //         tmp[i] = rand::thread_rng().gen();
    //     }
    //     inputs.push(tmp);
    // }
    // let inputs = gen_sha3mt_witness::<T>(&inputs, log_leaf_size);

    let mut c = Circuit::<T>::read_with_meta(
        &std::format!("resources/sha256mt/SHA256_64_merkle_{}_circuit.txt", log_leaf_size),
        &std::format!("resources/sha256mt/SHA256_64_merkle_{}_meta.txt", log_leaf_size));
    let inputs = c.layers[0].values.clone();
    // println!("wtf");
    // let mut c = rustzk::app::sha3mt::build_circuit(log_leaf_size);
    c.eval(&inputs);


    // let log_len = log_leaf_size + 9 + 1;
    let log_len = c.layers[0].bit_len;
    let n_ldt_iter = get_n_ldt_iter(log_len, T::SIZE * 8, 1 << LOG_RHO);
    let n_ldt_iter = 33;
    println!("LDT num of queries: {}", n_ldt_iter);

    let mut c_v = c.clone();
    let outputs = c.layers.last().unwrap().values.clone();
    // for i in 0..outputs.len() {
    //     assert_eq!(outputs[i], T::from(0));
    // }

    let elapsed = start.elapsed();
    let in_ms = elapsed.as_micros();
    println!("sha256 merkle tree generated in {} ms", in_ms as f64 / 1000.0);

    let start = Instant::now();
    let log_l = log_len + LOG_RHO;
    let mut coset = Coset::<T>::init(log_l, LOG_RHO, LOG_LAYER_RATIO);
    coset.precompute(n_ldt_iter, MPC_VPD);

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
        if MPC_VPD {
            prover.prove_plain(n_ldt_iter, &coset_p);
            // prover.prove_mpcvpd(&coset_p, n_ldt_iter, N_PARTY, N_MPC_ITER);
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
            if MPC_VPD {
                verifier.verify_plain(&outputs, N_PARTY, tape_len, n_ldt_iter, &coset_v)
                // verifier.verify_mpcvpd(&outputs, &coset_v, N_PARTY, tape_len, n_ldt_iter, N_MPC_ITER)
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
        .open(std::format!("log/zk_sha2mt_{}.log", type_name::<T>()))
        .unwrap();
    statistic(&mut file, log_leaf_size, n_ldt_iter, N_MPC_ITER);
}


fn test_zk_sha2mt<T: Field>() {
    println!("Field: {}", type_name::<T>());
    println!("RS code rate: 1/{}", 1 << LOG_RHO);
    println!("Coset ratio: {}", 1 << LOG_LAYER_RATIO);
    println!("MPC num of parties: {}", N_PARTY);
    println!("--------------------------------------------------------------------------");

    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(false)
        .truncate(true)
        .open(std::format!("log/zk_sha2mt_{}.log", type_name::<T>()))
        .unwrap();
    writeln!(file, "Field: {}", type_name::<T>()).unwrap();
    writeln!(file, "RS code rate: 1/{}", 1 << LOG_RHO).unwrap();
    writeln!(file, "Coset ratio: {}", 1 << LOG_LAYER_RATIO).unwrap();
    writeln!(file, "MPC num of parties: {}", N_PARTY).unwrap();

    for log_leaf_size in 0..10 {
        println!("num of leaves: {0}", 1 << log_leaf_size);
        zk_sha2mt_w_vpd::<T>(log_leaf_size);
    }
}

#[test]
fn test_zk_sha2mt_all() {
    // println!("");
    // test_sha2mt::<M61p2>();
    println!("");
    test_zk_sha2mt::<M61p2>();
    // test_zk_sha2mt::<M127p2>();
    // test_zk_sha2mt::<P192>();
    // test_zk_sha2mt::<GF2p32>();
    // test_zk_sha2mt::<GF2p64>();
    // test_zk_sha2mt::<GF2p128>();
    // test_zk_sha2mt::<GF2p256>();
}
