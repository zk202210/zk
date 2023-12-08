use std::time::Instant;
use std::thread;
use std::sync::Arc;
use std::fs::OpenOptions;
use std::io::Write;

use rustzk::matmul::*;
use rustzk::field::*;
use rustzk::channel::*;
use rustzk::statistic::*;
use rustzk::tape::*;
use rustzk::util::*;

// const N_PARTY: usize = 16;
// const N_LDT_ITER: usize = 33;
// const LOG_RHO: usize = 5;
// const LOG_LAYER_RATIO: usize = 3;

const N_PARTY: usize = 16;
const LOG_RHO: usize = 4;
// const N_LDT_ITER: usize = 65;
const LOG_LAYER_RATIO: usize = 1;
const MPC_VPD: bool = false;
const N_MPC_ITER: usize = 1;

// input + (round * 2 + 2) * depth + 3 * depth + (depth + 1)
// const tape_len: usize = (1 << (2 * log_len + 1)) + log_len * 2 + 2 + 3 + 2;
//


fn type_name<'a, T>() -> &'a str {
    let s = std::any::type_name::<T>();
    let parts: Vec<&str> = s.rsplit(":").collect();
    parts[0]
}


fn zk_matmul_w_vpd<T: Field>(log_len: usize) {
    stats_reset_all();
    let n_ldt_iter = get_n_ldt_iter(log_len * 2 + 1, T::SIZE * 8, 1 << LOG_RHO);
    println!("LDT num of queries: {}", n_ldt_iter);

    let n_ldt_iter = 5;
    // let N_MPC_ITER = 1;

    let tape_len: usize = log_len * 2 + 2 + 3 + 2;

    let mseed = Seed::random();

    let mut lc = LocalZKChannel::<T>::new(&mseed, tape_len, N_PARTY);
    let mut lc0 = PLocalZKChannel::<T>::new(&mut lc);
    let mut lc1 = PLocalZKChannel::<T>::new(&mut lc);

    let start = Instant::now();

    let mut inputs: Vec<T> = vec![];
    let mut outputs: Vec<T> = vec![];

    for _ in 0..(1 << (2 * log_len + 1)) {
        inputs.push(T::random());
    }

    //   h  l       h  l      h  l
    // C[i][j] += A[k][j] * B[k][i]
    let offset: usize = 1 << (2 * log_len);
    for i in 0..(1 << log_len) {
        for j in 0..(1 << log_len) {
            let mut prod = T::from(0);
            for k in 0..(1 << log_len) {
                prod += inputs[(k << log_len) + j] * inputs[(k << log_len) + i + offset];
            }
            outputs.push(prod);
        }
    }
    let outputs_v = outputs.to_vec();
    let elapsed = start.elapsed();
    let in_ms = elapsed.as_micros();
    println!("matrices generated in {} ms", in_ms as f64 / 1000.0);

    let start = Instant::now();
    let log_l = log_len * 2 + 1 + LOG_RHO;
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
            log_len, &inputs, &outputs, lc);
        if MPC_VPD {
            prover.prove_w_mpcvpd(&coset_p, n_ldt_iter, N_PARTY, N_MPC_ITER);
        } else {
            prover.prove_w_vpdv2(n_ldt_iter, &coset_p);
            // prover.prove_w_vpd(n_ldt_iter, &coset_p);
        }
    });

    let coset_v = coset.clone();
    let v_thread = thread::spawn(move || {
        let lc = &mut *lc1;
        let mut verifier = Verifier::<T, LocalZKChannel<T>>::new(
            log_len, lc);
        let res =
            if MPC_VPD {
                verifier.zk_verify_w_mpcvpd(&outputs_v, &coset_v, N_PARTY, tape_len, n_ldt_iter, N_MPC_ITER)
            } else {
                // verifier.zk_verify_w_vpd(&outputs_v, N_PARTY, tape_len, n_ldt_iter, &coset_v)
                verifier.zk_verify_w_vpdv2(&outputs_v, N_PARTY, tape_len, n_ldt_iter, &coset_v)
            };
        assert!(res);
    });

    v_thread.join().unwrap();
    p_thread.join().unwrap();
    lc.close();

    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .open(std::format!("log/zk_matmul_{}.log", type_name::<T>()))
        .unwrap();
    statistic(&mut file, log_len, n_ldt_iter, N_MPC_ITER);
}


fn test_zk_matmul<T: Field>() {
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
        .open(std::format!("log/zk_matmul_{}.log", type_name::<T>()))
        .unwrap();
    writeln!(file, "Field: {}", type_name::<T>()).unwrap();
    writeln!(file, "RS code rate: 1/{}", 1 << LOG_RHO).unwrap();
    writeln!(file, "Coset ratio: {}", 1 << LOG_LAYER_RATIO).unwrap();
    writeln!(file, "MPC num of parties: {}", N_PARTY).unwrap();

    for log_len in 4..11 {
        println!("matrix size: {0} x {0}", 1 << log_len);
        println!("log inputs: {0}", log_len * 2 + 1);
        zk_matmul_w_vpd::<T>(log_len);
    }
}

#[test]
fn test_zk_matmul_all() {
    println!("");
    // test_zk_matmul::<M61p2>();
    test_zk_matmul::<M127p2>();
    // test_zk_matmul::<P192>();
    // test_zk_matmul::<GF2p16>();
    // test_zk_matmul::<GF2p32>();
    // test_zk_matmul::<GF2p64>();
    // test_zk_matmul::<GF2p128>();
    // test_zk_matmul::<GF2p192>();
    // test_zk_matmul::<GF2p256>();
    // test_zk_matmul::<P256>();
}
