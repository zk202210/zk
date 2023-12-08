// Op Statistics
// Debug mode only
// #![cfg(debug_assertions)]

use std::sync::atomic::{AtomicUsize, Ordering};
use procfs::process::Process;
use std::io::Write;
use std::fs::File;

pub static MUL_COUNT: AtomicUsize = AtomicUsize::new(0);
pub static INV_COUNT: AtomicUsize = AtomicUsize::new(0);
pub static ADD_COUNT: AtomicUsize = AtomicUsize::new(0);
pub static SUB_COUNT: AtomicUsize = AtomicUsize::new(0);
pub static HASH_COUNT: AtomicUsize = AtomicUsize::new(0);


macro_rules! new_stat {
    ($U:tt, $get: tt, $inc: tt, $rst: tt) => {
        pub static $U: AtomicUsize = AtomicUsize::new(0);

        pub fn $get() -> usize {
            $U.load(Ordering::SeqCst)
        }

        pub fn $inc(x: usize) {
            $U.fetch_add(x, Ordering::SeqCst);
        }

        pub fn $rst() {
            $U.store(0, Ordering::SeqCst);
        }
    }
}

new_stat!(VPD_PRV, vpd_prv_get, vpd_prv_inc, vpd_prv_rst);
new_stat!(GKR_PRV, gkr_prv_get, gkr_prv_inc, gkr_prv_rst);
new_stat!(MPC_PRV, mpc_prv_get, mpc_prv_inc, mpc_prv_rst);
new_stat!(FLP_PRV, flp_prv_get, flp_prv_inc, flp_prv_rst);
new_stat!(CHANNEL0, channel0_get, channel0_inc, channel0_rst);
new_stat!(CHANNEL1, channel1_get, channel1_inc, channel1_rst);
new_stat!(VPD_VRF, vpd_vrf_get, vpd_vrf_inc, vpd_vrf_rst);
new_stat!(GKR_VRF, gkr_vrf_get, gkr_vrf_inc, gkr_vrf_rst);
new_stat!(ZKGKR_VRF, zkgkr_vrf_get, zkgkr_vrf_inc, zkgkr_vrf_rst);
new_stat!(MPC_VRF, mpc_vrf_get, mpc_vrf_inc, mpc_vrf_rst);
new_stat!(VPD_BW, vpd_bw_get, vpd_bw_inc, vpd_bw_rst);
new_stat!(GKR_BW, gkr_bw_get, gkr_bw_inc, gkr_bw_rst);
new_stat!(MPC_BW, mpc_bw_get, mpc_bw_inc, mpc_bw_rst);
new_stat!(VPD_BW_HASH, vpd_bw_hash_get, vpd_bw_hash_inc, vpd_bw_hash_rst);
new_stat!(GKR_BW_HASH, gkr_bw_hash_get, gkr_bw_hash_inc, gkr_bw_hash_rst);
new_stat!(MPC_BW_HASH, mpc_bw_hash_get, mpc_bw_hash_inc, mpc_bw_hash_rst);
new_stat!(DEBUG_FLAG, debug_flag_get, debug_flag_inc, debug_flag_rst);
new_stat!(VIRGO_FLAG, virgo_flag_get, virgo_flag_inc, virgo_flag_rst);

pub fn stats_reset_all() {
    vpd_prv_rst();
    gkr_prv_rst();
    mpc_prv_rst();
    flp_prv_rst();
    channel0_rst();
    channel1_rst();
    vpd_vrf_rst();
    gkr_vrf_rst();
    zkgkr_vrf_rst();
    mpc_vrf_rst();
    vpd_bw_rst();
    gkr_bw_rst();
    mpc_bw_rst();
    debug_flag_rst();
    virgo_flag_rst();
    mpc_bw_hash_rst();
    gkr_bw_hash_rst();
    vpd_bw_hash_rst();
}

pub fn add_count() -> usize {
    ADD_COUNT.load(Ordering::SeqCst)
}

pub fn sub_count() -> usize {
    SUB_COUNT.load(Ordering::SeqCst)
}

pub fn mul_count() -> usize {
    MUL_COUNT.load(Ordering::SeqCst)
}

pub fn inv_count() -> usize {
    INV_COUNT.load(Ordering::SeqCst)
}

pub fn hash_count() -> usize {
    HASH_COUNT.load(Ordering::SeqCst)
}


pub fn mul_call() {
    MUL_COUNT.fetch_add(1, Ordering::SeqCst);
}

pub fn inv_call() {
    INV_COUNT.fetch_add(1, Ordering::SeqCst);
}

pub fn add_call() {
    ADD_COUNT.fetch_add(1, Ordering::SeqCst);
}

pub fn sub_call() {
    SUB_COUNT.fetch_add(1, Ordering::SeqCst);
}

pub fn hash_call() {
    HASH_COUNT.fetch_add(1, Ordering::SeqCst);
}


pub fn mul_reset() {
    MUL_COUNT.store(0, Ordering::SeqCst);
}

pub fn inv_reset() {
    INV_COUNT.store(0, Ordering::SeqCst);
}

pub fn add_reset() {
    ADD_COUNT.store(0, Ordering::SeqCst);
}

pub fn sub_reset() {
    SUB_COUNT.store(0, Ordering::SeqCst);
}

pub fn hash_reset() {
    HASH_COUNT.store(0, Ordering::SeqCst);
}

pub fn procinfo() -> String {
    let p = Process::myself().unwrap();
    let info = p.statm().unwrap();
    format!("VSZ: {:.3} GB, RSS {:.3} GB",
            (info.size as f64) * 4.0 / 1024.0 / 1024.0,
            (info.data as f64) * 4.0 / 1024.0 / 1024.0)
}

pub fn statistic(file: &mut File, log_len: usize, n_ldt_iter: usize, n_mpc_iter: usize) {
    let vpd_prv = vpd_prv_get() as f64 / 1000.0;
    let gkr_prv = gkr_prv_get() as f64 / 1000.0;
    let mpc_prv = mpc_prv_get() as f64 / 1000.0;
    let vpd_vrf = vpd_vrf_get() as f64 / 1000.0;
    let gkr_vrf = gkr_vrf_get() as f64 / 1000.0;
    let zkgkr_vrf = zkgkr_vrf_get() as f64 / 1000.0;
    let mpc_vrf = mpc_vrf_get() as f64 / 1000.0;
    let pi_vpd = (vpd_bw_get() + gkr_bw_get()) as f64 / 1024.0;
    let pi_mpc = mpc_bw_get() as f64 / 1024.0;
    let hash_vpd = vpd_bw_hash_get() + gkr_bw_hash_get();
    let hash_mpc = mpc_bw_hash_get();
    println!("Prv vpd            {:>11.3} ms", vpd_prv);
    println!("Prv gkr            {:>11.3} ms", gkr_prv);
    println!("Prv mpc (per iter) {:>11.3} ms", mpc_prv);
    println!("");
    println!("Vrf vpd w/o gkr    {:>11.3} ms", vpd_vrf);
    println!("Vrf vpd gkr        {:>11.3} ms", gkr_vrf);
    println!("Vrf zkgkr          {:>11.3} ms", zkgkr_vrf);
    println!("Vrf mpc (per iter) {:>11.3} ms", mpc_vrf);
    println!("");
    println!("Prf vpd w/o gkr    {:>11.2} kB", vpd_bw_get() as f64 / 1024.0);
    println!("Prf vpd gkr        {:>11.2} kB", gkr_bw_get() as f64 / 1024.0);
    println!("Prf mpc (per iter) {:>11.2} kB", pi_mpc);
    println!("Prf vpd hash            {:>11}", hash_vpd);
    println!("Prf mpc hash (per iter) {:>11}", hash_mpc);
    println!("--------------------------------------------------------------------------");

    let n_mpc_iter = n_mpc_iter as f64;
    let prv = vpd_prv + gkr_prv + mpc_prv * n_mpc_iter;
    let vrf = vpd_vrf + gkr_vrf + zkgkr_vrf + mpc_vrf * n_mpc_iter;
    let prf = pi_vpd + pi_mpc * n_mpc_iter;
    writeln!(file, "{:>2}\t{:>2}\t{:>11.3}\t{:>11.3}\t{:>11.3}\t{:>11.3}\t{:>11.3}\t{:>11.3}\t{:>11.3}\t{:>11.3}\t{:>11.3}\t{:>11.3}\t{:>11.2}\t{:>11.2}",
             log_len,
             n_ldt_iter,
             prv,
             vrf,
             prf,
             vpd_prv,
             gkr_prv,
             mpc_prv,
             vpd_vrf,
             gkr_vrf,
             zkgkr_vrf,
             mpc_vrf,
             pi_vpd,
             pi_mpc).unwrap()
}
