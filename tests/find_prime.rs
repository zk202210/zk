use gmp_mpfr_sys::gmp;
use core::mem::MaybeUninit;

fn new_int() -> gmp::mpz_t {
    unsafe {
        let mut z = MaybeUninit::uninit();
        gmp::mpz_init(z.as_mut_ptr());
        z.assume_init()
    }
}


fn find_prime_with_subgroup(log_max: usize, log_size: usize) {
    unsafe {
        let mut a = new_int();
        let mut b = new_int();

        gmp::mpz_set_ui(&mut a, 1);
        gmp::mpz_mul_2exp(&mut a, &a, log_size as u64);
        gmp::mpz_set_ui(&mut b, 1);
        gmp::mpz_mul_2exp(&mut b, &b, (log_max - log_size) as u64);
        gmp::mpz_sub_ui(&mut b, &b, 1);
        let mut i = 1;
        loop {
            let mut p = new_int();
            gmp::mpz_set_ui(&mut p, 1);

            // p = 2^xx * b + 1
            gmp::mpz_addmul(&mut p, &a, &b);

            let test = gmp::mpz_probab_prime_p(&p, 100);

            if test != 0 {
                println!("possible {} prime: 2**{} - {:>3} * (2**{}) + 1", test, log_max, i, log_size);
                gmp::mpz_clear(&mut p);
                break;
            }
            gmp::mpz_clear(&mut p);

            i += 2;
            // b -= 2
            gmp::mpz_sub_ui(&mut b, &b, 2);
        }
        gmp::mpz_clear(&mut a);
        gmp::mpz_clear(&mut b);
    }
}

#[test]
fn find_prime() {
    println!("");
    for i in 30..65 {
        find_prime_with_subgroup(256, i);
    }
}
