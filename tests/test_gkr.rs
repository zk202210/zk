#[cfg(test)]
mod tests {
    use std::thread;
    use rustzk::gkr::*;
    use rustzk::circuit::*;
    use rustzk::field::*;
    use rustzk::channel::*;

    #[test]
    fn plain_gkr() {
        // type T = GF2p16;
        type T = M61p2;
        let mut lc = LocalChannel::<T>::new();
        let mut lc0 = PLocalChannel::<T>::new(&mut lc);
        let mut lc1 = PLocalChannel::<T>::new(&mut lc);
        let mut c = Circuit::<T>::random(10, 4, 10);
        let mut cv = c.clone();

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

        let mut outputs: Vec<T> = vec![];
        let output_layer = c.layers.last().unwrap();
        for i in 0..(1 << (output_layer.bit_len)) {
            outputs.push(output_layer.values[i]);
        }

        let p_thread = thread::spawn(move || {
            let lc = &mut *lc0;
            let mut prover = Prover::<T, LocalChannel<T>>::new(&mut c, lc);
            prover.prove();
        });

        let v_thread = thread::spawn(move || {
            let lc = &mut *lc1;
            let mut verifier = Verifier::<T, LocalChannel<T>>::new(&mut cv, lc);
            let res = verifier.verify(&outputs);
            assert!(res);
        });

        v_thread.join().unwrap();
        p_thread.join().unwrap();
        lc.close();
    }
}
