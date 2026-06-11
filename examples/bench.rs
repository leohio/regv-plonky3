//! Timing harness: keygen/encrypt once, then prove and verify batches of
//! encryption proofs at production-like parameters.
//!
//! ```sh
//! cargo run --release --example bench            # n = 1024
//! cargo run --release --example bench -- 2048    # n = 2048
//! ```

use std::time::Instant;

use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use regev_plonky3::*;

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .map(|s| s.parse().expect("n must be a power of two"))
        .unwrap_or(1024);
    let params = RegevParams { n, eta: 2 };

    let mut rng = SmallRng::seed_from_u64(0xbeef);
    let (pk, sk) = keygen(&mut rng, &params);
    let config = default_config();

    println!(
        "Regev/plonky3: n = {n}, q = {} (BabyBear), CBD(2) noise, ternary r",
        RegevParams::q()
    );

    let t = Instant::now();
    let iters = 1000;
    for _ in 0..iters {
        let m: Vec<u8> = (0..n).map(|_| rng.random_range(0..=1)).collect();
        let _ = encrypt(&mut rng, &params, &pk, &m);
    }
    println!("encrypt: {:?}/op", t.elapsed() / iters);

    for batch in [1usize, 8, 32] {
        let mut cts = Vec::new();
        let mut wits = Vec::new();
        let mut msgs = Vec::new();
        for _ in 0..batch {
            let m: Vec<u8> = (0..n).map(|_| rng.random_range(0..=1)).collect();
            let (ct, w) = encrypt(&mut rng, &params, &pk, &m);
            msgs.push(m);
            cts.push(ct);
            wits.push(w);
        }

        let t = Instant::now();
        let proof = prove_encryptions(&config, &params, &pk, &cts, &wits);
        let prove_time = t.elapsed();

        let t = Instant::now();
        verify_encryptions(&config, &params, &pk, &cts, &proof).expect("verify");
        let verify_time = t.elapsed();

        let size = postcard::to_allocvec(&proof).unwrap().len();
        println!(
            "batch {batch:>3}: prove {prove_time:>9.2?}  ({:>9.2?}/ct)   verify {verify_time:>9.2?}   proof {:>8} bytes ({} B/ct)",
            prove_time / batch as u32,
            size,
            size / batch
        );

        for (m, ct) in msgs.iter().zip(&cts) {
            assert_eq!(&decrypt(&params, &sk, ct), m);
        }
    }

    // Zero-knowledge variant.
    {
        let zk = zk_config();
        let batch = 8usize;
        let mut cts = Vec::new();
        let mut wits = Vec::new();
        for _ in 0..batch {
            let m: Vec<u8> = (0..n).map(|_| rng.random_range(0..=1)).collect();
            let (ct, w) = encrypt(&mut rng, &params, &pk, &m);
            cts.push(ct);
            wits.push(w);
        }
        let t = Instant::now();
        let proof = prove_encryptions(&zk, &params, &pk, &cts, &wits);
        let prove_time = t.elapsed();
        let t = Instant::now();
        verify_encryptions(&zk, &params, &pk, &cts, &proof).expect("zk verify");
        let verify_time = t.elapsed();
        let size = postcard::to_allocvec(&proof).unwrap().len();
        println!(
            "zk batch {batch}: prove {prove_time:>9.2?}  ({:>9.2?}/ct)   verify {verify_time:>9.2?}   proof {:>8} bytes ({} B/ct)",
            prove_time / batch as u32,
            size,
            size / batch
        );
    }

    // Encryption proof bundled with a 32-bit plaintext range proof
    // (value in [0, 1_000_000)).
    {
        let spec = RangeSpec {
            value_bits: 20,
            bound: 1_000_000,
        };
        let batch = 8usize;
        let mut cts = Vec::new();
        let mut wits = Vec::new();
        for k in 0..batch {
            // encode a value in the low `value_bits` bits of the message
            let value = (k as u64 * 111_111) % spec.bound;
            let mut m = vec![0u8; n];
            for (i, slot) in m.iter_mut().enumerate().take(spec.value_bits) {
                *slot = ((value >> i) & 1) as u8;
            }
            let (ct, w) = encrypt(&mut rng, &params, &pk, &m);
            cts.push(ct);
            wits.push(w);
        }
        let t = Instant::now();
        let proof = prove_encryptions_with_range(&config, &params, &pk, &cts, &wits, spec);
        let prove_time = t.elapsed();
        let t = Instant::now();
        verify_encryptions_with_range(&config, &params, &pk, &cts, &proof, spec)
            .expect("range verify");
        let verify_time = t.elapsed();
        let size = postcard::to_allocvec(&proof).unwrap().len();
        println!(
            "range b {batch}: prove {prove_time:>9.2?}  ({:>9.2?}/ct)   verify {verify_time:>9.2?}   proof {:>8} bytes ({} B/ct)",
            prove_time / batch as u32,
            size,
            size / batch
        );
    }
}
