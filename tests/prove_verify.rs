//! End-to-end tests: encrypt → prove → verify, plus soundness checks that
//! tampered statements and tampered proofs are rejected.

use p3_field::PrimeCharacteristicRing;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use regev_plonky3::*;

const TEST_PARAMS: RegevParams = RegevParams { n: 256, eta: 2 };

fn setup(
    seed: u64,
    batch: usize,
) -> (
    PublicKey,
    SecretKey,
    Vec<Vec<u8>>,
    Vec<Ciphertext>,
    Vec<EncryptionWitness>,
) {
    let mut rng = SmallRng::seed_from_u64(seed);
    let (pk, sk) = keygen(&mut rng, &TEST_PARAMS);
    let mut msgs = Vec::new();
    let mut cts = Vec::new();
    let mut wits = Vec::new();
    for _ in 0..batch {
        let m: Vec<u8> = (0..TEST_PARAMS.n).map(|_| rng.random_range(0..=1)).collect();
        let (ct, w) = encrypt(&mut rng, &TEST_PARAMS, &pk, &m);
        msgs.push(m);
        cts.push(ct);
        wits.push(w);
    }
    (pk, sk, msgs, cts, wits)
}

#[test]
fn prove_verify_single() {
    let (pk, sk, msgs, cts, wits) = setup(1, 1);
    assert_eq!(decrypt(&TEST_PARAMS, &sk, &cts[0]), msgs[0]);

    let config = test_config();
    let proof = prove_encryptions(&config, &TEST_PARAMS, &pk, &cts, &wits);
    verify_encryptions(&config, &TEST_PARAMS, &pk, &cts, &proof).expect("honest proof verifies");
}

#[test]
fn prove_verify_batch() {
    let (pk, _, _, cts, wits) = setup(2, 4);
    let config = test_config();
    let proof = prove_encryptions(&config, &TEST_PARAMS, &pk, &cts, &wits);
    verify_encryptions(&config, &TEST_PARAMS, &pk, &cts, &proof).expect("batch proof verifies");

    // Proof size for the curious (postcard, unoptimized).
    let bytes = postcard::to_allocvec(&proof).unwrap();
    println!("batch of 4: proof size = {} bytes", bytes.len());
}

#[test]
fn rejects_wrong_ciphertext() {
    let (pk, _, _, cts, wits) = setup(3, 1);
    let config = test_config();
    let proof = prove_encryptions(&config, &TEST_PARAMS, &pk, &cts, &wits);

    // Same proof, different claimed ciphertext. The ciphertext is part of
    // the public values, so the transcript already diverges and the inner
    // STARK rejects (before the evaluation-binding check even runs).
    let mut other = cts.clone();
    other[0].c2[7] += Val::ONE;
    assert!(verify_encryptions(&config, &TEST_PARAMS, &pk, &other, &proof).is_err());
}

/// The critical soundness path for the evaluation binding: a malicious
/// prover commits trace columns for a *different* (validly formed)
/// ciphertext while presenting the real one as the statement. The inner
/// STARK is fully self-consistent in that case — only the comparison of the
/// published evaluations `c1(z), c2(z)` against the claimed statement
/// catches it.
#[test]
fn rejects_malicious_prover_with_forged_columns() {
    use p3_batch_stark::common::ProverData;
    use p3_batch_stark::prover::StarkInstance;
    
    use p3_lookup::LookupAir;
    use regev_plonky3::air::{generate_trace, public_values, RegevEncAir};

    let (pk, _, _, cts, _) = setup(7, 1);
    let mut rng = SmallRng::seed_from_u64(77);

    // A forged-but-valid encryption of a different message.
    let m2: Vec<u8> = (0..TEST_PARAMS.n).map(|_| rng.random_range(0..=1)).collect();
    let (forged_ct, forged_wit) = encrypt(&mut rng, &TEST_PARAMS, &pk, &m2);

    let config = test_config();
    let air = RegevEncAir::new(
        TEST_PARAMS.n,
        Val::from_u32(RegevParams::delta()),
    );
    // Trace for the forged ciphertext, public values for the real one.
    let trace = generate_trace(&pk, &forged_ct, &forged_wit);
    let instances = vec![StarkInstance {
        air: &air,
        trace: &trace,
        public_values: public_values(&pk, &cts[0]),
        lookups: air.clone().get_lookups(),
    }];
    let prover_data = ProverData::from_instances(&config, &instances);
    let proof = regev_plonky3::stark::prove_batch(&config, &instances, &prover_data);

    let err = verify_encryptions(&config, &TEST_PARAMS, &pk, &cts, &proof).unwrap_err();
    match err {
        RegevVerifyError::StatementMismatch { poly, .. } => {
            assert!(poly == "c1" || poly == "c2");
        }
        e => panic!("expected StatementMismatch, got {e}"),
    }
}

#[test]
fn rejects_wrong_public_key() {
    let (pk, _, _, cts, wits) = setup(4, 1);
    let config = test_config();
    let proof = prove_encryptions(&config, &TEST_PARAMS, &pk, &cts, &wits);

    let mut other_pk = pk.clone();
    other_pk.b[0] += Val::ONE;
    // Changing the public key changes the transcript (public values), so the
    // inner STARK itself must fail.
    assert!(verify_encryptions(&config, &TEST_PARAMS, &other_pk, &cts, &proof).is_err());
}

#[test]
fn rejects_tampered_published_evaluation() {
    let (pk, _, _, cts, wits) = setup(5, 1);
    let config = test_config();
    let mut proof = prove_encryptions(&config, &TEST_PARAMS, &pk, &cts, &wits);

    // Tamper with the published evaluation of c1. Both the transcript and the
    // first-row constraint depend on it, so verification must fail.
    proof.global_lookup_data[0][2].expected_cumulated += Challenge::ONE;
    assert!(verify_encryptions(&config, &TEST_PARAMS, &pk, &cts, &proof).is_err());
}

#[test]
fn rejects_swapped_instances() {
    // A batch proof for (ct0, ct1) must not verify as (ct1, ct0).
    let (pk, _, _, cts, wits) = setup(6, 2);
    let config = test_config();
    let proof = prove_encryptions(&config, &TEST_PARAMS, &pk, &cts, &wits);

    let swapped = vec![cts[1].clone(), cts[0].clone()];
    assert!(verify_encryptions(&config, &TEST_PARAMS, &pk, &swapped, &proof).is_err());
}

#[test]
fn prove_verify_zero_knowledge_config() {
    use rand::rngs::SmallRng as ZkRng;
    let (pk, _, _, cts, wits) = setup(8, 2);
    // Same masking seed on both sides is NOT required — the verifier never
    // uses the prover's rng — but the config constructor takes one.
    let config = regev_plonky3::config::zk_config_seeded(ZkRng::seed_from_u64(99));
    let proof = prove_encryptions(&config, &TEST_PARAMS, &pk, &cts, &wits);
    assert!(proof.commitments.random.is_some(), "zk proof carries random commitment");
    verify_encryptions(&config, &TEST_PARAMS, &pk, &cts, &proof).expect("zk proof verifies");
}
