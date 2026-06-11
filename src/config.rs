//! Concrete Plonky3 configurations for the Regev encryption proof.
//!
//! Default: BabyBear + Poseidon2 transcript + two-adic FRI. BabyBear is both
//! the ciphertext modulus and the proof field, and (unlike Mersenne31) is
//! NTT-friendly, so the *encryptor* also gets fast negacyclic multiplication
//! natively.

use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::Field;
use p3_fri::{FriParameters, TwoAdicFriPcs};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::StarkConfig;
use rand::rngs::SmallRng;
use rand::SeedableRng;

pub type Val = BabyBear;
pub type Challenge = BinomialExtensionField<Val, 4>;

type Perm = Poseidon2BabyBear<16>;
type Hash = PaddingFreeSponge<Perm, 16, 8, 8>;
type Compress = TruncatedPermutation<Perm, 2, 8, 16>;
type ValMmcs =
    MerkleTreeMmcs<<Val as Field>::Packing, <Val as Field>::Packing, Hash, Compress, 2, 8>;
type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
pub type Challenger = DuplexChallenger<Val, Perm, 16, 8>;
type Dft = Radix2DitParallel<Val>;
pub type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
pub type RegevStarkConfig = StarkConfig<Pcs, Challenge, Challenger>;

/// Production-leaning config: blowup 2, 84 queries, 16-bit grinding
/// (~100-bit conjectured security with the quartic extension).
pub fn default_config() -> RegevStarkConfig {
    make_config(FriParametersSpec {
        log_blowup: 1,
        log_final_poly_len: 0,
        max_log_arity: 4,
        num_queries: 84,
        commit_proof_of_work_bits: 8,
        query_proof_of_work_bits: 16,
    })
}

/// Cheap config for tests: few queries, tiny grinding. NOT secure.
pub fn test_config() -> RegevStarkConfig {
    make_config(FriParametersSpec {
        log_blowup: 1,
        log_final_poly_len: 0,
        max_log_arity: 1,
        num_queries: 8,
        commit_proof_of_work_bits: 1,
        query_proof_of_work_bits: 1,
    })
}

struct FriParametersSpec {
    log_blowup: usize,
    log_final_poly_len: usize,
    max_log_arity: usize,
    num_queries: usize,
    commit_proof_of_work_bits: usize,
    query_proof_of_work_bits: usize,
}

fn make_config(spec: FriParametersSpec) -> RegevStarkConfig {
    // Fixed-seed RNG: the Poseidon2 round constants are public parameters.
    let mut rng = SmallRng::seed_from_u64(0x5245_4745_56u64); // "REGEV"
    let perm = Perm::new_from_rng_128(&mut rng);
    let hash = Hash::new(perm.clone());
    let compress = Compress::new(perm.clone());
    let val_mmcs = ValMmcs::new(hash, compress, 0);
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
    let fri_params = FriParameters {
        log_blowup: spec.log_blowup,
        log_final_poly_len: spec.log_final_poly_len,
        max_log_arity: spec.max_log_arity,
        num_queries: spec.num_queries,
        commit_proof_of_work_bits: spec.commit_proof_of_work_bits,
        query_proof_of_work_bits: spec.query_proof_of_work_bits,
        mmcs: challenge_mmcs,
    };
    let pcs = Pcs::new(Dft::default(), val_mmcs, fri_params);
    let challenger = Challenger::new(perm);
    RegevStarkConfig::new(pcs, challenger)
}
