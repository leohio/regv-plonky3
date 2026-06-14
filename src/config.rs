//! Concrete Plonky3 configurations for the Regev encryption proof.
//!
//! Default: BabyBear + Poseidon2 transcript + two-adic FRI. BabyBear is both
//! the ciphertext modulus and the proof field, and (unlike Mersenne31) is
//! NTT-friendly, so the *encryptor* also gets fast negacyclic multiplication
//! natively.

use p3_baby_bear::{default_babybear_poseidon2_16, BabyBear, Poseidon2BabyBear};
use p3_challenger::{DuplexChallenger, HashChallenger, SerializingChallenger32};
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::Field;
use p3_fri::{FriParameters, HidingFriPcs, TwoAdicFriPcs};
use p3_keccak::{Keccak256Hash, KeccakF};
use p3_merkle_tree::{MerkleTreeHidingMmcs, MerkleTreeMmcs};
use p3_symmetric::{
    CompressionFunctionFromHasher, PaddingFreeSponge, SerializingHasher, TruncatedPermutation,
};
use p3_uni_stark::StarkConfig;
use rand::rngs::SmallRng;
use rand::SeedableRng;

pub type Val = BabyBear;
pub type Challenge = BinomialExtensionField<Val, 4>;

pub type Perm = Poseidon2BabyBear<16>;

/// The canonical Poseidon2-BabyBear permutation used everywhere in this crate
/// (transcript, Merkle hashing). It uses Plonky3's published compile-time
/// round constants, **not** an RNG, so it is byte-identical on every target.
///
/// Building the permutation from `rand::rngs::SmallRng` instead would fork the
/// transcript across pointer widths (32-bit wasm vs 64-bit native); see
/// [`crate::portability`].
pub fn canonical_perm() -> Perm {
    default_babybear_poseidon2_16()
}
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

// --- Zero-knowledge variant (hiding FRI) ---------------------------------
//
// The witness (encryption randomness, noise, message) is secret, so for
// deployments where the proof itself must not leak trace information, use
// this config: `HidingFriPcs` masks committed polynomials with randomness
// and adds a random codeword to the FRI batch (statistical zk).

type ByteHash = Keccak256Hash;
type U64Hash = PaddingFreeSponge<KeccakF, 25, 17, 4>;
type ZkFieldHash = SerializingHasher<U64Hash>;
type ZkCompress = CompressionFunctionFromHasher<U64Hash, 2, 4>;
type ZkValMmcs = MerkleTreeHidingMmcs<
    [Val; p3_keccak::VECTOR_LEN],
    [u64; p3_keccak::VECTOR_LEN],
    ZkFieldHash,
    ZkCompress,
    SmallRng,
    2,
    4,
    4,
>;
type ZkChallengeMmcs = ExtensionMmcs<Val, Challenge, ZkValMmcs>;
pub type ZkChallenger = SerializingChallenger32<Val, HashChallenger<u8, ByteHash, 32>>;
pub type ZkPcs = HidingFriPcs<Val, Dft, ZkValMmcs, ZkChallengeMmcs, SmallRng>;
pub type RegevZkStarkConfig = StarkConfig<ZkPcs, Challenge, ZkChallenger>;

/// Zero-knowledge config (hiding commitments + randomized FRI batch).
/// Roughly 2x the proving cost of [`default_config`].
pub fn zk_config() -> RegevZkStarkConfig {
    zk_config_seeded(SmallRng::from_rng(&mut rand::rng()))
}

/// ZK config with caller-supplied masking randomness (for reproducible
/// tests; use [`zk_config`] in production).
pub fn zk_config_seeded(rng: SmallRng) -> RegevZkStarkConfig {
    let u64_hash = U64Hash::new(KeccakF {});
    let field_hash = ZkFieldHash::new(u64_hash);
    let compress = ZkCompress::new(u64_hash);
    let val_mmcs = ZkValMmcs::new(field_hash, compress, 0, rng.clone());
    let challenge_mmcs = ZkChallengeMmcs::new(val_mmcs.clone());
    let fri_params = FriParameters {
        log_blowup: 2,
        log_final_poly_len: 0,
        max_log_arity: 1,
        num_queries: 42,
        commit_proof_of_work_bits: 8,
        query_proof_of_work_bits: 16,
        mmcs: challenge_mmcs,
    };
    let pcs = ZkPcs::new(Dft::default(), val_mmcs, fri_params, 4, rng);
    let challenger = ZkChallenger::from_hasher(vec![], ByteHash {});
    RegevZkStarkConfig::new(pcs, challenger)
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
    // Poseidon2 permutation drives the Fiat-Shamir transcript, so its round
    // constants MUST be identical on every target. We use Plonky3's canonical
    // published BabyBear constants (`default_babybear_poseidon2_16`): they are
    // compile-time tables, not derived from any RNG.
    //
    // NOTE: do NOT build the permutation from `rand::rngs::SmallRng` here.
    // `SmallRng` is pointer-width-dependent (Xoshiro128++ on 32-bit wasm32,
    // Xoshiro256++ on 64-bit native), so a fixed seed yields *different* round
    // constants per target — which silently forks the transcript and makes a
    // wasm-generated proof fail native verification (FRI `InvalidPowWitness`)
    // and vice-versa. See `portability::transcript_digest` for the regression
    // guard and `tests/portability.rs` for the golden cross-target values.
    let perm = canonical_perm();
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
