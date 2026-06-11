//! # regev-plonky3
//!
//! Plonky3-friendly Regev (Ring-LWE) encryption with STARK proofs of correct
//! encryption, built around three design decisions:
//!
//! 1. **Native modulus.** The ciphertext modulus `q` *is* the proof-field
//!    prime (BabyBear, `q = 2^31 - 2^27 + 1`). Every ring operation is a
//!    native field operation; no limb decomposition, no reduction gadgets.
//!    LWE tolerates a large `q` — security comes from the dimension and the
//!    noise — and BabyBear's two-adicity also gives the encryptor a fast
//!    negacyclic NTT for free.
//!
//! 2. **Linear-size multiplication proof.** The products `a·r`, `b·r` in
//!    `Z_q[x]/(x^n+1)` are never computed in-circuit. Instead, a random
//!    point `z` is sampled (Fiat-Shamir) *after* the witness is committed,
//!    every polynomial is evaluated at `z` with an O(n) Horner running-sum
//!    column, and the identities
//!    `c1(z) = a(z)·r(z) + e1(z) − (z^n+1)·k1(z)`,
//!    `c2(z) = b(z)·r(z) + e2(z) + Δ·m(z) − (z^n+1)·k2(z)`
//!    are enforced at that single point (Schwartz-Zippel). The two-phase
//!    commitment uses the same multi-round challenge plumbing that Plonky3's
//!    logUp lookups use ([`gadget::EvalGadget`] implements the evaluation
//!    argument as a `LookupGadget`).
//!
//! 3. **Cheap smallness checks.** Ternary `r` and CBD(η=2) noise are
//!    enforced by degree-3 vanishing constraints (`x(x-1)(x+1)`,
//!    `x(x-1)(x-2)` on the CBD halves) — for ranges this small that is
//!    strictly cheaper than a logUp lookup, which would cost an extension
//!    field running-sum column per check.
//!
//! Batching: each ciphertext is one `p3-batch-stark` instance; the whole
//! batch shares one commitment per phase and a single FRI opening, so fixed
//! costs amortize across transactions.
//!
//! ## Privacy
//!
//! Only evaluations of *public* polynomials (`a, b, c1, c2`) are published.
//! The witness evaluations `r(z), e1(z), e2(z), m(z), k1(z), k2(z)` stay in
//! the (committed) permutation trace and enter the ring identities as
//! in-circuit constraints — publishing them would leak ~128 bits of
//! information per polynomial and allow dictionary attacks on low-entropy
//! messages. The default config is succinct but not zero-knowledge (plain
//! FRI); use [`config::zk_config`] (hiding FRI commitments + randomized
//! quotient chunks) when the proof itself must leak nothing about the
//! witness.
//!
//! ## Quick start
//!
//! ```rust
//! use regev_plonky3::*;
//! use rand::{rngs::SmallRng, RngExt, SeedableRng};
//!
//! let params = RegevParams { n: 256, eta: 2 }; // use N1024/N2048 in production
//! let mut rng = SmallRng::seed_from_u64(42);
//! let (pk, sk) = keygen(&mut rng, &params);
//!
//! let m: Vec<u8> = (0..params.n).map(|_| rng.random_range(0..=1)).collect();
//! let (ct, witness) = encrypt(&mut rng, &params, &pk, &m);
//! assert_eq!(decrypt(&params, &sk, &ct), m);
//!
//! let config = test_config();
//! let proof = prove_encryptions(&config, &params, &pk, &[ct.clone()], &[witness]);
//! verify_encryptions(&config, &params, &pk, &[ct], &proof).unwrap();
//! ```

extern crate alloc;

pub mod air;
pub mod config;
pub mod gadget;
pub mod ntt;
pub mod params;
pub mod prove;
pub mod regev;
pub mod stark;

pub use config::{default_config, test_config, zk_config, Challenge, RegevStarkConfig, RegevZkStarkConfig, Val};
pub use params::RegevParams;
pub use prove::{prove_encryptions, verify_encryptions, RegevProof, RegevVerifyError};
pub use regev::{
    decrypt, encrypt, keygen, Ciphertext, EncryptionWitness, PublicKey, SecretKey,
};
