//! Parameter sets for the Regev (Ring-LWE) scheme.
//!
//! The ciphertext modulus `q` is **fixed to the proof-system field prime**
//! (BabyBear: `q = 2^31 - 2^27 + 1 = 2013265921`). This is the key design
//! decision that makes the encryption proof cheap: every `mod q` operation in
//! the scheme is a native field operation in the STARK.
//!
//! # Security
//!
//! Ring-LWE hardness is governed by the ring dimension `n`, the modulus `q`,
//! and the noise width. With `q ≈ 2^31` and centered-binomial noise
//! `CBD(η=2)` (σ = 1), `n = 1024` lands at roughly 2^100 classical
//! core-SVP security — usable, but below the 128-bit target. `n = 2048`
//! gives a comfortable margin. **Run the lattice-estimator on your final
//! parameter choice** (see README) before deploying.

use p3_baby_bear::BabyBear;
use p3_field::PrimeField32;

/// Parameters for the Regev encryption scheme over `Z_q[x]/(x^n + 1)`,
/// with `q` equal to the BabyBear prime.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RegevParams {
    /// Ring dimension; must be a power of two (it is also the STARK trace
    /// height, so it must satisfy any minimum imposed by the FRI config).
    pub n: usize,
    /// Centered binomial noise parameter; the implementation is specialised
    /// to η = 2 (noise in `[-2, 2]`, decomposed as `u - v` with
    /// `u, v ∈ {0, 1, 2}` so smallness is a degree-3 constraint).
    pub eta: usize,
}

impl RegevParams {
    /// `n = 1024` as in the design sketch. Re-estimate security before use;
    /// roughly 2^100 classical with ternary secrets and CBD(2) noise.
    pub const N1024: Self = Self { n: 1024, eta: 2 };

    /// `n = 2048`, comfortable 128-bit-plus margin at `q ≈ 2^31`.
    pub const N2048: Self = Self { n: 2048, eta: 2 };

    /// The ciphertext modulus `q` (the BabyBear prime).
    pub const fn q() -> u32 {
        BabyBear::ORDER_U32
    }

    /// Message scaling factor `Δ = floor(q / 2)` for binary messages.
    pub const fn delta() -> u32 {
        BabyBear::ORDER_U32 / 2
    }

    /// log2 of the ring dimension.
    pub fn log_n(&self) -> usize {
        debug_assert!(self.n.is_power_of_two());
        self.n.trailing_zeros() as usize
    }
}
