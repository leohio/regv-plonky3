//! Radix-2 NTT over a two-adic prime field, used for fast polynomial
//! multiplication in `Z_q[x]/(x^n + 1)` where `q` is the proof-system field
//! prime (BabyBear by default).
//!
//! Because the ciphertext modulus equals the STARK field modulus, all ring
//! arithmetic here is *native* field arithmetic — no limb decomposition or
//! modular-reduction gadgets are ever needed, in or out of the circuit.

use p3_field::TwoAdicField;

/// In-place iterative Cooley-Tukey NTT (decimation in time).
///
/// `values.len()` must be a power of two and at most `2^F::TWO_ADICITY`.
/// If `inverse` is true, computes the inverse transform (including the `1/N`
/// scaling).
pub fn ntt<F: TwoAdicField>(values: &mut [F], inverse: bool) {
    let n = values.len();
    assert!(n.is_power_of_two(), "NTT size must be a power of two");
    let log_n = n.trailing_zeros() as usize;
    assert!(log_n <= F::TWO_ADICITY, "NTT size exceeds field two-adicity");

    // Bit-reversal permutation.
    for i in 0..n {
        let j = i.reverse_bits() >> (usize::BITS as usize - log_n);
        if i < j {
            values.swap(i, j);
        }
    }

    let root = if inverse {
        F::two_adic_generator(log_n).inverse()
    } else {
        F::two_adic_generator(log_n)
    };

    let mut len = 2;
    while len <= n {
        // w_len is a primitive `len`-th root of unity.
        let w_len = root.exp_u64((n / len) as u64);
        for start in (0..n).step_by(len) {
            let mut w = F::ONE;
            for k in start..start + len / 2 {
                let u = values[k];
                let v = values[k + len / 2] * w;
                values[k] = u + v;
                values[k + len / 2] = u - v;
                w *= w_len;
            }
        }
        len <<= 1;
    }

    if inverse {
        let n_inv = F::from_usize(n).inverse();
        for v in values.iter_mut() {
            *v *= n_inv;
        }
    }
}

/// Full product of two polynomials with coefficients `a` and `b` (degree
/// `< n` each), returned as `2n` coefficients (the top coefficient is zero).
///
/// Uses a size-`2n` cyclic NTT on zero-padded inputs, which is exact because
/// `deg(a*b) <= 2n - 2 < 2n`.
pub fn full_poly_mul<F: TwoAdicField>(a: &[F], b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    let size = 2 * n;

    let mut fa = vec![F::ZERO; size];
    let mut fb = vec![F::ZERO; size];
    fa[..n].copy_from_slice(a);
    fb[..n].copy_from_slice(b);

    ntt(&mut fa, false);
    ntt(&mut fb, false);
    for (x, y) in fa.iter_mut().zip(fb.iter()) {
        *x *= *y;
    }
    ntt(&mut fa, true);
    fa
}

/// Schoolbook negacyclic product, for testing the NTT path.
#[cfg(test)]
pub fn negacyclic_mul_naive<F: p3_field::Field>(a: &[F], b: &[F]) -> Vec<F> {
    let n = a.len();
    let mut out = vec![F::ZERO; n];
    for i in 0..n {
        for j in 0..n {
            let prod = a[i] * b[j];
            if i + j < n {
                out[i + j] += prod;
            } else {
                out[i + j - n] -= prod;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use p3_baby_bear::BabyBear;
    use p3_field::PrimeCharacteristicRing;
    use rand::rngs::SmallRng;
    use rand::{RngExt, SeedableRng};

    #[test]
    fn ntt_roundtrip() {
        let mut rng = SmallRng::seed_from_u64(0);
        let original: Vec<BabyBear> = (0..256).map(|_| rng.random()).collect();
        let mut values = original.clone();
        ntt(&mut values, false);
        ntt(&mut values, true);
        assert_eq!(values, original);
    }

    #[test]
    fn full_mul_matches_naive_negacyclic() {
        let mut rng = SmallRng::seed_from_u64(1);
        let n = 64;
        let a: Vec<BabyBear> = (0..n).map(|_| rng.random()).collect();
        let b: Vec<BabyBear> = (0..n).map(|_| rng.random()).collect();
        let full = full_poly_mul(&a, &b);
        let expected = negacyclic_mul_naive(&a, &b);
        let reduced: Vec<BabyBear> = (0..n).map(|i| full[i] - full[n + i]).collect();
        assert_eq!(reduced, expected);
        assert_eq!(full[2 * n - 1], BabyBear::ZERO);
    }
}
