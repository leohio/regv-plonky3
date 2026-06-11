//! The AIR proving that a Regev ciphertext is well formed:
//!
//! ```text
//! c1 = a·r + e1                 in Z_q[x]/(x^n + 1)
//! c2 = b·r + e2 + Δ·m           in Z_q[x]/(x^n + 1)
//! r  ternary, e1, e2 ∈ CBD(2), m binary
//! ```
//!
//! One trace = one ciphertext, with one row per polynomial coefficient
//! (height = ring dimension `n`). Multiple ciphertexts batch as multiple
//! instances of `p3-batch-stark`, sharing commitments and FRI.
//!
//! # Layout
//!
//! Main (base field) columns:
//!
//! | col | meaning                                  | constraint            |
//! |-----|------------------------------------------|-----------------------|
//! | a, b, c1, c2 | the public statement            | eval exposed at `z`   |
//! | r            | encryption randomness           | `r(r-1)(r+1) = 0`     |
//! | e1u, e1v     | CBD halves of `e1 = e1u - e1v`  | `x(x-1)(x-2) = 0`     |
//! | e2u, e2v     | CBD halves of `e2`              | `x(x-1)(x-2) = 0`     |
//! | m            | message bits                    | `m(m-1) = 0`          |
//! | k1, k2       | quotients by `x^n + 1`          | free                  |
//!
//! Smallness is enforced by these degree-≤3 vanishing constraints — for
//! ranges this tiny they are strictly cheaper than logUp lookups (zero
//! auxiliary columns and no batch inversions; a lookup would cost one
//! extension-field running-sum column *per* range check).
//!
//! Permutation (extension field) columns: one Horner running-evaluation
//! column per polynomial (see [`crate::gadget`]), at the challenge `z`
//! sampled after the main trace is committed.
//!
//! The multiplication is verified at the single random point `z`
//! (Schwartz-Zippel), via two first-row constraints over the running
//! evaluations:
//!
//! ```text
//! A·R + E1          − C1 − (z^n + 1)·K1 = 0
//! B·R + E2 + Δ·M    − C2 − (z^n + 1)·K2 = 0
//! ```
//!
//! Only `A, B, C1, C2` (evaluations of *public* polynomials) are published;
//! the witness evaluations stay inside the permutation trace, so the proof
//! leaks nothing about `r`, the noise, or the message beyond the STARK
//! itself.



use p3_air::symbolic::{BaseEntry, SymbolicVariable};
use p3_air::{Air, BaseAir, ExtensionBuilder, PermutationAirBuilder, WindowAccess};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_lookup::lookup_traits::{Kind, Lookup};
use p3_lookup::LookupAir;
use p3_matrix::dense::RowMajorMatrix;


use crate::regev::{centered_to_field, Ciphertext, EncryptionWitness, PublicKey};

// Main trace column indices.
pub const COL_A: usize = 0;
pub const COL_B: usize = 1;
pub const COL_C1: usize = 2;
pub const COL_C2: usize = 3;
pub const COL_R: usize = 4;
pub const COL_E1U: usize = 5;
pub const COL_E1V: usize = 6;
pub const COL_E2U: usize = 7;
pub const COL_E2V: usize = 8;
pub const COL_M: usize = 9;
pub const COL_K1: usize = 10;
pub const COL_K2: usize = 11;
pub const NUM_COLS: usize = 12;

// Auxiliary (permutation) column indices: one evaluation argument each.
// The first four are exposed; the verifier recomputes them from the
// statement. The rest stay hidden.
pub const AUX_A: usize = 0;
pub const AUX_B: usize = 1;
pub const AUX_C1: usize = 2;
pub const AUX_C2: usize = 3;
pub const AUX_R: usize = 4;
pub const AUX_E1: usize = 5;
pub const AUX_E2: usize = 6;
pub const AUX_M: usize = 7;
pub const AUX_K1: usize = 8;
pub const AUX_K2: usize = 9;
pub const NUM_AUX: usize = 10;

/// Interaction names for the exposed evaluation arguments, in aux-column
/// order. The verifier matches published values positionally; the names are
/// labels for debugging.
pub const EXPOSED_EVALS: [&str; 4] = ["eval:a", "eval:b", "eval:c1", "eval:c2"];

/// AIR for one Regev encryption. `n` is the ring dimension (= trace height)
/// and `delta` the message scaling `Δ`.
#[derive(Clone, Debug)]
pub struct RegevEncAir<F> {
    pub n: usize,
    pub delta: F,
}

impl<F: Field> RegevEncAir<F> {
    pub fn new(n: usize, delta: F) -> Self {
        assert!(n.is_power_of_two());
        Self { n, delta }
    }

    fn main_var(col: usize) -> SymbolicVariable<F> {
        SymbolicVariable::new(BaseEntry::Main { offset: 0 }, col)
    }
}

impl<F: Field> BaseAir<F> for RegevEncAir<F> {
    fn width(&self) -> usize {
        NUM_COLS
    }

    /// `[a || b || c1 || c2]`, binding the statement into the transcript
    /// before the evaluation challenge `z` is sampled.
    fn num_public_values(&self) -> usize {
        4 * self.n
    }

    /// No constraint reads the next row of the *main* trace (only the
    /// permutation trace steps rows), so main columns open at `zeta` only.
    fn main_next_row_columns(&self) -> Vec<usize> {
        vec![]
    }

    fn max_constraint_degree(&self) -> Option<usize> {
        Some(3)
    }
}

impl<F: Field> LookupAir<F> for RegevEncAir<F> {
    fn get_lookups(&mut self) -> Vec<Lookup<F>> {
        let col = |c: usize| -> Vec<Vec<_>> { vec![vec![Self::main_var(c).into()]] };
        let e1 = vec![vec![Self::main_var(COL_E1U) - Self::main_var(COL_E1V)]];
        let e2 = vec![vec![Self::main_var(COL_E2U) - Self::main_var(COL_E2V)]];

        let exposed = |name: &str| Kind::Global(name.to_string());

        // Aux column order must match the AUX_* constants.
        let specs: Vec<(Kind, Vec<Vec<_>>)> = vec![
            (exposed(EXPOSED_EVALS[0]), col(COL_A)),
            (exposed(EXPOSED_EVALS[1]), col(COL_B)),
            (exposed(EXPOSED_EVALS[2]), col(COL_C1)),
            (exposed(EXPOSED_EVALS[3]), col(COL_C2)),
            (Kind::Local, col(COL_R)),
            (Kind::Local, e1),
            (Kind::Local, e2),
            (Kind::Local, col(COL_M)),
            (Kind::Local, col(COL_K1)),
            (Kind::Local, col(COL_K2)),
        ];

        specs
            .into_iter()
            .enumerate()
            .map(|(aux, (kind, element_exprs))| Lookup {
                kind,
                element_exprs,
                multiplicities_exprs: vec![],
                columns: vec![aux],
            })
            .collect()
    }
}

impl<AB> Air<AB> for RegevEncAir<AB::F>
where
    AB: PermutationAirBuilder,
{
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local: &[AB::Var] = main.current_slice();

        // --- Smallness, all degree <= 3, no lookups needed ---------------

        // r ∈ {-1, 0, 1}
        let r: AB::Expr = local[COL_R].into();
        builder.assert_zero(r.clone() * (r.clone() - AB::Expr::ONE) * (r + AB::Expr::ONE));

        // CBD halves ∈ {0, 1, 2}
        for col in [COL_E1U, COL_E1V, COL_E2U, COL_E2V] {
            let x: AB::Expr = local[col].into();
            builder.assert_zero(
                x.clone() * (x.clone() - AB::Expr::ONE) * (x - AB::Expr::TWO),
            );
        }

        // m ∈ {0, 1}
        let m: AB::Expr = local[COL_M].into();
        builder.assert_bool(m);

        // --- Ring identities at the random point z -----------------------
        //
        // The Horner constraints themselves are emitted by the EvalGadget;
        // here we consume the running evaluations on the first row, where
        // perm[AUX_P][0] = P(z).

        let perm = builder.permutation();
        let s = |aux: usize| -> AB::ExprEF {
            perm.current(aux)
                .expect("permutation trace too narrow")
                .into()
        };

        let z: AB::ExprEF = builder.permutation_randomness()[0].into();

        // z^n + 1 by repeated squaring (n is a power of two).
        let mut zn = z;
        for _ in 0..self.n.trailing_zeros() {
            zn = zn.clone() * zn;
        }
        let zn1 = zn + AB::ExprEF::ONE;

        // A·R + E1 − C1 − (z^n+1)·K1 = 0
        let eq1 = s(AUX_A) * s(AUX_R) + s(AUX_E1) - s(AUX_C1) - zn1.clone() * s(AUX_K1);
        builder.when_first_row().assert_zero_ext(eq1);

        // B·R + E2 + Δ·M − C2 − (z^n+1)·K2 = 0
        let eq2 = s(AUX_B) * s(AUX_R) + s(AUX_E2) + s(AUX_M) * AB::Expr::from(self.delta)
            - s(AUX_C2)
            - zn1 * s(AUX_K2);
        builder.when_first_row().assert_zero_ext(eq2);
    }
}

/// Build the main trace for one ciphertext.
pub fn generate_trace(
    pk: &PublicKey,
    ct: &Ciphertext,
    witness: &EncryptionWitness,
) -> RowMajorMatrix<crate::regev::F> {
    use crate::regev::F;
    let n = pk.a.len();
    assert_eq!(witness.r.len(), n);

    let mut values = F::zero_vec(n * NUM_COLS);
    for i in 0..n {
        let row = &mut values[i * NUM_COLS..(i + 1) * NUM_COLS];
        row[COL_A] = pk.a[i];
        row[COL_B] = pk.b[i];
        row[COL_C1] = ct.c1[i];
        row[COL_C2] = ct.c2[i];
        row[COL_R] = centered_to_field(witness.r[i] as i64);
        row[COL_E1U] = F::from_u8(witness.e1u[i]);
        row[COL_E1V] = F::from_u8(witness.e1v[i]);
        row[COL_E2U] = F::from_u8(witness.e2u[i]);
        row[COL_E2V] = F::from_u8(witness.e2v[i]);
        row[COL_M] = F::from_u8(witness.m[i]);
        row[COL_K1] = witness.k1[i];
        row[COL_K2] = witness.k2[i];
    }
    RowMajorMatrix::new(values, NUM_COLS)
}

/// The public values for one instance: `[a || b || c1 || c2]`.
pub fn public_values(pk: &PublicKey, ct: &Ciphertext) -> Vec<crate::regev::F> {
    let mut pv = Vec::with_capacity(4 * pk.a.len());
    pv.extend_from_slice(&pk.a);
    pv.extend_from_slice(&pk.b);
    pv.extend_from_slice(&ct.c1);
    pv.extend_from_slice(&ct.c2);
    pv
}
