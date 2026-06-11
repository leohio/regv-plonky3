# regev-plonky3

Plonky3-friendly **Regev (Ring-LWE) encryption** with batched **STARK proofs of
correct encryption**. The whole design is organised around making the proof
circuit as small as possible:

| Design decision | Effect |
|---|---|
| ciphertext modulus `q` **=** proof-field prime (BabyBear, `2³¹−2²⁷+1`) | every ring operation is a native field op — no limb decomposition, no non-native reduction gadgets, ever |
| polynomial products verified by **random-point evaluation** (Schwartz–Zippel) | `c₁ = a·r + e₁` costs O(n) constraints instead of O(n²) (schoolbook) or O(n log n) (in-circuit NTT) |
| smallness via **degree-3 vanishing constraints** | ternary `r` and CBD(η=2) noise checked with zero auxiliary columns |
| one ciphertext = one `p3-batch-stark` instance | one commitment per phase + one FRI opening for the whole batch; fixed costs amortise |

## Performance

Apple Silicon, single proof batch, `n = 1024`, blowup 2, 84 queries, 16-bit
grinding (~100-bit conjectured FRI security):

```
encrypt:    252 µs/op
batch   1:  prove  9.2 ms              verify 12.9 ms   proof 226 KB
batch   8:  prove 43.2 ms (5.4 ms/ct)  verify 19.7 ms   proof  52 KB/ct
batch  32:  prove 163 ms  (5.1 ms/ct)  verify 46.7 ms   proof  33 KB/ct

zk batch 8: prove 228 ms (28.5 ms/ct)  verify  4.1 ms   proof  69 KB/ct
```

Reproduce with `cargo run --release --example bench` (pass `2048` for the
larger ring).

## The statement

For each ciphertext `(c₁, c₂)` under public key `(a, b)`, the proof shows
knowledge of `r, e₁, e₂, m` such that, in `Z_q[x]/(xⁿ+1)`:

```
c₁ = a·r + e₁                  r  ∈ {-1,0,1}ⁿ   (ternary)
c₂ = b·r + e₂ + Δ·m            e₁, e₂ ∈ [-2,2]ⁿ  (CBD η=2)
                               m  ∈ {0,1}ⁿ,  Δ = ⌊q/2⌋
```

## How the multiplication check works

The trace has one row per coefficient (height `n`); columns hold the
coefficients of `a, b, c₁, c₂, r, e1u, e1v, e2u, e2v, m, k₁, k₂` where
`k₁, k₂` are the quotients by `xⁿ+1`:

```
a·r + e₁        = c₁ + (xⁿ+1)·k₁     (over Z_q[x], degree < 2n−1)
b·r + e₂ + Δ·m  = c₂ + (xⁿ+1)·k₂
```

After the main trace is committed, Fiat–Shamir yields a challenge
`z ∈ F_(q⁴)`. A second (extension-field) trace of **Horner running sums**
`s[i] = p[i] + z·s[i+1]` evaluates every polynomial at `z` in O(n)
constraints, and the two identities are enforced **at the single point `z`**
by first-row constraints:

```
A·R + E₁       − C₁ − (zⁿ+1)·K₁ = 0
B·R + E₂ + Δ·M − C₂ − (zⁿ+1)·K₂ = 0
```

Soundness error is `< 2n/|F_(q⁴)| ≈ 2⁻¹¹³` per identity (Schwartz–Zippel).

The two-phase commit machinery is exactly Plonky3's logUp plumbing:
[`EvalGadget`](src/gadget.rs) implements the evaluation argument as a
`p3_lookup::LookupGadget`, and [`src/stark/`](src/stark) is a ~30-line-diff
vendor of `p3-batch-stark`'s prover/verifier that swaps the gadget and
samples one shared challenge `z`.

### Privacy

Only the evaluations of the **public** polynomials `a(z), b(z), c₁(z), c₂(z)`
are published (and cross-checked by the verifier against the actual
statement — that comparison is what binds the committed columns to the
claimed ciphertext; see `rejects_malicious_prover_with_forged_columns`).
The witness evaluations `r(z), e₁(z), e₂(z), m(z), k₁(z), k₂(z)` never leave
the committed permutation trace: publishing them would leak ~124 bits per
polynomial and allow dictionary attacks on low-entropy messages.

The default config is succinct but **not zero-knowledge** (plain FRI leaks
negligible-but-nonzero information about the trace). `zk_config()` runs the
same pipeline over Plonky3's `HidingFriPcs` (hiding Merkle commitments,
masked polynomials, randomized quotient chunks) for statistical zk. Note
hiding FRI requires `log_blowup >= 2` — with blowup 1 the masked polynomials
have rate 1 and verification fails.

### Why no logUp lookups for the smallness checks?

For ranges as small as `{-1,0,1}` / `{0,1,2}`, a vanishing polynomial
`x(x−1)(x+1)` is a single degree-3 constraint with **zero** extra columns.
A logUp lookup would cost an extension-field running-sum column (4 base
columns) per check plus batch inversions. Lookups only win for wide ranges
(e.g. 8/16-bit limbs of a balance range proof) — the gadget layer is already
in place if you want to add them.

## Usage

```rust
use regev_plonky3::*;
use rand::{rngs::SmallRng, RngExt, SeedableRng};

let params = RegevParams::N1024;
let mut rng = SmallRng::seed_from_u64(42);
let (pk, sk) = keygen(&mut rng, &params);

// Encrypt a 1024-bit message; the witness is the proof input.
let m: Vec<u8> = (0..params.n).map(|_| rng.random_range(0..=1)).collect();
let (ct, witness) = encrypt(&mut rng, &params, &pk, &m);
assert_eq!(decrypt(&params, &sk, &ct), m);

// Prove + verify (batch as many ciphertexts as you like).
let config = default_config();
let proof = prove_encryptions(&config, &params, &pk, &[ct.clone()], &[witness]);
verify_encryptions(&config, &params, &pk, &[ct], &proof).unwrap();
```

## Security notes — read before deploying

- **Lattice security.** `q ≈ 2³¹` is much larger than e.g. Kyber's 3329, and
  noise is narrow (σ = 1), so the dimension must carry the hardness.
  Ballpark (HE-standard-style extrapolation): `n = 1024` with ternary secrets
  and CBD(2) lands around **2¹⁰⁰ classical** — usable but below a 128-bit
  target; `n = 2048` has comfortable margin. **Run the
  [lattice-estimator](https://github.com/malb/lattice-estimator)** on your
  final `(n, q, σ)` before deploying. The proof system is agnostic to `n`
  (any power of two ≥ the FRI minimum).
- **Knowledge soundness vs. exact noise distribution.** The proof shows the
  noise is *in range* `[-2, 2]`, not that it is CBD-distributed — the
  standard relaxation for verifiable encryption; it does not affect IND-CPA
  of honestly generated ciphertexts and bounds decryption noise for proven
  ones.
- **FRI parameters.** `default_config()` targets ~100-bit conjectured
  soundness (blowup 2, 84 queries, 16-bit grinding); tune to taste.
  `test_config()` is *insecure* and for tests only.
- **No key-correctness proof.** The statement takes `(a, b)` as given; if
  your application needs `b = a·s + e` proven, the same machinery applies
  (one more instance with `s, e` as witness).

## Crate layout

```
src/ntt.rs     radix-2 NTT, full 2n-product with negacyclic quotient
src/regev.rs   keygen / encrypt / decrypt + witness (k₁, k₂, CBD halves)
src/params.rs  parameter sets (q is pinned to BabyBear)
src/gadget.rs  EvalGadget: Horner evaluation argument as a LookupGadget
src/air.rs     RegevEncAir: columns, smallness constraints, ring identities
src/stark/     vendored p3-batch-stark prover/verifier (gadget swapped)
src/prove.rs   prove_encryptions / verify_encryptions (statement binding)
src/config.rs  BabyBear + Poseidon2 + FRI configs
```

## Future work

- Poseidon2-based hiding MMCS for a faster zk prover (current `zk_config`
  uses Keccak, ~5x slower hashing than the non-zk Poseidon2 config).
- Observe `(a, b)` once per batch instead of per instance (verifier is
  currently dominated by transcript hashing of public values).
- Mersenne31 + Circle STARK backend (`q = 2³¹−1`; needs a Karatsuba/Toom
  encryptor since M31 has no two-adic NTT).
- logUp range lookups for wide ranges (balance range proofs on `m`).
