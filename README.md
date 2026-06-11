# regev-plonky3

A Plonky3-friendly **Regev (Ring-LWE) public-key encryption** scheme with
**value-level additive homomorphism**, batched **STARK proofs of correct
encryption**, an optional **plaintext range proof**, and a **transfer proof**
(`before = after + delta` over encrypted balances). The whole design is
organised around making the proof circuit as small as possible.

The proof circuit (the AIR / constraint system) lives entirely in this repo
([`src/air.rs`](src/air.rs)); [`src/stark/`](src/stark) is a ~30-line-diff
vendor of `p3-batch-stark`'s prover/verifier, and everything heavy (FRI,
Merkle commitments, Poseidon2, the field) comes from the published
[Plonky3](https://github.com/Plonky3/Plonky3) crates.

| Design decision | Effect |
|---|---|
| ciphertext modulus `q` **=** proof-field prime (BabyBear, `2³¹−2²⁷+1`) | every ring operation is a native field op — no limb decomposition, no non-native reduction gadgets, ever |
| polynomial products verified by **random-point evaluation** (Schwartz–Zippel) | `c₁ = a·r + e₁` costs O(n) constraints instead of O(n²) (schoolbook) or O(n log n) (in-circuit NTT) |
| smallness via **degree-3 vanishing constraints** | ternary `r` and CBD(η=2) noise checked with zero auxiliary columns |
| one ciphertext = one `p3-batch-stark` instance | one commitment per phase + one FRI opening for the whole batch; fixed costs amortise |
| plaintext range proof entirely in the **witness trace** | works identically under plain and hiding (zk) FRI; no preprocessed commitment to reconstruct |
| plaintext modulus `t = 2^plain_bits` with digit decoding (`Δ = ⌊q/t⌋`) | ciphertext addition is **value-level additive**: `decrypt_value(Enc(A) ⊞ Enc(B)) = A + B` |
| transfer AIR with a **ripple-carry column** | proves `before = after + delta` over three encrypted values in one instance, all degree ≤ 2 |

---

## Quick start

```rust
use regev_plonky3::*;
use rand::{rngs::SmallRng, RngExt, SeedableRng};

let params = RegevParams::N1024;          // n = 1024, CBD(η=2) noise
let mut rng = SmallRng::seed_from_u64(42);
let (pk, sk) = keygen(&mut rng, &params);

// Encrypt a 1024-bit message. `witness` is the secret input to the proof.
let m: Vec<u8> = (0..params.n).map(|_| rng.random_range(0..=1)).collect();
let (ct, witness) = encrypt(&mut rng, &params, &pk, &m);
assert_eq!(decrypt(&params, &sk, &ct), m);

// Prove + verify correct encryption (batch as many ciphertexts as you like).
let config = default_config();
let proof = prove_encryptions(&config, &params, &pk, &[ct.clone()], &[witness]);
verify_encryptions(&config, &params, &pk, &[ct], &proof).unwrap();
```

Run the test suite and the benchmark:

```sh
cargo test
cargo run --release --example bench          # n = 1024
cargo run --release --example bench -- 2048  # n = 2048
```

---

## Benchmarks

Apple Silicon (M-series), `release` profile, `n = 1024`, `default_config`
(blowup 2, 84 queries, 16-bit grinding ≈ 100-bit conjectured FRI security).
Reproduce with `cargo run --release --example bench`.

| Workload | Prove (total) | Prove / ct | Verify | Proof size | Proof / ct |
|---|---:|---:|---:|---:|---:|
| encrypt (no proof) | 241 µs | — | — | — | — |
| batch of 1 | 11.3 ms | 11.3 ms | 11.7 ms | 226 KB | 226 KB |
| batch of 8 | 42.4 ms | 5.3 ms | 18.4 ms | 415 KB | 52 KB |
| batch of 32 | 163 ms | 5.1 ms | 43.4 ms | 1.06 MB | 33 KB |
| zk, batch of 8 (`zk_config`) | 230 ms | 28.7 ms | 4.4 ms | 551 KB | 69 KB |
| + range proof, batch of 8 | 45.4 ms | 5.7 ms | 20.3 ms | 433 KB | 54 KB |
| transfer (3 cts + conservation), batch of 8 | 58.0 ms | 7.3 ms/tx | 30.8 ms | 722 KB | 90 KB/tx |

Takeaways:

- **Batching is essential.** FRI has a large fixed cost; at batch 32 the
  per-ciphertext prove time drops to ~5 ms and the proof to ~33 KB/ct.
- **The range proof is cheap**: ~1 ms/ct and ~2 KB/ct on top of the
  encryption proof (it adds 5 witness columns, all degree-≤2 constraints).
- **ZK costs ~2× to prove** (hiding commitments + randomized quotient chunks,
  and a Keccak-based MMCS) but verifies fast.
- **A transfer is cheaper than 3 separate encryption proofs** (7.3 ms vs
  ~16 ms): the three ciphertexts share one instance, one set of `a, b`
  columns and one evaluation challenge.

---

## What is proven

For each ciphertext `(c₁, c₂)` under public key `(a, b)`, the proof shows
knowledge of `r, e₁, e₂, m` such that, in `Z_q[x]/(xⁿ+1)`:

```
c₁ = a·r + e₁                  r  ∈ {-1,0,1}ⁿ   (ternary)
c₂ = b·r + e₂ + Δ·m            e₁, e₂ ∈ [-2,2]ⁿ  (CBD η=2)
                               m  ∈ {0,1}ⁿ,  Δ = ⌊q/t⌋, t = 2^plain_bits
```

With a range proof attached, the statement additionally asserts that the
integer encoded by the low `value_bits` message coefficients lies in
`[0, bound)` — **without revealing it** (see below).

### How the multiplication check works

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
`p3_lookup::LookupGadget`, and [`src/stark/`](src/stark) swaps that gadget in
and samples one shared challenge `z`.

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
hiding FRI requires `log_blowup ≥ 2` — with blowup 1 the masked polynomials
have rate 1 and verification fails.

---

## Additive homomorphism (value level)

Messages are scaled by `Δ = ⌊q/t⌋` with plaintext modulus `t = 2^plain_bits`
(default `t = 256`), and decryption decodes each coefficient to a **digit**
in `[0, t)` rather than a single bit. Values are encoded as little-endian
bits across coefficients; because radix-2 weights are linear, ciphertext
addition is exact integer addition with **no carry logic needed**:

```text
Σ (aᵢ + bᵢ) · 2^i  =  Σ aᵢ·2^i + Σ bᵢ·2^i
```

```rust
use regev_plonky3::*;
use rand::{rngs::SmallRng, SeedableRng};

let params = RegevParams::N1024;
let mut rng = SmallRng::seed_from_u64(0);
let (pk, sk) = keygen(&mut rng, &params);

let (ct_a, _) = encrypt(&mut rng, &params, &pk, &encode_value_message(5, params.n));
let (ct_b, _) = encrypt(&mut rng, &params, &pk, &encode_value_message(3, params.n));

let ct_sum = add_ciphertexts(&ct_a, &ct_b);
assert_eq!(decrypt_value(&params, &sk, &ct_sum), 8);   // not 5 XOR 3 = 6!
```

Budgets (default `plain_bits = 8`):

- **digits**: each coefficient digit must stay below `t`, so up to `t − 1 =
  255` stacked additions of binary-encoded values;
- **noise**: total noise must stay below `Δ/2 ≈ 2^22`; fresh-ciphertext noise
  is a few hundred, so thousands of additions — digits bind first.

`plain_bits` trades addition depth against noise margin and is a public
parameter (it does not affect lattice security, only correctness).

Combined with the encryption proofs this gives **verified homomorphic
sums**: verify the proofs of the input ciphertexts, then add them publicly —
by linearity the result is guaranteed to encrypt the sum of the proven
values (see the `homomorphic_sum_of_proven_ciphertexts` test).

---

## Transfer proofs: `before = after + delta`

For confidential-balance flows you also want the *opposite direction*: given
three independently encrypted values, prove in zero-knowledge that they
satisfy a conservation law. `prove_transfers` does this in a single STARK
instance per transfer:

> all three of `before`, `delta`, `after` are well-formed encryptions under
> `(a, b)`, **and** `before = after + delta` as n-bit integers.

Since `after` and `delta` are committed *bit vectors* (hence non-negative)
and the addition is exact, this simultaneously gives **no-underflow**:
`delta ≤ before` is implied. A third party can verify that a balance update
is conserving without learning any of the three values.

```rust,ignore
let (before, w_b) = encrypt(&mut rng, &params, &pk, &encode_value_message(100, n));
let (delta,  w_d) = encrypt(&mut rng, &params, &pk, &encode_value_message(42, n));
let (after,  w_a) = encrypt(&mut rng, &params, &pk, &encode_value_message(58, n));

let t = Transfer { before, delta, after };
let w = TransferWitness { before: w_b, delta: w_d, after: w_a };

let proof = prove_transfers(&config, &params, &pk, &[t.clone()], &[w]);
verify_transfers(&config, &params, &pk, &[t], &proof)?;
```

### How it works

The three message-bit columns are wired through a **ripple-carry adder**,
one bit per row, with a single extra carry column `c`:

```text
after[i] + delta[i] + c[i] = before[i] + 2·c[i+1]    (transition rows)
c[0] = 0                                             (first row)
after[n-1] + delta[n-1] + c[n-1] = before[n-1]       (last row ⇒ carry-out 0)
c[i] ∈ {0, 1}
```

All degree ≤ 2; the zero carry-out makes the equation hold over the
integers, not mod `2^n`. The instance has 33 main columns (shared `a, b` +
10 per ciphertext + carry) and 26 permutation columns (Horner evaluations;
`a, b, c1, c2`×3 exposed, witness evaluations hidden), so one transfer is
cheaper than three separate encryption proofs.

To additionally cap the amount (e.g. `delta < 2^16`), run a range-proof
instance on the same `delta` ciphertext — both proofs bind to the identical
public ciphertext, so accepting both yields conservation **and** the cap
(see the `transfer_composes_with_delta_range_proof` test).

---

## Plaintext range proofs

A confidential-balance application wants to prove not just that a ciphertext
is well formed, but that the *encrypted value* is in range (e.g. a transfer
amount is non-negative and below some cap) — all without revealing it. Attach
a [`RangeSpec`] and the AIR proves, for each ciphertext,

```
value = Σ_{i < value_bits} m[i] · 2^i   ∈   [0, bound)
```

where `value` is the little-endian integer encoded by the low `value_bits`
coefficients of the message. Higher message bits, if any, are unconstrained
by the range proof and can carry other payload.

### Usage

```rust
use regev_plonky3::*;
use rand::{rngs::SmallRng, RngExt, SeedableRng};

let params = RegevParams::N1024;
let mut rng = SmallRng::seed_from_u64(1);
let (pk, sk) = keygen(&mut rng, &params);

// Encode a value (e.g. a balance) in the low 20 bits of the message.
let value: u64 = 123_456;
let mut m = vec![0u8; params.n];
for i in 0..20 { m[i] = ((value >> i) & 1) as u8; }
let (ct, witness) = encrypt(&mut rng, &params, &pk, &m);

// Prove: well-formed encryption AND value ∈ [0, 1_000_000).
let spec = RangeSpec { value_bits: 20, bound: 1_000_000 };
let config = default_config();
let proof = prove_encryptions_with_range(&config, &params, &pk, &[ct.clone()], &[witness], spec);

// The verifier supplies the range it wants to enforce.
verify_encryptions_with_range(&config, &params, &pk, &[ct], &proof, spec).unwrap();
```

### How it works

The proof uses the standard complement technique: the prover supplies
complement bits `d[i]` with `Σ d[i] 2^i = bound − 1 − value`, and the AIR
enforces

```
Σ_{i < K} (m[i] + d[i]) · 2^i = bound − 1     (K = value_bits)
```

Since `value` and `d_value` are each a sum of `K` boolean bits — hence in
`[0, 2^K)` — and they sum to `bound − 1`, it follows that
`value ≤ bound − 1 < bound`.

Everything is materialised in **witness columns** (no preprocessed
commitment, deliberately — a preprocessed weights column would be *salted* by
the hiding PCS in zk mode and thus impossible for the verifier to
reconstruct). The five extra columns are:

- `flag[i]`: active indicator, `1` for the first `K` rows then `0`, pinned by
  `flag[0]=1`, a non-increasing constraint, and a count accumulator `cnt`
  with `cnt[0]=K` (this count is what stops a malicious prover from shrinking
  the window to hide value in dropped high bits);
- `w[i]`: the weight `2^i` while active, `0` afterwards, via `w[0]=1` and
  `w[i+1] = 2·w[i]·flag[i+1]`;
- `d[i]`: complement bits, forced to `0` where inactive;
- `acc`, `cnt`: suffix-sum accumulators so the bound and count checks are
  single first-row constraints.

All range constraints are degree ≤ 2, so the AIR's overall degree stays at 3
(set by the ternary/CBD smallness checks). The range is **private**: `value`
is never exposed; only the bound and `value_bits` are public, and the verifier
supplies the ones it wants to check — a proof for any other `(bound,
value_bits)` simply fails its constraints.

### Parameter limits

To rule out modular wraparound in `value + (bound−1−value) = bound−1` over
BabyBear, we require `2^(value_bits + 1) ≤ q`, i.e. **`value_bits ≤ 29`**
(`RangeSpec::MAX_VALUE_BITS`), and `1 ≤ bound ≤ 2^value_bits`. That covers
values up to ~536 million. For larger ranges (e.g. 64-bit balances),
decompose the value into limbs and run one range proof per limb — the
machinery generalises directly (future work).

### Why not logUp lookups for the smallness / range checks?

For ranges as small as `{-1,0,1}` / `{0,1,2}`, a vanishing polynomial
`x(x−1)(x+1)` is a single degree-3 constraint with **zero** extra columns.
A logUp lookup would cost an extension-field running-sum column (4 base
columns) per check plus batch inversions. Lookups only win for *wide* ranges
(e.g. 8/16-bit limbs) — the gadget layer is already in place if you want to
add them.

---

## API summary

| Function | Purpose |
|---|---|
| `keygen` / `encrypt` / `decrypt` | the Ring-LWE scheme |
| `prove_encryptions` / `verify_encryptions` | batched proof of correct encryption |
| `prove_encryptions_with_range` / `verify_encryptions_with_range` | same, plus a plaintext range proof |
| `prove_transfers` / `verify_transfers` | batched proof of `before = after + delta` over encrypted values |
| `add_ciphertexts` / `decrypt_value` / `encode_value_message` | value-level homomorphic addition |
| `default_config` | succinct, non-zk (plain FRI) |
| `zk_config` | statistical zero-knowledge (hiding FRI) |
| `test_config` | **insecure**, tiny parameters, for tests only |
| `RegevParams::N1024` / `N2048` | parameter sets |
| `RangeSpec { value_bits, bound }` | range-proof parameters |

All prove/verify functions are generic over the STARK config, so the same
code runs under `default_config` or `zk_config`.

---

## Security notes — read before deploying

- **Lattice security.** `q ≈ 2³¹` is much larger than e.g. Kyber's 3329, and
  the noise is narrow (σ = 1), so the dimension must carry the hardness.
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
  soundness; tune to taste. `test_config()` is *insecure* and for tests only.
- **No key-correctness proof.** The statement takes `(a, b)` as given; if
  your application needs `b = a·s + e` proven, the same machinery applies
  (one more instance with `s, e` as witness).

---

## Crate layout

```
src/ntt.rs     radix-2 NTT, full 2n-product with negacyclic quotient
src/regev.rs   keygen / encrypt / decrypt + witness (k₁, k₂, CBD halves)
src/params.rs  parameter sets (q is pinned to BabyBear)
src/gadget.rs  EvalGadget: Horner evaluation argument as a LookupGadget
src/air.rs     RegevEncAir: columns, smallness, ring identities, range proof
src/stark/     vendored p3-batch-stark prover/verifier (gadget swapped)
src/transfer.rs transfer AIR (ripple-carry conservation) + prove/verify
src/prove.rs   prove/verify wrappers (+ statement binding, range variants)
src/config.rs  BabyBear + Poseidon2 + (plain / hiding) FRI configs
```

Tests: [`tests/prove_verify.rs`](tests/prove_verify.rs) (end-to-end + soundness,
including a malicious-prover forgery and range-proof boundary/soundness cases)
and [`tests/zk_debug.rs`](tests/zk_debug.rs) (hiding-FRI isolation tests).

---

## Future work

- Limb-decomposed range proofs for values wider than 29 bits (64-bit balances).
- Fold the delta range cap directly into the transfer AIR (currently a
  second, composed instance).
- Poseidon2-based hiding MMCS for a faster zk prover (current `zk_config`
  uses Keccak, ~5× slower hashing than the non-zk Poseidon2 config).
- Observe `(a, b)` once per batch instead of per instance (the verifier is
  currently dominated by transcript hashing of public values).
- Mersenne31 + Circle STARK backend (`q = 2³¹−1`; needs a Karatsuba/Toom
  encryptor since M31 has no two-adic NTT).
