# fibonacci-rust

Huge Fibonacci numbers via **fast doubling + 3-prime NTT**, with Rayon parallelism.
A Rust port of `fibonacciV2` (Java), with its bugs fixed and verified bit-for-bit
against Python's exact big integers up to F(10⁷) (2,089,877 digits).

## Quick start

```bash
cargo build --release
./target/release/fibonacci-rust 1000000 -a 5 -w 20,10000
./target/release/fibonacci-rust 1000000 -p -s fib.txt
```

CLI mirrors the Java version:

| Flag | Meaning |
|------|---------|
| `n` | index of the Fibonacci number |
| `-a N`, `--average N` | run N times, print mean / 95% CI / stddev / skewness |
| `-w [RUNS,N]`, `--warmup [RUNS,N]` | warmup runs first (bare flag = `50,10000`) |
| `-p`, `--print` | print the full decimal result |
| `-s FILE`, `--save FILE` | save the decimal result to a file |

## Benchmarks

Machine: AMD Ryzen 5 3600X, 6 cores / 12 threads, Windows. Release build
(`opt-level=3, lto=true, codegen-units=1`, `target-cpu=native`, PGO — see the
optimization log for the exact recipe), warmed runs. Java:
`fibonacci-1.0-SNAPSHOT.jar` on JDK 21.0.8, same machine.

| n | decimal digits | Rust compute | Rust decimal | Java compute |
|---|---------------|--------------|--------------|--------------|
| 10⁴ | 2,090 | 0.1 ms | — | — |
| 10⁵ | 20,899 | 1.9 ms | ~1 ms | ~5 ms steady-state |
| 10⁶ | 208,988 | 9.0 ms | ~15 ms | **CRASHES** (`add: overflow`) |
| 10⁷ | 2,089,877 | ~70 ms | ~0.45 s | **CRASHES** (same bug) |
| 5·10⁷ | ~10,449,382 | 0.75 s | — | **CRASHES** (same bug) |
| 10⁸ | 20,898,764 | 1.77 s | ~14.3 s | **CRASHES** (same bug) |

(Compute figures are means over 3–5 warmed runs. End-to-end at 10⁸:
~16.0 s vs ~20.0 s before optimization (1.25×); at 10⁷: 0.83 → 0.52 s
(1.6×); at 10⁶: 48 → 24 ms (2.0×). Decimal conversion dominates at scale
and survived every attack below — see log. Original pre-optimization compute
for reference: 28.7 ms / 245 ms / 2.31 s / 4.53 s for
10⁶ / 10⁷ / 5·10⁷ / 10⁸.)

### Optimization log (measured, release build, same machine)

Method: warmed release runs (`-a 3–5 -w …`), one change at a time, `cargo test`
plus full-output SHA256/modular verification after every change. No privileged
profiler was available (Windows, unelevated), so phases were timed with
temporary `Instant` guards (reverted afterwards) plus complexity analysis.

Phase profile at F(10⁷): compute 0.22 s = forward NTT 0.18 + inverse NTT
0.26 + pointwise/Garner/alloc ≈ 0.06 (thread-seconds — NTT ≈ 85%, allocation
churn ≈ 2%, so pooling/allocator swaps were skipped on evidence); decimal
0.55 s = top-level `div_rem` 0.26 serial + deeper splits + `10^k` pow 0.15.

1. **Integer Barrett `mulmod` (kept): 1.3–1.6× compute.** The butterfly's
   `f64`-quotient estimate chained an `f64` multiply plus a `cvttsd2si`
   float→int conversion per operation. Replaced with `q̂ = (x·μ) >> 64`,
   `μ = ⌈2⁶⁴/p⌉` — pure integer, same ±1 single-fix-up contract, validated by
   the 20k×3 u128-oracle test plus end-to-end hashes. Compute: 28.7→20.9 ms
   (10⁶), 245→159 ms (10⁷), 2.31→1.77 s (5·10⁷). Decimal unchanged (it is
   dominated by num-bigint division, untouched).
2. **`pow10` via our NTT core (reverted — measured worse).** Binary
   exponentiation through our multiply: 184 ms vs num-bigint `pow` 48 ms for
   10^1044938 in isolation; end-to-end decimal regressed 0.58→0.88 s.
   Head-to-head showed our single multiply (126 ms @100k limbs) trailing
   num-bigint's Karatsuba (94 ms) — our 3-prime NTT constants only pay off at
   larger sizes. Reverted in full; lesson preserved here instead of in code.
 3. **Not pursued after Round 1, with reasons:** `ibig::to_string` as decimal backend
    (846 ms vs our parallel 621 ms @10⁷ — its Display doesn't beat parallel
    D&C); Barrett division for the top `div_rem` (≈3 mults + Newton-μ ≈ 10
    mults ≈ 780 ms vs current 250 ms — needs a faster multiply first, see 2);
    GMP/`rug` (no C toolchain on this machine — no vcpkg/MSYS/MinGW);
    cache-blocked NTT and buffer pooling (profile says bandwidth-bound with
    negligible alloc share — diminishing returns, stated plainly).

### Round 2 log (same method; light load budget: dev at ≤1M, verdicts at 10M)

4. **`target-cpu=native` (kept): ~1.7× compute, zero source risk.** One-line
   `.cargo/config.toml`; generic x86-64 was starving AVX2/BMI2 everywhere
   (our u128 Barrett, num-bigint's u64 kernels). Compute @10⁷: 159→94 ms;
   decimal unchanged (division doesn't vectorize). Binary is now
   machine-specific — deliberate.
5. **PGO (kept): ~6–8% compute.** `llvm-tools` + profile-generate → train
   (`1M -a2` + `10M -s NUL`) → merge → profile-use. First attempt showed
   +9% but overlapping CIs, so the final verdict came from an interleaved
   A/B on finished code (alternating binaries, PGO won 8/8 pairs: 1M
   9.6 vs 10.2 ms, 10M 74.5 vs 80.8 ms). Zero source risk. Recipe:
   ```powershell
   $env:RUSTFLAGS = "-C target-cpu=native -C profile-generate=$env:TEMP\pgo-data"
   cargo build --release
   .\target\release\fibonacci-rust.exe 1000000 -a 2
   .\target\release\fibonacci-rust.exe 10000000 -s NUL
   & "<toolchain>\lib\rustlib\x86_64-pc-windows-msvc\bin\llvm-profdata.exe" merge -o merged.profdata $env:TEMP\pgo-data
   $env:RUSTFLAGS = "-C target-cpu=native -C profile-use=$env:TEMP\merged.profdata"
   cargo build --release
   ```
   (Plain rebuilds without `RUSTFLAGS` silently drop PGO — re-run the recipe
   after source changes. Profile data lives outside the repo; regenerate, don't
   commit it.)
6. **Tuning sweep (all null, reverted):** `RAYON_NUM_THREADS` 12 vs 6
   (91 vs 96 ms — keep default), `PARALLEL_NTT_MIN_N` 4096 vs 16384
   (no difference — join overhead negligible as predicted),
   `TO_STRING_THRESHOLD_BITS` 15k/30k/60k (flat — keep 30k).
7. **u32 NTT buffers (kept): ~1.2× compute.** Working arrays and root tables
   `Vec<u64>` → `Vec<u32>` (values `< 2³²` by the laziness invariant;
   arithmetic stays `u64`), halving streaming traffic — plus `fill` becomes a
   `memcpy`. At 10⁸ the top transform's working set drops from 32 MB (spills
   from L3) to 16 MB (fits), which is why 10⁸ compute improved
   disproportionately (3.7 → 1.7 s). Compute @10⁷: 94→79 ms. Hash-identical
   output.
8. **Shift-split decimal division (reverted — no gain).** Idea: `10^k` divisor
   = `2^k·5^k`, peel the power of two with a free shift+mask, BZ-divide by
   the ~30%-smaller `5^k`. Measured 0.52→0.52 s — the O(n) shift/mask/or
   passes and extra `5^k` pow eat the division saving. Reverted in full.
9. **mimalloc (kept): ~10% decimal, tighter variance.** The parallel decimal
   phase allocates heavily (immutable-BigUint style: fresh quotient/remainder
   per `div_rem` across 12 threads). Decimal @10⁷: ~0.54 → ~0.49 s across
   4 samples, and the run-to-run band tightened visibly. Compute unchanged,
   as predicted (it barely allocates). Two-line change
   (`mimalloc = "0.1"` + `#[global_allocator]`); MSVC `cl.exe` from the
   installed VS 2022 toolchain builds the C sources fine.

### Round 3 log: custom Barrett division (built, verified, reverted)

10. **Barrett `div_rem` with Newton reciprocal (reverted — measured 3–4×
    worse).** Motivated by fresh evidence (our NTT multiply had become ~5×
    num-bigint's at 166k limbs), this was a full implementation, not a
    sketch: new `src/div.rs` (~400 lines) with word-doubling Newton reciprocal
    (u128-division seed), verify-by-property reciprocal finalization, and a
    bulletproof correction-loop Barrett core (exact from ANY estimate —
    proven terminating, verified by differential tests against constructed
    `u = q·d + r` oracles plus a 104k-digit end-to-end conversion oracle).
    Debugging it also fixed two real bugs along the way (an off-by-`b` slice
    in reciprocal finalization, caught by concrete example; an iteration-count
    shortfall for small `n`, caught by the tests' correction-count asserts).
    The finished, correct implementation then lost on the clock at every size
    (1M: 0.020→0.087 s; 5M: ~0.2→0.85 s; 10M: 0.52→1.88 s decimal) and was
    reverted in full — nothing of it remains in the tree (see git history).
    Two root causes, both measured:
    - *Per-node divisor explosion:* halves differ per node, so the tree built
      ~50 distinct Newton reciprocals per run instead of ~8 (measured with a
      counter) — canonical per-depth halves plus `OnceLock` dedup fixed the
      count but not the outcome.
    - *Nested-parallelism contention collapse:* Barrett mults parallelize via
      nested rayon joins *inside* an already-parallel tree; at 10⁷ the
      conversion hung past 600 s that way, while serial-inner mults finished
      (slowly). Parallelize at exactly one level — the outermost — or pay
      orders of magnitude, not percent.
    Standing lesson: their Burnikel-Ziegler division is ~3 mults at these
    sizes while ours needs ~6 (Newton-μ amortization fails with one division
    per divisor). beating it needs truncated (middle-product) multiplies,
    which is a project of its own — explicitly out of scope.

Correctness: full decimal strings hashed against Python (`hashlib.sha256`):
F(10⁵), F(2·10⁵) (also vs Java where it runs), F(10⁶), F(10⁷) — all identical.
F(10⁸) verified by exact digit count (20,898,764) plus F(10⁸) mod five
independent moduli (1,000,000,007 / 1,000,000,009 / 998,244,353 / 167,772,161 /
10¹⁸+3) computed by an independent Python fast-doubling implementation — all
match (false-pass odds ≈ 2⁻²⁰⁰).
Plus 15 `cargo test` unit tests (NTT round-trips, `mulmod` vs `u128` oracle
over the full `< 2³²` input range, schoolbook vs NTT cross-checks, known
F(100)/F(1000) values, doubling identities).

## Bugs found in fibonacciV2 (all fixed here)

**1. Signed-overflow in `NTT.mulmod` corrupts results, then crashes (the big one).**
`NTT.java:67` documents "0 ≤ a, b < p → a·b < 2⁶⁰", but the forward butterfly
legitimately holds *lazily-reduced* sums up to 2³²−1 (single-subtract instead
of full `mod p`), and these grow with transform size — measured residue
`3,690,333,181` at NTT size 1024. Squaring such a residue exceeds 2⁶³−1, and
Java's `long x = a * b` silently wraps negative. The Barrett quotient is then
taken of the wrapped value, which is congruent to `x − 2⁶⁴`, not `x (mod p)` —
the residue is garbage. Proof, all measured on the pristine checkout:
- F(300000), iteration at NTT size 1024: `B2` (b²) differs from a schoolbook
  cross-check in **all 1026 limbs from limb 0**, while `A2`/`AB` from the same
  call are exact — one wrapped prime track poisons the whole Garner CRT.
- Downstream the phantom high limbs trip the buffer checks: F(300000) dies
  with `copy: dest too small 16384 < 16386`, F(500000)/F(1000000) with
  `add: overflow 32768/65536`. Whether you get a crash or silently wrong
  digits depends on data-dependent headroom — F(200000)/F(400000) happen to
  be correct, F(300000)+ is not computable at all.
- Fix: pure-`u64` Barrett valid for the full `a, b < 2³²` contract
  (`a·b < 2⁶⁴`, quotient error `< 1e−4`, single correction provably suffices),
  plus a unit test over the whole range with a `u128` oracle — the exact
  regime Java never tested.

**2. Garner CRT reduction is one conditional short.**
`garnerCRT` reduces `(r2 − r1) mod P2` with a single `if (diff < 0) diff += p2`,
but `r1 < P1 ≈ 6·P2`, so the difference can stay negative (measured 37–60% of
coefficients). It only *usually* survives via the signed-`mulmod` rescue
branch. Fix: exact `rem_euclid` (O(1), always right).

**3. `Unsafe`-poked `BigInteger` aliasing.** `MutableBigInt.toBigIntegerFast`
hand-builds a `BigInteger` via `sun.misc.Unsafe` field offsets — brittle
across JDK releases. Fix: our limb layout (`Vec<u32>` LE) already matches
`num-bigint`'s, so conversion is one safe clone. No `unsafe` anywhere in
this crate.

**4. Buffer over-allocation + parallel-pool oversubscription.** Java sizes one
global `bufSize ≈ 4×` the estimate and caps every iteration buffer at it
(which is what turns bug 1's phantom limbs into a crash), and runs *two*
`ForkJoinPool`s (NTT + decimal) that oversubscribe the machine. Fix: every
intermediate is sized exactly (`next_pow2(2·max_len) + 2`, slack proven in
comments), values are simply owned (no pool, no `getLongArrayNoClear` garbage
class), and one global Rayon pool serves NTT + decimal.

**5. Doc inconsistencies.** README says P2/`TO_STRING_THRESHOLD` differ from
the code (50k bits vs the tuned 30k); one `AGENTS.md` copy lists P2 as
`1004535809` while the code uses `167772161`. Fix: single constants with the
tuned values, documented once. Also: past NTT size 2²³ (P1's limit, ≈F(190M))
Java would silently wrap — we return `FibError::TooLarge`.

## How it works — simple version

Computing F(n) digit-by-digit would take n steps. Instead the program uses
two classic speedups:

1. **Fast doubling (O(log n) steps).** Two identities let you jump from
   `(F(k), F(k+1))` straight to `(F(2k), F(2k+1))` with just a few big
   multiplications. Walking the ~24 bits of n = 10⁷ means ~24 doubling steps
   instead of ten million additions. Each step needs `a²`, `b²`, `a·b`.
2. **NTT multiplication (O(m log m) per multiply).** Big-number multiplication
   by convolution: transform both numbers into "frequency space" (Number
   Theoretic Transform — an exact, integer, finite-field FFT), multiply
   pointwise, transform back. One transform pair is reused for all three
   products (`a²`, `b²`, `a·b`), and three small prime moduli run on three
   thread groups, combined exactly at the end (Garner CRT). Small numbers
   (< 2048 bits) use plain schoolbook multiplication — NTT overhead isn't
   worth it there.

Printing 2M digits is its own problem (naive repeated division is
quadratic), so decimal conversion recursively splits the number in half
(`high · 10^half + low`) and converts both halves **in parallel**.

## How it works — detailed version

### 1. Fast-doubling loop (`src/fib.rs`)

From `F(2k) = F(k)·(2·F(k+1) − F(k))` and `F(2k+1) = F(k)² + F(k+1)²`: keep the
pair `(a, b) = (F(m), F(m+1))`. Seed `m` from the top 4 bits of `n` via a
16-entry table (skips the first iterations), then for each remaining bit from
MSB to LSB: double `(a,b) → (F(2m), F(2m+1))`; if the bit is 1, advance
`(a,b) → (b, a+b)`. Ends with `a = F(n)`. Cost: `⌊log₂ n⌋ − 3` doublings.

### 2. Limb representation (`src/bigint.rs`)

`BigInt { mag: Vec<u32>, len: usize }`: unsigned 32-bit limbs, little-endian
(`mag[0]` = least significant), `len` = used limbs, `mag[len..]` always zero.
Add/sub/shift are textbook carry/borrow loops. Capacities: schoolbook
products allocate `len + 2` spare limbs (proven room for the fused shift+add
of the doubling step); NTT products allocate `n + 2` where `n` is the
transform size (room for the CRT carry, which extends ≤ 2 limbs).

### 3. NTT over three primes (`src/ntt.rs`)

A length-`m` product is the convolution `c[i] = Σⱼ a[j]·b[i−j]`, done as
`NTT⁻¹(NTT(a)·NTT(b))` in `O(n log n)`, `n` = next power of two ≥ needed
length. One prime can't hold the full coefficient (`Σ` up to ~2⁸⁷ at max
size), so we transform mod three NTT-friendly primes
`998244353 (2²³)`, `167772161 (2²⁵)`, `469762049 (2²⁶)` (all primitive root 3;
P1 caps exact sizes at `n ≤ 2²³`) and reconstruct each coefficient with
Garner's CRT. Details that matter:
- **Barrett `mulmod`**: `q ≈ (a·b)/p` via precomputed `1.0/p` double, one
  conditional fix-up. Replaces `div` (~20–40 cycles) with multiply+adjust
  (~5). Valid over the full `a,b < 2³²` contract (quotient error `< 1e−4`;
  both over/under-estimate self-correct — proven in comments).
- **Layered roots**: `roots[1]=1`, `roots[k+j]` laid out so each butterfly
  stage reads sequentially (prefetcher-friendly) and the per-butterfly
  `w·wlen` update vanishes. Tables shared via `Arc` (clone = pointer copy).
- **Lazy reduction**: butterfly sums do one conditional subtract, so values
  roam up to 2³²−1 yet stay *congruent* — and `mulmod` is exact over that
  whole range, so exactness is preserved end to end (this is the invariant
  Java violated via signed overflow).
- **Serial/parallel cutoff**: `n < 4096` runs serially (Rayon overhead loses
  below that); above, the 3 primes split via nested `rayon::join`, and inside
  each track the two forward + three inverse transforms split again — up to
  ~12 threads busy, mirroring the Java fork structure but on one pool.

### 4. Fused 3-product doubling (`BigInt::fib_double`)

The doubling needs `A=a²`, `B=b²`, `C=a·b` from one `(fa, fb)` transform pair
per prime (3 forward transforms saved per step vs separate multiplies), then
`F(2k+1) = A+B`, `F(2k) = 2C−A`. Ownership does the buffer management Java's
`Workspace` did by hand: each track *moves* its `Vec<u64>`s into the Rayon
closure and returns the residues — no pool, no `NoClear` garbage class, no
use-after-release possible. Small operands (`< 64` limbs) take the
schoolbook path (`fib_double_small`).

### 5. Garner CRT (`BigInt::garner`)

Per coefficient: `v1 = r1`; `v2 = (r2−r1)·P1⁻¹ mod P2` (exact Euclidean
reduction — bug-2 fix); `t = (v1+P1·v2) mod P3` (double-Barrett, double
corrected); `v3 = (r3−t)·(P1P2)⁻¹ mod P3`; accumulate
`v1 + P1·v2 + P1P2·v3` into base-2³² words with 64-bit carry. Precomputed:
`P1·P2`, both modular inverses, and the `P1P2` lo/hi 32-bit split.

### 6. Parallel decimal (`src/decimal.rs`)

`digits ≈ bits·log₁₀(2)+1`; split `n = high·10^half + low` with one big
`div_rem`, recurse both halves with `rayon::join`, zero-pad `low` to `half`
digits. Below 30,000 bits, plain `to_string()`. `10^half` cached. The
`BigUint` bridge is one clone (identical limb layout).

### 7. CLI + stats (`src/main.rs`)

`clap` CLI (same flags as Java), NTT pre-warm before timing (like the Java
static block), warmup loop, then mean / sample-stddev / t-based 95% CI
(df 1–30 table, 1.96 beyond — matches Commons-Math `TDistribution`) /
moment skewness. Errors: `thiserror` enum (`FibError::TooLarge`) + `anyhow`
context in `main`; capacity breaches are `assert!`s (programmer invariants,
proven in comments — the footgun-checked choice).

## Layout

```text
fibonacci-rust/
├── Cargo.toml          (clap, rayon, num-bigint, thiserror, anyhow; release: lto, cg-units=1)
└── src/
    ├── main.rs         (CLI, warmup, stats, orchestration)
    ├── fib.rs          (fast-doubling loop + table tests + identity tests)
    ├── bigint.rs       (limbs, add/sub/shift, schoolbook, garner, fib_double)
    ├── ntt.rs          (primes, Barrett mulmod, layered roots, DIT transforms)
    └── decimal.rs      (parallel base conversion)
```

No `unsafe`, no build script, no workspace (single binary — split only when a
real second crate boundary appears).
