# fibonacci-rust — technical documentation

Algorithms, optimization history, and port notes. For installation and usage,
see [README.md](README.md).

## How it works — simple version

Computing F(n) digit-by-digit would take n steps. Instead the program uses
three ideas:

1. **Fast doubling (O(log n) steps).** Two identities jump from
   `(F(k), F(k+1))` straight to `(F(2k), F(2k+1))` with just a few big
   multiplications. Walking the 27 bits of n = 10⁸ means ~27 doubling steps
   instead of a hundred million additions. Each step needs `a²`, `b²`, `a·b`.
2. **NTT multiplication (O(m log m) per multiply).** Big-number multiplication
   by convolution: transform both numbers into "frequency space" (Number
   Theoretic Transform — an exact, integer, finite-field FFT), multiply
   pointwise, transform back. One transform pair is reused for all three
   products, and three small prime moduli run on three thread groups,
   combined exactly at the end (Garner CRT). Small numbers use plain
   schoolbook multiplication.
3. **Decimal limbs (no conversion step).** Limbs are base-10⁹ digits, so
   printing is just zero-padded formatting of each limb — the old
   divide-and-conquer conversion (14 s at 10⁸) is gone entirely.

## How it works — detailed version

### 1. Fast-doubling loop (`src/fib.rs`)

From `F(2k) = F(k)·(2·F(k+1) − F(k))` and `F(2k+1) = F(k)² + F(k+1)²`: keep the
pair `(a, b) = (F(m), F(m+1))`. Seed `m` from the top 4 bits of `n` via a
16-entry table (skips the first iterations), then for each remaining bit from
MSB to LSB: double `(a,b) → (F(2m), F(2m+1))`; if the bit is 1, advance
`(a,b) → (b, a+b)`. Ends with `a = F(n)`. Cost: `⌊log₂ n⌋ − 3` doublings.

### 2. Limb representation (`src/bigint.rs`)

`BigInt { mag: Vec<u32>, len: usize }`: base-10⁹ limbs, little-endian
(`mag[0]` = least-significant 9 digits), `len` = used limbs. Add/sub/double
are textbook carry/borrow loops over decimal digits. Capacities: schoolbook
products allocate `len + 2` spare limbs (room for the fused double+add
of the doubling step); NTT products allocate `n + 4` where `n` is the
transform size (room for the CRT carry drain, ≤ 2 limbs).

### 3. NTT over three primes (`src/ntt.rs`)

A length-`m` product is the convolution `c[i] = Σⱼ a[j]·b[i−j]`, done as
`NTT⁻¹(NTT(a)·NTT(b))` in `O(n log n)`, `n` = next power of two ≥ needed
length. One prime can't hold the full coefficient, so we transform mod three
NTT-friendly primes `998244353 (2²³)`, `167772161 (2²⁵)`, `469762049 (2²⁶)`
(all primitive root 3; P1 caps exact sizes at `n ≤ 2²³`) and reconstruct each
coefficient with Garner's CRT. Details that matter:
- **Integer Barrett `mulmod`**: `q̂ = (x·μ) >> 64` with `μ = ⌈2⁶⁴/p⌉`, pure
  integer (no float conversion on the hot loop), one conditional fix-up.
  Valid over the full `a,b < 2³²` contract — proven in comments.
- **Layered roots**: `roots[1]=1`, `roots[k+j]` laid out so each butterfly
  stage reads sequentially (prefetcher-friendly) and the per-butterfly
  `w·wlen` update vanishes. Tables shared via `Arc` (clone = pointer copy).
- **Lazy reduction**: butterfly sums do one conditional subtract, so values
  roam up to 2³²−1 yet stay *congruent* — and `mulmod` is exact over that
  whole range, so exactness is preserved end to end.
- **Serial/parallel cutoff**: `n < 4096` runs serially (Rayon overhead loses
  below that); above, the 3 primes split via nested `rayon::join`, and inside
  each track the two forward + three inverse transforms split again — up to
  ~12 threads busy on one global pool.

### 4. Fused 3-product doubling (`BigInt::fib_double`)

The doubling needs `A=a²`, `B=b²`, `C=a·b` from one `(fa, fb)` transform pair
per prime (3 forward transforms saved per step vs separate multiplies), then
`F(2k+1) = A+B`, `F(2k) = 2C−A`. Ownership does the buffer management: each
track *moves* its `Vec<u32>`s into the Rayon closure and returns the
residues — no pool, no shared mutable state. Small operands (`< 64` limbs)
take the schoolbook path (`fib_double_small`).

### 5. Garner CRT (`BigInt::garner`)

Per coefficient: `v1 = r1`; `v2 = (r2−r1)·P1⁻¹ mod P2` (exact Euclidean
reduction); `t = (v1+P1·v2) mod P3` (double-Barrett, double corrected);
`v3 = (r3−t)·(P1P2)⁻¹ mod P3`; then the full value `v1 + P1·v2 + P1P2·v3`
(+ carry, exact in `u128`) emits one base-10⁹ digit via Barrett divmod,
carrying the quotient. Max coefficient at top size is `n·(10⁹)² ≈ 2^81.7`,
far below CRT capacity `≈ 2^86.3`. Precomputed: `P1·P2`, both modular
inverses, and the `÷10⁹` Barrett constant.

### 6. Decimal output (`src/decimal.rs`)

There is no conversion step: limbs already are 9-digit groups, so output is
top limb plain plus every lower limb zero-padded to exactly 9 digits
(manual digit loop, no `format!` overhead). Lower-limb blocks format in
parallel via `par_chunks` and concatenate in order.

### 7. CLI + stats (`src/main.rs`)

`clap` CLI, NTT pre-warm before timing, then mean / sample-stddev / t-based
95% CI (df 1–30 table, 1.96 beyond) / moment skewness. Errors: `thiserror`
enum (`FibError::TooLarge`) + `anyhow` context in `main`; capacity breaches
are `assert!`s (programmer invariants, proven in comments).

## Optimization history

Method: warmed release runs (`-a 3–5`), one change at a time, `cargo test`
plus full-output SHA256 verification after every change. No privileged
profiler was available (Windows, unelevated), so phases were timed with
temporary `Instant` guards (reverted afterwards) plus complexity analysis.

Phase profile at F(10⁷) back when decimal conversion still existed:
compute 0.22 s = forward NTT 0.18 + inverse NTT 0.26 + pointwise/Garner/alloc
≈ 0.06 (thread-seconds — NTT ≈ 85%, allocation churn ≈ 2%); decimal 0.55 s =
top-level `div_rem` 0.26 serial + deeper splits + `10^k` pow 0.15. Every
change below was kept or reverted purely on measured numbers.

### Round 1

1. **Integer Barrett `mulmod` (kept): 1.3–1.6× compute.** The butterfly's
   `f64`-quotient estimate chained an `f64` multiply plus a `cvttsd2si`
   float→int conversion per operation. Replaced with `q̂ = (x·μ) >> 64`,
   `μ = ⌈2⁶⁴/p⌉` — pure integer, same ±1 single-fix-up contract, validated by
   the 20k×3 u128-oracle test plus end-to-end hashes.
2. **`pow10` via our NTT core (reverted — measured worse).** Binary
   exponentiation through our multiply: 184 ms vs num-bigint `pow` 48 ms for
   10^1044938 in isolation; end-to-end decimal regressed 0.58→0.88 s.
   Head-to-head showed our single multiply (126 ms @100k limbs) trailing
   num-bigint's Karatsuba (94 ms) — our 3-prime NTT constants only pay off at
   larger sizes. Reverted in full.
3. **Not pursued after Round 1, with reasons:** `ibig::to_string` as decimal
   backend (846 ms vs our parallel 621 ms @10⁷); Barrett division for the top
   `div_rem` (≈3 mults + Newton-μ ≈ 10 mults ≈ 780 ms vs 250 ms — needs a
   faster multiply first); GMP/`rug` (no C toolchain — no vcpkg/MSYS/MinGW);
   cache-blocked NTT and buffer pooling (bandwidth-bound profile, negligible
   alloc share).

### Round 2 (light load budget: dev at ≤1M, verdicts at 10M)

4. **`target-cpu=native` (kept): ~1.7× compute, zero source risk.** One-line
   `.cargo/config.toml`; generic x86-64 was starving AVX2/BMI2 everywhere.
   Binary is machine-specific — deliberate.
5. **PGO (kept): ~6–8% compute.** `llvm-tools` + profile-generate → train
   → merge → profile-use. The final verdict came from an interleaved A/B on
   finished code (alternating binaries, PGO won 8/8 pairs). Zero source risk.
   See [Reproducing benchmark builds](#reproducing-benchmark-builds).
6. **Tuning sweep (all null, reverted):** `RAYON_NUM_THREADS` 12 vs 6,
   `PARALLEL_NTT_MIN_N` 4096 vs 16384 (join overhead negligible as predicted),
   decimal split threshold 15k/30k/60k (flat).
7. **u32 NTT buffers (kept): ~1.2× compute.** Working arrays and root tables
   `Vec<u64>` → `Vec<u32>` (values `< 2³²` by the laziness invariant),
   halving streaming traffic. At 10⁸ the top transform's working set drops
   from 32 MB (spills from L3) to 16 MB (fits). Hash-identical output.
8. **Shift-split decimal division (reverted — no gain).** Idea: `10^k`
   divisor = `2^k·5^k`, peel the power of two with a free shift+mask,
   divide by the ~30%-smaller `5^k`. Measured 0.52→0.52 s — the O(n)
   shift/mask/or passes and extra pow eat the saving. Reverted in full.
9. **mimalloc (kept): ~10% decimal, tighter variance.** The parallel decimal
   phase allocated heavily across 12 threads. Two-line change
   (`mimalloc = "0.1"` + `#[global_allocator]`); MSVC `cl.exe` from the
   installed VS 2022 toolchain builds the C sources fine.

### Round 3: custom Barrett division (built, verified, reverted)

10. **Barrett `div_rem` with Newton reciprocal (reverted — measured 3–4×
    worse).** Motivated by fresh evidence (our NTT multiply had become ~5×
    num-bigint's at 166k limbs): new `src/div.rs` (~400 lines) with
    word-doubling Newton reciprocal (u128-division seed), verify-by-property
    reciprocal finalization, and a bulletproof correction-loop Barrett core
    (exact from ANY estimate — proven terminating, verified by differential
    tests against constructed `u = q·d + r` oracles plus a 104k-digit
    end-to-end conversion oracle). Debugging it fixed two real bugs along the
    way (an off-by-`b` slice in reciprocal finalization; an iteration-count
    shortfall for small `n`). The finished, correct implementation then lost
    on the clock at every size (1M: 0.020→0.087 s; 10M: 0.52→1.88 s decimal)
    and was reverted in full. Two root causes, both measured:
    - *Per-node divisor explosion:* halves differ per node, so the tree built
      ~50 distinct Newton reciprocals per run instead of ~8 — canonical
      per-depth halves plus `OnceLock` dedup fixed the count but not the
      outcome.
    - *Nested-parallelism contention collapse:* Barrett mults parallelize via
      nested rayon joins *inside* an already-parallel tree; at 10⁷ the
      conversion hung past 600 s that way, while serial-inner mults finished
      (slowly). Parallelize at exactly one level — the outermost — or pay
      orders of magnitude, not percent.
    Standing lesson: their Burnikel-Ziegler division is ~3 mults at these
    sizes while ours needs ~6 (Newton-μ amortization fails with one division
    per divisor). Beating it needs truncated (middle-product) multiplies.

### Round 4: decimal limbs end-to-end (kept — ~10× end-to-end at 10⁸)

11. **Base-10⁹ limbs everywhere (kept): conversion deleted.** Instead of
    fighting division throughput, every limb became a 9-digit group: NTT
    inputs only need values `< 2³²` (10⁹ qualifies), max convolution
    coefficients (`n·(10⁹)² ≈ 2^81.7` at top size) stay far below CRT
    capacity (`≈ 2^86.3`), and Garner emits decimal digits directly (one
    Barrett divmod per coefficient, ~20 cycles). Printing is then
    zero-padded formatting — O(n), parallel over limb blocks. Rewrote
    add/sub/double/`from_u64`/schoolbook (u128 accumulators) for decimal
    carry; `num-bigint`/`num-integer` left `[dependencies]` entirely.
    Compute slowed ~17% (normalization overhead), decimal went 14.3 s → 7 ms
    @10⁸. Verified: 17 tests + F(10⁶)/F(10⁷)/F(10⁸) hashes identical.

## Port notes: bugs found in fibonacciV2 (all fixed here)

**1. Signed-overflow in `NTT.mulmod` corrupts results, then crashes.**
`NTT.java` documents "0 ≤ a, b < p → a·b < 2⁶⁰", but the forward butterfly
legitimately holds *lazily-reduced* sums up to 2³²−1, and these grow with
transform size — measured residue `3,690,333,181` at NTT size 1024. Squaring
such a residue exceeds 2⁶³−1, and Java's `long x = a * b` silently wraps
negative; the Barrett quotient is then taken of the wrapped value. Proof,
all measured on the pristine checkout: at F(300000), NTT size 1024, `B2`
(b²) differs from a schoolbook cross-check in **all 1026 limbs from limb 0**
while `A2`/`AB` from the same call are exact — one wrapped prime track
poisons the whole Garner CRT. Downstream the phantom limbs trip the buffer
checks (`copy: dest too small`, `add: overflow`); F(200000)/F(400000) happen
to be correct, F(300000)+ is not computable at all. Fix here: pure-`u64`
Barrett valid for the full `a, b < 2³²` contract, plus a unit test over the
whole range with a `u128` oracle.

**2. Garner CRT reduction is one conditional short.** `(r2 − r1) mod P2`
with a single `if (diff < 0) diff += p2` fails because `r1 < P1 ≈ 6·P2`
(measured 37–60% of coefficients stay negative). Fix: exact `rem_euclid`.

**3. `Unsafe`-poked `BigInteger` aliasing.** Hand-built `BigInteger` via
`sun.misc.Unsafe` field offsets — brittle across JDK releases. Moot here:
no `Unsafe` anywhere, and decimal output needs no conversion at all.

**4. Buffer over-allocation + parallel-pool oversubscription.** One global
`bufSize ≈ 4×` the estimate caps every buffer (which turns bug 1's phantom
limbs into a crash), and *two* `ForkJoinPool`s oversubscribe the machine.
Fix: exactly-sized owned buffers, one global Rayon pool.

**5. Doc inconsistencies.** Primes and thresholds differ between README and
code; past NTT size 2²³ Java would silently wrap — here it's `FibError`.

## Reproducing benchmark builds

The committed `.cargo/config.toml` already sets `target-cpu=native`. The
table numbers additionally use PGO (not committed — regenerate it):

```powershell
$env:RUSTFLAGS = "-C target-cpu=native -C profile-generate=$env:TEMP\pgo-data"
cargo build --release
.\target\release\fibonacci-rust.exe 1000000 -a 2
.\target\release\fibonacci-rust.exe 10000000 -s NUL
& "<toolchain>\lib\rustlib\x86_64-pc-windows-msvc\bin\llvm-profdata.exe" merge -o merged.profdata $env:TEMP\pgo-data
$env:RUSTFLAGS = "-C target-cpu=native -C profile-use=$env:TEMP\merged.profdata"
cargo build --release
```

Plain rebuilds without `RUSTFLAGS` silently drop PGO — re-run the recipe
after source changes. (`llvm-tools` rustup component required for merge.)

## Design notes

- **Ownership instead of pooling.** Values move into Rayon closures and back;
  dropping frees deterministically. No `Rc<RefCell<_>>`, no buffer pool, no
  uninitialized-read bug class.
- **Errors vs panics.** Fallible runtime conditions (`TooLarge`) are
  `Result`s via `thiserror`; capacity breaches are `assert!`s on proven
  programmer invariants; `.unwrap()` appears only in tests.
- **No `unsafe` in the crate.** The hottest loop is safe indexed access;
  bounds checks are structured so LLVM elides them.
