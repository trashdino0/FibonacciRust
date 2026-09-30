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

### Round 5: transform arithmetic (one kept threshold, three kills)

12. **Schoolbook threshold 64 → 768 limbs (kept).** A head-to-head sweep
    (schoolbook vs Karatsuba vs NTT, 64–8192 limbs) showed NTT's 9-transform
    constant loses below ~700 limbs: at 512 schoolbook wins 1.5×, at 1024 NTT
    is back ahead 1.4×. One-line change; small-n end-to-end improved
    (10⁵: 2.1 → 1.5 ms, 10⁶: 11.7 → 10.7 ms native).
13. **Karatsuba (built, tested, deleted).** Correct implementation with
    differential tests — but the sweep above left it a useful band only
    ~700–900 limbs wide (everywhere else schoolbook or NTT wins outright).
    ~150 lines for ~0.05 ms end-to-end: deleted, threshold raise keeps nearly
    all the gain.
14. **4× loop unroll (reverted — no gain).** Extracted `butterfly()` + unrolled
    the stage loop for ILP on the mulhi latency chain: 92.8 vs ~85 ms
    (overlapping noise, if anything worse). LLVM had it covered.
15. **Split-radix NTT (designed, spiked, killed by measurement).** The math
    checks out (−25% mults, verified formulas) and a first implementation
    even passed differential tests on first run — but measured **1.85×
    slower** (16.3 vs 8.8 ms @262144 points): gather/scatter traffic plus
    ~30k freelist operations per transform swamped the arithmetic saving.
    Also killed on paper first: complex-packing (needs a conjugation
    symmetry 𝔽ₚ doesn't have — proved impossible, not just slow) and
    Stockham/blocking (roofline measured 4.5 GB/s vs ~40 GB/s: CPU-bound,
    traffic tricks can't help).
    Standing lesson: the remaining transform cost is ALU-throughput on scalar
    integer ops (asm confirms: branchless BMI2, no vectorization). The only
    unconquered lever is explicit SIMD (~1.3× possible, high complexity +
    new dep) — deferred, poor ROI/risk at this point.

### Round 6: Montgomery reduction (units proved, production switch killed by A/B)

16. **Montgomery `R = 2³²` domain (built, hash-verified, reverted).** Units
    (`montmul`: two `u64` mults vs Barrett's two `u64` + `u128`-high;
    `mont_in`/`mont_out_scaled`; Hensel `PP` + `R²` consts) passed a
    3000-random × 3-prime differential test vs an independent `R⁻¹` oracle +
    roundtrips — after fixing two test bugs of ours (the `p·PP ≡ −1` check
    must mask to the low 32 bits; edge values must respect the `< p`
    contract). Full integration (Montgomery twiddle caches,
    `ntt_forward/inverse_mont` with strict `< p` discipline, entry as
    `REDC(x·R²)` needing no pre-reduction, exit fused with the `1/n`
    scaling, pointwise `montmul`) was hash-exact end to end (F(10M)
    `DEE6…`, F(100M) `098A…`). Verdict by back-to-back stash A/B at the
    target: Barrett 1.569 s vs Montgomery 1.609 s mean @100M (95% CIs fully
    overlapping; 10M likewise flat) — **no win, reverted**. Why: ~⅓ cheaper
    butterfly mults save ~3%, but the two added `O(n)` entry/exit streaming
    passes cost ~3–4% traffic — net zero. Standing lesson, sharpening Round
    5: Barrett's *lazy* reduction (raw `u32` limbs straight into the
    butterflies, zero entry passes) is a structural traffic advantage — any
    challenger must beat it *including* conversion cost, not just per-mult
    ALU. Resurrection pointer: units + transforms deleted here live in git
    history should entry/exit ever fuse into existing passes.

### Round 7: profile-guided build + killing the bit-reversal passes (both kept)

17. **PGO (kept, ~+4%).** Routine builds never used it (no `profile-use`
    anywhere, no `.profdata` in tree). Pipeline: instrumented
    (`-C profile-generate`), train (10M ×2 + 100M ×1), `llvm-profdata
    merge`, rebuild (`-C profile-use`). Interleaved 10M: 72–74 ms vs
    76–77 ms base. Cautionary data point: the first PGO 100M read 2.02 s
    vs 1.57 s non-PGO — pure heat soak after the 40 s instrumented
    training run, not a regression (interleaving proved it). Lesson holds:
    never compare across sessions without re-measuring. Profile goes stale
    on hot-loop edits — regenerate after code settles, not after every
    spike (all code gates below run non-PGO for fairness).
18. **DIF forward + `nobitrev` inverse (kept, ~−18% @100M).** Analysis:
    bit-reversal of streaming data costs a full random redistribution
    (~128 B/element with line amplification) — moving it (permuted
    pointwise writes) or re-emitting it (inverse-DIF, whose bit-reversed
    output would corrupt downstream positional BigInt arithmetic) both
    fail; only absorbing the permutation *into the stage order* wins.
    `ntt_forward_dif` (natural in, bit-reversed out, no `bit_rev`) pairs
    with stageless-entry `ntt_inverse_nobitrev`; pointwise between them is
    order-preserving. DIF needs strict `< 2p` discipline (sums of two
    stored values — DIT's lazy `< 2³²` doesn't transfer): entry reduces via
    one sequential `mulmod(·, 1)` pass, sums double-subtract. Same layered
    twiddle table, same `rt[len+j]` indexing, loop top-down. Proven by two
    new tests (full-`u32`/edge identity + DIF-path-matches-DIT-path), then
    hash-exact end to end. 10M flat (cache-resident: bitrev is cheap in L3
    there), 100M ABBA: 1.66 s → 1.37 s.     DIT forward + bit-reversed inverse
    stay as the test-oracle path.

### Round 8: radix-4 (correct first try, killed — heavier passes lose)

19. **Radix-4 DIF/DIT (built, tested, deleted).** One radix-4 stage per two
    radix-2 stages: same 4 mults per 4 elements, half the streaming passes,
    dedicated triple-twiddle tables, directional `i` (`−i` for inverse —
    derived, not guessed). Both new tests green on first run (identity +
    full-path match incl. odd-size fallback). Verdict, both orders, both
    regimes: 10M 93–99 ms vs 77–79 ms DIF (−25%, L3-resident), 100M 1.70 s
    vs 1.42 s DIF (−20%, streaming). Half the passes yet slower everywhere:
    post-DIF the bottleneck is the butterfly itself, and radix-4 lengthens
    its dependent chains (triple-subtract sums, `i`-mult gating y1/y3),
    doubles twiddle traffic, runs 8 concurrent streams in outer stages
    (vs 4), and odd sizes stay radix-2 anyway (dilution). Same family as
    the split-radix kill: pass count is not the cost — pass *weight* is.
    Deleted per precedent; resurrection pointer in git history.

### Round 9: exit fusion killed, entry fusion killed on paper

20. **Fusing `1/n` into the last inverse stage (built, measured, reverted).**
    Same mult count, one less streaming pass — modeled as a ~2% win.
    Measured 10M: 122 ms vs 75 ms (**+65%**), base binary 74 ms same
    session (not heat). Cause: the `if last` branch plus extra live values
    (`ni`, staged outputs) perturb the hot loop's codegen in *every* stage
    — spills/uop-cache pressure — and serialize two dependent mults into
    the longest-stride stage; the standalone exit loop was perfectly
    prefetchable streaming. Reverted; 75 ms recovered immediately.
21. **Entry-into-`fill_u32` fusion (killed without building).** Same pattern
    (fuse mults into an existing loop), same ~2% ceiling, same codegen
    fragility just demonstrated. Standing lesson from both: the     butterfly
    loops are finely balanced — micro-traffic fusions that touch hot loops
    are guilty until proven innocent, and ~2% ceilings don't pay for the
    risk.
22. **PGO regenerated on settled code: neutral, routine stays non-PGO.**
    Two fresh profiles on the DIF code (thin + rich training) both measured
    ~0% across 6+ interleaved 10M pairs (75–82 ms either way), vs +4% on the
    pre-DIF shape — PGO's payoff is code-shape-dependent, not a constant.
    Verdict: keep the procedure documented, build routine non-PGO
    (simpler, reproducible). Final state this session: plain
    `cargo build --release` + `target-cpu=native`, 10M ≈ 75 ms,
    100M ≈ 1.4 s (cool box; heat-soaked reads up to ~1.7 s — see
    `learned.md` §8).

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
