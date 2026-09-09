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
(`opt-level=3, lto=true, codegen-units=1`), warmed NTT tables, compute-only
unless noted. Java: `fibonacci-1.0-SNAPSHOT.jar` on JDK 21.0.8, same machine.

| n | decimal digits | Rust compute | Rust decimal | Java compute |
|---|---------------|--------------|--------------|--------------|
| 10⁴ | 2,090 | 0.3 ms | — | — |
| 10⁵ | 20,899 | 6.3 ms | ~1 ms | ~5 ms steady-state |
| 10⁶ | 208,988 | 28.7 ms | 19 ms | **CRASHES** (`add: overflow`) |
| 10⁷ | 2,089,877 | 219 ms | 506 ms | **CRASHES** (same bug) |
| 5·10⁷ | ~10,449,382 | ~2.3 s | — | **CRASHES** (same bug) |
| 10⁸ | 20,898,764 | 4.53 s | 15.5 s | **CRASHES** (same bug) |

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
