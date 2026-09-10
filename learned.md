# What this project taught me (fibonacci-rust → prime calculator)

Everything below was learned the hard way on huge-Fibonacci computation and
applies directly to bignum-heavy work like primality proving. Each lesson
includes *why*, not just *what* — the measurement or bug behind it.

## 1. Measure first, one change at a time, revert losers

- **Baseline before touching anything**, release build, warmed runs, with
  statistics (we used mean / 95% CI / stddev / skewness via `-a N`). A single
  number is not a measurement.
- **One change + re-measure per step.** Several of our "obvious wins" lost on
  the clock (NTT-`pow10`: 0.58 → 0.88 s; shift-split division: 0.52 → 0.52 s;
  Barrett division: 3–4× *worse*). Every one was reverted in full. If you
  batch changes, you can't attribute anything.
- **Interleaved A/B for close calls.** Alternating binaries (PGO vs native,
  3+ pairs each) beat comparing two noisy means. PGO won 11/11 pairs at
  ~6–8% — believable only because of the pairing.
- **Distrust overlapping confidence intervals.** An early PGO reading showed
  +9% with overlapping CIs; the honest verdict was "undecided, re-test on
  final code" — which later confirmed it properly.
- **Keep a decision log with numbers** (ours lives in `docs.md`): what was
  tried, before/after, kept or reverted and why. Negative results are the
  most valuable entries — they stop you retrying dead ends.

## 2. Profile before optimizing (and know your profilers' limits)

- No privileged profiler was available (Windows, unelevated — samply
  requires admin), so phases were timed with temporary `Instant` guards,
  reverted afterwards. Ugly but honest; state the method alongside numbers.
- The profile said: NTT ≈ 85% of compute, allocation ≈ 2% — which immediately
  killed pooling/allocator ideas *for compute* (but mimalloc later earned its
  place on the allocation-heavy decimal side — profile each phase separately).
- Read your dependencies' vendored sources (`~/.cargo/registry/src/...`).
  Discovering num-bigint's Burnikel-Ziegler division + u64 limbs + Toom-3
  killed two planned optimizations before a line was written.

## 3. Parallelize at exactly one level

The single most expensive lesson. Barrett-division mults parallelized via
nested rayon joins *inside* an already-parallel divide-and-conquer tree:
- At small sizes it looked fine (even 2× faster — misleading!).
- At 10⁷ the conversion **hung past 600 s**; serial-inner mults finished.
- Nested fine-grained joins contend destructively (stealing storms, parked
  threads) while the outer tree already saturates all cores.

Rule: rayon at the outermost parallel level only; everything beneath it runs
serially. If tempted otherwise, measure at the largest size first — contention
collapse is superlinear and invisible in small tests.

## 4. Big-integer arithmetic that transfers directly

Primality work (Miller–Rabin, sieves with modular arithmetic, primorials)
lives and dies on modular multiplication — everything below applies as-is:

- **Integer Barrett reduction** (`q̂ = (x·μ) >> 64`, `μ = ⌈2⁶⁴/p⌉`): replaces
  `% p` and float quotient estimates. Provably within ±1 for `x < 2⁶⁴`, one
  fix-up, no `div` instruction, no float→int conversion latency (dropping the
  `f64` version was a measured 1.6× on butterflies). Validate against a
  `u128` oracle over the FULL input range, including `2³²−1` × `2³²−1`.
- **Know your overflow regime exactly.** The Java original did `long x = a*b`
  with residues that legitimately reached 3.69×10⁹ — silent signed wrap past
  2⁶³, garbage residues, then crashes-or-wrong-digits depending on luck. In
  Rust: `u64` products `< 2⁶⁴` are exact; `u128` for the quotient math;
  `debug_assert!` every contract (`a, b < 2³²`).
- **Newton reciprocal + verify-by-property finalization.** For repeated
  division by one modulus: Newton-double from a u128-division seed
  (~60 good bits saves ~5 iterations), then finalize by the *exact defining
  property* (`μ·d ≤ b²ⁿ < (μ+1)·d`) rather than trusting convergence math.
  Correction loops must terminate structurally (monotone progress toward a
  bounded target) so they return exact answers from ANY starting estimate —
  Newton quality then affects only speed, which tests assert (`corrections ≤ 4`).
- **Three-modulus NTT + Garner CRT** for exact large multiplication: three
  `< 2³⁰` primes (pairwise products fit `u64`), layered root tables, lazy
  single-subtract reduction (values roam `< 2³²` but stay *congruent* — prove
  the invariant, don't just hope). Bottleneck prime caps exact sizes
  (ours: 2²³); return errors past it, never wrap.
- **Single subtract, don't `%`, in butterflies**; precompute roots once and
  share via `Arc`; gate rayon below ~4096 points (measured join-overhead
  cutoff, not guessed).
- **Base-10⁹ limbs** when decimal output dominates: NTT only needs inputs
  `< 2³²`, and output formatting becomes trivial. Cost: decimal carry in
  add/sub/mul (u128 accumulators) + per-multiply normalization. Check the
  CRT capacity inequality (`n·base² < P₁·P₂·P₃`) before committing to a base.
- **Karatsuba/Toom threshold thinking:** our 3-prime NTT constants lose to
  Karatsuba below ~100k limbs — know your crossover with a head-to-head
  measurement, don't assume asymptotically-faster always wins.

## 5. Division is the boss (plan around it)

- Division throughput decides base-conversion cost, not multiplication.
  Burnikel-Ziegler-class division runs ~3 mults; a Barrett scheme needs
  ~6 (Newton-μ amortization fails with one division per divisor — only pays
  with many divisions per modulus).
- `10^k = 2^k·5^k` splitting (free shift + smaller divide) sounds great and
  measured exactly zero — the O(n) shift/mask passes eat the saving. Measure,
  don't admire.
- If output is decimal and huge, changing the *limb base* beats changing the
  *division algorithm* (see §4, last bullet).

## 6. Correctness infrastructure (non-negotiable, cheap)

- **Differential tests with constructed oracles:** `u = q·d + r` built
  explicitly ⇒ division must return exactly `(q, r)` — no external oracle
  needed, covers boundaries (q/r = 0/max) by construction.
- **Property tests:** reciprocal verified by `μ·d ≤ b²ⁿ < (μ+1)·d`, exactly.
- **End-to-end hashes:** full outputs SHA256-compared against an independent
  implementation (Python) at every scale after every kept change. Three
  mid-build failures were all caught this way (two stale test oracles, one
  real off-by-limb bug) — plus one reversed array. Never ship a timing number
  without a hash beside it.
- **Concrete tiny examples beat abstract reasoning.** A `b = 10` decimal
  analogue of the Newton math exposed an off-by-`b` slice instantly after
  hours of correct-looking algebra failed to explain a test failure.
- **Tests must discriminate:** the Newton suite caught a real iteration-count
  bug (28M finalize corrections) — a test that can't fail is decoration.

## 7. Rust specifics that paid off

- Owned values moving through rayon closures replace buffer pools entirely —
  no `Rc<RefCell<_>>`, no use-after-release class, deterministic free.
- `thiserror` for library errors (matchable), `anyhow` + context at the
  binary edge; `assert!` for proven programmer invariants,
  `debug_assert!` for expensive ones; `.unwrap()` only in tests.
- `u128` arithmetic compiles to a few integer ops (fine); u128 *division*
  is a libcall (avoid on hot paths — Barrett-chunk it).
- `#[global_allocator] mimalloc`: ~10% on allocation-heavy parallel phases
  with tighter variance; needs a C compiler at build time (VS `cl.exe`
  worked — MSYS2/MinGW was *not* needed, unlike GMP).
- `target-cpu=native` was the single biggest cheap win (~1.7× compute):
  generic x86-64 starves AVX2/BMI2 everywhere including dependency code.
  Binary becomes machine-specific — say so in the README.
- PGO recipe that works on stable Windows: `llvm-tools` component →
  profile-generate → train on representative inputs → `llvm-profdata merge`
  → profile-use. ~6–8%, zero source risk — but re-run after source changes
  (stale profiles silently apply) and note it clearly (plain rebuilds drop it).
- `LazyLock` + `Mutex<HashMap>` + `OnceLock` for exactly-once cached
  computation shared across threads (reciprocals, powers, root tables).
- `cargo clippy --all-targets` + `cargo fmt --check` + `cargo test` as the
  gate before every commit; `cargo doc --no-deps` catches broken intra-doc
  links.
- `BigUint::new` takes `Vec<u32>` even though internals are u64 — check
  public digit APIs before assuming layout parity.

## 8. Benchmark hygiene on a noisy box

- Sustained all-core runs throttle: identical binaries varied 25%+ across a
  long session (3.70 vs 4.61 s). Interleave, repeat, report bands — never
  single samples for close calls.
- **Kill hygiene:** timed-out benchmark processes may survive the timeout and
  starve everything after them (this masqueraded as a 600 s "hang" twice).
  After any kill: `Get-Process <name>` and terminate strays before trusting
  further numbers.
- **Pipe buffering lies:** `2>&1 | Select-Object` output is lost on kill;
  redirect to a file (`> out.txt 2>&1`) so partial output survives. Related:
  `print!` without newline stays buffered to files — don't read absence of a
  line as proof a phase never started.
- Keep heavy runs rare by policy (dev ≤ small, verdicts at medium, one final
  large validation) — and say the policy in the docs so future-you obeys it.

## 9. What I'd do first on the prime calculator

1. Same harness day one: `-a N` stats, hash-vs-known-answer checks
   (known prime tables / sympy as oracle), clippy+fmt+test gate.
2. Miller–Rabin on integer-Barrett `mulmod` (§4) — no NTT needed until
   numbers pass ~2000 bits; Karatsuba/Toom threshold measured, not guessed.
3. Sieve work is cache work: segmented + wheel + bit-packed odds, benchmark
   at L3-overflow sizes specifically.
4. If huge exact products appear: the 3-prime NTT + Garner here ports almost
   verbatim (only the output base changes).
5. If decimal output of huge primorials matters: base-10⁹ limbs from the
   start (§4, last bullet) — don't build a conversion pipeline you'll delete.
6. Write `learned.md` for it too — future-you will thank present-you.
