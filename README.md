# fibonacci-rust

Compute huge Fibonacci numbers — exact, fast, with full decimal output.
`F(10⁸)` (20,898,764 digits) in about 2 seconds on a 6-core desktop.

Fast doubling keeps the step count at `O(log n)`; NTT-based multiplication
keeps each step at `O(m log m)`; base-10⁹ limbs make decimal output essentially
free. See [docs.md](docs.md) for algorithms, benchmarks history, and port notes.

## Features

- Exact `F(n)` for `n` up to ~1.8×10⁸, verified bit-for-bit against Python
- Full decimal output at any size (no separate conversion bottleneck)
- Statistical benchmarking built in (mean, 95% CI, stddev, skewness)
- Single static binary, no runtime dependencies, no `unsafe`

## Performance

AMD Ryzen 5 3600X (6C/12T), warmed runs, release+PGO build:

| n | digits | compute | decimal | total |
|---|--------|---------|---------|-------|
| 10⁶ | 208,988 | 10 ms | ~0 ms | ~11 ms |
| 10⁷ | 2,089,877 | 86 ms | ~1 ms | ~0.09 s |
| 5·10⁷ | ~10,449,382 | 0.85 s | ~5 ms | ~0.85 s |
| 10⁸ | 20,898,764 | 2.00 s | ~7 ms | ~2.0 s |

## Installation

Prerequisites: a recent stable Rust toolchain (≥ 1.80), plus a C compiler —
[VS Build Tools](https://visualstudio.microsoft.com/downloads/) with MSVC on
Windows, `cc`/`gcc` on Linux/macOS (needed to build the bundled `mimalloc`).

```bash
git clone https://github.com/trashdino0/FibonacciRust.git
cd FibonacciRust
cargo build --release
```

The default build already sets `target-cpu=native` (see `.cargo/config.toml`),
so the binary is tuned for the build machine. For the last few percent, an
optional PGO pass is documented in [docs.md](docs.md#reproducing-benchmark-builds).

## Usage

```bash
# Time the computation only
./target/release/fibonacci-rust 10000000

# Average 5 runs with statistics
./target/release/fibonacci-rust 10000000 -a 5

# Print the full number, or save it to a file (includes decimal output time)
./target/release/fibonacci-rust 1000000 -p
./target/release/fibonacci-rust 100000000 -s fib100m.txt
```

| Flag | Meaning |
|------|---------|
| `n` | index of the Fibonacci number to compute |
| `-a N`, `--average N` | run N times, print mean / 95% CI / stddev / skewness |
| `-p`, `--print` | print the full decimal result |
| `-s FILE`, `--save FILE` | save the decimal result to a file |

## How it works (summary)

1. **Fast doubling** jumps `(F(k), F(k+1)) → (F(2k), F(2k+1))` per bit of `n`
   (~27 big-multiplication steps for `n = 10⁸` instead of a hundred million
   additions).
2. **NTT multiplication** does each big multiply by exact integer convolution
   over three prime moduli, fused so one transform pair yields all three
   products each step needs (`a²`, `b²`, `a·b`).
3. **Decimal limbs** (base 10⁹) make printing zero-padded formatting — there
   is no conversion step.

Full detail: [docs.md](docs.md#how-it-works--detailed-version).

## Correctness

- Full decimal outputs hashed against Python's exact integers: F(10⁵),
  F(10⁶), F(10⁷), F(10⁸) — identical (`AEF6…`, `DEE6…`, `098AC40F…`).
- 17 `cargo test` unit tests (NTT round-trips, modular-arithmetic oracles,
  known values, doubling identities).
- Inputs past the exact NTT size limit return an error instead of garbage.

## Limitations

- Exact for `n` up to ~1.8×10⁸ (NTT transform cap 2²³); larger `n` is refused
  with an error, not silently miscomputed.
- Memory scales with output size (a few hundred MB peak at `n = 10⁸`).
- The binary is CPU-specific (`target-cpu=native`); rebuild without it for
  portable binaries.

## Layout

```text
├── Cargo.toml          (clap, rayon, thiserror, anyhow, mimalloc)
└── src/
    ├── main.rs         (CLI, stats, orchestration)
    ├── fib.rs          (fast-doubling loop + tests)
    ├── bigint.rs       (limbs, add/sub/double, schoolbook, garner, fib_double)
    ├── ntt.rs          (primes, Barrett mulmod, layered roots, transforms)
    └── decimal.rs      (parallel 9-digit formatting)
```

## License

MIT (see `license` in `Cargo.toml`).
