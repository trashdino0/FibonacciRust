//! `fibonacci-rust` — huge Fibonacci numbers via fast doubling + 3-prime NTT.
//!
//! CLI mirrors `fibonacciV2` (`HugeFibonacciMutable`):
//! ```text
//! fibonacci-rust [OPTIONS] <N>
//!   -a, --average <N>     run N times and print statistics
//!   -w, --warmup [R,N]    warmup RUNS iterations of F(N) first (default 50,10000)
//!   -p, --print           print the full decimal result
//!   -s, --save <FILE>     save the decimal result to a file
//! ```

mod bigint;
mod decimal;
mod fib;
mod ntt;

use anyhow::{Context, Result};
use clap::Parser;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(
    name = "fibonacci-rust",
    version,
    about = "Huge Fibonacci via fast doubling + NTT"
)]
struct Args {
    /// Index of the Fibonacci number to calculate.
    n: u64,

    /// Number of runs for averaging (default = 1).
    #[arg(short = 'a', long = "average", default_value_t = 1)]
    runs: u32,

    /// Warmup: RUNS,N (default when bare flag: 50,10000).
    #[arg(short = 'w', long = "warmup", num_args = 0..=1, default_missing_value = "50,10000")]
    warmup: Option<String>,

    /// Print F(n) to the console after calculation.
    #[arg(short = 'p', long = "print")]
    print: bool,

    /// Save F(n) as decimal text to the given file.
    #[arg(short = 's', long = "save")]
    save: Option<PathBuf>,
}

/// Two-sided 95% t critical values, df = 1..=30, then normal 1.96.
/// (Matches Java's `TDistribution.inverseCumulativeProbability(0.975)`.)
fn t_crit_95(df: usize) -> f64 {
    const T: [f64; 30] = [
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160,
        2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056,
        2.052, 2.048, 2.045, 2.042,
    ];
    if df == 0 {
        f64::NAN
    } else if df <= 30 {
        T[df - 1]
    } else {
        1.96
    }
}

fn parse_warmup(s: &str) -> Result<(u32, u64)> {
    let mut it = s.split(',');
    let runs: u32 = it
        .next()
        .unwrap_or("50")
        .trim()
        .parse()
        .context("warmup RUNS must be an integer")?;
    let n: u64 = it
        .next()
        .unwrap_or("10000")
        .trim()
        .parse()
        .context("warmup N must be an integer")?;
    Ok((runs, n))
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Pre-warm NTT root tables (like the Java static block) so the timed
    // runs don't pay table-construction costs.
    ntt::prewarm();

    if let Some(ref w) = args.warmup {
        let (warmup_runs, warmup_n) = parse_warmup(w)?;
        println!("Warming up: {warmup_runs} x F({warmup_n})");
        for _ in 0..warmup_runs {
            let _ = fib::compute_fib(warmup_n)?;
        }
        println!("Warmup complete.");
    }

    println!("Calculating F({}) with Parallel NTT (Rust)...", args.n);

    let mut times: Vec<f64> = Vec::with_capacity(args.runs.max(1) as usize);
    let mut last: Option<bigint::BigInt> = None;
    for x in 0..args.runs.max(1) {
        let start = Instant::now();
        let result = fib::compute_fib(args.n)?;
        let secs = start.elapsed().as_secs_f64();
        times.push(secs);
        last = Some(result);
        if args.runs > 1 {
            println!("  Run {}: {:.4} s", x + 1, secs);
        }
    }

    let mean = times.iter().sum::<f64>() / times.len() as f64;
    if times.len() > 1 {
        let n = times.len() as f64;
        let var = times.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / (n - 1.0);
        let sd = var.sqrt();
        let crit = t_crit_95(times.len() - 1);
        let moe = crit * (sd / n.sqrt());
        // Skewness (moment estimator, like Commons Math).
        let skew = if sd > 0.0 {
            times.iter().map(|t| ((t - mean) / sd).powi(3)).sum::<f64>() / n
        } else {
            0.0
        };
        println!("\n--- Statistical Analysis ---");
        println!("Runs:     {}", times.len());
        println!("Mean:     {mean:.4} s");
        println!("95% CI:   [{:.4} s, {:.4} s]", mean - moe, mean + moe);
        println!("StdDev:   {sd:.4} s");
        println!("Skewness: {skew:.4}");
    } else {
        println!("F({}) calculated in {:.4} s", args.n, mean);
    }

    if let Some(v) = last {
        if args.print || args.save.is_some() {
            println!("Digits:   ~{} (estimate)", decimal::digit_count(&v));
            print!("Converting to decimal...");
            let cs = Instant::now();
            let s = decimal::to_decimal_string(&v);
            println!(" ({:.3} s)", cs.elapsed().as_secs_f64());
            println!("Digits:   {}", s.len());

            if args.print {
                println!("{s}");
            }
            if let Some(path) = args.save {
                fs::write(&path, &s)
                    .with_context(|| format!("failed to write {}", path.display()))?;
                println!("Saved to {}", path.display());
            }
        }
    }

    Ok(())
}
