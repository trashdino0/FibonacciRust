//! Three-prime Number-Theoretic Transform (NTT), ported from `NTT.java`.
//!
//! Multiplication is convolution via NTT over three NTT-friendly primes,
//! each < 2^30 so pairwise products fit in a `u64`:
//!
//! | Prime | Value     | Factorisation   | Max transform |
//! |-------|-----------|-----------------|---------------|
//! | P1    | 998244353 | 119 * 2^23 + 1  | 2^23          |
//! | P2    | 167772161 | 5 * 2^25 + 1    | 2^25          |
//! | P3    | 469762049 | 7 * 2^26 + 1    | 2^26          |
//!
//! All three share primitive root 3. The bottleneck is P1 (2^23), so results
//! are exact for NTT sizes `n <= 2^23`, i.e. roughly `F(n)` with `n <= ~190M`.
//! Callers must enforce that limit (see [`MAX_LOG_N`]); Java silently wraps
//! past it, we return an error instead.
//!
//! Optimisations carried over from the Java version, then improved:
//! * Integer Barrett reduction replaces every `% p` in the butterfly
//!   (hardware `div` is ~20-40 cycles, this is a multiply-high plus fix-up).
//!   Java used a `double`-based quotient estimate; we use a pure-integer one
//!   (`q̂ = (x·μ) >> 64`, `μ = ⌈2^64/p⌉`) to avoid the `f64`→`u64` conversion
//!   latency chain on the hottest loop.
//! * Layered root-of-unity tables: the inner loop reads
//!   `roots[len], roots[len+1], ...` sequentially (prefetcher-friendly) and
//!   the per-butterfly `w = w * wlen % p` update disappears entirely.
//! * Tables are built once and shared via `Arc` (Java used a
//!   `ConcurrentHashMap`); cloning the `Arc` is pointer-sized, the 64 MiB
//!   tables themselves are never copied.

/// Bottleneck prime P1 supports transforms up to 2^23.
pub const MAX_LOG_N: u32 = 23;

pub const P1: u64 = 998_244_353;
pub const P2: u64 = 167_772_161;
pub const P3: u64 = 469_762_049;

const G: u64 = 3;

pub const PRIMES: [u64; 3] = [P1, P2, P3];

/// `⌈2^64 / p⌉` multipliers for integer Barrett reduction (one per prime).
const fn barrett_mu(p: u64) -> u64 {
    (1u128 << 64).div_ceil(p as u128) as u64
}

pub const MU_P1: u64 = barrett_mu(P1);
pub const MU_P2: u64 = barrett_mu(P2);
pub const MU_P3: u64 = barrett_mu(P3);

const MU_PRIMES: [u64; 3] = [MU_P1, MU_P2, MU_P3];

/// `1.0 / P3`, kept for the one double-Barrett remainder left in Garner's
/// `partial % P3` (`bigint.rs`), which is exact and off the hot path.
pub const INV_P3: f64 = 1.0 / P3 as f64;

use std::collections::HashMap;
use std::sync::{Arc, LazyLock, RwLock};

static FWD_CACHE: LazyLock<RwLock<HashMap<u32, Arc<Vec<u64>>>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));
static INV_CACHE: LazyLock<RwLock<HashMap<u32, Arc<Vec<u64>>>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));
static NINV_CACHE: LazyLock<RwLock<HashMap<u32, u64>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));

#[inline(always)]
fn cache_key(log_n: u32, pi: usize) -> u32 {
    log_n * 3 + pi as u32
}

/// `(a * b) mod p` via integer Barrett reduction.
///
/// Preconditions: `a, b < 2^32` (NOT `< p` — the forward butterfly legitimately
/// holds lazily-reduced sums up to `2^32 - 1`; see below), so `a*b < 2^64`
/// fits in a `u64` with no 128-bit division. With `μ = ⌈2^64/p⌉`, the estimate
/// `q̂ = ⌊x·μ/2^64⌋` is within ±1 of `⌊x/p⌋` for all `x < 2^64`
/// (Möller–Granlund "Improved division by invariant integers"), so a single
/// conditional add/subtract corrects it. Pure integer arithmetic — no `f64`
/// quotient, no `cvttsd2si` latency chain — while staying exact, which is also
/// what fixes the `NTT.java` signed-`long`-overflow bug by construction
/// (see the write-up).
///
/// Laziness argument (why sums `< 2^32` are fine): every butterfly output is
/// congruent mod `p` to the true value (single subtract preserves congruence,
/// differences are exact), and `mulmod` is exact for any inputs `< 2^32`, so
/// by induction every stage output — and the final convolution — is exact.
#[inline(always)]
pub fn mulmod(a: u64, b: u64, p: u64, mu: u64) -> u64 {
    debug_assert!(
        a < (1u64 << 32) && b < (1u64 << 32),
        "mulmod inputs out of range"
    );
    let x = a * b;
    let q = ((x as u128 * mu as u128) >> 64) as u64;
    // `q·p` can exceed `u64` by up to `p` when `x` is near `2^64`, so compare
    // in `u128` — a single widening multiply, still no division.
    let qp = q as u128 * p as u128;
    let x128 = x as u128;
    if qp <= x128 {
        let r = (x128 - qp) as u64;
        if r >= p {
            r - p
        } else {
            r
        }
    } else {
        // `q` overestimated by (at most) 1: `qp - x <= p`, so `p - (qp - x)`
        // is the true remainder with no underflow.
        let d = qp - x128;
        debug_assert!(d <= p as u128, "Barrett estimate off by more than 1");
        (p as u128 - d) as u64
    }
}

pub fn powmod(mut base: u64, mut exp: u64, modu: u64) -> u64 {
    let mut res = 1u64;
    base %= modu;
    while exp > 0 {
        if exp & 1 == 1 {
            res = (res * base) % modu;
        }
        base = (base * base) % modu;
        exp >>= 1;
    }
    res
}

pub fn mod_inverse(n: u64, modu: u64) -> u64 {
    powmod(n, modu - 2, modu)
}

/// Layered root table for a DIT NTT of size `2^log_n`.
///
/// Layout (`i` in `[1, 2^log_n)`): `roots[1] = 1`, and for each
/// `k = 2, 4, 8, ..., n/2`: `roots[k] = roots[k/2]`,
/// `roots[k+j] = roots[k+j-1] * e` for `j = 1..k`, where
/// `e = g^((p-1)/(2k))` is the primitive `(2k)`-th root of unity.
/// Inverse tables use `g^-1` as generator.
fn build_roots(log_n: u32, pi: usize, inverse: bool) -> Vec<u64> {
    let n = 1usize << log_n;
    let p = PRIMES[pi];
    let mut g = G;
    if inverse {
        g = powmod(g, p - 2, p);
    }
    let mu = MU_PRIMES[pi];

    let mut rt = vec![0u64; n.max(2)];
    rt[1] = 1;
    let mut k = 2usize;
    while k < n {
        let e = powmod(g, (p - 1) / (2 * k as u64), p);
        for i in k..2 * k {
            rt[i] = if i & 1 == 0 {
                rt[i >> 1]
            } else {
                mulmod(rt[i - 1], e, p, mu)
            };
        }
        k <<= 1;
    }
    rt
}

/// Forward root table for size `2^log_n`, prime `pi` (shared `Arc`, cheap to clone).
pub fn fwd_roots(log_n: u32, pi: usize) -> Arc<Vec<u64>> {
    let key = cache_key(log_n, pi);
    if let Some(hit) = FWD_CACHE.read().expect("lock").get(&key) {
        return Arc::clone(hit);
    }
    let mut cache = FWD_CACHE.write().expect("lock");
    Arc::clone(
        cache
            .entry(key)
            .or_insert_with(|| Arc::new(build_roots(log_n, pi, false))),
    )
}

/// Inverse root table (built from `g^-1`).
pub fn inv_roots(log_n: u32, pi: usize) -> Arc<Vec<u64>> {
    let key = cache_key(log_n, pi);
    if let Some(hit) = INV_CACHE.read().expect("lock").get(&key) {
        return Arc::clone(hit);
    }
    let mut cache = INV_CACHE.write().expect("lock");
    Arc::clone(
        cache
            .entry(key)
            .or_insert_with(|| Arc::new(build_roots(log_n, pi, true))),
    )
}

/// `n^-1 mod p`, cached.
pub fn n_inv(log_n: u32, pi: usize) -> u64 {
    let key = cache_key(log_n, pi);
    if let Some(&hit) = NINV_CACHE.read().expect("lock").get(&key) {
        return hit;
    }
    let mut cache = NINV_CACHE.write().expect("lock");
    *cache
        .entry(key)
        .or_insert_with(|| powmod(1u64 << log_n, PRIMES[pi] - 2, PRIMES[pi]))
}

/// Pre-warm root tables for common sizes (mirrors the Java static block that
/// warmed `logN 1..=20` for all 3 primes synchronously at startup).
pub fn prewarm() {
    for log_n in 1..=20u32 {
        for pi in 0..3 {
            let _ = fwd_roots(log_n, pi);
            let _ = inv_roots(log_n, pi);
            let _ = n_inv(log_n, pi);
        }
    }
}

/// Bit-reversal permutation, amortised O(n) via the binary-counter trick.
fn bit_rev(a: &mut [u64]) {
    let n = a.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
}

/// Forward in-place NTT. `data` must have length exactly `2^log_n`.
pub fn ntt_forward(data: &mut [u64], log_n: u32, pi: usize) {
    if log_n == 0 {
        return;
    }
    debug_assert_eq!(data.len(), 1usize << log_n);
    let p = PRIMES[pi];
    let mu = MU_PRIMES[pi];
    let rt = fwd_roots(log_n, pi);
    let n = 1usize << log_n;

    bit_rev(data);

    let mut len = 1usize;
    while len < n {
        let two = len << 1;
        let mut i = 0usize;
        while i < n {
            for j in 0..len {
                let u = data[i + j];
                let v = mulmod(data[i + j + len], rt[len + j], p, mu);
                let sum = u + v;
                data[i + j] = if sum >= p { sum - p } else { sum };
                data[i + j + len] = if u < v { u + p - v } else { u - v };
            }
            i += two;
        }
        len <<= 1;
    }
}

/// Inverse in-place NTT (includes `1/n` scaling).
pub fn ntt_inverse(data: &mut [u64], log_n: u32, pi: usize) {
    if log_n == 0 {
        return;
    }
    debug_assert_eq!(data.len(), 1usize << log_n);
    let p = PRIMES[pi];
    let mu = MU_PRIMES[pi];
    let rt = inv_roots(log_n, pi);
    let ni = n_inv(log_n, pi);
    let n = 1usize << log_n;

    bit_rev(data);

    let mut len = 1usize;
    while len < n {
        let two = len << 1;
        let mut i = 0usize;
        while i < n {
            for j in 0..len {
                let u = data[i + j];
                let v = mulmod(data[i + j + len], rt[len + j], p, mu);
                let sum = u + v;
                data[i + j] = if sum >= p { sum - p } else { sum };
                data[i + j + len] = if u < v { u + p - v } else { u - v };
            }
            i += two;
        }
        len <<= 1;
    }

    for v in data.iter_mut() {
        *v = mulmod(*v, ni, p, mu);
    }
}

/// `ceil(log2(x))` for `x >= 1`: smallest `l` with `2^l >= x`.
#[inline]
pub fn ceil_log2(x: usize) -> u32 {
    assert!(x >= 1, "ceil_log2(0) is undefined");
    usize::BITS - (x - 1).leading_zeros()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mulmod_matches_remainder() {
        // Spot-check across all three primes over the FULL contract range
        // `a, b < 2^32` (forward butterflies hold lazily-reduced sums that
        // big — testing only `< p` would miss the Java overflow regime).
        // Oracle is exact u128 arithmetic.
        for (pi, &p) in PRIMES.iter().enumerate() {
            let mu = MU_PRIMES[pi];
            let mut x = 1u64;
            for _ in 0..20000 {
                // Wrapping LCG (deliberate overflow — hence wrapping ops).
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let a = x & 0xFFFF_FFFF;
                let b = (x >> 11) & 0xFFFF_FFFF;
                let want = ((a as u128 * b as u128) % p as u128) as u64;
                assert_eq!(mulmod(a, b, p, mu), want, "p={p} a={a} b={b}");
            }
            // Edges: max product, near-2^63 wrap boundary, zeros and ones.
            for (a, b) in [
                (0xFFFF_FFFF, 0xFFFF_FFFF),
                (p - 1, p - 1),
                (0xFFFF_FFFF, p - 1),
                (0, p - 1),
                (1, 1),
                (3_037_000_499, 3_037_000_499), // just under sqrt(2^63)
                (3_037_000_500, 3_037_000_500), // just over sqrt(2^63)
            ] {
                let want = ((a as u128 * b as u128) % p as u128) as u64;
                assert_eq!(mulmod(a, b, p, mu), want, "p={p} a={a} b={b}");
            }
        }
    }

    #[test]
    fn table_forward_then_inverse_is_identity() {
        for log_n in [1u32, 2, 4, 8, 12] {
            for (pi, &p) in PRIMES.iter().enumerate() {
                let n = 1usize << log_n;
                let mut v: Vec<u64> = (0..n as u64).map(|i| (i * 7 + 1) % p).collect();
                let orig = v.clone();
                ntt_forward(&mut v, log_n, pi);
                assert_ne!(v, orig, "forward should change data");
                ntt_inverse(&mut v, log_n, pi);
                assert_eq!(v, orig, "log_n={log_n} pi={pi}");
            }
        }
    }

    #[test]
    fn ceil_log2_cases() {
        let cases = [
            (1usize, 0u32),
            (2, 1),
            (3, 2),
            (4, 2),
            (5, 3),
            (1024, 10),
            (1025, 11),
        ];
        for (input, expected) in cases {
            assert_eq!(ceil_log2(input), expected, "input={input}");
        }
    }
}
