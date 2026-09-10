//! Little-endian limb big-int + NTT multiplication + Garner CRT.
//!
//! Limbs are base-10^9 digits: printing is then just zero-padded formatting,
//! so binary-to-decimal conversion vanishes. NTT inputs only need values
//! `< 2^32` (10^9 qualifies), and max convolution coefficients stay far
//! below CRT capacity — see `garner`.
//!
//! Differences versus the Java original it was ported from:
//! * No buffer pool: owned values move in and out of Rayon closures, so
//!   every buffer is zeroed or filled before it is read, by construction.
//! * No `Unsafe`: our limb layout already matches `num-bigint`'s.
//! * Oversize inputs are an error, not silent garbage (see [`FibError`]).

use crate::ntt;
use thiserror::Error;

#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum FibError {
    #[error("n={0} needs NTT size 2^{1}, exceeding the exact limit 2^23 (P1 bottleneck, ~F(190M)); refusing instead of returning garbage")]
    TooLarge(u64, u32),
}

/// Below this many limbs (64 limbs ~= 2048 bits) schoolbook O(n^2) beats NTT.
pub const NTT_THRESHOLD: usize = 64;
/// Below this NTT size, Rayon overhead exceeds the parallel gain.
pub const PARALLEL_NTT_MIN_N: usize = 4096;

// Precomputed CRT constants for the 3-modulus Garner.
fn p1p2() -> u64 {
    ntt::P1 * ntt::P2
}
fn p1_inv_p2() -> u64 {
    ntt::mod_inverse(ntt::P1, ntt::P2)
}
fn p12_inv_p3() -> u64 {
    ntt::mod_inverse(p1p2() % ntt::P3, ntt::P3)
}

/// Limb base: each limb holds 9 decimal digits (`< 1_000_000_000`).
pub const LIMB_BASE: u64 = 1_000_000_000;

/// Unsigned little-endian big-int: `mag[0]` holds the least-significant
/// 9 digits, `len` counts used limbs (`mag[len..]` is always zero).
#[derive(Debug, Clone)]
pub struct BigInt {
    mag: Vec<u32>,
    len: usize,
}

/// Barrett multiplier for `÷ 10^9`: `⌈2^64 / 10^9⌉`.
const MU_1E9: u64 = (u64::MAX / LIMB_BASE) + 1;

/// Exact `(quo, rem)` of `x ÷ 10^9` for `x < 2^64`, integer-Barrett style
/// (one multiply-high + fix-ups, no hardware `div`): `q̂` lands within ±1.
#[inline(always)]
fn divmod_1e9_u64(x: u64) -> (u64, u32) {
    let q = ((x as u128 * MU_1E9 as u128) >> 64) as u64;
    let qp = q as u128 * LIMB_BASE as u128;
    let x128 = x as u128;
    if qp <= x128 {
        let r = (x128 - qp) as u64;
        if r >= LIMB_BASE {
            (q + 1, (r - LIMB_BASE) as u32)
        } else {
            (q, r as u32)
        }
    } else {
        let d = qp - x128;
        debug_assert!(d <= LIMB_BASE as u128, "Barrett-1e9 estimate off by > 1");
        (q - 1, (LIMB_BASE as u128 - d) as u32)
    }
}

/// Exact `(quo, rem)` of `x ÷ 10^9` for `x < 2^87` (covers every Garner
/// coefficient and schoolbook word): Knuth 1-word division folding three
/// base-2^32 digits through the Barrett step above.
fn divmod_1e9(x: u128) -> (u64, u32) {
    debug_assert!(x < (1u128 << 87));
    let d2 = (x >> 64) as u64; // < 2^23, contributes no quotient digit
    let d1 = ((x >> 32) & 0xFFFF_FFFF) as u64;
    let d0 = (x & 0xFFFF_FFFF) as u64;
    let (q1, r1) = divmod_1e9_u64((d2 << 32) | d1);
    let (q0, r0) = divmod_1e9_u64((r1 as u64) << 32 | d0);
    ((q1 << 32) | q0, r0)
}

impl BigInt {
    /// Zero with room for `capacity` limbs.
    pub fn zeros(capacity: usize) -> Self {
        Self {
            mag: vec![0u32; capacity.max(1)],
            len: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_zero(&self) -> bool {
        self.len == 0
    }

    /// Store a `u64` as base-10^9 digits (≤ 3 limbs).
    pub fn from_u64(mut value: u64, capacity: usize) -> Self {
        let mut out = Self::zeros(capacity.max(3));
        if value == 0 {
            return out;
        }
        let mut i = 0;
        while value > 0 {
            assert!(i < out.mag.len(), "from_u64: dest too small");
            out.mag[i] = (value % LIMB_BASE) as u32;
            value /= LIMB_BASE;
            i += 1;
        }
        out.len = i;
        out
    }

    fn trim(&mut self) {
        while self.len > 0 && self.mag[self.len - 1] == 0 {
            self.len -= 1;
        }
    }

    /// `self += other`. Caller guarantees capacity; overflow is a programmer
    /// bug, so we assert rather than wrap.
    pub fn add_assign(&mut self, other: &Self) {
        let mut carry = 0u64;
        let max = self.len.max(other.len);
        let mut i = 0usize;
        while i < max || carry != 0 {
            assert!(i < self.mag.len(), "add: overflow {}", self.mag.len());
            let mut sum = carry;
            if i < self.len {
                sum += self.mag[i] as u64;
            }
            if i < other.len {
                sum += other.mag[i] as u64;
            }
            // sum < 2·10^9 + 1: a single conditional subtract normalizes.
            if sum >= LIMB_BASE {
                sum -= LIMB_BASE;
                carry = 1;
            } else {
                carry = 0;
            }
            self.mag[i] = sum as u32;
            if i >= self.len {
                self.len = i + 1;
            }
            i += 1;
        }
    }

    /// `self -= other`. Requires `self >= other`.
    pub fn sub_assign(&mut self, other: &Self) {
        debug_assert!(self.ge_mag(other), "sub_assign requires self >= other");
        let mut borrow = 0i64;
        for i in 0..self.len {
            let mut diff = self.mag[i] as i64 - borrow;
            if i < other.len {
                diff -= other.mag[i] as i64;
            }
            if diff < 0 {
                diff += LIMB_BASE as i64;
                borrow = 1;
            } else {
                borrow = 0;
            }
            self.mag[i] = diff as u32;
        }
        debug_assert_eq!(borrow, 0, "sub_assign underflowed: self < other");
        self.trim();
    }

    fn ge_mag(&self, other: &Self) -> bool {
        if self.len != other.len {
            return self.len > other.len;
        }
        for i in (0..self.len).rev() {
            if self.mag[i] != other.mag[i] {
                return self.mag[i] > other.mag[i];
            }
        }
        true
    }

    /// `self *= 2` (limbs are decimal, so this is doubling, not a shift).
    pub fn mul2(&mut self) {
        if self.len == 0 {
            return;
        }
        assert!(self.len < self.mag.len(), "mul2: overflow");
        let mut carry = 0u64;
        for i in 0..self.len {
            let v = self.mag[i] as u64 * 2 + carry;
            if v >= LIMB_BASE {
                self.mag[i] = (v - LIMB_BASE) as u32;
                carry = 1;
            } else {
                self.mag[i] = v as u32;
                carry = 0;
            }
        }
        if carry != 0 {
            self.mag[self.len] = carry as u32;
            self.len += 1;
        }
    }

    /// Read-only limb view (base-10^9 digits, top nonzero).
    pub fn limbs(&self) -> &[u32] {
        &self.mag[..self.len]
    }

    /// Exact decimal digit count.
    pub fn decimal_digits(&self) -> usize {
        if self.len == 0 {
            return 1;
        }
        let mut top = self.mag[self.len - 1];
        let mut digits = 0u32;
        while top > 0 {
            digits += 1;
            top /= 10;
        }
        (self.len - 1) * 9 + digits as usize
    }

    // ------------------------------------------------------------------
    // Multiplication
    // ------------------------------------------------------------------

    /// Schoolbook O(n^2) multiply. Wide `u128` row accumulators stay exact
    /// (each holds `< 64·(10^9)^2 < 2^76` at threshold sizes); one normalize
    /// pass at the end, keeping 2 spare limbs for what follows the multiply.
    fn mul_schoolbook(a: &Self, b: &Self) -> Self {
        if a.len == 0 || b.len == 0 {
            return Self::zeros(1);
        }
        let base = LIMB_BASE as u128;
        let mut acc = vec![0u128; a.len + b.len];
        for (i, &ai) in a.mag[..a.len].iter().enumerate() {
            for (j, &bj) in b.mag[..b.len].iter().enumerate() {
                acc[i + j] += ai as u128 * bj as u128;
            }
        }
        let mut mag = vec![0u32; a.len + b.len + 2];
        let mut carry = 0u128;
        for (i, &v) in acc.iter().enumerate() {
            let cur = v + carry;
            mag[i] = (cur % base) as u32;
            carry = cur / base;
        }
        let mut out = Self {
            mag,
            len: acc.len(),
        };
        while carry > 0 {
            assert!(out.len < out.mag.len(), "schoolbook carry overflow");
            out.mag[out.len] = (carry % base) as u32;
            carry /= base;
            out.len += 1;
        }
        out.trim();
        out
    }

    /// Copy used limbs into a zeroed `dst` (one `memcpy` — same element type).
    fn fill_u32(src: &Self, dst: &mut [u32]) {
        dst[..src.len].copy_from_slice(&src.mag[..src.len]);
    }

    /// Garner CRT: combine residues `r1/r2/r3` (length `n`) into base-10^9
    /// limbs. `v1`/`v2`/`v3` reconstruction is base-independent integer math
    /// (identical to before); only the digit assembly changed: each exact
    /// coefficient (`< 2^87`) plus carry emits one decimal limb via Barrett
    /// divmod, carrying the quotient. Max coefficient at top size is
    /// `n·(10^9)^2 ≈ 2^81.7` (strict `< 10^9` inputs), far below CRT
    /// capacity `P1·P2·P3 ≈ 2^86.3` — the exactness argument is unchanged.
    fn garner(r1: &[u32], r2: &[u32], r3: &[u32], n: usize) -> Self {
        let p1 = ntt::P1;
        let p2 = ntt::P2;
        let p3 = ntt::P3;
        let inv_p3 = ntt::INV_P3;
        let p1_inv_p2 = p1_inv_p2();
        let p12_inv_p3 = p12_inv_p3();
        let p1p2 = p1p2();

        let mut mag = vec![0u32; n + 4];
        let mut carry: u64 = 0;
        for i in 0..n {
            let v1 = r1[i] as u64;
            // `(r2 - r1) mod P2` needs the full reduction: r1 ranges over
            // P1 ≈ 6x P2, so a single conditional add can leave it negative.
            let diff2 = (r2[i] as i64 - v1 as i64).rem_euclid(p2 as i64) as u64;
            let v2 = ntt::mulmod(diff2, p1_inv_p2, p2, ntt::MU_P2);

            let partial = v1 + p1 * v2;
            // `partial % p3` via double-Barrett, corrected twice.
            let q = (partial as f64 * inv_p3) as u64;
            let mut t = partial as i64 - (q as i64) * (p3 as i64);
            if t < 0 {
                t += p3 as i64;
            } else if t as u64 >= p3 {
                t -= p3 as i64;
            }
            if t < 0 {
                t += p3 as i64;
            } else if t as u64 >= p3 {
                t -= p3 as i64;
            }

            let mut diff3 = r3[i] as i64 - t;
            if diff3 < 0 {
                diff3 += p3 as i64;
            }
            let v3 = ntt::mulmod(diff3 as u64, p12_inv_p3, p3, ntt::MU_P3);
            // Full exact coefficient + incoming carry (< 2^87, exact in
            // u128), then one base-10^9 digit out.
            let cur = partial as u128 + p1p2 as u128 * v3 as u128 + carry as u128;
            let (q, r) = divmod_1e9(cur);
            mag[i] = r;
            carry = q;
        }

        let mut ci = n;
        while carry > 0 {
            assert!(ci < mag.len(), "garner carry overflow at {ci}");
            let (q, r) = divmod_1e9(carry as u128);
            mag[ci] = r;
            carry = q;
            ci += 1;
        }

        let mut out = Self {
            len: ci.max(n),
            mag,
        };
        out.trim();
        out
    }

    /// Single-product NTT multiply (test-only; production fuses all three
    /// products into one transform pair per prime in `fib_double`).
    #[cfg(test)]
    fn mul_ntt(a: &Self, b: &Self, log_n: u32, n: usize) -> Self {
        let build = |pi: usize| -> Vec<u32> {
            let mut fa = vec![0u32; n];
            let mut fb = vec![0u32; n];
            Self::fill_u32(a, &mut fa);
            Self::fill_u32(b, &mut fb);
            ntt::ntt_forward(&mut fa, log_n, pi);
            ntt::ntt_forward(&mut fb, log_n, pi);
            let p = ntt::PRIMES[pi];
            let mu = [ntt::MU_P1, ntt::MU_P2, ntt::MU_P3][pi];
            let mut r = vec![0u32; n];
            for i in 0..n {
                r[i] = ntt::mulmod(fa[i] as u64, fb[i] as u64, p, mu) as u32;
            }
            ntt::ntt_inverse(&mut r, log_n, pi);
            r
        };

        let (r1, r2, r3) = if n < PARALLEL_NTT_MIN_N {
            (build(0), build(1), build(2))
        } else {
            let (x, (y, z)) = rayon::join(|| build(0), || rayon::join(|| build(1), || build(2)));
            (x, y, z)
        };
        Self::garner(&r1, &r2, &r3, n)
    }

    // ------------------------------------------------------------------
    // Fast-doubling core: (a, b) = (F(k), F(k+1)) -> (F(2k), F(2k+1))
    // ------------------------------------------------------------------

    /// One prime track of the 3-product double (`a²`, `b²`, `a·b` from one
    /// forward-transform pair). Takes buffers by value, returns residues.
    fn double_track(
        mut fa: Vec<u32>,
        mut fb: Vec<u32>,
        log_n: u32,
        pi: usize,
    ) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
        let p = ntt::PRIMES[pi];
        let mu = [ntt::MU_P1, ntt::MU_P2, ntt::MU_P3][pi];
        let n = fa.len();

        // Fork `forward fb` alongside inline `forward fa`.
        if n >= PARALLEL_NTT_MIN_N {
            rayon::join(
                || ntt::ntt_forward(&mut fa, log_n, pi),
                || ntt::ntt_forward(&mut fb, log_n, pi),
            );
        } else {
            ntt::ntt_forward(&mut fa, log_n, pi);
            ntt::ntt_forward(&mut fb, log_n, pi);
        }

        let mut ra = vec![0u32; n];
        let mut rb = vec![0u32; n];
        let mut rab = vec![0u32; n];
        for i in 0..n {
            let ai = fa[i] as u64;
            let bi = fb[i] as u64;
            ra[i] = ntt::mulmod(ai, ai, p, mu) as u32;
            rb[i] = ntt::mulmod(bi, bi, p, mu) as u32;
            rab[i] = ntt::mulmod(ai, bi, p, mu) as u32;
        }

        // Two forked inverses, one inline.
        if n >= PARALLEL_NTT_MIN_N {
            let ((), ()) = rayon::join(
                || ntt::ntt_inverse(&mut ra, log_n, pi),
                || {
                    let ((), ()) = rayon::join(
                        || ntt::ntt_inverse(&mut rb, log_n, pi),
                        || ntt::ntt_inverse(&mut rab, log_n, pi),
                    );
                },
            );
        } else {
            ntt::ntt_inverse(&mut ra, log_n, pi);
            ntt::ntt_inverse(&mut rb, log_n, pi);
            ntt::ntt_inverse(&mut rab, log_n, pi);
        }
        (ra, rb, rab)
    }

    /// Doubling step reusing one forward-transform pair to get all three
    /// products: `F(2k+1) = A + B`, `F(2k) = 2C - A`
    /// with `A = a^2`, `B = b^2`, `C = a*b`.
    pub fn fib_double(a: &Self, b: &Self) -> Result<(Self, Self), FibError> {
        let max_len = a.len.max(b.len).max(1);
        let raw = 2 * max_len;
        let log_n = ntt::ceil_log2(raw);
        if log_n > ntt::MAX_LOG_N {
            return Err(FibError::TooLarge(u64::MAX, log_n));
        }
        let n = 1usize << log_n;

        // Fresh zeroed buffers per track.
        let mut fa0 = vec![0u32; n];
        let mut fb0 = vec![0u32; n];
        Self::fill_u32(a, &mut fa0);
        Self::fill_u32(b, &mut fb0);
        // One buffer pair per prime.
        let (fa1, fa2) = (fa0.clone(), fa0.clone());
        let (fb1, fb2) = (fb0.clone(), fb0.clone());

        let ((a2_0, b2_0, ab_0), ((a2_1, b2_1, ab_1), (a2_2, b2_2, ab_2))) =
            if n < PARALLEL_NTT_MIN_N {
                (
                    Self::double_track(fa0, fb0, log_n, 0),
                    (
                        Self::double_track(fa1, fb1, log_n, 1),
                        Self::double_track(fa2, fb2, log_n, 2),
                    ),
                )
            } else {
                rayon::join(
                    || Self::double_track(fa0, fb0, log_n, 0),
                    || {
                        rayon::join(
                            || Self::double_track(fa1, fb1, log_n, 1),
                            || Self::double_track(fa2, fb2, log_n, 2),
                        )
                    },
                )
            };

        let a2 = Self::garner(&a2_0, &a2_1, &a2_2, n);
        let b2 = Self::garner(&b2_0, &b2_1, &b2_2, n);
        let ab = Self::garner(&ab_0, &ab_1, &ab_2, n);

        // NLL lets the borrow end before the move; no clone needed.
        // Garner leaves cap `n + 2` at `len <= n + 1`, so the doubling
        // (+1 limb max) and the add (+1 max) always fit.
        let mut out_a = ab;
        out_a.mul2();
        // `2ab >= a^2` for k >= 1 (`b > a`), so no underflow.
        out_a.sub_assign(&a2);
        let mut out_b = a2;
        out_b.add_assign(&b2);
        Ok((out_a, out_b))
    }

    /// Schoolbook doubling step for small operands.
    pub fn fib_double_small(a: &Self, b: &Self) -> (Self, Self) {
        let a2 = Self::mul_schoolbook(a, a);
        let b2 = Self::mul_schoolbook(b, b);
        let mut out_b = a2.clone();
        out_b.add_assign(&b2);
        let mut ab = Self::mul_schoolbook(a, b);
        ab.mul2();
        let mut out_a = ab;
        out_a.sub_assign(&a2);
        (out_a, out_b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small(v: u64) -> BigInt {
        BigInt::from_u64(v, 4)
    }

    #[test]
    fn divmod_1e9_matches_operators() {
        // Full u64 range edges plus pseudo-random u128 values < 2^87
        // (every Garner coefficient and schoolbook word lives there).
        let mut cases = vec![
            0u128,
            1,
            999_999_999,
            1_000_000_000,
            1_000_000_001,
            u64::MAX as u128,
            (1u128 << 87) - 1,
            (1u128 << 86) + 12345,
        ];
        let mut x = 0x1234_5678_9ABCu64;
        for _ in 0..2000 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            cases.push(((x as u128) << 23) ^ (x as u128).wrapping_mul(0x9E3779B9));
        }
        for &v in &cases {
            let v = v % (1u128 << 87);
            let (q, r) = divmod_1e9(v);
            assert_eq!(q as u128, v / 1_000_000_000, "quo of {v}");
            assert_eq!(r as u128, v % 1_000_000_000, "rem of {v}");
        }
    }

    #[test]
    fn add_carry_chain() {
        let mut a = small(999_999_999);
        a.add_assign(&small(1));
        assert_eq!(a.limbs(), &[0, 1]);
    }

    #[test]
    fn sub_borrow_chain() {
        let mut a = BigInt::from_u64(1_000_000_000, 4);
        a.sub_assign(&small(1));
        assert_eq!(a.limbs(), &[999_999_999]);
    }

    #[test]
    fn mul2_carries_out() {
        let mut a = small(600_000_000);
        a.mul2();
        assert_eq!(a.limbs(), &[200_000_000, 1]);
    }

    #[test]
    fn schoolbook_matches_u128() {
        // Deterministic pseudo-random pairs verified against native u128 math
        // (expected value folded back into base-10^9 limbs).
        let mut x = 0x1234_5678_9ABCu64;
        for _ in 0..300 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let y = x ^ 0x9E3779B97F4A7C15;
            let (a, b) = (small(x), small(y));
            let got = BigInt::mul_schoolbook(&a, &b);
            let mut want = (x as u128) * (y as u128);
            let mut limbs = Vec::new();
            while want > 0 {
                limbs.push((want % 1_000_000_000) as u32);
                want /= 1_000_000_000;
            }
            assert_eq!(got.limbs(), &limbs, "x={x} y={y}");
        }
    }

    #[test]
    fn ntt_multiply_matches_schoolbook() {
        // Sizes straddling the NTT threshold, NTT forced via mul_ntt.
        // Inputs are normalized (< 10^9 per limb) like all production values.
        for limbs in [64usize, 100, 300] {
            let mut a = BigInt::zeros(limbs + 2);
            let mut b = BigInt::zeros(limbs + 2);
            for i in 0..limbs {
                a.mag[i] = (i as u32).wrapping_mul(2654435761).wrapping_add(1) % 1_000_000_000;
                b.mag[i] = (i as u32).wrapping_mul(40503).wrapping_add(7) % 1_000_000_000;
            }
            if a.mag[limbs - 1] == 0 {
                a.mag[limbs - 1] = 1;
            }
            if b.mag[limbs - 1] == 0 {
                b.mag[limbs - 1] = 1;
            }
            a.len = limbs;
            b.len = limbs;
            a.trim();
            b.trim();
            let want = BigInt::mul_schoolbook(&a, &b);
            let raw = a.len + b.len;
            let log_n = ntt::ceil_log2(raw);
            let got = BigInt::mul_ntt(&a, &b, log_n, 1 << log_n);
            assert_eq!(got.limbs(), want.limbs(), "limbs={limbs}");
        }
    }

    #[test]
    fn fib_double_matches_identity() {
        // Doubling identity: F(2k)=F(k)(2F(k+1)-F(k)), F(2k+1)=F(k)^2+F(k+1)^2.
        // F(10)=55, F(11)=89 -> F(20)=6765, F(21)=10946.
        let (a, b) = BigInt::fib_double_small(&small(55), &small(89));
        assert_eq!(a.limbs(), small(6765).limbs());
        assert_eq!(b.limbs(), small(10946).limbs());

        // Same via the NTT path.
        let (a, b) = BigInt::fib_double(&small(55), &small(89)).unwrap();
        assert_eq!(a.limbs(), small(6765).limbs());
        assert_eq!(b.limbs(), small(10946).limbs());
    }
}
