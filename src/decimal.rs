//! Parallel decimal conversion, ported from `ToStringTask`.
//!
//! A 2M-digit number cannot be stringified by repeated `divmod(10)` — that is
//! O(digits^2). Instead we split at the decimal midpoint:
//!
//! ```text
//! half    = digits / 2
//! divisor = 10^half
//! (high, low) = n.div_rem(divisor)   // one big division
//! string  = to_string(high) + zero_padded_to_half(to_string(low))
//! ```
//!
//! Both halves recurse **in parallel** via `rayon::join`.
//! Below the threshold bits we fall through to `to_string()` directly.

use crate::bigint::BigInt;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::Zero;
use std::collections::HashMap;
use std::sync::{LazyLock, RwLock};

/// Below this bit length, convert directly without splitting.
/// Single tuned threshold (30k bits ≈ 9k digits).
pub const TO_STRING_THRESHOLD_BITS: usize = 30_000;

static POW10_CACHE: LazyLock<RwLock<HashMap<u32, BigUint>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));

fn pow10(exp: u32) -> BigUint {
    if let Some(hit) = POW10_CACHE.read().expect("lock").get(&exp) {
        return hit.clone();
    }
    let mut cache = POW10_CACHE.write().expect("lock");
    cache
        .entry(exp)
        .or_insert_with(|| BigUint::from(10u32).pow(exp))
        .clone()
}

/// `num-bigint` uses `u32` LE digits like our limbs: one safe clone.
pub fn to_biguint(v: &BigInt) -> BigUint {
    if v.is_zero() {
        BigUint::zero()
    } else {
        BigUint::new(v.limbs().to_vec())
    }
}

fn par_to_string(n: &BigUint) -> String {
    if n.bits() < TO_STRING_THRESHOLD_BITS as u64 {
        return n.to_string();
    }
    // digits ~= bits * log10(2).
    let digits = (n.bits() as usize * 30_103) / 100_000 + 1;
    let half = digits / 2;
    let divisor = pow10(half as u32);
    let (high, low) = n.div_rem(&divisor);
    let (left, right) = rayon::join(|| par_to_string(&high), || par_to_string(&low));
    let right = if right.len() < half {
        let mut padded = String::with_capacity(half);
        padded.push_str(&"0".repeat(half - right.len()));
        padded.push_str(&right);
        padded
    } else {
        right
    };
    let mut out = String::with_capacity(left.len() + right.len());
    out.push_str(&left);
    out.push_str(&right);
    out
}

/// Decimal string for a computed Fibonacci number (parallel for huge sizes).
pub fn to_decimal_string(v: &BigInt) -> String {
    par_to_string(&to_biguint(v))
}

/// Fast digit count without building the string.
pub fn digit_count(v: &BigInt) -> usize {
    if v.is_zero() {
        return 1;
    }
    (v.bit_length() as f64 * std::f64::consts::LOG10_2).floor() as usize + 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_values_match_to_string() {
        for n in [0u64, 1, 2, 10, 100, 1000] {
            let v = crate::fib::compute_fib(n).unwrap();
            assert_eq!(to_decimal_string(&v), to_biguint(&v).to_string(), "n={n}");
        }
    }

    #[test]
    fn padding_is_correct() {
        let v = crate::fib::compute_fib(100).unwrap();
        assert_eq!(to_decimal_string(&v), "354224848179261915075");
        assert_eq!(digit_count(&v), 21);
    }
}
