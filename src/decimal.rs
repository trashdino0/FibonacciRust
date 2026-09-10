//! Decimal output: limbs already ARE base-10^9 digits, so printing is just
//! zero-padded formatting — no conversion, no division tree, no BigUint.

use crate::bigint::BigInt;
use rayon::prelude::*;

/// Below this many lower limbs, format serially (rayon overhead isn't worth
/// it for a few hundred groups).
const PAR_CHUNK_MIN_LIMBS: usize = 256;

/// Append exactly 9 digits with leading zeros (manual loop: no `format!`
/// overhead per limb, no `unsafe`, always valid ASCII).
fn push_9(out: &mut String, mut w: u32) {
    let mut buf = [b'0'; 9];
    for i in (0..9).rev() {
        buf[i] = b'0' + (w % 10) as u8;
        w /= 10;
    }
    out.push_str(std::str::from_utf8(&buf).expect("digits are ASCII"));
}

/// Decimal string: top limb plain, every lower limb exactly 9 digits.
/// Lower-limb blocks format in parallel; order is preserved on concat.
pub fn to_decimal_string(v: &BigInt) -> String {
    if v.is_zero() {
        return String::from("0");
    }
    let limbs = v.limbs();
    let top = limbs[limbs.len() - 1].to_string();
    let rest = &limbs[..limbs.len() - 1];
    if rest.is_empty() {
        return top;
    }
    let mut out = String::with_capacity(top.len() + rest.len() * 9);
    out.push_str(&top);
    if rest.len() < PAR_CHUNK_MIN_LIMBS {
        for &w in rest.iter().rev() {
            push_9(&mut out, w);
        }
        return out;
    }
    let n_threads = rayon::current_num_threads().max(1);
    let chunk = rest.len().div_ceil(n_threads);
    let parts: Vec<String> = rest
        .par_chunks(chunk)
        .map(|c| {
            let mut s = String::with_capacity(c.len() * 9);
            for &w in c.iter().rev() {
                push_9(&mut s, w);
            }
            s
        })
        .collect();
    for p in parts.iter().rev() {
        out.push_str(p);
    }
    out
}

/// Fast digit count without building the string.
pub fn digit_count(v: &BigInt) -> usize {
    v.decimal_digits()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_values_match_known_strings() {
        for (n, want) in [
            (0u64, "0"),
            (1, "1"),
            (2, "1"),
            (10, "55"),
            (100, "354224848179261915075"),
        ] {
            let v = crate::fib::compute_fib(n).unwrap();
            assert_eq!(to_decimal_string(&v), want, "n={n}");
        }
    }

    #[test]
    fn padding_is_correct() {
        let v = crate::fib::compute_fib(100).unwrap();
        assert_eq!(to_decimal_string(&v), "354224848179261915075");
        assert_eq!(digit_count(&v), 21);
    }

    #[test]
    fn parallel_chunks_reassemble() {
        // F(20000) has 465 limbs (parallel path): length matches the digit
        // count, every non-top 9-char group parses below 10^9, and the top
        // group has no leading zeros.
        let v = crate::fib::compute_fib(20_000).unwrap();
        let s = to_decimal_string(&v);
        assert_eq!(s.len(), digit_count(&v));
        let limbs = v.limbs();
        let top_len = limbs[limbs.len() - 1].to_string().len();
        assert_eq!(s.len(), top_len + 9 * (limbs.len() - 1));
        assert!(!s.starts_with('0'));
        for chunk in s.as_bytes()[top_len..].chunks(9) {
            assert_eq!(chunk.len(), 9);
            let g = std::str::from_utf8(chunk).unwrap().parse::<u32>().unwrap();
            assert!(g < 1_000_000_000);
        }
    }
}
