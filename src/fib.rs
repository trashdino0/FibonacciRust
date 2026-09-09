//! Fast-doubling Fibonacci driver, ported from `HugeFibonacciMutable.computeFib`.
//!
//! Identities (O(log n) doublings):
//! * `F(2k)   = F(k) * (2*F(k+1) - F(k))`
//! * `F(2k+1) = F(k)^2 + F(k+1)^2`
//!
//! The loop walks the bits of `n` from MSB to LSB, seeded from a 4-bit lookup
//! table so the first iterations are skipped.

use crate::bigint::{BigInt, FibError, NTT_THRESHOLD};

/// `FIB_BASE[i] = F(i)` for `i = 0..=16` — covers any 4-bit seed `0..=15`
/// plus its successor.
const FIB_BASE: [u64; 17] = [
    0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
];

/// `log2(phi)` — `F(n)` has about `n * LOG2_PHI` bits.
const LOG2_PHI: f64 = 0.694_241_913_630_617_8;

/// Compute `F(target)`.
///
/// Small operands use the schoolbook doubling step, large ones the 3-product
/// NTT step. Returns [`FibError::TooLarge`] past the exact NTT limit instead
/// of Java's silent wrap-around.
pub fn compute_fib(target: u64) -> Result<BigInt, FibError> {
    if target < FIB_BASE.len() as u64 {
        return Ok(BigInt::from_u64(FIB_BASE[target as usize], 2));
    }

    // Tight capacity story (a Java bug fixed here): Java estimated
    // `estLimbs = n*0.022 + 32` then over-allocated `bufSize ~= 4x` that
    // estimate "for headroom". We size every intermediate exactly
    // (`next_pow2(2*max_len) + 2`) inside the doubling steps, so no global
    // buffer — and no 2.4x profligate overallocation — is needed at all.
    let _estimated_bits = (target as f64 * LOG2_PHI).ceil() as usize + 64;

    let msb = 63 - target.leading_zeros();
    let prefix_bits = (msb + 1).min(4);
    let prefix_val = (target >> (msb + 1 - prefix_bits)) as usize;

    let mut a = BigInt::from_u64(FIB_BASE[prefix_val], 2);
    let mut b = BigInt::from_u64(FIB_BASE[prefix_val + 1], 2);

    for i in (0..=(msb - prefix_bits)).rev() {
        let (na, nb) = if a.len() < NTT_THRESHOLD && b.len() < NTT_THRESHOLD {
            BigInt::fib_double_small(&a, &b)
        } else {
            BigInt::fib_double(&a, &b).map_err(|e| match e {
                FibError::TooLarge(_, log_n) => FibError::TooLarge(target, log_n),
            })?
        };
        // Bit set: advance (F(k), F(k+1)) -> (F(k+1), F(k)+F(k+1)).
        // Reuse `na`'s allocation for the sum — no fresh alloc.
        if (target >> i) & 1 == 1 {
            let mut next_b = na;
            next_b.add_assign(&nb);
            a = nb;
            b = next_b;
        } else {
            a = na;
            b = nb;
        }
    }
    Ok(a)
}

#[cfg(test)]
mod tests {
    use super::*;

    // First 31 Fibonacci numbers — table-driven, one test, failure names the case.
    const SMALL: [u64; 31] = [
        0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765,
        10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040,
    ];

    #[test]
    fn first_31_match_table() {
        for (n, &want) in SMALL.iter().enumerate() {
            let got = compute_fib(n as u64).unwrap();
            let mut limbs = Vec::new();
            let mut w = want;
            while w > 0 {
                limbs.push(w as u32);
                w >>= 32;
            }
            assert_eq!(got.limbs(), &limbs, "F({n})");
        }
    }

    #[test]
    fn f100_matches_known_value() {
        let got = compute_fib(100).unwrap();
        let want: Vec<u32> = num_bigint::BigUint::parse_bytes(b"354224848179261915075", 10)
            .unwrap()
            .to_u32_digits();
        assert_eq!(got.limbs(), &want);
    }

    #[test]
    fn f1000_has_209_digits() {
        let got = compute_fib(1000).unwrap();
        let s = crate::decimal::to_decimal_string(&got);
        assert_eq!(s.len(), 209, "F(1000) digit count");
        assert_eq!(
            s,
            "43466557686937456435688527675040625802564660517371780402481729089536555417949051890403879840079255169295922593080322634775209689623239873322471161642996440906533187938298969649928516003704476137795166849228875"
        );
    }

    #[test]
    fn doubling_identity_holds() {
        // For several k: fib_double(F(k),F(k+1)) == (F(2k),F(2k+1)).
        for k in [2u64, 3, 5, 50, 500, 5000] {
            let (fk, fk1) = (compute_fib(k).unwrap(), compute_fib(k + 1).unwrap());
            let (f2k, f2k1) = if fk.len() < NTT_THRESHOLD && fk1.len() < NTT_THRESHOLD {
                BigInt::fib_double_small(&fk, &fk1)
            } else {
                BigInt::fib_double(&fk, &fk1).unwrap()
            };
            assert_eq!(f2k.limbs(), compute_fib(2 * k).unwrap().limbs(), "F(2*{k})");
            assert_eq!(
                f2k1.limbs(),
                compute_fib(2 * k + 1).unwrap().limbs(),
                "F(2*{k}+1)"
            );
        }
    }
}
