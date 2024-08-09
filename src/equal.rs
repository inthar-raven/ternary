use std::ops::{Add, AddAssign};

use nalgebra::{ArrayStorage, RowSVector};
use num_traits::Pow;

use crate::{
    interval::{Dyad, JiRatio},
    ji::odd_limit,
    ji_ratio::RawJiRatio,
    monzo::Monzo,
    primes::{SMALL_PRIMES, SMALL_PRIMES_COUNT},
};

/// In regular temperament theory, a [*val*](https://en.xen.wiki/w/Val) is a covector,
/// i.e. an element of the [dual module](https://en.wikipedia.org/wiki/Dual_module) of
/// interval vectors.
/// (the space of linear maps from the abelian group of intervals to `\mathbb{Z}`).
#[derive(Debug, Copy, Clone, Hash, PartialEq)]
pub struct Val(RowSVector<i32, SMALL_PRIMES_COUNT>);

impl Val {
    /// The zero map.
    pub const ZERO: Self = Self(
        nalgebra::RowSVector::<i32, SMALL_PRIMES_COUNT>::from_array_storage(ArrayStorage(
            [[0i32; 1]; SMALL_PRIMES_COUNT],
        )),
    );
    /// Unwrap the nalgebra row vector.
    pub fn into_inner(&self) -> RowSVector<i32, SMALL_PRIMES_COUNT> {
        self.0
    }
    /// Get a monzo from a slice. Panics if `slice.len() > SMALL_PRIMES_COUNT`.
    pub fn from_slice(slice: &[i32]) -> Self {
        Self(
            nalgebra::RowSVector::<i32, SMALL_PRIMES_COUNT>::from_row_slice(
                &[slice, &vec![0i32; SMALL_PRIMES_COUNT - slice.len()]].concat(),
            ),
        )
    }
    /// Find how many steps a given val maps a monzo (a prime-factorized JI ratio).
    pub fn evaluate(&self, monzo: Monzo) -> i32 {
        (self.0 * monzo.into_inner())[0]
    }
}

impl Add for Val {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Val(self.0 + rhs.0)
    }
}

impl AddAssign for Val {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}

/// The [generalized patent val](https://en.xen.wiki/w/Generalized_patent_val) of x-edo where x is any positive real.
pub fn gpval(edo: f64) -> Val {
    Val(RowSVector::from_iterator(SMALL_PRIMES.into_iter().map(
        |p| {
            f64::round(edo * f64::log(p as f64, std::f64::consts::E) / std::f64::consts::LN_2)
                as i32
        },
    )))
}

/// The best approximation of `ratio` in steps of `ed`ed<`equave`>.
pub fn direct_approx<J: JiRatio>(ratio: J, ed: f64, equave: J) -> i32 {
    f64::round(
        ed * (f64::log(ratio.numer() as f64, std::f64::consts::E)
            - f64::log(ratio.denom() as f64, std::f64::consts::E))
            / (f64::log(equave.numer() as f64, std::f64::consts::E)
                - f64::log(equave.denom() as f64, std::f64::consts::E)),
    ) as i32
}

/// `steps`\`ed`<`equave`> in cents.
pub fn steps_as_cents(steps: i32, ed: f64, equave: RawJiRatio) -> f64 {
    (steps as f64) / ed * equave.cents()
}

/// All integer ed`equave` tunings for `step_sig` scales below `ed_bound`.
pub fn ed_tunings_for_ternary(
    step_sig: &[usize],
    equave: RawJiRatio,
    ed_bound: i32,
    aber_lower: f64,
    aber_upper: f64,
) -> Vec<Vec<i32>> {
    (3..ed_bound)
        .into_iter()
        .flat_map(|l| (2..l).flat_map(move |m| (1..m).map(move |s| vec![l, m, s])))
        .filter(|v| {
            let ed: i32 =
                v[0] * step_sig[0] as i32 + v[1] * step_sig[1] as i32 + v[2] * step_sig[2] as i32;
            let aber_size = steps_as_cents(v[2], ed as f64, equave);
            v[0] * step_sig[0] as i32 + v[1] * step_sig[1] as i32 + v[2] * step_sig[2] as i32
                <= ed_bound
                && aber_lower <= aber_size
                && aber_size <= aber_upper
        })
        .collect()
}

/// Relative error of the patent val mapping for `monzo`.
pub fn relative_error(monzo: Monzo, edo: f64) -> f64 {
    let val = gpval(edo);
    let steps = val.evaluate(monzo);
    return (steps_as_cents(steps, edo, RawJiRatio::OCTAVE) - monzo.cents()) * edo / 1200.0;
}

/// L^1 error on a specified odd limit.
pub fn odd_limit_l1_error(odd: u64, edo: f64) -> f64 {
    odd_limit(odd)
        .into_iter()
        .filter(|&r| r * r < RawJiRatio::OCTAVE)
        .map(|r| Monzo::try_from_ratio(r).unwrap())
        .map(|monzo| relative_error(monzo, edo).abs())
        .sum()
}

/// L^1 error on a specified odd limit.
pub fn odd_limit_l2_error(odd: u64, edo: f64) -> f64 {
    odd_limit(odd)
        .into_iter()
        .filter(|&r| r * r < RawJiRatio::OCTAVE)
        .map(|r| Monzo::try_from_ratio(r).unwrap())
        .map(|monzo| relative_error(monzo, edo).pow(2))
        .sum::<f64>()
        .sqrt()
}

#[macro_export]
/// Creates a `Val`.
macro_rules! val {
    () => (
        $crate::equal::Val::ZERO
    );
    ($elem:expr; $n:expr) => (
        $crate::equal::Val(nalgebra::RowSVector::<i32, {crate::primes::SMALL_PRIMES_COUNT}>::from_row_slice(&[$elem; crate::primes::SMALL_PRIMES_COUNT]))
    );
    ($($x:expr),+ $(,)?) => (
        $crate::equal::Val::from_slice(&[$($x),+])
    );
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    #[allow(unused)]
    use super::*;
    use crate::monzo;

    #[test]
    fn test_rel_error() {
        let five_to_four = monzo![-2, 0, 1];
        let rel_error_in_12edo = relative_error(five_to_four, 12.0);
        assert!(rel_error_in_12edo - 0.1369 < 0.005);
        let rel_error_in_31edo = relative_error(five_to_four, 31.0);
        println!("{}", rel_error_in_31edo);
        assert!(rel_error_in_31edo - 0.020 < 0.005);
    }
    /*
    #[test]
    fn test_error_on_odd_limit() {
        println!("unweighted L^1 error:\n");
        for i in 1..=15 { // odd limits 1, ..., 31
            // test 46edo
            println!("{}-odd-limit:\t31\t{:.2}\t34\t{:.2}", 2*i+1, odd_limit_l1_error(2*i+1, 31.0), odd_limit_l1_error(2*i+1, 34.0));
        }
        println!("unweighted L^1 error:\n");
        for i in 1..=15 { // odd limits 1, ..., 31
            // test 46edo
            println!("{}-odd-limit:\t41\t{:.2}\t46\t{:.2}", 2*i+1, odd_limit_l1_error(2*i+1, 41.0), odd_limit_l1_error(2*i+1, 46.0));
        }
        println!();
        for i in 1..=15 { // odd limits 1, ..., 31
            // test 46edo
            println!("{}-odd-limit:\t53\t{:.2}\t58\t{:.2}", 2*i+1, odd_limit_l1_error(2*i+1, 53.0), odd_limit_l1_error(2*i+1, 58.0));
        }
        println!();
    }
    */
    #[test]
    fn test_ternary_ed_tunings() {
        let blackdye_tunings =
            ed_tunings_for_ternary(&[5, 2, 3], RawJiRatio::OCTAVE, 53, 20.0, 60.0);
        assert_eq!(
            BTreeSet::from_iter(blackdye_tunings.into_iter()),
            BTreeSet::from_iter(
                vec![
                    vec![3, 2, 1],
                    vec![4, 2, 1],
                    vec![4, 3, 1],
                    vec![5, 2, 1],
                    vec![5, 3, 1],
                    vec![5, 4, 1],
                    vec![6, 2, 1],
                    vec![6, 3, 1],
                    vec![6, 4, 1],
                    vec![6, 3, 2],
                    vec![7, 2, 1],
                    vec![6, 5, 1],
                    vec![6, 4, 2],
                    vec![7, 3, 1],
                    vec![6, 5, 2],
                    vec![7, 4, 1],
                    vec![7, 3, 2],
                    vec![8, 2, 1],
                    vec![7, 5, 1],
                    vec![7, 4, 2],
                    vec![8, 3, 1],
                    vec![7, 6, 1],
                    vec![7, 5, 2],
                    vec![7, 4, 1],
                    vec![8, 4, 1],
                    vec![8, 3, 2],
                    vec![9, 2, 1],
                    vec![7, 6, 2],
                    vec![8, 5, 1],
                ]
                .into_iter()
            )
        );
    }

    #[test]
    fn test_val_macro() {
        let twelve_edo_5_lim = val![12, 19, 28];
        let syntonic_comma = monzo![-4, 4, -1];
        assert_eq!(twelve_edo_5_lim.evaluate(syntonic_comma), 0);
    }

    #[test]
    fn test_gpval() {
        let val_12 = gpval(12.0);
        assert_eq!(val_12.0[0], 12);
        assert_eq!(val_12.0[1], 19);
        assert_eq!(val_12.0[2], 28);
        assert_eq!(val_12.0[3], 34);
        let val_311 = gpval(311.0);
        assert_eq!(val_311.0[0], 311);
        assert_eq!(val_311.0[1], 493);
        assert_eq!(val_311.0[2], 722);
        assert_eq!(val_311.0[3], 873);
        let val_sixteen_point_nine = gpval(16.9);
        assert_eq!(val_sixteen_point_nine.0[0], 17);
        assert_eq!(val_sixteen_point_nine.0[1], 27);
        assert_eq!(val_sixteen_point_nine.0[2], 39);
    }
}
