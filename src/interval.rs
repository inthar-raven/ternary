use std::cmp::Ordering;

/// Trait for any type representing concrete interval sizes.
/// The trait provides the abstraction of a ℤ-module (free module over integers) under interval stacking.
/// Operations treat intervals as elements that can be combined additively.
pub trait Dyad: Copy + Eq + PartialEq + std::fmt::Debug + Send + Sync {
    /// The result of stacking two intervals.
    fn stack(self, rhs: Self) -> Self
    where
        Self: Sized;
    /// The result of stacking `rhs` down from `self`.
    fn unstack(self, rhs: Self) -> Self
    where
        Self: Sized;
    /// The logarithmic inverse of an interval.
    fn log_inv(self) -> Self
    where
        Self: Sized;
    /// The size of an interval in nepers.
    fn ln(self) -> f64;
    /// The cents size of an interval.
    fn cents(self) -> f64;
    /// The unison.
    fn unison() -> Self
    where
        Self: Sized;
    /// Stacking n copies of an interval where n is an integer.
    fn pow(self, n: i32) -> Self;
    /// The logarithmic size of `self` modulo `modulo` (e.g., octave reduction).
    /// Panics if `modulo == unison()`.
    fn rd(self, modulo: Self) -> Self
    where
        Self: Sized,
    {
        if modulo == Self::unison() {
            panic!("division by zero (log division by unison)")
        } else {
            let modulo = if modulo.cmp_dyad(&Self::unison()) == Ordering::Less {
                modulo.log_inv()
            } else {
                modulo
            };
            let mut ret = self; // Cloning is not necessary since we made `RawJiRatio` `Copy`.
            match ret.cmp_dyad(&Self::unison()) {
                Ordering::Less => {
                    while ret.cmp_dyad(&Self::unison()) == Ordering::Less {
                        ret = ret.stack(modulo);
                    }
                    ret
                }
                _ => {
                    while ret.cmp_dyad(&modulo) != Ordering::Less {
                        ret = ret.unstack(modulo);
                    }
                    ret
                }
            }
        }
    }
    /// Comparison for dyad sizes (based on logarithmic magnitude in cents).
    fn cmp_dyad(self, other: &Self) -> Ordering {
        self.cents().total_cmp(&other.cents())
    }
}

/// Trait for types representing JI ratios (just intonation ratios).
pub trait JiRatio: Dyad + Copy + Sync + Send {
    /// The numerator of a JI ratio
    fn numer(&self) -> u64;
    /// The denominator of a JI ratio
    fn denom(&self) -> u64;
}
