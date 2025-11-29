use std::cmp::Ordering;
use std::fmt;
use std::ops::{Div, DivAssign, Mul, MulAssign};

use num_traits::{CheckedDiv, CheckedMul};

use crate::helpers::gcd;
use crate::interval::Dyad;

// ERRORS

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
/// Error type for attempts to construct a RawJiRatio from a non-positive ratio.
pub struct IllegalJiRatio {
    // Show what the attempt was
    numer: u64,
    denom: u64,
}

impl fmt::Display for IllegalJiRatio {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "tried to create invalid JI ratio {}/{}",
            self.numer, self.denom
        )
    }
}

impl std::error::Error for IllegalJiRatio {}

/// Error type for invalid outputs of JI interval arithmetic.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum BadJiArith {
    /// logarithmic division by 0
    LogDivByUnison,
}

impl std::fmt::Display for BadJiArith {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            &Self::LogDivByUnison => {
                write!(f, "tried to log divide a JI ratio by the unison")
            }
        }
    }
}

impl std::error::Error for BadJiArith {}

// STRUCTS

/// A quasi-primitive wrapper class for a JI ratio. Implements `Copy`.
/// Will eventually be deprecated for most uses, but may be used for quick computations.
#[derive(Debug, Default, Clone, Copy)]
#[repr(C)]
pub struct RawJiRatio {
    numer: u64,
    denom: u64,
}

impl PartialEq for RawJiRatio {
    fn eq(&self, other: &Self) -> bool {
        self.numer * other.denom == other.numer * self.denom
    }
}

impl Eq for RawJiRatio {}

impl PartialOrd for RawJiRatio {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl std::iter::Product for RawJiRatio {
    fn product<I: Iterator<Item = RawJiRatio>>(iter: I) -> Self {
        iter.fold(RawJiRatio::UNISON, |x, y| x * y)
    }
}

impl Ord for RawJiRatio {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.numer * other.denom).cmp(&(other.numer * self.denom))
    }
}

impl RawJiRatio {
    /// The Pythagorean major pentatonic scale.
    pub const PYTH_5: [RawJiRatio; 5] = [
        RawJiRatio { numer: 9, denom: 8 },
        RawJiRatio {
            numer: 81,
            denom: 64,
        },
        RawJiRatio::PYTH_5TH,
        RawJiRatio {
            numer: 27,
            denom: 16,
        },
        RawJiRatio::OCTAVE,
    ];

    /// The Pythagorean Ionian mode.
    pub const PYTH_7: [RawJiRatio; 7] = [
        RawJiRatio { numer: 9, denom: 8 },
        RawJiRatio {
            numer: 81,
            denom: 64,
        },
        RawJiRatio::PYTH_4TH,
        RawJiRatio::PYTH_5TH,
        RawJiRatio {
            numer: 27,
            denom: 16,
        },
        RawJiRatio {
            numer: 243,
            denom: 128,
        },
        RawJiRatio::OCTAVE,
    ];

    /// The JI scale 12:14:16:18:21:24.
    pub const TAS_5: [RawJiRatio; 5] = [
        RawJiRatio { numer: 7, denom: 6 },
        RawJiRatio::PYTH_4TH,
        RawJiRatio::PYTH_5TH,
        RawJiRatio { numer: 7, denom: 4 },
        RawJiRatio::OCTAVE,
    ];

    /// The septal diasem scale.
    pub const TAS_9: [RawJiRatio; 9] = [
        RawJiRatio { numer: 9, denom: 8 },
        RawJiRatio { numer: 7, denom: 6 },
        RawJiRatio {
            numer: 21,
            denom: 16,
        },
        RawJiRatio::PYTH_4TH,
        RawJiRatio::PYTH_5TH,
        RawJiRatio {
            numer: 14,
            denom: 9,
        },
        RawJiRatio { numer: 7, denom: 4 },
        RawJiRatio {
            numer: 16,
            denom: 9,
        },
        RawJiRatio::OCTAVE,
    ];

    /// The pental blackdye scale.
    pub const BLACKDYE: [RawJiRatio; 10] = [
        RawJiRatio {
            numer: 10,
            denom: 9,
        },
        RawJiRatio { numer: 9, denom: 8 },
        RawJiRatio { numer: 6, denom: 5 },
        RawJiRatio::PYTH_4TH,
        RawJiRatio {
            numer: 27,
            denom: 20,
        },
        RawJiRatio::PYTH_5TH,
        RawJiRatio { numer: 8, denom: 5 },
        RawJiRatio {
            numer: 16,
            denom: 9,
        },
        RawJiRatio { numer: 9, denom: 5 },
        RawJiRatio::OCTAVE,
    ];

    /// The pental Zarlino scale, aka Ptolemy's intense diatonic.
    pub const ZARLINO: [RawJiRatio; 7] = [
        RawJiRatio { numer: 9, denom: 8 },
        RawJiRatio { numer: 5, denom: 4 },
        RawJiRatio::PYTH_4TH,
        RawJiRatio::PYTH_5TH,
        RawJiRatio { numer: 5, denom: 3 },
        RawJiRatio {
            numer: 15,
            denom: 8,
        },
        RawJiRatio::OCTAVE,
    ];

    /// The septal version of Zarlino, modifying Pythagorean intervals by 64/63 rather than flattening them by 81/80.
    pub const ARCHYLINO: [RawJiRatio; 7] = [
        RawJiRatio { numer: 9, denom: 8 },
        RawJiRatio { numer: 9, denom: 7 },
        RawJiRatio::PYTH_4TH,
        RawJiRatio::PYTH_5TH,
        RawJiRatio {
            numer: 12,
            denom: 7,
        },
        RawJiRatio {
            numer: 27,
            denom: 14,
        },
        RawJiRatio::OCTAVE,
    ];

    /// Zhea Erose's Eurybia scale.
    pub const EURYBIA: [RawJiRatio; 12] = [
        RawJiRatio {
            numer: 23,
            denom: 22,
        },
        RawJiRatio {
            numer: 25,
            denom: 22,
        },
        RawJiRatio {
            numer: 13,
            denom: 11,
        },
        RawJiRatio {
            numer: 14,
            denom: 11,
        },
        RawJiRatio {
            numer: 15,
            denom: 11,
        },
        RawJiRatio {
            numer: 31,
            denom: 22,
        },
        RawJiRatio::PYTH_5TH,
        RawJiRatio {
            numer: 35,
            denom: 22,
        },
        RawJiRatio {
            numer: 37,
            denom: 22,
        },
        RawJiRatio {
            numer: 39,
            denom: 22,
        },
        RawJiRatio {
            numer: 21,
            denom: 11,
        },
        RawJiRatio::OCTAVE,
    ];

    pub fn checked_product<I: Iterator<Item = RawJiRatio>>(iter: &mut I) -> Option<Self> {
        iter.try_fold(RawJiRatio::UNISON, |x, y| x.checked_mul(&y))
    }

    pub fn checked_pow(&self, n: i32) -> Option<Self> {
        if n >= 0 {
            let mut result = Some(Self::UNISON);
            for _ in 0..n {
                result = result?.checked_mul(self);
                result?;
            }
            result
        } else {
            let mut result = Some(Self::UNISON);
            for _ in 0..-n {
                result = result?.checked_div(self);
                result?;
            }
            result
        }
    }

    /// Check the input for validity before creating a new `RawJiRatio`.
    #[inline(always)]
    pub fn try_new(numer: u64, denom: u64) -> Result<RawJiRatio, IllegalJiRatio> {
        if (denom == 0) || (numer == 0) {
            Err(IllegalJiRatio { numer, denom })
        } else {
            let d = gcd(numer, denom);
            Ok(RawJiRatio {
                numer: numer / d,
                denom: denom / d,
            })
        }
    }
    /// The numerator of a `RawJiRatio`.
    #[inline(always)]
    pub fn numer(&self) -> u64 {
        self.numer
    }
    /// The denominator of a `RawJiRatio`.
    #[inline(always)]
    pub fn denom(&self) -> u64 {
        self.denom
    }
    /// The reciprocal of a `RawJiRatio`.
    #[inline(always)]
    pub fn reciprocal(&self) -> RawJiRatio {
        RawJiRatio {
            numer: self.denom,
            denom: self.numer,
        }
    }
    /// 1/1, the unison.
    pub const UNISON: RawJiRatio = RawJiRatio { numer: 1, denom: 1 };
    /// 2/1, the octave.
    pub const OCTAVE: RawJiRatio = RawJiRatio { numer: 2, denom: 1 };
    /// 3/2, the Pythagorean perfect fifth.
    pub const PYTH_5TH: RawJiRatio = RawJiRatio { numer: 3, denom: 2 };
    /// 4/3, the Pythagorean perfect fourth.
    pub const PYTH_4TH: RawJiRatio = RawJiRatio { numer: 4, denom: 3 };
    /// 3/2, the tritave.
    pub const TRITAVE: RawJiRatio = RawJiRatio { numer: 3, denom: 1 };
    /// 5/4, the pental major third.
    pub const PENTAL_MAJ3: RawJiRatio = RawJiRatio { numer: 5, denom: 4 };
    /// 7/4, the harmonic seventh.
    pub const SEPTIMAL_MIN7: RawJiRatio = RawJiRatio { numer: 7, denom: 4 };
    /// 11/8, the harmonic half-sharp fourth.
    pub const SEMIAUGMENTED_4TH: RawJiRatio = RawJiRatio {
        numer: 11,
        denom: 8,
    };

    /// Get the nth harmonic.
    pub fn harm(n: u64) -> Result<Self, IllegalJiRatio> {
        if n == 0 {
            Err(IllegalJiRatio { numer: 0, denom: 1 })
        } else {
            // SAFETY: gcd(n, 1) == 1
            Ok(RawJiRatio { numer: n, denom: 1 })
        }
    }

    /// Reduce the ratio r modulo the logarithmic absolute value of m.
    pub fn checked_rd(self, equave: Self) -> Result<Self, BadJiArith> {
        if equave.numer() == equave.denom() {
            Err(BadJiArith::LogDivByUnison)
        } else {
            let equave = if equave.numer() < equave.denom() {
                equave.reciprocal()
            } else {
                equave
            };
            let mut ret = self; // Cloning is not necessary since we made `RawJiRatio` `Copy`.
            if ret >= RawJiRatio::UNISON {
                while ret >= equave {
                    ret /= equave;
                }
                Ok(ret)
            } else {
                while ret < RawJiRatio::UNISON {
                    ret *= equave;
                }
                Ok(ret)
            }
        }
    }
    /// Logarithmic absolute value of a JI ratio (if m > n, m/n; otherwise n/m).
    pub fn magnitude(self) -> Self {
        if self < Self::UNISON {
            self.reciprocal()
        } else {
            self
        }
    }
}

impl Dyad for RawJiRatio {
    fn stack(self, rhs: Self) -> Self {
        self * rhs
    }
    fn unstack(self, rhs: Self) -> Self {
        self / rhs
    }
    fn log_inv(self) -> Self {
        self.reciprocal()
    }
    fn cents(self) -> f64 {
        ((self.numer as f64) / (self.denom as f64)).log2() * 1200.0
    }
    fn ln(self) -> f64 {
        ((self.numer as f64) / (self.denom as f64)).ln()
    }
    fn unison() -> Self {
        Self::UNISON
    }
    fn pow(self, n: i32) -> Self {
        (0..n).fold(Self::UNISON, |acc, _| self * acc)
    }

    fn rd(self, modulo: Self) -> Self {
        if modulo == RawJiRatio::UNISON {
            panic!("division by zero (log division by unison)")
        } else {
            let modulo = if modulo.numer() < modulo.denom() {
                modulo.reciprocal()
            } else {
                modulo
            };
            let mut ret = self;
            if ret >= Self::UNISON {
                while ret >= modulo {
                    ret /= modulo;
                }
                ret
            } else {
                while ret < Self::UNISON {
                    ret *= modulo;
                }
                ret
            }
        }
    }
}

impl fmt::Display for RawJiRatio {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}", self.numer, self.denom)
    }
}

impl Mul for RawJiRatio {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let d = gcd(self.numer * other.numer, self.denom * other.denom);
        RawJiRatio {
            numer: self.numer * other.numer / d,
            denom: self.denom * other.denom / d,
        }
    }
}
impl MulAssign for RawJiRatio {
    fn mul_assign(&mut self, other: Self) {
        let d = gcd(self.numer * other.numer, self.denom * other.denom);
        self.numer *= other.numer;
        self.numer /= d;
        self.denom *= other.denom;
        self.denom /= d;
    }
}

impl CheckedMul for RawJiRatio {
    fn checked_mul(&self, other: &Self) -> Option<Self> {
        let d = gcd(self.numer * other.numer, self.denom * other.denom);
        Some(RawJiRatio {
            numer: self.numer.checked_mul(other.numer)? / d,
            denom: self.denom.checked_mul(other.denom)? / d,
        })
    }
}

impl Div for RawJiRatio {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let d = gcd(self.numer * other.denom, self.denom * other.numer);
        RawJiRatio {
            numer: self.numer * other.denom / d,
            denom: self.denom * other.numer / d,
        }
    }
}

impl DivAssign for RawJiRatio {
    fn div_assign(&mut self, other: Self) {
        let d = gcd(self.numer * other.denom, self.denom * other.numer);
        self.numer *= other.denom;
        self.numer /= d;
        self.denom *= other.numer;
        self.denom /= d;
    }
}

impl CheckedDiv for RawJiRatio {
    fn checked_div(&self, other: &Self) -> Option<Self> {
        let d = gcd(self.numer * other.denom, self.denom * other.numer);
        Some(RawJiRatio {
            numer: self.numer.checked_mul(other.denom)? / d,
            denom: self.denom.checked_mul(other.numer)? / d,
        })
    }
}
