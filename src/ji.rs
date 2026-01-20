//! Just Intonation scale operations and tuning solvers.
//!
//! This module provides tools for working with JI (Just Intonation) scales:
//! analyzing their structure, finding tunings for step signatures, and
//! constructing scales from generators or harmonic series segments.
//!
//! # Key Concepts
//!
//! - **Odd limit**: The set of JI intervals with odd numerator and denominator
//!   up to a given limit.
//! - **Cumulative form**: A scale represented as intervals from the tonic
//!   (e.g., `[9/8, 5/4, 4/3, 3/2, 5/3, 15/8, 2/1]`).
//! - **Step form**: A scale represented as consecutive intervals
//!   (e.g., `[9/8, 10/9, 16/15, 9/8, 10/9, 9/8, 16/15]`).
//! - **Constant structure (CS)**: A scale where each interval class is unique — no
//!   interval appears as both an n-step and an m-step (n ≠ m).
//! - **Interleaved scale**: A scale formed by duplicating a "strand" and
//!   offsetting copies by notes of an "offset chord", where for every pair of
//!   strands one note of one strand always occurs strictly between two notes of
//!   the other.
//!
//! # Examples
//!
//! ```
//! use ternary::ji::{mode, step_form, cumulative_form, is_cs_ji_scale};
//! use ternary::ji_ratio::RawJiRatio;
//!
//! // Convert between cumulative and step forms
//! let zarlino = vec![
//!     RawJiRatio::try_new(9, 8).unwrap(),   // 9/8 from tonic
//!     RawJiRatio::try_new(5, 4).unwrap(),   // 5/4 from tonic
//!     RawJiRatio::try_new(4, 3).unwrap(),   // etc.
//!     RawJiRatio::try_new(3, 2).unwrap(),
//!     RawJiRatio::try_new(5, 3).unwrap(),
//!     RawJiRatio::try_new(15, 8).unwrap(),
//!     RawJiRatio::OCTAVE,
//! ];
//!
//! let steps = step_form(&zarlino);
//! assert_eq!(steps.len(), 7);  // 7 steps in the scale
//!
//! // Get a different mode (rotation)
//! let dorian = mode(&zarlino, 1);
//! assert_eq!(dorian.len(), 7);
//!
//! // Check if the scale is a constant structure
//! assert!(is_cs_ji_scale(&zarlino));
//! ```

use itertools::Itertools;
use num_integer::{gcd, lcm};
use std::collections::BTreeSet;
use std::ops::ControlFlow;

use crate::helpers::{ScaleError, is_sorted_strictly_desc, pairs};
use crate::interval::{Dyad, JiRatio};
use crate::ji_ratio::{BadJiArith, RawJiRatio};
use crate::monzo::Monzo;
use crate::simple_steps::SIMPLE_STEPS;
use crate::words::{CountVector, rotate};

/// Given a list of odd numbers, return the octave-reduced intervals in the corresponding odd-limit,
/// not including the unison.
pub fn specified_odd_limit(odds: &[u32]) -> Vec<RawJiRatio> {
    let odds: Vec<u32> = odds
        .iter()
        .filter(|x| **x > 0 && **x % 2 == 1)
        .copied()
        .collect(); // filter out invalid input
    pairs(&odds, &odds)
        .into_iter()
        .map(|(m, n)| {
            RawJiRatio::try_new(m, n)
                .expect("should have removed all non-positive ints")
                .rd(RawJiRatio::OCTAVE)
        })
        .filter(|ratio| *ratio > RawJiRatio::UNISON)
        .collect()
}

/// Returns the octave-reduced intervals of a specified odd limit, not including the unison.
/// Generates all ratios with odd numerator and denominator up to the limit.
///
/// # Examples
///
/// ```
/// use ternary::ji::odd_limit;
/// use ternary::ji_ratio::RawJiRatio;
/// use itertools::Itertools;
///
/// // The 5-odd-limit contains ratios like 3/2, 5/4, 5/3, etc.
/// let five_limit: Vec<_> = odd_limit(5).into_iter().sorted().collect();
///
/// assert_eq!(five_limit, vec![
///     RawJiRatio::try_new(6, 5).unwrap(),
///     RawJiRatio::try_new(5, 4).unwrap(),
///     RawJiRatio::try_new(4, 3).unwrap(),
///     RawJiRatio::try_new(3, 2).unwrap(),
///     RawJiRatio::try_new(8, 5).unwrap(),
///     RawJiRatio::try_new(5, 3).unwrap(),
/// ]);
///
/// // Higher odd-limits contain more intervals
/// let nine_limit = odd_limit(9);
/// assert!(nine_limit.len() > five_limit.len());
/// ```
pub fn odd_limit(limit: u32) -> Vec<RawJiRatio> {
    let odds = (0..=(limit - 1) / 2).map(|i| 2 * i + 1).collect::<Vec<_>>();
    pairs(&odds, &odds)
        .into_iter()
        .map(|(m, n)| {
            RawJiRatio::try_new(m, n)
                .expect("should have removed all non-positive ints")
                .rd(RawJiRatio::OCTAVE)
        })
        .sorted_unstable()
        .dedup()
        .filter(|ratio| *ratio > RawJiRatio::UNISON)
        .collect()
}

/// Solutions to a step signature (with decreasing step sizes).
/// Steps are required to be between `cents_lower_bound` and `cents_upper_bound`.
/// All but the smallest step are required to be in SIMPLE_STEPS.
/// If `allow_neg_aber` is `true`, then the smallest step size is allowed to be negative.
pub fn solve_step_sig_simple_steps(
    step_sig: &[usize],
    equave: Monzo,
    cents_lower_bound: f64,
    cents_upper_bound: f64,
    allow_neg_aber: bool,
) -> Vec<Vec<Monzo>> {
    let small_steps: Vec<_> = SIMPLE_STEPS
        .into_iter()
        .filter(|monzo| monzo.cents() > cents_lower_bound && monzo.cents() < cents_upper_bound)
        .collect();
    let prod = (0..step_sig.len() - 1)
        .map(|_| small_steps.to_vec())
        .multi_cartesian_product()
        .collect::<Vec<_>>();
    let mut result = vec![];
    for steps in prod {
        if is_sorted_strictly_desc(&steps) {
            let multiplied_steps = steps
                .iter()
                .copied()
                .enumerate()
                .map(|(i, v)| v * (step_sig[i] as i32));
            let sum = multiplied_steps.into_iter().sum();
            let residue = equave - sum;
            if residue.is_divisible_by(step_sig[step_sig.len() - 1] as i32)
                && (allow_neg_aber || residue.is_positive())
            {
                let smallest_step = residue / (step_sig[step_sig.len() - 1] as i32);
                let mut soln = steps;

                if smallest_step < soln[soln.len() - 1] {
                    // Check if the last step is actually the smallest to validate the solution.
                    soln.push(smallest_step);
                    result.push(soln);
                }
            }
        }
    }
    result
}

/// Multiset of `subword_length`-step intervals in a JI scale.
/// Counts the frequency of each interval subtype.
///
/// # Examples
///
/// ```
/// use ternary::ji::spectrum;
/// use ternary::ji_ratio::RawJiRatio;
///
/// // Pythagorean pentatonic: 1/1 9/8 81/64 3/2 27/16 2/1
/// let pentatonic = vec![
///     RawJiRatio::try_new(9, 8).unwrap(),
///     RawJiRatio::try_new(81, 64).unwrap(),
///     RawJiRatio::try_new(3, 2).unwrap(),
///     RawJiRatio::try_new(27, 16).unwrap(),
///     RawJiRatio::OCTAVE,
/// ];
///
/// // Get the 1-step spectrum (seconds)
/// let seconds = spectrum(&pentatonic, 1);
/// // Pentatonic has exactly two step sizes: 9/8 (whole tone) and 32/27 (minor third)
/// assert_eq!(seconds.keys_count(), 2);
/// ```
pub fn spectrum(scale: &[RawJiRatio], subword_length: usize) -> CountVector<RawJiRatio> {
    let mut result = std::collections::BTreeMap::new();
    let intervals = (0..scale.len()).map(|degree| {
        (scale[(degree + subword_length) % scale.len()] / scale[degree]).rd(RawJiRatio::OCTAVE)
    });

    for interval in intervals {
        if let Some(update_this) = result.get_mut(&interval) {
            *update_this += 1;
        } else {
            result.insert(interval, 1);
        }
    }
    CountVector::from_btree_map(result)
}

/// Display a JI scale as a list of pitches from the tonic.
pub fn disp_ji_scale(scale: &[RawJiRatio]) -> String {
    let mut ret: String = String::from("");
    for item in scale {
        ret.push_str(&format!("{item}"));
        ret.push(' ');
    }
    ret
}

/// Display a JI scale as a JI chord (ratio of two or more integers).
/// Computes the LCM of denominators to express as a single chord.
pub fn disp_ji_scale_as_enumerated_chord(scale: &[RawJiRatio]) -> String {
    let mut ell_cee_emm: u32 = 1;
    for item in scale {
        ell_cee_emm = lcm(ell_cee_emm, item.denom());
    }
    let mut ret: String = ell_cee_emm.to_string();
    ret.push(':');
    for i in 0..scale.len() {
        assert_eq!(ell_cee_emm % (scale[i].denom()), 0);
        let multiply_by = ell_cee_emm / (scale[i].denom());
        ret.push_str(&(scale[i].numer() * multiply_by).to_string());
        if i < scale.len() - 1 {
            ret.push(':');
        }
    }
    ret
}

/// Get a specific mode of a JI scale in cumulative form.
/// Rotates by the given degree and returns the intervals from the new root.
pub fn mode(scale: &[RawJiRatio], degree: usize) -> Vec<RawJiRatio> {
    let steps: Vec<_> = step_form(scale);
    // rotate() already does degree % scale.len()
    let steps_rotated = rotate(&steps, degree);
    cumulative_form(&steps_rotated)
}

/// Convert a cumulative form into a step form.
/// Each step is the interval from one note to the next in the scale.
pub fn step_form(cumul_form: &[RawJiRatio]) -> Vec<RawJiRatio> {
    [
        &[cumul_form[0]],
        (0..cumul_form.len() - 1)
            .map(|i| cumul_form[i + 1] / cumul_form[i])
            .collect::<Vec<_>>()
            .as_slice(),
    ]
    .concat()
}

/// Convert a step form into a cumulative form.
/// Each note is the product of all steps up to that point.
pub fn cumulative_form(step_form: &[RawJiRatio]) -> Vec<RawJiRatio> {
    step_form
        .iter()
        .scan(RawJiRatio::UNISON, |acc, &step| {
            *acc *= step;
            Some(*acc)
        })
        .collect()
}

/// All modes of a JI scale written in cumulative form.
pub fn ji_scale_modes(scale: &[RawJiRatio]) -> Vec<Vec<RawJiRatio>> {
    (0..scale.len()).map(|degree| mode(scale, degree)).collect()
}

/// Returns the harmonic series mode `mode_num`: mode_num:...:(2*mode_num).
/// Includes the octave duplication.
///
/// # Examples
///
/// ```
/// use ternary::ji::harmonic_mode;
/// use ternary::ji_ratio::RawJiRatio;
///
/// // Harmonic mode 4 is the 4:5:6:7:8 chord
/// let mode_4 = harmonic_mode(4).unwrap();
/// assert_eq!(mode_4.len(), 4);  // 5/4, 6/4, 7/4, 8/4
/// assert_eq!(mode_4[0], RawJiRatio::try_new(5, 4).unwrap());
/// assert_eq!(mode_4[3], RawJiRatio::OCTAVE);
///
/// // Harmonic mode 8 is the 8:9:10:11:12:13:14:15:16 chord
/// let mode_8 = harmonic_mode(8).unwrap();
/// assert_eq!(mode_8.len(), 8);
/// ```
pub fn harmonic_mode(mode_num: u32) -> Result<Vec<RawJiRatio>, ScaleError> {
    if mode_num < 1 {
        Err(ScaleError::CannotMakeScale)
    } else {
        Ok((mode_num + 1..=(2 * mode_num))
            .map(|x| RawJiRatio::try_new(x, mode_num).expect("`numer` should be > `denom` here"))
            .collect())
    }
}

/// Returns the harmonic series mode `mode_num` without octave: mode_num:...:(2*mode_num - 1).
pub fn harmonic_mode_no_oct(mode_num: u32) -> Result<Vec<RawJiRatio>, ScaleError> {
    if mode_num < 1 {
        Err(ScaleError::CannotMakeScale)
    } else {
        Ok((mode_num + 1..=(2 * mode_num - 1))
            .map(|x| RawJiRatio::try_new(x, mode_num).expect("`numer` should be > `denom`"))
            .collect())
    }
}

/// Compute if `offset` is valid for an interleaved scale with strand `strand`.
/// Validates that all strand rotations differ from the offset by intervals in the strand.
pub fn is_valid_offset(strand: &[RawJiRatio], offset: RawJiRatio) -> bool {
    let _ = (1..=(strand.len() - 1)).try_for_each(|i| {
        let i_steps = spectrum(strand, i);
        let min_i_step = *i_steps
            .into_inner()
            .first_key_value()
            .expect("`i_steps` should not be empty")
            .0;
        let max_i_step = *i_steps
            .into_inner()
            .last_key_value()
            .expect("`i_steps` should not be empty")
            .0;
        if offset >= min_i_step && offset <= max_i_step {
            // If `offset` falls in [`min_i_step`, `max_i_step`], it's not valid for an interleaved scale with `strand`.
            ControlFlow::Break(false)
        } else {
            ControlFlow::Continue(())
        }
    });
    true
}

/// Compute if `offset_chord` is valid for an interleaved scale with strand `strand`.
/// The offset chord is given as a Vec of offsets from the unison.
/// Validates that all pairwise differences of offset notes are intervals in the strand.
pub fn is_valid_offset_chord(strand: &[RawJiRatio], offset_chord: &[RawJiRatio]) -> bool {
    let mut offsets_with_unison: Vec<RawJiRatio> = offset_chord.to_vec();
    offsets_with_unison.push(RawJiRatio::UNISON);
    // Check validity of every interval in the offset_chord.
    offsets_with_unison.sort_by(|a, b| (*a).cmp(b));
    for (i, ratio1) in offsets_with_unison.iter().copied().enumerate() {
        for (j, ratio2) in offsets_with_unison.iter().copied().enumerate() {
            if i != j && !is_valid_offset(strand, ratio2 / ratio1) {
                return false;
            }
        }
    }
    true
}

/// Compute the set of valid offset chords for an interleaved scale on harmonic series mode `m`.
pub fn valid_offset_chords(
    strand: &[RawJiRatio],
    m: u32,
) -> Result<Vec<Vec<RawJiRatio>>, ScaleError> {
    let mode = harmonic_mode_no_oct(m)?; // harmonic mode m
    Ok(crate::helpers::powerset(&mode)
        .into_iter()
        .filter(|x| is_valid_offset_chord(strand, x))
        .collect())
}

/// Compute the set of *maximal* valid offset chords for an interleaved scale on harmonic series mode `m`.
pub fn maximal_valid_offset_chords(
    strand: &[RawJiRatio],
    m: u32,
) -> Result<Vec<Vec<RawJiRatio>>, ScaleError> {
    let mode = harmonic_mode_no_oct(m)?; // harmonic mode m
    let power_set = crate::helpers::powerset(&mode);
    let valid_offset_chords: Vec<_> = power_set
        .into_iter()
        .filter(|x| is_valid_offset_chord(strand, x))
        .collect();
    Ok(valid_offset_chords
        .iter()
        .filter(|x| crate::helpers::is_maximal_in(x, &valid_offset_chords))
        .cloned()
        .collect())
}

/// Compute an interleaved scale with the given strand and offset chord, if it is valid.
/// A *strand* is a smaller scale that is duplicated and copies of it offset by the notes of the offset chord.
/// The resulting scale made of interleaved strands is *interleaved* if all the strands are actually interleaved.
/// Returns error if the combination is not valid.
pub fn interleaved_scale_ji(
    strand: &[RawJiRatio],
    offset_chord: &[RawJiRatio],
) -> Result<Vec<RawJiRatio>, ScaleError> {
    if is_valid_offset_chord(strand, offset_chord) {
        let mut vec: Vec<RawJiRatio> = vec![];
        for offset in offset_chord {
            for note in strand {
                vec.push((*offset * *note).rd(RawJiRatio::OCTAVE));
            }
            vec.push((*offset).rd(RawJiRatio::OCTAVE));
        }
        for note in strand {
            vec.push((*note).rd(RawJiRatio::OCTAVE));
        }
        vec.sort();
        vec.dedup();
        vec.push(RawJiRatio::OCTAVE);
        Ok(vec.drain(1..).collect()) // Remove 1/1
    } else {
        Err(ScaleError::NotInterleavable)
    }
}

/// Given `arr` a periodic JI scale given in JI ratios from the tonic, is `arr` a CS (constant structure)?
/// Assumes `arr[0]` = the 1-step from the tonic, ..., `arr[arr.len() - 1]` = the equave;
/// arr.len() == the scale size.
///
/// A constant structure (CS) is a scale where no interval appears in two different
/// interval classes. For example, 3/2 cannot be both a 4-step and a 5-step.
///
/// # Examples
///
/// ```
/// use ternary::ji::is_cs_ji_scale;
/// use ternary::ji_ratio::RawJiRatio;
///
/// // Pythagorean major scale is CS
/// let pyth_major = vec![
///     RawJiRatio::try_new(9, 8).unwrap(),
///     RawJiRatio::try_new(81, 64).unwrap(),
///     RawJiRatio::try_new(4, 3).unwrap(),
///     RawJiRatio::try_new(3, 2).unwrap(),
///     RawJiRatio::try_new(27, 16).unwrap(),
///     RawJiRatio::try_new(243, 128).unwrap(),
///     RawJiRatio::OCTAVE,
/// ];
/// assert!(is_cs_ji_scale(&pyth_major));
/// ```
pub fn is_cs_ji_scale(arr: &[RawJiRatio]) -> bool {
    let n = arr.len();
    let mut interval_classes = vec![vec![RawJiRatio::UNISON; n]; n - 1];
    // interval_classes[i] is the set of (i+1)-steps in the scale. Get 1-steps, ..., (n-1)-steps.
    for i in 1..=(n - 1) {
        // i is the increment.
        for j in 0..n {
            // j is the 0-indexed degree.
            let unreduced_interval = if i + j >= n {
                let equave = arr[n - 1];
                arr[(i + j) % n] * equave
            } else {
                arr[i + j]
            };
            // Unstack by arr[j] so we have the interval on the j-degree.
            interval_classes[i - 1][j] = unreduced_interval / arr[j];
        }
    }
    // Check for pairwise intersections between step classes.
    // Range the first class over 1-steps, ..., n/2-steps and range the second over all classes with more steps than the first class.
    // (Watch out for off-by-1 errors!)
    for i in 0..(n / 2) {
        // This loop makes at most (n-1)(n-2)/2 comparisons between sets.
        let unique_i_plus_1_steps: BTreeSet<RawJiRatio> =
            interval_classes[i].iter().cloned().collect();
        for class in interval_classes.iter().take(n - 1).skip(i + 1) {
            let unique_j_plus_1_steps: BTreeSet<RawJiRatio> = class.iter().cloned().collect();
            if !unique_i_plus_1_steps.is_disjoint(&unique_j_plus_1_steps) {
                // If two different classes have a non-empty intersection, return false.
                return false;
            }
        }
    }
    true
}

/// Given a generator sequence `gs`,
/// return an `n`-note generator sequence scale with equave `equave`
/// formed by stacking and reducing `n - 1` intervals of `gs` in turn.
///
/// # Examples
///
/// ```
/// use ternary::ji::gs_scale;
/// use ternary::ji_ratio::RawJiRatio;
///
/// // Build a pentatonic scale by stacking 3/2 fifths
/// let generators = [RawJiRatio::try_new(3, 2).unwrap()];
/// let pentatonic = gs_scale(&generators, 5, RawJiRatio::OCTAVE).unwrap();
/// assert_eq!(pentatonic.len(), 5);
///
/// // Alternating generators create more complex scales
/// let gens = [
///     RawJiRatio::try_new(7, 6).unwrap(),
///     RawJiRatio::try_new(8, 7).unwrap(),
/// ];
/// let scale = gs_scale(&gens, 5, RawJiRatio::OCTAVE).unwrap();
/// assert_eq!(scale.len(), 5);
/// ```
pub fn gs_scale(
    gs: &[RawJiRatio],
    n: usize,
    equave: RawJiRatio,
) -> Result<Vec<RawJiRatio>, Box<dyn std::error::Error>> {
    if gs.is_empty() || n == 0 || equave == RawJiRatio::UNISON {
        Err(Box::new(ScaleError::CannotMakeScale))
    } else if equave == RawJiRatio::UNISON {
        Err(Box::new(BadJiArith::LogDivByUnison))
    } else {
        let equave = equave.magnitude(); // Take the equave's magnitude
        let mut result = vec![];
        let mut gs_cycle = gs.iter().cycle();
        let mut accumulator = RawJiRatio::UNISON;
        for _ in 0..n - 1 {
            accumulator = (accumulator * *gs_cycle.next().expect("`gs.len() > 0` in this branch, thus `gs.into_iter().cycle()` is infinite and can never run out") ).rd(equave);
            result.push(accumulator);
        }
        result.push(equave);
        result.sort();
        result.dedup(); // The stacking may have resulted in duplicate notes; TODO: notify the user of duplicates.
        Ok(result)
    }
}

/// Return a well-formed necklace of `gener_class`-steps in a given JI scale.
/// A necklace is a rotationally-invariant way to arrange generators.
pub fn well_formed_necklace_in_ji_scale(
    scale: &[RawJiRatio],
    gener_class: usize,
) -> Result<Vec<RawJiRatio>, ScaleError> {
    if gcd(scale.len() as i32, gener_class as i32) == 1 {
        Ok((0..(scale.len()))
            .map(|k| {
                (scale[(gener_class * (k + 1)) % scale.len()]
                    / scale[(gener_class * k) % scale.len()])
                .rd(RawJiRatio::OCTAVE)
            })
            .collect()) // For every `k` (including the last one), get the `k`th stacked `gener_class`-step.
    } else {
        Err(ScaleError::NonCoprimeGenError)
    }
}

#[cfg(test)]
mod tests {
    #[allow(unused)]
    use super::*;
    use crate::ji_ratio::RawJiRatio;

    use crate::monzo;

    #[test]
    fn test_ji_spectrum() {
        use crate::ji::spectrum;
        use crate::ji_ratio::RawJiRatio;
        // Pythagorean pentatonic: 1/1 9/8 81/64 3/2 27/16 2/1
        let pentatonic = vec![
            RawJiRatio::try_new(9, 8).unwrap(),
            RawJiRatio::try_new(81, 64).unwrap(),
            RawJiRatio::try_new(3, 2).unwrap(),
            RawJiRatio::try_new(27, 16).unwrap(),
            RawJiRatio::OCTAVE,
        ];
        // Get the 1-step spectrum (seconds)
        let seconds = spectrum(&pentatonic, 1);
        // Pentatonic has exactly two step sizes: 9/8 (whole tone) and 32/27 (minor third)
        assert_eq!(seconds.keys_count(), 2);
    }

    #[test]
    fn test_ji_scale_modes() {
        let diasem_modes = ji_scale_modes(&RawJiRatio::TAS_9);
        let correct_diasem_modes = vec![
            vec![
                RawJiRatio::try_new(9, 8).unwrap(),
                RawJiRatio::try_new(7, 6).unwrap(),
                RawJiRatio::try_new(21, 16).unwrap(),
                RawJiRatio::try_new(4, 3).unwrap(),
                RawJiRatio::try_new(3, 2).unwrap(),
                RawJiRatio::try_new(14, 9).unwrap(),
                RawJiRatio::try_new(7, 4).unwrap(),
                RawJiRatio::try_new(16, 9).unwrap(),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::try_new(28, 27).unwrap(),
                RawJiRatio::try_new(7, 6).unwrap(),
                RawJiRatio::try_new(32, 27).unwrap(),
                RawJiRatio::try_new(4, 3).unwrap(),
                RawJiRatio::try_new(112, 81).unwrap(),
                RawJiRatio::try_new(14, 9).unwrap(),
                RawJiRatio::try_new(128, 81).unwrap(),
                RawJiRatio::try_new(16, 9).unwrap(),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::try_new(9, 8).unwrap(),
                RawJiRatio::try_new(8, 7).unwrap(),
                RawJiRatio::try_new(9, 7).unwrap(),
                RawJiRatio::try_new(4, 3).unwrap(),
                RawJiRatio::try_new(3, 2).unwrap(),
                RawJiRatio::try_new(32, 21).unwrap(),
                RawJiRatio::try_new(12, 7).unwrap(),
                RawJiRatio::try_new(27, 14).unwrap(),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::try_new(64, 63).unwrap(),
                RawJiRatio::try_new(8, 7).unwrap(),
                RawJiRatio::try_new(32, 27).unwrap(),
                RawJiRatio::try_new(4, 3).unwrap(),
                RawJiRatio::try_new(256, 189).unwrap(),
                RawJiRatio::try_new(32, 21).unwrap(),
                RawJiRatio::try_new(12, 7).unwrap(),
                RawJiRatio::try_new(16, 9).unwrap(),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::try_new(9, 8).unwrap(),
                RawJiRatio::try_new(7, 6).unwrap(),
                RawJiRatio::try_new(21, 16).unwrap(),
                RawJiRatio::try_new(4, 3).unwrap(),
                RawJiRatio::try_new(3, 2).unwrap(),
                RawJiRatio::try_new(27, 16).unwrap(),
                RawJiRatio::try_new(7, 4).unwrap(),
                RawJiRatio::try_new(63, 32).unwrap(),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::try_new(28, 27).unwrap(),
                RawJiRatio::try_new(7, 6).unwrap(),
                RawJiRatio::try_new(32, 27).unwrap(),
                RawJiRatio::try_new(4, 3).unwrap(),
                RawJiRatio::try_new(3, 2).unwrap(),
                RawJiRatio::try_new(14, 9).unwrap(),
                RawJiRatio::try_new(7, 4).unwrap(),
                RawJiRatio::try_new(16, 9).unwrap(),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::try_new(9, 8).unwrap(),
                RawJiRatio::try_new(8, 7).unwrap(),
                RawJiRatio::try_new(9, 7).unwrap(),
                RawJiRatio::try_new(81, 56).unwrap(),
                RawJiRatio::try_new(3, 2).unwrap(),
                RawJiRatio::try_new(27, 16).unwrap(),
                RawJiRatio::try_new(12, 7).unwrap(),
                RawJiRatio::try_new(27, 14).unwrap(),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::try_new(64, 63).unwrap(),
                RawJiRatio::try_new(8, 7).unwrap(),
                RawJiRatio::try_new(9, 7).unwrap(),
                RawJiRatio::try_new(4, 3).unwrap(),
                RawJiRatio::try_new(3, 2).unwrap(),
                RawJiRatio::try_new(32, 21).unwrap(),
                RawJiRatio::try_new(12, 7).unwrap(),
                RawJiRatio::try_new(16, 9).unwrap(),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::try_new(9, 8).unwrap(),
                RawJiRatio::try_new(81, 64).unwrap(),
                RawJiRatio::try_new(21, 16).unwrap(),
                RawJiRatio::try_new(189, 128).unwrap(),
                RawJiRatio::try_new(3, 2).unwrap(),
                RawJiRatio::try_new(27, 16).unwrap(),
                RawJiRatio::try_new(7, 4).unwrap(),
                RawJiRatio::try_new(63, 32).unwrap(),
                RawJiRatio::OCTAVE,
            ],
        ];
        assert_eq!(diasem_modes, correct_diasem_modes);
    }
    #[test]
    fn test_gs() {
        let tas_scales: Vec<_> = (1..=8)
            .map(|i| {
                gs_scale(
                    &[
                        RawJiRatio::try_new(7, 6).unwrap(),
                        RawJiRatio::try_new(8, 7).unwrap(),
                    ],
                    i,
                    RawJiRatio::OCTAVE,
                )
                .unwrap()
            })
            .collect();
        assert_eq!(
            tas_scales,
            vec![
                vec![RawJiRatio::OCTAVE],
                vec![RawJiRatio::try_new(7, 6).unwrap(), RawJiRatio::OCTAVE],
                vec![
                    RawJiRatio::try_new(7, 6).unwrap(),
                    RawJiRatio::try_new(4, 3).unwrap(),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::try_new(7, 6).unwrap(),
                    RawJiRatio::try_new(4, 3).unwrap(),
                    RawJiRatio::try_new(14, 9).unwrap(),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::try_new(7, 6).unwrap(),
                    RawJiRatio::try_new(4, 3).unwrap(),
                    RawJiRatio::try_new(14, 9).unwrap(),
                    RawJiRatio::try_new(16, 9).unwrap(),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::try_new(28, 27).unwrap(),
                    RawJiRatio::try_new(7, 6).unwrap(),
                    RawJiRatio::try_new(4, 3).unwrap(),
                    RawJiRatio::try_new(14, 9).unwrap(),
                    RawJiRatio::try_new(16, 9).unwrap(),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::try_new(28, 27).unwrap(),
                    RawJiRatio::try_new(7, 6).unwrap(),
                    RawJiRatio::try_new(32, 27).unwrap(),
                    RawJiRatio::try_new(4, 3).unwrap(),
                    RawJiRatio::try_new(14, 9).unwrap(),
                    RawJiRatio::try_new(16, 9).unwrap(),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::try_new(28, 27).unwrap(),
                    RawJiRatio::try_new(7, 6).unwrap(),
                    RawJiRatio::try_new(32, 27).unwrap(),
                    RawJiRatio::try_new(4, 3).unwrap(),
                    RawJiRatio::try_new(112, 81).unwrap(),
                    RawJiRatio::try_new(14, 9).unwrap(),
                    RawJiRatio::try_new(16, 9).unwrap(),
                    RawJiRatio::OCTAVE,
                ],
            ]
        );

        let zil_gs = [
            RawJiRatio::try_new(7, 4).unwrap(),
            RawJiRatio::try_new(12, 7).unwrap(),
            RawJiRatio::try_new(7, 4).unwrap(),
            RawJiRatio::try_new(12, 7).unwrap(),
            RawJiRatio::try_new(7, 4).unwrap(),
            RawJiRatio::try_new(12, 7).unwrap(),
            RawJiRatio::try_new(7, 4).unwrap(),
            RawJiRatio::try_new(320, 189).unwrap(),
            RawJiRatio::try_new(7, 4).unwrap(),
            RawJiRatio::try_new(12, 7).unwrap(),
        ];
        let zil_24 = gs_scale(&zil_gs, 24, RawJiRatio::OCTAVE).unwrap();
        let correct_zil_24 = vec![
            RawJiRatio::try_new(525, 512).unwrap(),
            RawJiRatio::try_new(135, 128).unwrap(),
            RawJiRatio::try_new(35, 32).unwrap(),
            RawJiRatio::try_new(9, 8).unwrap(),
            RawJiRatio::try_new(4725, 4096).unwrap(),
            RawJiRatio::try_new(75, 64).unwrap(),
            RawJiRatio::try_new(315, 256).unwrap(),
            RawJiRatio::try_new(5, 4).unwrap(),
            RawJiRatio::try_new(21, 16).unwrap(),
            RawJiRatio::try_new(675, 512).unwrap(),
            RawJiRatio::try_new(2835, 2048).unwrap(),
            RawJiRatio::try_new(45, 32).unwrap(),
            RawJiRatio::try_new(189, 128).unwrap(),
            RawJiRatio::try_new(3, 2).unwrap(),
            RawJiRatio::try_new(1575, 1024).unwrap(),
            RawJiRatio::try_new(405, 256).unwrap(),
            RawJiRatio::try_new(105, 64).unwrap(),
            RawJiRatio::try_new(27, 16).unwrap(),
            RawJiRatio::try_new(7, 4).unwrap(),
            RawJiRatio::try_new(225, 128).unwrap(),
            RawJiRatio::try_new(945, 512).unwrap(),
            RawJiRatio::try_new(15, 8).unwrap(),
            RawJiRatio::try_new(63, 32).unwrap(),
            RawJiRatio::OCTAVE,
        ];
        assert_eq!(zil_24, correct_zil_24);
    }

    #[test]
    fn test_cs() {
        assert!(is_cs_ji_scale(&RawJiRatio::PYTH_5));
        assert!(is_cs_ji_scale(&RawJiRatio::PYTH_7));
        assert!(is_cs_ji_scale(&RawJiRatio::ZARLINO));
        assert!(is_cs_ji_scale(&RawJiRatio::TAS_5));
        assert!(is_cs_ji_scale(&RawJiRatio::TAS_9));
        assert!(is_cs_ji_scale(&RawJiRatio::BLACKDYE));
    }

    #[test]
    fn test_solve_81_odd_lim() {
        let diatonic_solns: Vec<Vec<Monzo>> =
            solve_step_sig_simple_steps(&[5, 2], Monzo::OCTAVE, 20.0, 300.0, false);
        assert_eq!(diatonic_solns, vec![vec![monzo![-3, 2], monzo![8, -5]]]);
        let blackdye_solns: Vec<Vec<Monzo>> =
            solve_step_sig_simple_steps(&[5, 2, 3], Monzo::OCTAVE, 20.0, 300.0, false);
        assert!(blackdye_solns.contains(&vec![
            monzo![1, -2, 1],
            monzo![4, -1, -1],
            monzo![-4, 4, -1]
        ]));
    }
}
