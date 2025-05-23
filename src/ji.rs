use itertools::Itertools;
use num_integer::{gcd, lcm};
use std::collections::BTreeSet;
use std::ops::ControlFlow;

use crate::helpers::{is_sorted_strictly_desc, pairs, ScaleError};
use crate::interval::Dyad;
use crate::ji_ratio::{BadJiArith, RawJiRatio};
use crate::monzo::Monzo;
use crate::words::{rotate, CountVector};

/// Given a list of odd numbers, return the octave-reduced intervals in the corresponding odd-limit,
/// not including the unison.
pub fn specified_odd_limit(odds: &[u64]) -> Vec<RawJiRatio> {
    let odds: Vec<u64> = odds
        .iter()
        .filter(|x| **x > 0 && **x % 2 == 1)
        .copied()
        .collect(); // filter out invalid input
    pairs(&odds, &odds.clone())
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
pub fn odd_limit(limit: u64) -> Vec<RawJiRatio> {
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
/// Steps are required to be smaller than 300c,
/// and all but the smallest step are required to be in the `SMALL_PRIMES`-prime-limited 81-odd-limit.
/// If `allow_neg_aber` is `true`, then the smallest step size is allowed to be negative.
pub fn solve_step_sig_81_odd_limit(
    step_sig: &[usize],
    equave: Monzo,
    allow_neg_aber: bool,
) -> Vec<Vec<Monzo>> {
    let cents_bound: f64 = 300.0;
    let small_steps: Vec<_> = Monzo::EIGHTY_ONE_ODD_LIMIT
        .into_iter()
        .filter(|monzo| monzo.cents() < cents_bound)
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
pub fn spectrum(scale: &[RawJiRatio], subword_length: usize) -> CountVector<RawJiRatio> {
    let mut result = std::collections::BTreeMap::new();
    for interval in (0..scale.len()).map(|degree| {
        (scale[(degree + subword_length) % scale.len()] / scale[degree]).rd(RawJiRatio::OCTAVE)
    }) {
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
        ret.push_str(&format!("{}", item));
        ret.push(' ');
    }
    ret
}

/// Display a JI scale as a JI chord (ratio of two or more integers).
pub fn disp_ji_scale_as_enumerated_chord(scale: &[RawJiRatio]) -> String {
    let mut ell_cee_emm: u64 = 1;
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
pub fn mode(scale: &[RawJiRatio], degree: usize) -> Vec<RawJiRatio> {
    let steps: Vec<_> = step_form(scale);
    // rotate() already does degree % scale.len()
    let steps_rotated = rotate(&steps, degree);
    cumulative_form(&steps_rotated)
}

/// Convert a cumulative form into a step form.
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
pub fn cumulative_form(step_form: &[RawJiRatio]) -> Vec<RawJiRatio> {
    step_form
        .iter()
        .scan(RawJiRatio::UNISON, |acc, &step| {
            *acc *= step;
            Some(*acc)
        })
        .collect()
}

/// Modes of a JI scale written in cumulative form.
pub fn ji_scale_modes(scale: &[RawJiRatio]) -> Vec<Vec<RawJiRatio>> {
    (0..scale.len()).map(|degree| mode(scale, degree)).collect()
}

/// Returns the JI chord `mode_num`:...:`2*mode_num`.
pub fn harmonic_mode(mode_num: u64) -> Result<Vec<RawJiRatio>, ScaleError> {
    if mode_num < 1 {
        Err(ScaleError::CannotMakeScale)
    } else {
        Ok((mode_num + 1..=(2 * mode_num))
            .map(|x| RawJiRatio::try_new(x, mode_num).expect("`numer` should be > `denom` here"))
            .collect())
    }
}

/// Returns the JI chord `mode_num`:...:`2*mode_num - 1`.
pub fn harmonic_mode_no_oct(mode_num: u64) -> Result<Vec<RawJiRatio>, ScaleError> {
    if mode_num < 1 {
        Err(ScaleError::CannotMakeScale)
    } else {
        Ok((mode_num + 1..=(2 * mode_num - 1))
            .map(|x| RawJiRatio::try_new(x, mode_num).expect("`numer` should be > `denom`"))
            .collect())
    }
}

/// Compute if `offset` is valid for an interleaved scale with strand `strand`.
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

/// Compute if `polyoffset` is valid for an interleaved scale with strand `strand`.
/// The polyoffset is given as a Vec of offsets from the unison.
pub fn is_valid_polyoffset(strand: &[RawJiRatio], polyoffset: &[RawJiRatio]) -> bool {
    let mut offsets_with_unison: Vec<RawJiRatio> = polyoffset.to_vec();
    offsets_with_unison.push(RawJiRatio::UNISON);
    // Check validity of every interval in the polyoffset.
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

/// Compute the set of valid polyoffsets for an interleaved scale on harmonic series mode `m`.
pub fn valid_polyoffsets(
    strand: &[RawJiRatio],
    m: u64,
) -> Result<Vec<Vec<RawJiRatio>>, ScaleError> {
    let mode = harmonic_mode_no_oct(m)?; // harmonic mode m
    Ok(crate::helpers::powerset(&mode)
        .into_iter()
        .filter(|x| is_valid_polyoffset(strand, x))
        .collect())
}

/// Compute the set of *maximal* valid polyoffsets for an interleaved scale on harmonic series mode `m`.
pub fn maximal_valid_polyoffsets(
    strand: &[RawJiRatio],
    m: u64,
) -> Result<Vec<Vec<RawJiRatio>>, ScaleError> {
    let mode = harmonic_mode_no_oct(m)?; // harmonic mode m
    let power_set = crate::helpers::powerset(&mode);
    let valid_polyoffsets: Vec<_> = power_set
        .into_iter()
        .filter(|x| is_valid_polyoffset(strand, x))
        .collect();
    Ok(valid_polyoffsets
        .clone()
        .into_iter()
        .filter(|x| crate::helpers::is_maximal_in(x.to_vec(), &valid_polyoffsets))
        .collect())
}

/// Compute an interleaved scale with the given strand and polyoffset, if it is valid.
/// A *strand* is a smaller scale that is duplicated and copies of it offset by the notes of the *polyoffset* chord.
/// The resulting scale made of interleaved strands is *interleaved* if all the strands are actually interleaved.
pub fn interleaved_scale_ji(
    strand: &[RawJiRatio],
    polyoffset: &[RawJiRatio],
) -> Result<Vec<RawJiRatio>, ScaleError> {
    if is_valid_polyoffset(strand, polyoffset) {
        let mut vec: Vec<RawJiRatio> = vec![];
        for offset in polyoffset {
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

/// Given `arr` a periodic JI scale given in JI ratios from the tonic, is `arr` CS?
/// Assumes arr[0] = the 1-step from the tonic, ..., arr[n-1] == the equave;
/// arr.len() == the scale size.
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

/// Return a well-formed necklace of `gen_class`-steps in a given JI scale.
pub fn well_formed_necklace_in_ji_scale(
    scale: &[RawJiRatio],
    gen_class: usize,
) -> Result<Vec<RawJiRatio>, ScaleError> {
    if gcd(scale.len() as i64, gen_class as i64) == 1 {
        Ok((0..(scale.len()))
            .map(|k| {
                (scale[(gen_class * (k + 1)) % scale.len()] / scale[(gen_class * k) % scale.len()])
                    .rd(RawJiRatio::OCTAVE)
            })
            .collect()) // For every `k` (including the last one), get the `k`th stacked `gen_class`-step.
    } else {
        Err(ScaleError::NonCoprimeGenError)
    }
}

#[cfg(test)]
mod tests {
    #[allow(unused)]
    use super::*;
    use crate::{ji::odd_limit, ji_ratio::RawJiRatio};

    use crate::monzo;

    #[test]
    fn test_ji_scale_modes() {
        let diasem_modes = ji_scale_modes(&RawJiRatio::TAS_9);
        let correct_diasem_modes = vec![
            vec![
                RawJiRatio::new(9, 8),
                RawJiRatio::new(7, 6),
                RawJiRatio::new(21, 16),
                RawJiRatio::new(4, 3),
                RawJiRatio::new(3, 2),
                RawJiRatio::new(14, 9),
                RawJiRatio::new(7, 4),
                RawJiRatio::new(16, 9),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::new(28, 27),
                RawJiRatio::new(7, 6),
                RawJiRatio::new(32, 27),
                RawJiRatio::new(4, 3),
                RawJiRatio::new(112, 81),
                RawJiRatio::new(14, 9),
                RawJiRatio::new(128, 81),
                RawJiRatio::new(16, 9),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::new(9, 8),
                RawJiRatio::new(8, 7),
                RawJiRatio::new(9, 7),
                RawJiRatio::new(4, 3),
                RawJiRatio::new(3, 2),
                RawJiRatio::new(32, 21),
                RawJiRatio::new(12, 7),
                RawJiRatio::new(27, 14),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::new(64, 63),
                RawJiRatio::new(8, 7),
                RawJiRatio::new(32, 27),
                RawJiRatio::new(4, 3),
                RawJiRatio::new(256, 189),
                RawJiRatio::new(32, 21),
                RawJiRatio::new(12, 7),
                RawJiRatio::new(16, 9),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::new(9, 8),
                RawJiRatio::new(7, 6),
                RawJiRatio::new(21, 16),
                RawJiRatio::new(4, 3),
                RawJiRatio::new(3, 2),
                RawJiRatio::new(27, 16),
                RawJiRatio::new(7, 4),
                RawJiRatio::new(63, 32),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::new(28, 27),
                RawJiRatio::new(7, 6),
                RawJiRatio::new(32, 27),
                RawJiRatio::new(4, 3),
                RawJiRatio::new(3, 2),
                RawJiRatio::new(14, 9),
                RawJiRatio::new(7, 4),
                RawJiRatio::new(16, 9),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::new(9, 8),
                RawJiRatio::new(8, 7),
                RawJiRatio::new(9, 7),
                RawJiRatio::new(81, 56),
                RawJiRatio::new(3, 2),
                RawJiRatio::new(27, 16),
                RawJiRatio::new(12, 7),
                RawJiRatio::new(27, 14),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::new(64, 63),
                RawJiRatio::new(8, 7),
                RawJiRatio::new(9, 7),
                RawJiRatio::new(4, 3),
                RawJiRatio::new(3, 2),
                RawJiRatio::new(32, 21),
                RawJiRatio::new(12, 7),
                RawJiRatio::new(16, 9),
                RawJiRatio::OCTAVE,
            ],
            vec![
                RawJiRatio::new(9, 8),
                RawJiRatio::new(81, 64),
                RawJiRatio::new(21, 16),
                RawJiRatio::new(189, 128),
                RawJiRatio::new(3, 2),
                RawJiRatio::new(27, 16),
                RawJiRatio::new(7, 4),
                RawJiRatio::new(63, 32),
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
                    &[RawJiRatio::new(7, 6), RawJiRatio::new(8, 7)],
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
                vec![RawJiRatio::new(7, 6), RawJiRatio::OCTAVE],
                vec![
                    RawJiRatio::new(7, 6),
                    RawJiRatio::new(4, 3),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::new(7, 6),
                    RawJiRatio::new(4, 3),
                    RawJiRatio::new(14, 9),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::new(7, 6),
                    RawJiRatio::new(4, 3),
                    RawJiRatio::new(14, 9),
                    RawJiRatio::new(16, 9),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::new(28, 27),
                    RawJiRatio::new(7, 6),
                    RawJiRatio::new(4, 3),
                    RawJiRatio::new(14, 9),
                    RawJiRatio::new(16, 9),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::new(28, 27),
                    RawJiRatio::new(7, 6),
                    RawJiRatio::new(32, 27),
                    RawJiRatio::new(4, 3),
                    RawJiRatio::new(14, 9),
                    RawJiRatio::new(16, 9),
                    RawJiRatio::OCTAVE,
                ],
                vec![
                    RawJiRatio::new(28, 27),
                    RawJiRatio::new(7, 6),
                    RawJiRatio::new(32, 27),
                    RawJiRatio::new(4, 3),
                    RawJiRatio::new(112, 81),
                    RawJiRatio::new(14, 9),
                    RawJiRatio::new(16, 9),
                    RawJiRatio::OCTAVE,
                ],
            ]
        );

        let zil_gs = [
            RawJiRatio::new(7, 4),
            RawJiRatio::new(12, 7),
            RawJiRatio::new(7, 4),
            RawJiRatio::new(12, 7),
            RawJiRatio::new(7, 4),
            RawJiRatio::new(12, 7),
            RawJiRatio::new(7, 4),
            RawJiRatio::new(320, 189),
            RawJiRatio::new(7, 4),
            RawJiRatio::new(12, 7),
        ];
        let zil_24 = gs_scale(&zil_gs, 24, RawJiRatio::OCTAVE).unwrap();
        let correct_zil_24 = vec![
            RawJiRatio::new(525, 512),
            RawJiRatio::new(135, 128),
            RawJiRatio::new(35, 32),
            RawJiRatio::new(9, 8),
            RawJiRatio::new(4725, 4096),
            RawJiRatio::new(75, 64),
            RawJiRatio::new(315, 256),
            RawJiRatio::new(5, 4),
            RawJiRatio::new(21, 16),
            RawJiRatio::new(675, 512),
            RawJiRatio::new(2835, 2048),
            RawJiRatio::new(45, 32),
            RawJiRatio::new(189, 128),
            RawJiRatio::new(3, 2),
            RawJiRatio::new(1575, 1024),
            RawJiRatio::new(405, 256),
            RawJiRatio::new(105, 64),
            RawJiRatio::new(27, 16),
            RawJiRatio::new(7, 4),
            RawJiRatio::new(225, 128),
            RawJiRatio::new(945, 512),
            RawJiRatio::new(15, 8),
            RawJiRatio::new(63, 32),
            RawJiRatio::OCTAVE,
        ];
        assert_eq!(zil_24, correct_zil_24);
    }
    #[test]
    fn test_odd_limit() {
        let mut eighty_one_odd_limit = odd_limit(81);
        eighty_one_odd_limit.sort();
        eighty_one_odd_limit.dedup();
        println!("{:?}", eighty_one_odd_limit);
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
            solve_step_sig_81_odd_limit(&[5, 2], Monzo::OCTAVE, false);
        assert_eq!(diatonic_solns, vec![vec![monzo![-3, 2], monzo![8, -5]]]);
        let blackdye_solns: Vec<Vec<Monzo>> =
            solve_step_sig_81_odd_limit(&[5, 2, 3], Monzo::OCTAVE, false);
        assert!(blackdye_solns.contains(&vec![
            monzo![1, -2, 1],
            monzo![4, -1, -1],
            monzo![-4, 4, -1]
        ]));
    }
}
