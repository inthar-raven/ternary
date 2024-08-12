use itertools::Itertools;
use serde::Serialize;
use std::cmp::{max, Ordering};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::hash::Hash;

use crate::utils::{gcd, modinv, ScaleError};

pub type Letter = usize;

/// Types such that when you put them into vectors, it makes sense to interpret the vectors as scales and take intervals from these scales.
pub trait Subtendable: Clone + Send + Sync + Sized {
    type Interval: Send;

    fn interval_from_slice(slice: &[Self]) -> Self::Interval;

    fn dyad_on_degree(scale: &[Self], degree: usize) -> Self::Interval {
        Self::interval_from_slice(&rotate(scale, degree))
    }
}

impl Subtendable for Letter {
    type Interval = CountVector<Letter>;
    fn interval_from_slice(slice: &[Self]) -> Self::Interval {
        CountVector::from_slice(slice)
    }
}

impl Subtendable for CountVector<Letter> {
    type Interval = CountVector<Letter>;
    fn interval_from_slice(slice: &[Self]) -> Self::Interval {
        slice
            .iter()
            .fold(CountVector::ZERO, |a, b| CountVector::add(&a, b))
    }
}

/// The [chirality](https://en.xen.wiki/w/Chirality) of a scale.
#[derive(Copy, Clone, Debug, Hash, PartialEq)]
pub enum Chirality {
    /// Lexicographically first mode is greater than that of reversed scale word.
    Left,
    /// Equal as circular word to reversed scale word.
    Achiral,
    /// Lexicographically first mode is less than that of reversed scale word.
    Right,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord, Serialize)]
/// Wrapper type for the free abelian group on T.
pub struct CountVector<T>(BTreeMap<T, i32>);

impl<T> CountVector<T> {
    /// Wraps BTreeMap<T, i32> in `CountVector<T>`.
    pub fn from_btree_map(m: BTreeMap<T, i32>) -> Self {
        CountVector(m)
    }
    /// Creates a zero count vector.
    pub const ZERO: Self = CountVector(BTreeMap::<T, i32>::new());

    /// Whether the `CountVector` is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The sum of the absolute values of the components.
    pub fn len(&self) -> usize {
        self.0.values().map(|v| v.unsigned_abs() as usize).sum()
    }
    /// The sum of two count vectors. Each component gets added.
    pub fn add(&self, w: &Self) -> Self
    where
        T: Ord + Clone,
    {
        let mut result = self.0.clone();
        for (key, value) in w.0.iter() {
            if let Some(update_this) = result.get_mut(key) {
                *update_this += value;
                if *update_this == 0 {
                    result.remove(key);
                }
            } else {
                result.insert((*key).clone(), *value);
            }
        }
        Self(result)
    }
    /// Additive inverse of a `CountVector`.
    pub fn neg(&self) -> Self
    where
        T: Ord + Clone + Send + Sync,
    {
        Self(self.0.iter().map(|(k, v)| (k.clone(), -v)).collect())
    }

    /// Multiply a CountVector by a scalar.
    pub fn scalar_mul(&self, lambda: i32) -> Self
    where
        T: Ord + Clone + Send + Sync,
    {
        Self(
            self.0
                .iter()
                .map(|(k, v)| (k.clone(), lambda * v))
                .collect(),
        )
    }

    /// Convert a slice to a CountVector.
    pub fn from_slice(slice: &[T]) -> Self
    where
        T: Ord + Clone,
    {
        let mut result = BTreeMap::new();
        for key in slice {
            if let Some(update_this) = result.get_mut(key) {
                *update_this += 1;
            } else {
                result.insert((*key).clone(), 1);
            }
        }
        Self(result)
    }
    /// Convert a `BTreeSet` to a `CountVector` (with every nonzero component equal to 1).
    pub fn from_btree_set(set: BTreeSet<T>) -> CountVector<T>
    where
        T: Ord + Send,
    {
        Self(
            set.into_iter()
                .map(|t| (t, 1))
                .collect::<BTreeMap<T, i32>>(),
        )
    }
    /// Convert a `BTreeSet` to a `CountVector` (with every nonzero component equal to 1).
    pub fn from_tuples(iter: impl Iterator<Item = (T, i32)>) -> CountVector<T>
    where
        T: Ord + Send,
    {
        Self::from_btree_map(BTreeMap::from_iter(iter))
    }
    /// Unwrap the underlying `BTreeMap`.
    pub fn into_inner(&self) -> BTreeMap<T, i32>
    where
        T: Clone,
    {
        self.clone().0
    }
    /// The multiset of `CountVector`s that occur as `subword_length`-steps in `scale`.
    pub fn spectrum(scale: &[T], subword_length: usize) -> CountVector<CountVector<T>>
    where
        T: Ord + Clone,
    {
        let mut result: BTreeMap<CountVector<T>, i32> = BTreeMap::new();
        for key in (0..scale.len())
            .map(|degree| CountVector::from_slice(&word_on_degree(scale, degree, subword_length)))
        {
            if let Some(update_this) = result.get_mut(&key) {
                *update_this += 1;
            } else {
                result.insert(key, 1);
            }
        }
        CountVector(result)
    }

    /// The set of `CountVector`s that occur as `subword_length`-steps in `scale`.
    pub fn distinct_spectrum(scale: &[T], subword_length: usize) -> HashSet<CountVector<T>>
    where
        T: Hash + Ord + Clone + Send + Sync,
    {
        (0..scale.len())
            .map(|degree| CountVector::from_slice(&word_on_degree(scale, degree, subword_length)))
            .collect::<HashSet<_>>()
    }
    /// Get the component for the specified key.
    pub fn get(&self, arg: &T) -> Option<&i32>
    where
        T: Ord,
    {
        self.0.get(arg)
    }
}

/// Treating `scale` as a circular string (that is, "scale[i] == scale[i % scale.len()]"),
/// take a slice of length `subword_length` from `degree`; assumes `subword_length` <= `scale`.len().
/// Reduce `degree` first.
pub fn word_on_degree<T>(scale: &[T], degree: usize, subword_length: usize) -> Vec<T>
where
    T: Clone,
{
    let rotated = rotate(scale, degree);
    if subword_length < scale.len() {
        Vec::<_>::from(&rotated[..subword_length])
    } else {
        let prefix_for_div = rotated
            .iter()
            .cloned()
            .cycle()
            .take(subword_length / scale.len())
            .collect::<Vec<_>>();
        let suffix_for_rem = &rotated[0..subword_length];
        [&prefix_for_div, suffix_for_rem].concat()
    }
}

/// The dyad on the specified degree of the scale as a `CountVector`.
pub fn dyad_on_degree<T>(scale: &[T], degree: usize, interval_class: usize) -> CountVector<T>
where
    T: Ord + Clone,
{
    CountVector::from_slice(&word_on_degree(scale, degree, interval_class))
}

#[allow(unused)]
fn strand_on_degree<T>(
    scale: &[T],
    degree: usize,
    interval_class: usize,
) -> Result<Vec<CountVector<T>>, ScaleError>
where
    T: Ord + Clone + Send + Sync,
{
    if scale.len() % interval_class == 0 {
        Ok((0..(scale.len() / interval_class))
            .map(|i| dyad_on_degree(scale, degree + i * interval_class, interval_class))
            .collect())
    } else {
        Err(ScaleError::NonDivisibleSubsetError)
    }
}

/// Return the step class of the offset of scale2 to the right relative to scale1 if strings are conjugate, otherwise None.
pub fn offset_string<T>(scale1: &[T], scale2: &[T]) -> Option<usize>
where
    T: Clone + PartialEq,
{
    if scale1.len() != scale2.len() {
        None
    } else {
        let length = scale1.len();
        for i in 0..length {
            if rotate(scale2, i) == scale1.to_vec() {
                // scale2 has been rotated i steps to the left.
                return Some(i);
            }
        }
        None
    }
}

/// Return the offset of vec2 to the right relative to vec1 if the Vecs are conjugate, otherwise None.
pub fn offset_vec<T>(vec1: &[T], vec2: &[T]) -> Option<usize>
where
    T: PartialEq + std::clone::Clone,
{
    if vec1.len() != vec2.len() {
        // if lengths are not equal, cannot be conjugate
        None
    } else {
        let len = vec1.len();
        for i in 0..len {
            let vec2_rotated = rotate(vec2, i); // vec2 has been rotated i steps to the left.
            if vec2_rotated == vec1 {
                return Some(i);
            }
        }
        None
    }
}

/// Computes the maximum variety of a scale word. (Not actual interval sizes)
pub fn maximum_variety<T>(scale: &[T]) -> usize
where
    T: Hash + Ord + Clone + Sync + Send,
{
    let mut result = 0;
    let floor_half: usize = scale.len() / 2;
    for subword_length in 1..(floor_half + 1) {
        let sizes = CountVector::distinct_spectrum(scale, subword_length);
        result = max(result, sizes.len()); // update result
    }
    result
}
/// Whether `scale` as a word is strict variety.
pub fn is_strict_variety<T>(scale: &[T]) -> bool
where
    T: Hash + Ord + Clone + Sync + Send,
{
    let mut prev: usize = 0;
    let floor_half: usize = scale.len() / 2;
    for subword_length in 1..=floor_half {
        let sizes: HashSet<CountVector<T>> = CountVector::distinct_spectrum(scale, subword_length);
        if prev == 0 {
            prev = sizes.len();
        } else if prev != sizes.len() {
            return false;
        }
    }
    true
}

/// Return the block_balance of `s` = `s`(x, y, z), defined as max { | |w|_{x_i} - |w'|_{x_i} | : x_i is a letter of `s` and k = len(w) = len(w') }.
pub fn block_balance<T>(scale: &[T]) -> usize
where
    T: Hash + Ord + Clone + PartialEq + Send + Sync,
{
    if scale.len() <= 1 {
        scale.len()
    } else {
        let mut max = 0;
        let floor_half: usize = scale.len() / 2;
        let distinct_letters = scale.iter().cloned().collect::<BTreeSet<T>>();
        for subword_length in 1..=floor_half {
            for letter in distinct_letters.iter() {
                let counts = CountVector::distinct_spectrum(scale, subword_length)
                    .iter()
                    .filter_map(|dyad| dyad.get(letter))
                    .copied()
                    .collect::<BTreeSet<_>>();
                let diff: i32 = *counts.last().expect("`counts` should be nonempty")
                    - *counts.first().expect("`counts` should be nonempty"); // this will always be >= 0 for a nonempty `BTreeSet`
                if diff > max {
                    max = diff;
                }
            }
        }
        max as usize
    }
}

/// Return the darkest mode of the MOS axby and the dark generator, using the Bresenham line algorithm.
/// We chose the darkest mode rather than the brightest because this is the mode with brightness == 0.
pub fn darkest_mos_mode_and_gen_bresenham(
    a: usize,
    b: usize,
) -> (Vec<Letter>, CountVector<Letter>) {
    let d = gcd(a as u64, b as u64) as usize;
    if d == 1 {
        let count_gen_steps = modinv(a as i64, a as i64 + b as i64)
                .expect("The dark generator is a (|L|⁻¹ mod |scale|)-step, since stacking it |L| times results in the s step (mod period).")
                as usize;
        let mut result_scale: Vec<usize> = vec![];
        let (mut current_x, mut current_y) = (0usize, 0usize); // Start from the (0, 0) and walk until the dark generator is reached; we now know how many steps to walk.
        while current_x < a || current_y < b {
            if a * (current_y) >= b * (current_x + 1) {
                // If going east (making a (1, 0) step) doesn't lead to going below the line y == b/a*x,
                current_x += 1; // append the x step and reflect the change in the plane vector.
                result_scale.push(0);
            } else {
                // Else, make a (0, 1) step.
                current_y += 1;
                result_scale.push(1);
            }
        }
        let result_gen = CountVector::from_slice(&result_scale[0..count_gen_steps]);
        (result_scale, result_gen)
    } else {
        let (prim_mos, gen) = darkest_mos_mode_and_gen_bresenham(a / d, b / d);
        (prim_mos.repeat(d), gen)
    }
}

/// Return the darkest mode of the MOS axby and the dark generator, using Bjorklund's algorithm.
/// We chose the darkest mode rather than the brightest because this is the mode with brightness == 0.
pub fn darkest_mos_mode_and_gen_bjorklund(
    a: usize,
    b: usize,
) -> (Vec<Letter>, CountVector<Letter>) {
    let d = gcd(a as u64, b as u64) as usize;
    if d == 1 {
        // These are the seed strings we build the brightest MOS word from.
        // The algorithm uses two subwords at each step, iteratively appending the
        // lexicographically second subword to the lexicographically first subword to ensure
        // that the lexicographically first mode is returned.
        // Since lex. order makes strict prefixes strictly "brighter",
        // we find the brightest mode and reverse it.
        let (mut first, mut second) = (vec![0], vec![1]);
        let (mut count_first, mut count_second) = (a, b); // axby(0) is the brightest mode of bxay with x and y swapped.
        while count_second != 1 {
            // Possibly after switching, are there more copies of `first` than `second`?
            // Then all the `second`s get appended to the first `count_second` copies of `first`s,
            // and the new `second`s are the remaining copies of `first`s.
            let old_first = first.clone();
            first.extend_from_slice(&second);
            if count_first > count_second {
                second = old_first;
                (count_first, count_second) = (count_second, count_first - count_second);
            }
            // Otherwise, there are strictly fewer `first`s than `second`s (as gcd(a, b) == 1),
            // and *all* the `first`s get modified, whereas `second` is unchanged since copies of it remain.
            // `count_first` is also unchanged.
            else {
                count_second -= count_first;
            }
            // At the current step we have `count_first` `first` substrings and `count_second` `second` substrings,
            // where we must guarantee that `first < second`.
            // Thus if `first > second`, then swap them and swap the count variables.
            // Do this step before checking the while condition; we know the desired lex. ordering holds for the first step,
            // and our stopping condition requires that `first < second` actually hold to really behave correctly.
            if first > second {
                (first, second) = (second, first);
                (count_first, count_second) = (count_second, count_first);
            }
        }
        // At the end, we have `count_first` `first`s and 1 `second`,
        // so first substitute 'z's for 'x's, and then
        // return (`first`)^`count_first` `second` (in standard mathematical word notation).
        // The dark generator is obtained by subtracting the Parikh vector of `first` form that of the MOS scale.
        let gen = CountVector(BTreeMap::from_iter(
            [a as i32, b as i32].into_iter().enumerate(),
        ))
        .add(&CountVector::from_slice(&first).neg()); // Do this before consuming `first` in the next line.
        let mut scale: Vec<usize> = first.repeat(count_first);
        scale.extend_from_slice(&second);
        // Reverse the scale word, since Bjorklund's algorithm has given us the brightest mode.
        scale.reverse();
        (scale, gen)
    } else {
        let (primitive_mos, gen) = darkest_mos_mode_and_gen_bjorklund(a / d, b / d);
        (primitive_mos.repeat(d), gen)
    }
}

/// The mode of the MOS aLbs with a given brightness (count of bright generators up from root).
pub fn mos_mode(a: usize, b: usize, brightness: usize) -> Result<Vec<Letter>, ScaleError> {
    if brightness >= a + b {
        Err(ScaleError::CannotMakeScale)
    } else {
        let (mos, dark_gen) = darkest_mos_mode_and_gen_bresenham(a, b);
        let dark_gen_step_count: usize = dark_gen.len();
        Ok(rotate(&mos, brightness * dark_gen_step_count))
    }
}

/// Rotate an array. Returns a Vec.
pub fn rotate<T: std::clone::Clone>(slice: &[T], degree: usize) -> Vec<T> {
    let degree = degree % slice.len();
    if degree == 0 {
        slice.to_vec()
    } else {
        [&slice[degree..slice.len()], &slice[0..degree]].concat()
    }
}

pub fn are_conjugate<T>(s1: &[T], s2: &[T]) -> bool
where
    T: Clone + Eq,
{
    (s1.len() == s2.len()) && { (0..s1.len()).any(|i| rotate(s1, i) == s2.to_vec()) }
}

/// [Letterwise substitution](https://en.xen.wiki/w/MOS_substitution) for scale words.
///
/// Note: This function does not fail even if the number of times `x` occurs in `template`
/// does not divide `filler.len()`.
pub fn subst<T>(template: &[T], x: &T, filler: &[T]) -> Vec<T>
where
    T: PartialEq + Clone,
{
    let mut ret = vec![];
    let mut i: usize = 0;
    if !filler.is_empty() {
        for letter in template {
            if *letter == *x {
                ret.push(filler[i % filler.len()].clone());
                i += 1;
            } else {
                ret.push((*letter).clone());
            }
        }
    } else {
        return template
            .iter()
            .filter(|letter| **letter != *x)
            .map(|y| (*y).clone())
            .collect();
    }
    ret
}

/// Return the collection of all MOS substitution scales `subst n0 x (n1 y n2 z)`
/// where the template MOS is assumed to have step signature `n0*0 (n1 + n2)*X` (`X` is the slot letter)
/// and the filling MOS has step signature `n1*1 n2*2`.
fn mos_substitution_scales_one_perm(n0: usize, n1: usize, n2: usize) -> Vec<Vec<Letter>> {
    let (template, _) = darkest_mos_mode_and_gen_bresenham(n0, n1 + n2);
    let (filler, gener) = darkest_mos_mode_and_gen_bresenham(n1, n2);
    let filler = filler.into_iter().map(|x| x + 1).collect::<Vec<_>>();
    let gener_size = gener.len();
    (0..(n1 + n2))
        .map(|i| {
            subst::<usize>(
                &template,
                &1usize,
                &rotate(&filler, (i * gener_size) % filler.len()),
            )
        })
        .collect()
}

/// The lexicographically brightest mode of a word (where the letters are in their usual order).
/// (Booth's algorithm)
pub fn least_mode(scale: &[Letter]) -> Vec<Letter> {
    rotate(scale, booth(scale))
}

/// The lexicographically brightest mode of a word (where the letters are in their usual order).
/// (Booth's algorithm)
pub fn booth(scale: &[Letter]) -> usize {
    let s = scale;
    let n = scale.len();
    let mut f = vec![usize::MAX; 2 * n];
    let mut k: usize = 0;
    for j in 1..2 * n {
        let mut i = f[j - k - 1];
        while i != usize::MAX && s[j % n] != s[k.wrapping_add(i).wrapping_add(1) % n] {
            if s[j % n] < s[k.wrapping_add(i).wrapping_add(1) % n] {
                k = j.wrapping_sub(i).wrapping_sub(1);
            }
            i = f[i];
        }
        if i == usize::MAX && s[j % n] != s[k.wrapping_add(i).wrapping_add(1) % n] {
            if s[j % n] < s[k.wrapping_add(i).wrapping_add(1) % n] {
                k = j;
            }
            f[j - k] = usize::MAX;
        } else {
            f[j - k] = i.wrapping_add(1);
        }
    }
    k
}

/// The set of all [MOS substitution](https://en.xen.wiki/w/User:Inthar/MOS_substitution) ternary scales.
pub fn mos_substitution_scales(sig: &[usize]) -> Vec<Vec<Letter>> {
    let (n0, n1, n2) = (sig[0], sig[1], sig[2]);

    // Only need 3 permutations of (0, 1, 2) for the MOS substitution patterns n0*_ (n1*_ n2*_)
    let redundant_list = [
        // n0L (n1m n2s)
        mos_substitution_scales_one_perm(n0, n1, n2),
        // n1m (n0L n2s)
        mos_substitution_scales_one_perm(n1, n2, n0)
            .into_iter()
            .map(|scale| scale.into_iter().map(|x| (x + 1) % 3).collect())
            .collect(),
        // n2s (n0L n1m)
        mos_substitution_scales_one_perm(n2, n0, n1)
            .into_iter()
            .map(|scale| {
                scale
                    .into_iter()
                    .map(|x| if x == 0 { 2 } else { (x - 1) % 3 })
                    .collect()
            })
            .collect(),
    ]
    .concat();
    // Canonicalize every scale and remove duplicates
    redundant_list
        .into_iter()
        .map(|scale| least_mode(&scale))
        .sorted()
        .dedup()
        .collect()
}

/// Return the number of distinct steps in `scale`.
pub fn step_variety(scale: &[Letter]) -> usize {
    scale.iter().collect::<BTreeSet<_>>().len()
}

/// subst but single letters.
pub fn replace(scale: &[Letter], from: Letter, to: Letter) -> Vec<Letter> {
    subst(scale, &from, &[to])
}

/// Delete all instances of one letter.
pub fn delete(scale: &[Letter], letter: Letter) -> Vec<Letter> {
    scale.iter().filter(|x| **x != letter).copied().collect()
}

/// If `scale` is ternary, return whether identifying L = m, m = s, and s = 0 results in a MOS.
/// Returns `false` if the scale is not ternary.
pub fn is_monotone_mos(scale: &[Letter]) -> bool {
    step_variety(scale) == 3
        && maximum_variety(&replace(scale, 1, 0)) == 2 // L = m
        && maximum_variety(&replace(scale, 2, 1)) == 2 // m = s
        && maximum_variety(&delete(scale, 2)) == 2 // s = 0
}

/// Check if the result of equating L = m is a MOS. Assumes the scale is ternary.
pub fn monotone_lm(scale: &[Letter]) -> bool {
    maximum_variety(&replace(scale, 1, 0)) == 2
}

/// Check if the result of equating m = s is a MOS. Assumes the scale is ternary.
pub fn monotone_ms(scale: &[Letter]) -> bool {
    maximum_variety(&replace(scale, 2, 1)) == 2
}

/// Check if the result of equating s = 0 is a MOS. Assumes the scale is ternary.
pub fn monotone_s0(scale: &[Letter]) -> bool {
    maximum_variety(&delete(scale, 2)) == 2
}

/// Check if pairiwse identifications of two of the step sizes always results in a MOS.
/// Returns `false` if the scale is not ternary.
pub fn is_pairwise_mos(scale: &[Letter]) -> bool {
    step_variety(scale) == 3
        && maximum_variety(&replace(scale, 1, 0)) == 2
        && maximum_variety(&replace(scale, 1, 2)) == 2
        && maximum_variety(&replace(scale, 2, 0)) == 2
}

/// The repeating portion of a slice.
pub fn period<T>(slice: &[T]) -> Vec<T>
where
    T: PartialEq + Clone + Send + Sync,
{
    let l = 1
        + (1..slice.len())
            .map(|i| slice.iter().take(i).cycle())
            .enumerate()
            .take_while(|(i, iter)| {
                slice.to_vec()
                    != iter
                        .clone()
                        .take((i + 1) * (slice.len() / (i + 1))) // Take one more than the index, then repeat
                        .cloned()
                        .collect::<Vec<T>>()
            })
            .collect::<Vec<_>>()
            .len();
    slice[..l].to_vec()
}

/// The repeating portion of `slice`,
/// but `slice` is only required to be equal to some prefix of the infinite repetition of the weak period.
/// For example, `[0, 1]` is a weak period but not a strong period of `[0, 1, 0, 1, 0]`.
pub fn weak_period<T>(slice: &[T]) -> Vec<T>
where
    T: PartialEq + Clone,
{
    let l = 1
        + (1..slice.len())
            .map(|i| slice.iter().take(i).cycle())
            .take_while(|iter| {
                slice.to_vec() != iter.clone().take(slice.len()).cloned().collect::<Vec<T>>()
            })
            .collect::<Vec<_>>()
            .len();
    slice[..l].to_vec()
}

/// The collection of rotations of a word. Contains redundant rotations if the word is not primitive.
pub fn rotations<T>(word: &[T]) -> Vec<Vec<T>>
where
    T: Clone + Eq,
{
    (0..word.len()).map(|i| rotate(word, i)).collect()
}

/// The chirality of a scale word.
pub fn chirality(word: &[Letter]) -> Chirality {
    let least_mode_word = least_mode(word);

    let word_rev: Vec<usize> = word.iter().cloned().rev().collect();
    let least_mode_word_rev = least_mode(&word_rev);

    match least_mode_word.cmp(&least_mode_word_rev) {
        Ordering::Less => Chirality::Right,
        Ordering::Equal => Chirality::Achiral,
        Ordering::Greater => Chirality::Left,
    }
}

#[cfg(test)]
mod tests {
    #[allow(unused)]
    use crate::utils::gcd;

    use super::*;

    #[test]
    fn test_booth() {
        let blackdye = [2, 0, 1, 0, 2, 0, 1, 0, 2, 0];
        assert_eq!(least_mode(&blackdye), vec![0, 1, 0, 2, 0, 1, 0, 2, 0, 2]);
    }
    #[test]
    fn test_word_on_degree() {
        let blackdye = [2, 0, 1, 0, 2, 0, 1, 0, 2, 0];
        let even_two_steps = (0..5)
            .map(|i| word_on_degree(&blackdye, 2 * i, 2))
            .collect::<Vec<_>>();
        assert_eq!(
            even_two_steps,
            vec![vec![2, 0], vec![1, 0], vec![2, 0], vec![1, 0], vec![2, 0],]
        );
        let odd_two_steps = (0..5)
            .map(|i| word_on_degree(&blackdye, 2 * i + 1, 2))
            .collect::<Vec<_>>();
        assert_eq!(
            odd_two_steps,
            vec![vec![0, 1], vec![0, 2], vec![0, 1], vec![0, 2], vec![0, 2],]
        );
    }
    #[test]
    fn test_spectrum() {
        let diachrome_5sc = [0, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 1];
        for i in 1..=11 {
            assert_eq!(
                // .len() of `CountVector` is the taxicab norm, not the number of keys
                CountVector::spectrum(&diachrome_5sc, i).into_inner().len(),
                CountVector::distinct_spectrum(&diachrome_5sc, i).len()
            );
        }
    }
    #[test]
    fn test_chirality() {
        let diasem_2sr = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let blackdye = [2, 0, 1, 0, 2, 0, 1, 0, 2, 0];
        let diasem_2sl = [0, 2, 0, 1, 0, 2, 0, 1, 0];

        assert_eq!(chirality(&diasem_2sr), Chirality::Right);
        assert_eq!(chirality(&blackdye), Chirality::Achiral);
        assert_eq!(chirality(&diasem_2sl), Chirality::Left);

        let diachrome_5sr = [0, 1, 2, 0, 2, 0, 2, 0, 1, 2, 0, 2];
        let diachrome_5sc = [0, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 1];
        let diachrome_5sl = [0, 2, 1, 0, 2, 0, 2, 0, 2, 1, 0, 2];

        assert_eq!(chirality(&diachrome_5sr), Chirality::Right);
        assert_eq!(chirality(&diachrome_5sc), Chirality::Achiral);
        assert_eq!(chirality(&diachrome_5sl), Chirality::Left);
    }

    #[test]
    fn test_period() {
        let word_012 = [0, 1, 2, 0, 1, 2, 0, 1, 2];
        let pentawood = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1];
        let diasem = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        assert_eq!(period(&word_012).len(), 3);
        assert_eq!(period(&pentawood).len(), 2);
        assert_eq!(period(&diasem).len(), 9);
    }
    #[test]
    fn test_weak_period() {
        let word_012 = [0, 1, 2, 0, 1, 2, 0, 1, 2];
        let pentawood = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1];
        let diasem = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        assert_eq!(weak_period(&word_012).len(), 3);
        assert_eq!(weak_period(&pentawood).len(), 2);
        assert_eq!(weak_period(&diasem).len(), 4);
    }
    #[test]
    fn test_maximum_variety() {
        assert_eq!(maximum_variety::<usize>(&[]), 0);
        assert_eq!(maximum_variety(&[0, 0, 0, 0]), 1); // check even and odd
        assert_eq!(maximum_variety(&[1, 1, 1, 1, 1]), 1);
        assert_eq!(maximum_variety(&[0, 0, 0, 1, 0, 0, 1]), 2); // diatonic has max variety 2
        assert_eq!(maximum_variety(&[1, 1, 1, 1, 2, 1, 2]), 3); // altered diatonic has max variety 3
        assert_eq!(maximum_variety(&[0, 1, 0, 2, 0, 1, 0, 2, 0]), 3); // diasem has max variety 3
        assert_eq!(maximum_variety(&[0, 1, 0, 2, 0, 1, 0, 2, 0, 1]), 4); // blackdye has max variety 4
                                                                         // MOS scales should be MV2.
        for a in 1usize..=10 {
            for b in 1usize..=10 {
                if gcd(a as u64, b as u64) == 1 {
                    for br in 0..(a + b) / (gcd(a as u64, b as u64) as usize) {
                        let mos = mos_mode(a, b, br);
                        assert_eq!(maximum_variety(&mos.unwrap()), 2);
                    }
                }
            }
        }
    }

    #[test]
    fn test_bjorklund_and_bresenham() {
        // Bjorklund and Bresenham should agree.
        for a in 1usize..=20 {
            for b in 1usize..=20 {
                if gcd(a as u64, b as u64) == 1 {
                    let mos_bjorklund = darkest_mos_mode_and_gen_bjorklund(a, b);
                    let mos_bresenham = darkest_mos_mode_and_gen_bresenham(a, b);
                    assert_eq!(mos_bjorklund, mos_bresenham);
                }
            }
        }
    }
    #[test]
    fn test_mos_block_balanced() {
        // MOS scales should have block balance 1.
        for a in 3..=20 {
            for b in 3..=20 {
                if gcd(a, b) == 1 {
                    for br in 0..(a + b) / gcd(a, b) {
                        let mos = mos_mode(a as usize, b as usize, br as usize)
                            .expect("`br` should be in range `0..(a + b)/gcd(a, b)`");
                        assert_eq!(block_balance(&mos), 1);
                    }
                }
            }
        }
    }
}
