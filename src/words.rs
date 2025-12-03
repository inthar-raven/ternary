use itertools::Itertools;
use serde::Serialize;
use std::cmp::{Ordering, max};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::hash::Hash;

use crate::helpers::{ScaleError, gcd, modinv};

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
#[derive(Copy, Clone, Debug, Hash, PartialEq, Serialize)]
pub enum Chirality {
    /// Lexicographically first mode is greater than that of reversed scale word.
    Left,
    /// Equal as circular word to reversed scale word.
    Achiral,
    /// Lexicographically first mode is less than that of reversed scale word.
    Right,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
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
            match result.get_mut(key) {
                Some(update_this) => {
                    *update_this += value;
                    if *update_this == 0 {
                        result.remove(key);
                    }
                }
                _ => {
                    result.insert((*key).clone(), *value);
                }
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
            match result.get_mut(key) {
                Some(update_this) => {
                    *update_this += 1;
                }
                _ => {
                    result.insert((*key).clone(), 1);
                }
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
            match result.get_mut(&key) {
                Some(update_this) => {
                    *update_this += 1;
                }
                _ => {
                    result.insert(key, 1);
                }
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

pub fn countvector_to_slice(v: CountVector<usize>) -> Vec<i32> {
    if v.is_empty() {
        vec![]
    } else if let Some(max_key_value) = v.into_inner().last_key_value() {
        let max = max_key_value.0;
        let mut result = vec![0; max + 1];
        for key in v.into_inner().keys() {
            if let Some(value) = v.get(key) {
                result[*key] = *value;
            } else {
                return vec![];
            }
        }
        result
    } else {
        vec![]
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
    if scale.len().is_multiple_of(interval_class) {
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

/// Return the block balance of `s` = `s`(x, y, z), defined as max { | |w|_{x_i} - |w'|_{x_i} | : x_i is a letter of `s` and k = len(w) = len(w') }.
pub fn block_balance<T>(scale: &[T]) -> usize
where
    T: Hash + Ord + Clone + PartialEq + Send + Sync,
{
    if scale.len() <= 1 {
        scale.len()
    } else {
        // Only need to check half of the step size classes.
        let floor_half: usize = scale.len() / 2;
        let distinct_letters = scale.iter().cloned().collect::<BTreeSet<T>>();
        let mb_max = (1..=floor_half)
            .flat_map(|subword_length| {
                distinct_letters.iter().map(move |letter| {
                    let counts = CountVector::distinct_spectrum(scale, subword_length)
                        .iter()
                        .filter_map(|dyad| dyad.get(letter))
                        .copied()
                        .collect::<BTreeSet<_>>();
                    // the differences to collect
                    *counts.last().expect("`counts` should be nonempty")
                        - *counts.first().expect("`counts` should be nonempty") // this will always be >= 0 for a nonempty `BTreeSet`
                })
            })
            .max();
        if let Some(max) = mb_max {
            max as usize
        } else {
            usize::MAX
        }
    }
}

/// Return the brightest mode of the MOS aLbs and the bright generator, using the Bresenham line algorithm.
/// The brightest mode is the lexicographically first rotation.
pub fn brightest_mos_mode_and_gen_bresenham(
    a: usize,
    b: usize,
) -> (Vec<Letter>, CountVector<Letter>) {
    let d = gcd(a as u64, b as u64) as usize;
    if d == 1 {
        let count_gener_steps = modinv(a as i64, a as i64 + b as i64)
                .expect("The bright generator is a (|L|⁻¹ mod |scale|)-step, since stacking it |L| times results in the s step (mod period).")
                as usize;
        let mut result_scale: Vec<usize> = vec![];
        // Start from the origin (0, 0)
        let (mut current_x, mut current_y) = (0usize, 0usize);
        while (current_x, current_y) != (a, b) {
            if a * (current_y + 1) <= b * current_x {
                // If going north (making a (0, 1) step) doesn't lead to going above the line y == b/a*x,
                current_y += 1; // append the y step and reflect the change in the plane vector.
                result_scale.push(1);
            } else {
                // Else, make a (1, 0) step.
                current_x += 1;
                result_scale.push(0);
            }
        }
        // Get the bright generator. We know how many steps and that this will give the perfect generator, not the
        // diminished one, since we just got the brightest mode.
        let result_gener = CountVector::from_slice(&result_scale[0..count_gener_steps]);
        (result_scale, result_gener)
    } else {
        let (prim_mos, gener) = brightest_mos_mode_and_gen_bresenham(a / d, b / d);
        (prim_mos.repeat(d), gener)
    }
}

/// Return the brightest mode of the MOS aLbs and the bright generator, using Bjorklund's algorithm.
/// The brightest mode is the lexicographically first rotation.
pub fn brightest_mos_mode_and_gen_bjorklund(
    a: usize,
    b: usize,
) -> (Vec<Letter>, CountVector<Letter>) {
    let d = gcd(a as u64, b as u64) as usize;
    if d == 1 {
        let count_gener_steps = modinv(a as i64, a as i64 + b as i64)
                .expect("The bright generator is a (|L|⁻¹ mod |scale|)-step, since stacking it |L| times results in the s step (mod period).")
                as usize;
        // These are the seed strings we build the brightest MOS word from.
        // The algorithm uses two subwords at each step, iteratively appending the
        // lexicographically second subword to the lexicographically first subword to ensure
        // that the lexicographically first mode is returned.
        let (mut first, mut second) = (vec![0], vec![1]);
        let (mut count_first, mut count_second) = (a, b); // aLbs(0) is the brightest mode of bsaL with L and s swapped.
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
        // so return (`first`)^`count_first` `second` (in standard mathematical word notation).
        let mut scale: Vec<usize> = first.repeat(count_first);
        scale.extend_from_slice(&second);
        // The bright generator is the first `count_gener_steps` of the scale.
        let gener = CountVector::from_slice(&scale[0..count_gener_steps]);
        (scale, gener)
    } else {
        let (primitive_mos, gener) = brightest_mos_mode_and_gen_bjorklund(a / d, b / d);
        (primitive_mos.repeat(d), gener)
    }
}

/// The mode of the MOS aLbs with a given brightness (count of bright generators up from root).
/// Brightness is taken modulo (a + b), so any non-negative value is valid.
/// Brightness 0 returns the darkest mode, brightness (a + b - 1) returns the brightest mode.
pub fn mos_mode(a: usize, b: usize, brightness: usize) -> Vec<Letter> {
    let scale_len = a + b;
    let brightness = brightness % scale_len;
    let (mos, bright_gen) = brightest_mos_mode_and_gen_bresenham(a, b);
    let bright_gen_step_count: usize = bright_gen.len();
    // Rotate backwards from brightest mode by (scale_len - 1 - brightness) bright generators
    // which is equivalent to rotating forward by brightness dark generators from darkest mode
    let steps_from_brightest = (scale_len - 1 - brightness) * bright_gen_step_count;
    rotate(&mos, steps_from_brightest)
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

/// Whether two slices with elements of type T are rotationally equivalent.
pub fn rotationally_equivalent<T>(s1: &[T], s2: &[T]) -> bool
where
    T: Clone + Eq,
{
    (s1.len() == s2.len()) && { (0..s1.len()).any(|i| rotate(s1, i) == s2.to_vec()) }
}

/// The lexicographically least mode of a word (where the letters are in their usual order).
pub fn least_mode(scale: &[Letter]) -> Vec<Letter> {
    rotate(scale, booth(scale))
}

/// The rotation required from the current word to the
/// lexicographically least mode of a word.
/// Booth's algorithm requires at most 3*n* comparisons and *n* storage locations where *n* is the input word's length.
/// See Booth, K. S. (1980). Lexicographically least circular substrings.
/// Information Processing Letters, 10(4-5), 240–242. doi:10.1016/0020-0190(80)90149-0
pub fn booth(scale: &[Letter]) -> usize {
    let scale_len = scale.len();
    // `failure_func` is the failure function of the least rotation; `usize::MAX` is used as a null value.
    // null indicates that the failure function does not point backwards in the string.
    // `usize::MAX` will behave the same way as -1 does, assuming wrapping unsigned addition
    let mut failure_func = vec![usize::MAX; 2 * scale_len];
    let mut least_rotation: usize = 0;
    // `scan_pos` loops over `scale` twice.
    for scan_pos in 1..2 * scale_len {
        let mut match_len = failure_func[scan_pos - least_rotation - 1];
        while match_len != usize::MAX
            && scale[scan_pos % scale_len]
                != scale[least_rotation.wrapping_add(match_len).wrapping_add(1) % scale_len]
        {
            // (1) If the scan_pos-th letter is less than s[(least_rotation + match_len + 1) % scale_len] then change least_rotation to scan_pos - match_len - 1,
            // in effect left-shifting the failure function and the input string.
            // This appropriately compensates for the new, shorter least substring.
            if scale[scan_pos % scale_len]
                < scale[least_rotation.wrapping_add(match_len).wrapping_add(1) % scale_len]
            {
                least_rotation = scan_pos.wrapping_sub(match_len).wrapping_sub(1);
            }
            match_len = failure_func[match_len];
        }
        if match_len == usize::MAX
            && scale[scan_pos % scale_len]
                != scale[least_rotation.wrapping_add(match_len).wrapping_add(1) % scale_len]
        {
            // See note (1) above.
            if scale[scan_pos % scale_len]
                < scale[least_rotation.wrapping_add(match_len).wrapping_add(1) % scale_len]
            {
                least_rotation = scan_pos;
            }
            failure_func[scan_pos - least_rotation] = usize::MAX;
        } else {
            failure_func[scan_pos - least_rotation] = match_len.wrapping_add(1);
        }
        // The induction hypothesis is that
        // at this point `failure_func[0..scan_pos - least_rotation]` is the failure function of `s[least_rotation..(least_rotation+scan_pos)%scale_len]`,
        // and `least_rotation` is the lexicographically least subword of the letters scanned so far.
    }
    least_rotation
}

/// [Letterwise substitution](https://en.xen.wiki/w/MOS_substitution) for scale words.
///
/// Note: This function does not fail even if the number of times `x` occurs in `template`
/// does not divide `filler.len()`.
pub fn subst(template: &[Letter], x: Letter, filler: &[Letter]) -> Vec<Letter> {
    let mut ret = vec![];
    let mut i: usize = 0;
    if !filler.is_empty() {
        for &letter in template {
            if letter == x {
                // Use the currently pointed-to letter of `filler` in place of `x`.
                ret.push(filler[i % filler.len()]);
                // Only update `i` when an `x` is replaced.
                i += 1;
            } else {
                ret.push(letter);
            }
        }
    } else {
        // If `filler` is empty, we return `template` but with all `x`s removed.
        return delete(template, x);
    }
    ret
}

/// Return the collection of all MOS substitution scales `subst n0 x (n1 y n2 z)`
/// where the template MOS is assumed to have step signature `n0*0 (n1 + n2)*X` (`X` is the slot letter)
/// and the filling MOS has step signature `n1*1 n2*2`.
fn mos_substitution_scales_one_perm(n0: usize, n1: usize, n2: usize) -> Vec<Vec<Letter>> {
    let (template, _) = brightest_mos_mode_and_gen_bresenham(n0, n1 + n2);
    let (filler, gener) = brightest_mos_mode_and_gen_bresenham(n1, n2);
    let filler = filler.into_iter().map(|x| x + 1).collect::<Vec<_>>();
    let gener_size = gener.len();
    (0..(n1 + n2))
        .map(|i| {
            subst(
                &template,
                1usize,
                &rotate(&filler, (i * gener_size) % filler.len()),
            )
        })
        .collect()
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

/// Whether `scale` is a MOS substitution scale subst at(bf1 cf2).
pub fn is_mos_subst(scale: &[Letter], t: Letter, f1: Letter, f2: Letter) -> bool {
    step_variety(scale) == 3 // Is it ternary?
        && maximum_variety(&delete(scale, t)) == 2 // Is the result of deleting t a MOS?
        && maximum_variety(&replace(scale, f1, f2)) == 2 // Is the result of identifying letters of the filling MOS a MOS?
}

/// Return the number of distinct steps in `scale`.
pub fn step_variety(scale: &[Letter]) -> usize {
    scale.iter().collect::<BTreeSet<_>>().len()
}

/// `subst()` but the filler is just one letter.
pub fn replace(scale: &[Letter], from: Letter, to: Letter) -> Vec<Letter> {
    subst(scale, from, &[to])
}

/// Delete all instances of one letter.
pub fn delete(scale: &[Letter], letter: Letter) -> Vec<Letter> {
    scale.iter().filter(|x| **x != letter).cloned().collect()
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
    use crate::helpers::gcd;

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
                        assert_eq!(maximum_variety(&mos), 2);
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
                    let mos_bjorklund = brightest_mos_mode_and_gen_bjorklund(a, b);
                    let mos_bresenham = brightest_mos_mode_and_gen_bresenham(a, b);
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
                        let mos = mos_mode(a as usize, b as usize, br as usize);
                        assert_eq!(block_balance(&mos), 1);
                    }
                }
            }
        }
    }
}
