// #![deny(warnings)]
pub mod comb;
#[macro_use]
pub mod equal;
pub mod guide;
pub mod interval;
pub mod ji;
pub mod ji_ratio;
#[macro_use]
pub mod monzo;
pub mod helpers;
pub mod primes;
pub mod vector;
#[macro_use]
pub mod subgroup_monzo;
pub mod words;

// use equal::is_in_tuning_range;
use itertools::Itertools;
use ji_ratio::RawJiRatio;
// use monzo::Monzo;
// use nalgebra::{Matrix3, Vector3};
// use subgroup_monzo::SubgroupMonzo;
use wasm_bindgen::prelude::*;
use words::{chirality, is_mos_subst};
use words::{Chirality, Letter};

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
    // Use `js_namespace` here to bind `console.log(..)` instead of just
    // `log(..)`
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

#[allow(unused_macros)]
macro_rules! console_log {
    // Note that this is using the `log` function imported above during
    // `bare_bones`
    ($($t:tt)*) => (log(&format_args!($($t)*).to_string()))
}

const STEP_LETTERS: [&str; 12] = [
    "",                                                     // 0
    "X",                                                    // 1
    "Ls",                                                   // 2
    "Lms",                                                  // 3
    "Lmns",                                                 // 4
    "HLmns",                                                // 5
    "HLmnst",                                               // 6
    "BHLmnst",                                              // 7
    "BHLmnstw",                                             // 8
    "BCHLmnstw",                                            // 9
    "BCHLmnpstw",                                           // 10
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz", // >= 11
];

use std::cmp::min;
use std::collections::HashSet;

use serde::Serialize;
use serde_wasm_bindgen::to_value;

use guide::guide_frames;
use guide::GuideFrame;
// use interval::{Dyad, JiRatio};
use words::{least_mode, maximum_variety, monotone_lm, monotone_ms, monotone_s0, CountVector};

// for the edo search
pub const EDO_BOUND: i32 = 111;
pub const S_LOWER_BOUND: f64 = 20.0;
pub const S_UPPER_BOUND: f64 = 200.0;

/// Compute the determinant of a 3x3 matrix formed by three row vectors.
/// Used to check if vectors form a unimodular basis (determinant ±1).
fn det3(row0: &[u8], row1: &[u8], row2: &[u8]) -> i16 {
    row0[0] as i16 * row1[1] as i16 * row2[2] as i16
        + row0[1] as i16 * row1[2] as i16 * row2[0] as i16
        + row0[2] as i16 * row1[0] as i16 * row2[1] as i16
        - row0[2] as i16 * row1[1] as i16 * row2[0] as i16
        - row0[1] as i16 * row1[0] as i16 * row2[2] as i16
        - row0[0] as i16 * row1[2] as i16 * row2[1] as i16
}

// A representation of a GuideFrame that should be WASM-readable
#[derive(Clone, Debug, Serialize)]
pub struct GuideResult {
    /// Either Guided GS or multiple interleaved Guided GSes
    /// `guided_gs` generates a guided generator sequence (detempered single-period MOS) subscale.
    /// The `JsValue` is an array of 3 numbers where each entry is the count of the corresp. step size.
    pub gs: Vec<Vec<u8>>,
    /// The aggregate generator
    pub aggregate: Vec<u8>,
    /// `polyoffset` is the set of intervals that each guided generator sequence chain is based on. Always includes the unison.
    /// The `JsValue` is an array of 3 numbers where each entry is the count of the corresp. step size.
    pub polyoffset: Vec<Vec<u8>>,
    /// complexity result
    /// The base GS chains in a multiple GS structure don't form interleaved scales. Instead they form a detempered copy of m-edo.
    pub multiplicity: u8,
    pub complexity: u8,
}

/// A representation of a Scale Profile. Doesn't include tunings.
#[derive(Debug, Serialize)]
pub struct ScaleProfile {
    /// brightest word
    word: String,
    /// unimodular basis for lattice is there is one
    lattice_basis: Option<Vec<Vec<u8>>>,
    /// chirality
    chirality: Chirality,
    /// brightest mode of reversed word
    reversed: String,
    /// lowest-complexity guide frame structure provided there is one
    structure: Option<GuideResult>,
    /// whether scale is L=m monotone MOS
    lm: bool,
    /// whether scale is m=s monotone MOS
    ms: bool,
    /// whether scale is s=0 monotone MOS
    s0: bool,
    /// whether scale is a subst aL(bmcs)
    subst_l_ms: bool,
    /// whether scale is a subst bm(aLcs)
    subst_m_ls: bool,
    /// whether scale is a subst cs(aLbm)
    subst_s_lm: bool,
    /// maximum variety of scale
    mv: u8,
}

#[derive(Debug, Serialize)]
pub struct SigResult {
    profiles: Vec<ScaleProfile>,
    ji_tunings: Vec<Vec<String>>,
    ed_tunings: Vec<Vec<String>>,
}

#[derive(Debug, Serialize)]
pub struct WordResult {
    profile: ScaleProfile,
    ji_tunings: Vec<Vec<String>>,
    ed_tunings: Vec<Vec<String>>,
}

fn string_to_numbers(word: &str) -> Vec<usize> {
    let mut result = vec![];
    let arity = word.chars().collect::<HashSet<_>>().len();
    for c in word.chars() {
        if let Some(letter) = STEP_LETTERS[min(arity, 12)].find(c) {
            result.push(letter);
        }
    }
    result
}

fn word_to_sig(input: &[usize]) -> Vec<usize> {
    let mut result = vec![0, 0, 0];
    for &i in input {
        if i < 3 {
            result[i] += 1;
        }
    }
    result
}

fn numbers_to_string(word: &[usize]) -> String {
    let mut result = "".to_string();
    let arity = word.iter().collect::<HashSet<_>>().len();
    for i in word {
        if *i <= arity {
            result.push(STEP_LETTERS[min(arity, 12)].chars().nth(*i).unwrap_or('?'));
        }
    }
    result
}

/// Convert a CountVector to a 3-element u8 vector for serialization
fn countvector_to_u8_vec(count_vector: &CountVector<usize>) -> Vec<u8> {
    let btreemap = count_vector.into_inner();
    vec![
        *btreemap.get(&0).unwrap_or(&0) as u8,
        *btreemap.get(&1).unwrap_or(&0) as u8,
        *btreemap.get(&2).unwrap_or(&0) as u8,
    ]
}

fn guide_frame_to_result(structure: &GuideFrame) -> GuideResult {
    let GuideFrame { gs, polyoffset } = structure;
    let aggregate_cv: CountVector<usize> = gs
        .iter()
        .fold(CountVector::<usize>::ZERO, |acc, v| acc.add(v));

    GuideResult {
        gs: gs.iter().map(countvector_to_u8_vec).collect(),
        aggregate: countvector_to_u8_vec(&aggregate_cv),
        polyoffset: polyoffset.iter().map(countvector_to_u8_vec).collect(),
        multiplicity: structure.multiplicity() as u8,
        complexity: structure.complexity() as u8,
    }
}

fn get_unimodular_basis(
    structures: &[GuideFrame],
    step_sig: &[u8],
) -> Option<(Vec<Vec<u8>>, GuideResult)> {
    /*
    if (structure["multiplicity"] === 1) {
      if (structure["polyoffset"].length === 1) {
        // Check for two unequal step vectors.
        outer: for (let i = 0; i < gs.length; ++i) {
          for (let j = i; j < gs.length; ++j) {
            if (!isEqual(gs[i], gs[j])) {
              g = [...Object.values(gs[i])];
              h = [...Object.values(gs[j])];
              break outer;
            }
          }
        }
      } else {
        g = [...Object.values(structure["aggregate"])];
        h = [...Object.values(structure["polyoffset"][1])];
      }
    } else {
      g = [...Object.values(structure["aggregate"])];
      h = [...Object.values(dyadOnDegree(
        scaleWord,
        scaleWord.length / structure["multiplicity"],
        scaleWord.length / structure["multiplicity"],
      ))];
    }
     */
    for structure in structures {
        if structure.multiplicity() == 1 {
            let result = guide_frame_to_result(structure);
            let gs = &result.gs;
            for i in 0..gs.len() {
                for j in i..gs.len() {
                    if det3(step_sig, &gs[i], &gs[j]).abs() == 1 {
                        return Some((vec![gs[i].clone(), gs[j].clone()], result));
                    }
                }
            }
            for v in &result.polyoffset {
                for w in gs {
                    if det3(step_sig, v, w).abs() == 1 {
                        return Some((vec![v.clone(), w.clone()], result));
                    }
                }
            }
        } else {
            // this branch handles multiplicity > 1 scales
            let result = guide_frame_to_result(structure);
            if let (Some(vec_for_gs_element), Some(vec_for_offset)) =
                (result.gs.first(), result.polyoffset.last())
            {
                if det3(step_sig, vec_for_gs_element, vec_for_offset).abs() == 1 {
                    return Some((
                        vec![vec_for_gs_element.clone(), vec_for_offset.clone()],
                        result,
                    ));
                }
            }
        }
    }
    None
}

/*

// Mutates `arr` a flat array matrix, swapping rows `i1` and `i2`.
function swapRows(arr, m, n, i1, i2) {
  if (i1 < 0 || i1 >= m) {
    throw new Error(
      `switchRows(): matrix index out of bounds! i1: ${i1} but m: ${m}`,
    );
  }
  if (i2 < 0 || i2 >= m) {
    throw new Error(
      `switchRows(): matrix index out of bounds! i2: ${i2} but m: ${m}`,
    );
  }
  for (let j = 0; j < n; j++) {
    const new1 = arr[n * i2 + j];
    const new2 = arr[n * i1 + j];
    [arr[n * i1 + j], arr[n * i2 + j]] = [new1, new2];
  }
}

// Scales row `i` of a matrix by `coeff`.
function multiplyRow(arr, m, n, i, coeff) {
  for (let j = 0; j < n; j++) {
    arr[n * i + j] *= coeff;
  }
}

// Does the operation "[row i2 of arr] += coeff * [row i1 of arr]".
function addMultipleOfFirstRowToSecond(arr, _, n, i1, i2, coeff) {
  for (let j = 0; j < n; j++) {
    arr[n * i2 + j] += coeff * arr[n * i1 + j];
  }
}

// Mutates `v` a flat Vec matrix in row major order, swapping rows `i1` and `i2`.
fn swap_rows(v: &mut [f64], m: usize, n: usize, i1: usize, i2: usize) {
    if i1 > m {
        console_log!("swap_rows(): matrix index i1 out of bounds");
        panic!();
    }
    if i2 > m {
        console_log!("swap_rows(): matrix index i2 out of bounds");
        panic!();
    }
    for j in 0..n {
        let new1 = v[n * i2 + j];
        let new2 = v[n * i1 + j];

        v[n * i1 + j] = new1;
        v[n * i2 + j] = new2;
        console_log!("{:?}", v);
    }
    console_log!("{:?}", v);
}

// Scales row `i` of `v` a flat Vec matrix in row major order by `coeff`.
fn multiply_row(v: &mut [f64], n: usize, i: usize, coeff: f64) {
    for j in 0..n {
        v[n * i + j] *= coeff;
    }
}

// Does the operation "[row i2 of arr] += coeff * [row i1 of arr]".
fn add_multiple_of_first_row_to_second(v: &mut [f64], n: usize, i1: usize, i2: usize, coeff: f64) {
    for j in 0..n {
        v[n * i2 + j] += coeff * v[n * i1 + j];
    }
}

// Use Gaussian elimination to solve a linear system `left` * x = `right`.
// Both `left` and `right` are assumed to have `m` rows;
// `left` has `n1` columns, and `right` has `n2`.
// The function returns the RHS after this process.
/*

// Use Gaussian elimination to solve a linear system `left` * x = `right`.
// Both `left` and `right` are assumed to have `m` rows;
// `left` has `n1` columns, and `right` has `n2`.
// The function returns the RHS after this process.
function gaussianElimination(left, right, m, n1, n2) {
  // Don't mutate the input
  let leftClone = structuredClone(left);
  let rightClone = structuredClone(right);
  // Iterating over columns, clear all entries below the pivot
  for (let j = 0; j < Math.min(m, n1); ++j) {
    let i = j + 1;
    // Try to make the (j, j) entry nonzero by swapping rows
    while (i < m && Math.abs(leftClone[n1 * j + j]) < Number.EPSILON) {
      swapRows(leftClone, m, n1, j, i);
      swapRows(rightClone, m, n2, j, i);
      i++;
    }
    // Return null to indicate a singular matrix
    if (Math.abs(leftClone[n1 * j + j]) < Number.EPSILON) {
      return null;
    }
    // Clear lower triangle{
    const pivot = leftClone[n1 * j + j];
    for (let i2 = j + 1; i2 < m; ++i2) {
      const target = leftClone[n1 * i2 + j];
      if (Math.abs(target) >= Number.EPSILON) {
        addMultipleOfFirstRowToSecond(leftClone, m, n1, j, i2, -target / pivot);
        addMultipleOfFirstRowToSecond(
          rightClone,
          m,
          n2,
          j,
          i2,
          -target / pivot,
        );
      }
    }
  }
  // Clear upper triangle
  for (let j = Math.min(m, n1) - 1; j >= 0; --j) {
    const pivot = leftClone[n1 * j + j];
    for (let i2 = 0; i2 < j; ++i2) {
      const target = leftClone[n1 * i2 + j];
      if (Math.abs(target) >= Number.EPSILON) {
        addMultipleOfFirstRowToSecond(leftClone, m, n1, j, i2, -target / pivot);
        addMultipleOfFirstRowToSecond(
          rightClone,
          m,
          n2,
          j,
          i2,
          -target / pivot,
        );
      }
    }
    // Scale rows so LHS gets 1 on diag
    // (Mutation can happen via another alias, though not via the `const` aliases we defined above.)
    multiplyRow(leftClone, m, n1, j, 1 / pivot);
    multiplyRow(rightClone, m, n2, j, 1 / pivot);
  }
  return rightClone;
}
 */
fn gaussian_elimination(
    left: &[f64],
    right: &[f64],
    m: usize,
    n1: usize,
    n2: usize,
) -> Option<Vec<f64>> {
    let mut left_clone = left.to_owned();
    let mut right_clone = right.to_owned();
    for j in 0..min(m, n1) {
        console_log!("{:?}", j);
        let mut i = j + 1;
        // Try to make the (j, j) entry nonzero by swapping rows
        console_log!("{:?}", i);
        while i < m && f64::abs(left_clone[n1 * j + j]) < f64::EPSILON {
            swap_rows(&mut left_clone, m, n1, j, i);
            swap_rows(&mut right_clone, m, n2, j, i);
            i += 1;
        }
        if i < m && f64::abs(left_clone[n1 * j + j]) < f64::EPSILON {
            console_log!("Matrix {:?} not invertible", left);
            return None; // singular matrix, no solution
        }
        // Clear lower triangle
        let pivot = left_clone[n1 * j + j];
        for i2 in (j + 1)..m {
            let target = left_clone[n1 * i2 + j];
            if f64::abs(target) >= f64::EPSILON {
                add_multiple_of_first_row_to_second(&mut left_clone, n2, j, i2, -target / pivot);
                add_multiple_of_first_row_to_second(&mut right_clone, n2, j, i2, -target / pivot);
            }
        }
    }
    // Clear upper triangle
    for j in (0..=(min(m, n1) - 1)).rev() {
        let pivot = left_clone[n1 * j + j];
        for i2 in 0..j {
            let target = left_clone[n1 * i2 + j];
            if f64::abs(target) >= f64::EPSILON {
                add_multiple_of_first_row_to_second(&mut left_clone, n1, j, i2, -target / pivot);
                add_multiple_of_first_row_to_second(&mut right_clone, n2, j, i2, -target / pivot);
            }
        }
        // Scale rows so LHS gets 1 on diag
        // (Mutation can happen via another alias, though not via the `const` aliases we defined above.)
        multiply_row(&mut left_clone, n1, j, 1.0 / pivot);
        multiply_row(&mut right_clone, n2, j, 1.0 / pivot);
    }
    Some(right_clone)
}
*/

pub fn word_to_profile(query: &[usize]) -> ScaleProfile {
    let brightest = numbers_to_string(&least_mode(query));
    let lm = monotone_lm(query);
    let ms = monotone_ms(query);
    let s0 = monotone_s0(query);
    let chirality = chirality(query);
    let reversed = least_mode(&query.iter().copied().rev().collect::<Vec<usize>>());
    let reversed = numbers_to_string(&reversed);
    let mv = maximum_variety(query) as u8;
    let step_sig = word_to_sig(query)
        .iter()
        .map(|x| *x as u8)
        .collect::<Vec<u8>>();
    let subst_l_ms = is_mos_subst(query, 0, 1, 2);
    let subst_m_ls = is_mos_subst(query, 1, 0, 2);
    let subst_s_lm = is_mos_subst(query, 2, 0, 1);
    if let Some(pair) = get_unimodular_basis(&guide_frames(query), &step_sig) {
        let (lattice_basis, structure) = pair;
        ScaleProfile {
            word: brightest,
            lattice_basis: Some(lattice_basis),
            // ploidacot: Ploidacot::try_get_ploidacot(query),
            chirality,
            reversed,
            structure: Some(structure),
            lm,
            ms,
            s0,
            subst_l_ms,
            subst_m_ls,
            subst_s_lm,
            mv,
        }
    } else {
        ScaleProfile {
            word: brightest,
            lattice_basis: None,
            // ploidacot: Ploidacot::try_get_ploidacot(query),
            chirality,
            reversed,
            structure: None,
            lm,
            ms,
            s0,
            subst_l_ms,
            subst_m_ls,
            subst_s_lm,
            mv,
        }
    }
}

#[wasm_bindgen]
pub fn word_result(query: String) -> Result<JsValue, JsValue> {
    let word_as_numbers = string_to_numbers(&query);
    let step_sig = word_to_sig(&word_as_numbers);

    Ok(to_value(&WordResult {
        profile: word_to_profile(&word_as_numbers),
        ji_tunings: sig_to_ji_tunings(&step_sig, RawJiRatio::OCTAVE),
        ed_tunings: sig_to_ed_tunings(
            &step_sig,
            RawJiRatio::OCTAVE,
            EDO_BOUND,
            S_LOWER_BOUND,
            S_UPPER_BOUND,
        ),
    })?)
}

#[wasm_bindgen]
pub fn word_to_brightest(query: String) -> String {
    let word_in_numbers = string_to_numbers(&query);
    let brightest = least_mode(&word_in_numbers);
    numbers_to_string(&brightest)
}

#[wasm_bindgen]
pub fn word_to_mv(query: String) -> u8 {
    let word_in_numbers = string_to_numbers(&query);
    maximum_variety(&word_in_numbers) as u8
}

#[allow(unused_variables)]
pub fn sig_to_ji_tunings(step_sig: &[usize], equave: RawJiRatio) -> Vec<Vec<String>> {
    /*
    let (a, b, c) = (step_sig[0] as i32, step_sig[1] as i32, step_sig[2] as i32);
    let ji_tunings =
        {
            let equave_monzo_result = Monzo::try_from_ratio(equave);
            if let Ok(equave_monzo) = equave_monzo_result {
                Monzo::EIGHTY_ONE_ODD_LIMIT
                    .iter()
                    .permutations(2)
                    .flat_map(|v| {
                        let step_sig = step_sig.into_iter().map(|x| *x as i32).collect::<Vec<_>>();
                        let (q, r) = (*v[0], *v[1]); // for each ordered pair of distinct (q, r)
                        let step_vectors = [(0..=a), (0..=b), (0..=c)]
                            .into_iter()
                            .multi_cartesian_product();
                        step_vectors.flat_map(|step_vector| {
                            let q_reduced = q.rd(equave_monzo);
                            if is_in_tuning_range(q_reduced.cents(), &step_sig, &step_vector, equave) {
                                let (x, y, z) = (step_vector[0], step_vector[1], step_vector[2]);
                                let step_vectors_2 = [(0..=a), (0..=b), (0..=c)]
                                    .into_iter()
                                    .multi_cartesian_product();
                                step_vectors_2.filter_map(|step_vector_2| {
                                    let r_reduced = r.rd(equave_monzo);
                                    if is_in_tuning_range(r_reduced.cents(), &step_sig, &step_vector_2, equave) {
                                        let (x_, y_, z_) = (step_vector_2[0], step_vector_2[1], step_vector_2[2]);
                                        let det = a*y*z_ + b*z*x_ + c*x*y_ - a*z*y_ - b*x*z_ - c*y*x_;
                                        if det == 1 || det == -1 {
                                            let one_zero_zero = [1.0, 0.0, 0.0];
                                            let zero_one_zero = [0.0, 1.0, 0.0];
                                            let zero_zero_one = [0.0, 0.0, 1.0];
                                            let m = [
                                                    a as f64, b as f64, c as f64,
                                                    x as f64, y as f64, z as f64,
                                                    x_ as f64, y_ as f64, z_ as f64,
                                                ];
                                            console_log!("Getting sol1");
                                            let sol1: Vec<f64> = gaussian_elimination(&m, &one_zero_zero, 3, 3, 1)
                                                .expect("`m` should be invertible");
                                            console_log!("Getting sol2");
                                            let sol2: Vec<f64> = gaussian_elimination(&m, &zero_one_zero, 3, 3, 1)
                                                .expect("`m` should be invertible");
                                            console_log!("Getting sol3");
                                            let sol3: Vec<f64> = gaussian_elimination(&m, &zero_zero_one, 3, 3, 1)
                                                .expect("`m` should be invertible");
                                            console_log!("Getting smonzo1");
                                            let smonzo1 = SubgroupMonzo::try_new(
                                                &[equave_monzo, q, r],
                                                &[f64::round(sol1[0]) as i32, f64::round(sol1[1]) as i32, f64::round(sol1[2]) as i32]);
                                            console_log!("Getting smonzo2");
                                            let smonzo2 = SubgroupMonzo::try_new(
                                                &[equave_monzo, q, r],
                                                &[f64::round(sol2[0]) as i32, f64::round(sol2[1]) as i32, f64::round(sol2[2]) as i32]);
                                            console_log!("Getting smonzo3");
                                            let smonzo3 = SubgroupMonzo::try_new(
                                                &[equave_monzo, q, r],
                                                &[f64::round(sol3[0]) as i32, f64::round(sol3[1]) as i32, f64::round(sol3[2]) as i32]);
                                            console_log!("Returning Some(...)");
                                            if let (Ok(sm1), Ok(sm2), Ok(sm3)) = (smonzo1, smonzo2, smonzo3) {
                                                Some(vec![sm1.to_monzo(), sm2.to_monzo(), sm3.to_monzo()])
                                            } else {
                                                None
                                            }
                                        } else {
                                            None
                                        }
                                        None
                                    } else {
                                        None
                                    }
                                })
                                .collect::<Vec::<Vec::<Monzo>>>()
                            } else {
                                vec![]
                            }
                        })
                        .collect::<Vec::<Vec::<Monzo>>>()
                    })
                    .collect::<Vec::<Vec::<Monzo>>>() // Each Vec<Monzo> is one solution
            } else {
                vec![]
            }
        };
    console_log!("Got all tunings in monzos");
    ji_tunings
        .into_iter()
        .map(|v| {
            v.into_iter()
                .map(|monzo| {
                    monzo.to_string()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
                                            */
    vec![] // Return nothing until I figure out how to make this more efficient
}

pub fn sig_to_ed_tunings(
    step_sig: &[usize],
    equave: RawJiRatio,
    ed_bound: i32,
    s_lower: f64,
    s_upper: f64,
) -> Vec<Vec<String>> {
    let ed_tunings =
        crate::equal::ed_tunings_for_ternary(step_sig, equave, ed_bound, s_lower, s_upper);
    let is_octave = equave.numer() == 2 && equave.denom() == 1;
    ed_tunings
        .into_iter()
        .map(|v| {
            let ed: i32 = v
                .iter()
                .enumerate()
                .map(|(i, steps)| step_sig[i] as i32 * steps)
                .sum();
            if is_octave {
                v.iter().map(|i| format!("{i}\\{ed}")).collect::<Vec<_>>()
            } else {
                v.iter()
                    .map(|i| format!("{i}\\{ed}<{}/{}>", equave.numer(), equave.denom()))
                    .collect::<Vec<_>>()
            }
        })
        .collect::<Vec<_>>()
}

#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn sig_result(
    query: Vec<u8>,
    lm: bool,
    ms: bool,
    s0: bool,
    ggs_len: u8,
    ggs_len_constraint: String,
    complexity: u8,
    complexity_constraint: String,
    mv: u8,
    mv_constraint: String,
    mos_subst: String,
    equave_num: u32,
    equave_den: u32,
    ed_bound: i32,
    s_lower: f64,
    s_upper: f64,
) -> Result<JsValue, JsValue> {
    let equave =
        RawJiRatio::try_new(equave_num as u64, equave_den as u64).unwrap_or(RawJiRatio::OCTAVE);
    let step_sig = query;
    let filtering_cond = |scale: &[Letter]| {
        (!lm || monotone_lm(scale))
            && (!ms || monotone_ms(scale))
            && (!s0 || monotone_s0(scale))
            && (match ggs_len {
                0 => true,
                l => {
                    let guide_frames = guide_frames(scale);
                    if ggs_len_constraint == "exactly" {
                        !guide_frames.is_empty() && guide_frames[0].gs.len() == l as usize
                    } else {
                        !guide_frames.is_empty() && guide_frames[0].gs.len() <= l as usize
                    }
                }
            })
            && (match mv {
                0 => true,
                mv => {
                    if mv_constraint == "exactly" {
                        maximum_variety(scale) == mv as usize
                    } else {
                        maximum_variety(scale) <= mv as usize
                    }
                }
            })
            && (match complexity {
                0 => true,
                c => {
                    let guide_frames = guide_frames(scale);
                    if complexity_constraint == "exactly" {
                        !guide_frames.is_empty() && guide_frames[0].complexity() == c as usize
                    } else {
                        !guide_frames.is_empty() && guide_frames[0].complexity() <= c as usize
                    }
                }
            })
    };
    let step_sig = step_sig.iter().map(|x| *x as usize).collect::<Vec<_>>();
    let scales = if mos_subst == "on" {
        words::mos_substitution_scales(&step_sig)
    } else {
        crate::comb::necklaces_fixed_content(&step_sig)
    }; // Now filter
    let scales = scales
        .into_iter()
        .filter(|scale| filtering_cond(scale))
        .collect::<Vec<_>>();
    Ok(to_value(&SigResult {
        profiles: scales
            .iter()
            .map(|scale| word_to_profile(scale))
            .sorted_by_key(|profile| {
                if let Some(guide) = &profile.structure {
                    guide.complexity
                } else {
                    u8::MAX
                }
            })
            .collect(),
        ji_tunings: sig_to_ji_tunings(&step_sig, equave),
        ed_tunings: sig_to_ed_tunings(&step_sig, equave, ed_bound, s_lower, s_upper),
    })?)
}
