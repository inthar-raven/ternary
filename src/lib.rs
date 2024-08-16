#![deny(warnings)]
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
pub mod words;

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}

const STEP_LETTERS: [&str; 12] = [
    "",                                                     // 0
    "X",                                                    // 1
    "Ls",                                                   // 2
    "LMs",                                                  // 3
    "LMns",                                                 // 4
    "HLMns",                                                // 5
    "HLMnst",                                               // 6
    "BHLMnst",                                              // 7
    "BHLMnstw",                                             // 8
    "BCHLMnstw",                                            // 9
    "BCHLMnpstw",                                           // 10
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz", // >= 11
];

use std::cmp::min;
use std::collections::HashSet;

use itertools::Itertools;
use serde::Serialize;
use serde_wasm_bindgen::to_value;

use guide::guide_structures;
use guide::GuideFrame;
use interval::JiRatio;
use words::{
    least_mode, maximum_variety, monotone_lm, monotone_ms, monotone_s0, CountVector, Letter,
};

// for the edo search
pub const EDO_BOUND: i32 = 53;
pub const S_LOWER_BOUND: f64 = 20.0;
pub const S_UPPER_BOUND: f64 = 200.0;

// A representation of a GuideFrame that should be WASM-readable
#[derive(Debug, Serialize)]
pub struct GuideResult {
    /// Either WFGS or multiple interleaved WFGSes
    /// `wfgs` generates a well-formed generator sequence (detempered single-period MOS) subscale.
    /// The `JsValue` is an array of 3 numbers where each entry is the count of the corresp. step size.
    pub gs: Vec<Vec<u8>>,
    /// The aggregate generator
    pub aggregate: Vec<u8>,
    /// `polyoffset` is the set of intervals that each well-formed generator sequence chain is based on. Always includes the unison.
    /// The `JsValue` is an array of 3 numbers where each entry is the count of the corresp. step size.
    pub polyoffset: Vec<Vec<u8>>,
    /// complexity result
    /// The base GS chains in a multiple GS structure don't form interleaved scales. Instead they form a detempered copy of m-edo.
    pub multiplicity: u8,
    pub complexity: u8,
}

// A representation of a Scale Profile. Doesn't include tunings.
#[derive(Debug, Serialize)]
pub struct ScaleProfile {
    word: String,
    structure: GuideResult,
    lm: bool,
    ms: bool,
    s0: bool,
    mv: u8,
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

fn numbers_to_string(word: &[usize]) -> String {
    let mut result = "".to_string();
    let arity = word.iter().collect::<HashSet<_>>().len();
    for i in word {
        if *i <= arity {
            result.push(STEP_LETTERS[min(arity, 12)].chars().nth(*i).unwrap());
        }
    }
    result
}

#[wasm_bindgen]
pub fn word_to_profile(query: String) -> Result<JsValue, JsValue> {
    let word_in_numbers = string_to_numbers(&query);
    let brightest = numbers_to_string(&least_mode(&word_in_numbers));
    if let Some(structure) = guide_structures(&string_to_numbers(&query)).first() {
        let GuideFrame {
            gs,
            polyoffset,
            multiplicity,
        } = structure;
        let structure = if *multiplicity == 1 {
            GuideResult {
                gs: gs
                    .iter()
                    .map(|cv| {
                        let btreemap = cv.into_inner();
                        vec![
                            *btreemap.get(&0).unwrap_or(&0) as u8,
                            *btreemap.get(&1).unwrap_or(&0) as u8,
                            *btreemap.get(&2).unwrap_or(&0) as u8,
                        ]
                    })
                    .collect(),
                aggregate: {
                    let cv: CountVector<usize> = gs
                        .iter()
                        .fold(CountVector::<usize>::ZERO, |acc, v| acc.add(v));
                    let btreemap = cv.into_inner();
                    vec![
                        *btreemap.get(&0).unwrap_or(&0) as u8,
                        *btreemap.get(&1).unwrap_or(&0) as u8,
                        *btreemap.get(&2).unwrap_or(&0) as u8,
                    ]
                },
                polyoffset: polyoffset
                    .iter()
                    .map(|cv| {
                        let btreemap = cv.into_inner();
                        vec![
                            *btreemap.get(&0).unwrap_or(&0) as u8,
                            *btreemap.get(&1).unwrap_or(&0) as u8,
                            *btreemap.get(&2).unwrap_or(&0) as u8,
                        ]
                    })
                    .collect(),
                multiplicity: 1,
                complexity: structure.complexity() as u8,
            }
        } else {
            GuideResult {
                gs: gs
                    .iter()
                    .map(|cv| {
                        let btreemap = cv.into_inner();
                        vec![
                            *btreemap.get(&0).unwrap_or(&0) as u8,
                            *btreemap.get(&1).unwrap_or(&0) as u8,
                            *btreemap.get(&2).unwrap_or(&0) as u8,
                        ]
                    })
                    .collect(),
                aggregate: {
                    let cv: CountVector<usize> = gs
                        .iter()
                        .fold(CountVector::<usize>::ZERO, |acc, v| acc.add(v));
                    let btreemap = cv.into_inner();
                    vec![
                        *btreemap.get(&0).unwrap_or(&0) as u8,
                        *btreemap.get(&1).unwrap_or(&0) as u8,
                        *btreemap.get(&2).unwrap_or(&0) as u8,
                    ]
                },
                polyoffset: vec![vec![0, 0, 0]],
                multiplicity: *multiplicity as u8,
                complexity: structure.complexity() as u8,
            }
        };
        let lm = monotone_lm(&word_in_numbers);
        let ms = monotone_ms(&word_in_numbers);
        let s0 = monotone_s0(&word_in_numbers);
        let mv = maximum_variety(&word_in_numbers) as u8;
        Ok(to_value(&ScaleProfile {
            word: brightest,
            structure,
            lm,
            ms,
            s0,
            mv,
        })?)
    } else {
        Err(JsValue::NULL)
    }
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

#[wasm_bindgen]
pub fn word_to_ji_tunings(query: String) -> Vec<JsValue> {
    let word_in_numbers = string_to_numbers(&query);
    let step_sig = {
        let mut r = vec![0; 3];
        for letter in word_in_numbers {
            if letter < 3 {
                r[letter] += 1;
            }
        }
        r
    };
    let ji_tunings =
        crate::ji::solve_step_sig_81_odd_limit(&step_sig, crate::monzo::Monzo::OCTAVE, false);
    ji_tunings
        .into_iter()
        .map(|v| {
            v.into_iter()
                .map(|monzo| format!("{}/{}", monzo.numer(), monzo.denom()))
                .collect::<Vec<_>>()
                .into()
        })
        .collect::<Vec<_>>()
}

#[wasm_bindgen]
pub fn word_to_ed_tunings(query: String) -> Vec<JsValue> {
    let word_in_numbers = string_to_numbers(&query);
    let step_sig = {
        let mut r = vec![0; 3];
        for letter in word_in_numbers {
            if letter < 3 {
                r[letter] += 1;
            }
        }
        r
    };
    let ed_tunings = crate::equal::ed_tunings_for_ternary(
        &step_sig,
        crate::ji_ratio::RawJiRatio::OCTAVE,
        EDO_BOUND,
        S_LOWER_BOUND,
        S_UPPER_BOUND,
    );
    ed_tunings
        .into_iter()
        .map(|v| {
            let edo: i32 = v
                .iter()
                .enumerate()
                .map(|(i, steps)| step_sig[i] as i32 * steps)
                .sum();
            v.iter()
                .map(|i| format!("{}\\{}", i, edo))
                .collect::<Vec<_>>()
                .into()
        })
        .collect::<Vec<_>>()
}

#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn sig_to_profiles(
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
) -> Vec<JsValue> {
    let step_sig = query;
    let filtering_cond = |scale: &[Letter]| {
        (!lm || monotone_lm(scale))
            && (!ms || monotone_ms(scale))
            && (!s0 || monotone_s0(scale))
            && (match ggs_len {
                0 => true,
                l => {
                    let guide_frames = guide_structures(scale);
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
                    let guide_frames = guide_structures(scale);
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
        .sorted_unstable_by_key(|scale| {
            if let Some(first) = guide_structures(scale).first() {
                first.complexity()
            } else {
                usize::MAX
            }
        })
        .collect::<Vec<_>>();
    let mut result = vec![];
    for scale in scales {
        result.push(word_to_profile(numbers_to_string(&scale)).unwrap_or(JsValue::NULL));
    }
    result
}
