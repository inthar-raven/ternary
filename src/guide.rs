use std::collections::BTreeSet;

use itertools::Itertools;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
    // Use `js_namespace` here to bind `console.log(..)` instead of just
    // `log(..)`
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

use crate::{
    helpers::gcd,
    primes::factorize,
    words::{
        CountVector, Letter, Subtendable, dyad_on_degree, offset_vec, rotate, rotations,
        weak_period_pattern, word_on_degree,
    },
};

// Given a necklace of stacked k-steps, where k is fixed,
// get all valid Guided GSes using k-steps on any rotation.
fn guided_gs_chains<T>(chain: &[T]) -> Vec<Vec<T>>
where
    T: core::fmt::Debug + PartialEq + Clone + Eq + Send,
{
    // println!("chain: {:?}", chain);
    let len = chain.len();
    rotations(chain)
        .into_iter()
        .filter(|list| {
            // println!("list: {:?}, {}", list, !(list[..len - 1].contains(&list[len - 1])));
            !(list[..len - 1].contains(&list[len - 1]))
        })
        .map(|chain| weak_period_pattern(&chain[..len - 1]))
        .collect::<Vec<_>>()
}
/*
    A guided generator sequence (Guided GS) is a generator sequence made of stacked `k`-steps,
    where `k` is fixed and `gcd(k, scale.len()) == 1`.
    The interval left over after stacking (which is a `k`-step) is different from all the others,
    and where the generators in the generator sequence are distinct from any non-`k`-step interval.
    There is no need to check the last condition for step vectors.
*/
fn step_class_guided_gs_list(step_class: usize, scale: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    guided_gs_chains(&stacked_step_class(step_class, scale))
}

pub fn stacked_step_class<T>(step_class: usize, scale: &[T]) -> Vec<T::Interval>
where
    T: Subtendable + std::fmt::Debug,
{
    // println!("scale: {:?}", scale);
    (0..scale.len())
        .map(|i| word_on_degree(scale, step_class * i, step_class))
        .map(|subword| <T as Subtendable>::interval_from_slice(&subword))
        .collect()
}

/// TODO: To get all Guided GSes, we actually need to iterate through all the step classes, not just ones <= len / 2.
/// All Guided GSes that generate a given abstract necklace.
/// Guided GS: generator sequence using fixed k-steps where gcd(k, scale.len()) == 1.
pub fn guided_gs_list(scale: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    let len = scale.len();
    (2..=len - 2) // Don't include 1-step GSes
        .filter(|&step_class| gcd(step_class as u32, len as u32) == 1)
        .flat_map(|step_class| step_class_guided_gs_list(step_class, scale))
        .collect()
}

/// All Guided GSes of length `gs_length` that generate a given abstract necklace.
pub fn guided_gs_list_of_len(gs_length: usize, scale: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    let scale_len = scale.len();
    (2..=scale_len / 2) // Don't include 1-step GSes
        .filter(|&step_class| gcd(step_class as u32, scale_len as u32) == 1)
        .flat_map(|step_class| step_class_guided_gs_list(step_class, scale))
        .filter(|vs| gs_length == vs.len())
        .collect()
}

#[allow(unused)]
fn step_class_guided_gs_list_for_subscale(
    step_class: usize,
    subscale: &[CountVector<usize>],
) -> Vec<Vec<CountVector<usize>>> {
    guided_gs_chains(&stacked_step_class(step_class, subscale))
}

/// All Guided GSes of a chain which is represented as `Vec<CountVector<usize>>` rather than `Vec<usize>`.
fn guided_gs_list_for_subscale(subscale: &[CountVector<usize>]) -> Vec<Vec<CountVector<usize>>> {
    if subscale.len() == 2 {
        vec![vec![subscale[0].clone()]]
    } else {
        let len = subscale.len();
        (2..=len / 2) // Don't include 1-step GSes
            .filter(|&step_class| gcd(step_class as u32, len as u32) == 1)
            .flat_map(|step_class| step_class_guided_gs_list_for_subscale(step_class, subscale))
            .collect()
    }
}

/// A guide frame structure for a scale word, consisting of a generator sequence together with a set of offsets or a multiplicity.
/// Multiplicity > 1 is a generalization of [even-regular MV3 scales](https://en.xen.wiki/w/Even-regular_MV3_scale);
/// scales of this type always have number of notes divisible by the multiplicity.
/// Consists of m copies of same generator sequence offset by `scale.len()` / m steps (where m is the multiplicity).
#[derive(Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct GuideFrame {
    /// Either Guided GS or multiple interleaved GSes (Guided GSes when considered individually).
    /// `gs` generates a guided generator sequence (detempered single-period_pattern MOS) subscale.
    pub gs: Vec<CountVector<usize>>,
    /// `offset_chord` is the set of intervals that each guided generator sequence chain is based on. Always includes `CountVector::ZERO`.
    pub offset_chord: Vec<CountVector<usize>>,
}

impl GuideFrame {
    /// Create a simple guide frame (single generator sequence).
    pub fn new_simple(gs: Vec<CountVector<usize>>) -> Self {
        Self {
            gs,
            offset_chord: vec![CountVector::ZERO],
        }
    }
    /// Create a multiple/interleaved guide frame (with offset_chord).
    pub fn new_multiple(
        gs: Vec<CountVector<usize>>,
        offset_chord: Vec<CountVector<usize>>,
    ) -> Self {
        Self { gs, offset_chord }
    }
    /// The complexity of a guide frame (size of generator sequence plus size of offset_chord).
    pub fn complexity(&self) -> usize {
        self.gs.len() * self.offset_chord.len()
    }
    /// The multiplicity (number of times generator sequence is repeated with offsets).
    pub fn multiplicity(&self) -> usize {
        self.offset_chord.len()
    }
    /// Try to get simple guide frames with k-step generators.
    pub fn try_simple(scale: &[usize], step_class: usize) -> Vec<Self> {
        if scale.is_empty() || gcd(scale.len() as u32, step_class as u32) != 1 {
            vec![]
        } else {
            step_class_guided_gs_list(step_class, scale)
                .into_iter()
                .map(|gs| Self {
                    gs,
                    offset_chord: vec![CountVector::ZERO],
                })
                .sorted()
                .dedup()
                .collect::<Vec<_>>()
        }
    }
    pub fn try_multiple(scale: &[usize], multiplicity: usize, step_class: usize) -> Vec<Self> {
        // The scale cannot be empty and its size must be divisible by `multiplicity`.
        if multiplicity == 1 || scale.is_empty() || !scale.len().is_multiple_of(multiplicity) {
            vec![]
        } else {
            let gcd_value = gcd(step_class as u32, scale.len() as u32) as usize;
            let coprime_part = scale.len() / gcd_value;
            if !coprime_part.is_multiple_of(multiplicity) {
                if gcd_value == multiplicity {
                    // println!("interleaved");
                    // It's an interleaved scale.
                    let subscales = (0..gcd_value)
                        .map(|degree| rotate(scale, degree))
                        .map(|rotation| {
                            stacked_step_class(gcd_value, &rotation[..scale.len()])
                                [..scale.len() / gcd_value]
                                .to_vec()
                        })
                        .collect::<Vec<_>>();
                    // since we checked that `scale` is nonempty, the following operation should be infallible
                    let subscale_on_root = subscales[0].clone();
                    // println!("subscale_on_root: {:?}", subscale_on_root);

                    // All of the subscales must be rotations of one another.
                    // `offset_vec()` returns a witness to rotational equivalence (an offset) if there is any;
                    // the offsets are combined to form the offset_chord.
                    // If it returns `None` for any subscale, the whole procedure fails.
                    let maybe_offsets = subscales
                        .into_iter()
                        .enumerate()
                        .map(|(i, subscale)| {
                            offset_vec(&subscale_on_root, &subscale).map(|offset| {
                                // `.map()` returns `None` if the previous result is `None`
                                // and functorially applies the closure to `Some`s.
                                CountVector::from_slice(&word_on_degree(
                                    scale,
                                    0,
                                    offset * gcd_value + i,
                                ))
                            })
                        })
                        // `.collect()` returns `None` if there is any `None` returned by `map`.
                        .collect::<Option<Vec<CountVector<usize>>>>();
                    // println!("subscale_on_root: {:?}", subscale_on_root);
                    // println!("maybe_offsets: {:?}", maybe_offsets);
                    if let Some(offsets) = maybe_offsets {
                        // sort list of offsets by step class
                        // If offset_chord is {0} use multiplicity 1
                        if offsets.len() == 1 {
                            guided_gs_list(scale)
                                .into_iter()
                                .map(|gs| Self {
                                    gs,
                                    offset_chord: offsets.to_owned(),
                                })
                                .sorted()
                                .dedup()
                                .collect::<Vec<_>>()
                        } else {
                            let offsets_sorted: Vec<CountVector<usize>> =
                                offsets.into_iter().sorted_by_key(|v| v.len()).collect();
                            guided_gs_list_for_subscale(&subscale_on_root)
                                .into_iter()
                                .map(|gs| Self {
                                    gs,
                                    offset_chord: offsets_sorted.to_owned(),
                                })
                                .sorted()
                                .dedup()
                                .collect::<Vec<_>>()
                        }
                    } else {
                        vec![]
                    }
                } else {
                    vec![]
                }
            } else {
                // println!("not interleaved");
                // stack at most this many k-steps
                let chain_length: usize = scale.len() / multiplicity;
                if chain_length == 1 {
                    vec![]
                } else {
                    // For every degree of `scale`, get stack of `gs_length_limit`-many `step_class`-steps on that degree.
                    // We'll need the enumeration index later to get the offset_chord components from each index.
                    let gener_chains_enumerated: Vec<(usize, Vec<CountVector<usize>>)> =
                        rotations(scale)
                            .into_iter()
                            .map(|mode| {
                                // get stacked `step_class`-steps on this rotation
                                stacked_step_class(step_class, &mode)[0..chain_length].to_vec()
                            })
                            .enumerate() // get rotation index for each mode
                            .filter(|(_, stack)| {
                                // Each stack is generated by a GS,
                                // but for the GS to be a guided GS, the last element must differ from all previous elements.
                                !stack[0..chain_length - 1].contains(&stack[chain_length - 1])
                            })
                            .collect();
                    let gses: Vec<Vec<CountVector<Letter>>> = gener_chains_enumerated
                        .iter()
                        // Take prefix of `gs_length_limit - 1` elements and get the GS it is generated by
                        .map(|(_, chain)| weak_period_pattern(&chain[0..chain_length - 1]))
                        .sorted()
                        .dedup()
                        .collect();
                    // println!("{:?}", gses);
                    gses.iter()
                        .map(|gs| {
                            (
                                gs,
                                gener_chains_enumerated
                                    .iter()
                                    // Check only the prefix of gs_length_limit - 1 elements,
                                    // because that's what the guided GS is based on.
                                    // The filtering keeps all degrees on which this particular GS (`gs`) occurs.
                                    .filter(|(_, gener_chain)| {
                                        weak_period_pattern(&gener_chain[..chain_length - 1]) == *gs
                                    })
                                    .map(|(i, _)| *i)
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .filter(|(_, polyoffset_indices)| {
                            // Filter by whether the number of gen chains each GS occurs on is equal to the multiplcity.
                            // We also need to check that all of the chains are disjoint.
                            let mut union_of_chains: Vec<_> = polyoffset_indices
                                .iter()
                                .flat_map(|first| {
                                    (0..chain_length)
                                        .map(|i| (first + i * step_class) % scale.len())
                                        .collect::<Vec<_>>()
                                })
                                .collect();
                            union_of_chains.sort();
                            union_of_chains.dedup();
                            // println!("union_of_chains: {:?}", union_of_chains);
                            let chains_are_disjoint = union_of_chains.len() == scale.len();
                            chains_are_disjoint && polyoffset_indices.len() == multiplicity
                        })
                        .map(|(gs, polyoffset_indices)| {
                            // println!("gs: {:?}", gs);
                            // println!("polyoffset_indices: {:?}", polyoffset_indices);
                            let first_deg = polyoffset_indices[0];
                            let offset_chord: Vec<CountVector<Letter>> = polyoffset_indices
                                .iter()
                                .map(|degree| {
                                    dyad_on_degree(
                                        &rotate(scale, first_deg),
                                        first_deg,
                                        degree - first_deg,
                                    )
                                })
                                .collect();
                            Self {
                                gs: gs.to_owned(),
                                offset_chord,
                            }
                        })
                        .collect()
                }
            }
        }
    }
    /// Get both Simple and Multiple guide frames.
    fn try_all_variants(scale: &[usize], step_class: usize) -> Vec<Self> {
        // Let's just do primes for now.
        let prime_factors: Vec<usize> = factorize(scale.len() as u32)
            .into_iter()
            .dedup()
            .map(|prime| prime as usize)
            .collect();
        let simple_guide_moses: Vec<GuideFrame> = Self::try_simple(scale, step_class);
        let multiple_guide_moses: Vec<GuideFrame> = if BTreeSet::from_iter(scale.iter()).len() > 1 {
            prime_factors
                .into_iter()
                .flat_map(|prime| Self::try_multiple(scale, prime, step_class))
                .collect()
        } else {
            vec![]
        };

        let mut guide_frames = [simple_guide_moses, multiple_guide_moses].concat();
        guide_frames.sort_by_key(GuideFrame::complexity);
        // println!("{:?}", guide_frames);
        guide_frames
    }
}

/// Return the collection of guide frames for the given scale word, sorted by complexity.
pub fn guide_frames(scale: &[usize]) -> Vec<GuideFrame> {
    (2..=scale.len() / 2) // steps subtended by generator used for the guided generator sequence
        .flat_map(|step_class| GuideFrame::try_all_variants(scale, step_class))
        .sorted_by_key(GuideFrame::complexity)
        .collect()
}
#[cfg(test)]
mod tests {
    #[allow(unused)]
    use std::collections::{BTreeMap, BTreeSet, HashSet};

    use crate::words::{CountVector, Letter};

    use super::*;
    /*
    #[test]
    fn test_ploidacot() {
        let pinedye = [0, 0, 1, 0, 1, 0, 0, 2];
        assert_eq!(
            // Pinedye is monocot
            Ploidacot::try_get_ploidacot(&pinedye),
            Some(Ploidacot {
                ploid: 1,
                cot: 1,
                shear: 0,
            })
        );
        let diasem = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        assert_eq!(
            // Diasem is 1-sheared dicot
            Ploidacot::try_get_ploidacot(&diasem),
            Some(Ploidacot {
                ploid: 1,
                cot: 2,
                shear: 1,
            })
        );
        let diamech_4sl: [usize; 11] = [1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2];
        assert_eq!(
            // diaslen is tricot
            Ploidacot::try_get_ploidacot(&diamech_4sl),
            Some(Ploidacot {
                ploid: 1,
                cot: 3,
                shear: 0,
            })
        );
        let diachrome_5sc = [0, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 1];
        assert_eq!(
            // Central diachrome is diploid monocot
            Ploidacot::try_get_ploidacot(&diachrome_5sc),
            Some(Ploidacot {
                ploid: 2,
                cot: 1,
                shear: 0,
            })
        );
        let blackdye: [usize; 10] = [0, 1, 0, 2, 0, 1, 0, 2, 0, 2];
        assert_eq!(
            // Blackdye is haploid monocot (ignore the contorsion caused by the offset)
            Ploidacot::try_get_ploidacot(&blackdye),
            Some(Ploidacot {
                ploid: 1,
                cot: 1,
                shear: 0,
            })
        );
    }
    */
    #[test]
    fn test_lllmllms() {
        let bad_scale: [usize; 8] = [0, 0, 0, 1, 0, 0, 1, 2];
        let complexity_2_gses = GuideFrame::try_multiple(&bad_scale, 2, 2);
        println!("{complexity_2_gses:?}");
        assert!(complexity_2_gses.is_empty());
    }

    #[test]
    fn test_blackdye() {
        let blackdye: [usize; 10] = [0, 1, 0, 2, 0, 1, 0, 2, 0, 2];
        let should_have_mult_2 = GuideFrame::try_multiple(&blackdye, 2, 4);
        assert!(!should_have_mult_2.is_empty());
    }
    #[test]
    fn test_fix_bug_for_4sr() {
        let diamech_4sr: [Letter; 11] = [0, 2, 0, 1, 0, 2, 0, 2, 0, 1, 2];

        assert!(guided_gs_list_of_len(3, &diamech_4sr).contains(&vec![
            CountVector::from_slice(&[0, 2]),
            CountVector::from_slice(&[0, 1]),
            CountVector::from_slice(&[0, 2]),
        ]));
        assert!(guided_gs_list(&diamech_4sr).contains(&vec![
            CountVector::from_slice(&[0, 2]),
            CountVector::from_slice(&[0, 1]),
            CountVector::from_slice(&[0, 2]),
        ]));
        let guide_frames = GuideFrame::try_simple(&diamech_4sr, 2);
        // println!("{:?}", guide_frames);
        assert!(guide_frames.contains(&GuideFrame::new_simple(vec![
            CountVector::from_slice(&[0, 2]),
            CountVector::from_slice(&[0, 1]),
            CountVector::from_slice(&[0, 2]),
        ])));
    }
    #[test]
    fn test_guided_gs_based_guide_frame() {
        let pinedye = [0, 0, 1, 0, 1, 0, 0, 2];
        let pinedye_guide_moses = guide_frames(&pinedye);
        assert!(pinedye_guide_moses.contains(&GuideFrame::new_simple(vec![
            CountVector::from_slice(&[0, 0, 2]),
            CountVector::from_slice(&[0, 0, 1]),
            CountVector::from_slice(&[0, 0, 1]),
        ])));

        let diasem = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let diasem_guide_moses = guide_frames(&diasem);
        assert!(diasem_guide_moses.contains(&GuideFrame::new_simple(vec![
            CountVector::from_slice(&[0, 1]),
            CountVector::from_slice(&[0, 2])
        ])));
        assert_eq!(
            GuideFrame::new_simple(vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2])
            ],)
            .complexity(),
            2
        );
        let blackdye: [usize; 10] = [0, 1, 0, 2, 0, 1, 0, 2, 0, 2];
        let blackdye_guide_moses = guide_frames(&blackdye);
        assert!(blackdye_guide_moses.contains(&GuideFrame::new_multiple(
            vec![CountVector::from_slice(&[0, 0, 1, 2]),],
            vec![CountVector::ZERO, CountVector::from_slice(&[0]),]
        )));

        let diamech_4sl: [usize; 11] = [1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2];
        let diamech_guide_moses = guide_frames(&diamech_4sl);
        assert!(diamech_guide_moses.contains(&GuideFrame::new_simple(vec![
            CountVector::from_slice(&[0, 2]),
            CountVector::from_slice(&[0, 2]),
            CountVector::from_slice(&[0, 1]),
        ],)));
        assert_eq!(
            GuideFrame::new_simple(vec![
                CountVector::from_slice(&[0, 2]),
                CountVector::from_slice(&[0, 2]),
                CountVector::from_slice(&[0, 1]),
            ])
            .complexity(),
            3
        );

        let diachrome_5sc = [0, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 1];
        let diachrome_guide_moses = guide_frames(&diachrome_5sc);
        assert!(
            diachrome_guide_moses.contains(&GuideFrame::new_multiple(
                vec![CountVector::from_slice(&[0, 0, 1, 2, 2]),],
                vec![
                    CountVector::ZERO,
                    CountVector::from_slice(&[0, 0, 0, 1, 2, 2]),
                ],
            )) || diachrome_guide_moses.contains(&GuideFrame::new_multiple(
                vec![CountVector::from_slice(&[0, 0, 1, 2, 2]),],
                vec![
                    CountVector::ZERO,
                    CountVector::from_slice(&[0, 0, 1, 2, 2, 2]),
                ],
            ))
        );
        assert_eq!(
            GuideFrame::new_multiple(
                vec![CountVector::from_slice(&[0, 0, 1, 2, 2])],
                vec![
                    CountVector::ZERO,
                    CountVector::from_slice(&[0, 0, 1, 2, 2, 2]),
                ],
            )
            .complexity(),
            2
        );

        assert_eq!(
            GuideFrame::new_multiple(
                vec![CountVector::from_slice(&[0, 0, 1, 2, 2])],
                vec![
                    CountVector::ZERO,
                    CountVector::from_slice(&[0, 0, 1, 2, 2, 2]),
                ],
            )
            .multiplicity(),
            2
        );
    }

    #[test]
    fn test_stacked_step_class() {
        let diasem: [usize; 9] = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let one_steps = stacked_step_class(1, &diasem);
        assert_eq! {
            one_steps,
            vec![
                CountVector::from_slice(&[0]),
                CountVector::from_slice(&[1]),
                CountVector::from_slice(&[0]),
                CountVector::from_slice(&[2]),
                CountVector::from_slice(&[0]),
                CountVector::from_slice(&[1]),
                CountVector::from_slice(&[0]),
                CountVector::from_slice(&[2]),
                CountVector::from_slice(&[0]),
            ]
        }
        let two_steps = stacked_step_class(2, &diasem);
        assert_eq! {
            two_steps,
            vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2]),
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2]),
                CountVector::from_slice(&[0, 0]),
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2]),
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2]),
            ]
        }
        let three_steps = stacked_step_class(3, &diasem);
        assert_eq! {
            three_steps,
            vec![
                CountVector::from_slice(&[0, 1, 0]),
                CountVector::from_slice(&[2, 0, 1]),
                CountVector::from_slice(&[0, 2, 0]),
                CountVector::from_slice(&[0, 1, 0]),
                CountVector::from_slice(&[2, 0, 1]),
                CountVector::from_slice(&[0, 2, 0]),
                CountVector::from_slice(&[0, 1, 0]),
                CountVector::from_slice(&[2, 0, 1]),
                CountVector::from_slice(&[0, 2, 0]),
            ]
        }
    }
    #[test]
    fn test_guided_gs_chains() {
        let diasem: [usize; 9] = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let two_steps = stacked_step_class(2, &diasem);
        let chains = guided_gs_chains(two_steps.as_slice());
        assert_eq!(
            chains,
            vec![vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2]),
            ]]
        );
        let four_steps = stacked_step_class(4, &diasem);
        let chains = guided_gs_chains(four_steps.as_slice());
        assert_eq!(
            BTreeSet::from_iter(chains.into_iter()),
            BTreeSet::from_iter(
                vec![
                    vec![
                        CountVector::from_slice(&[1, 0, 2, 0]),
                        CountVector::from_slice(&[1, 0, 2, 0]),
                        CountVector::from_slice(&[0, 1, 0, 2]),
                        CountVector::from_slice(&[0, 1, 0, 2]),
                        CountVector::from_slice(&[0, 0, 1, 0]),
                    ],
                    vec![
                        CountVector::from_slice(&[2, 0, 1, 0]),
                        CountVector::from_slice(&[2, 0, 0, 1]),
                        CountVector::from_slice(&[0, 2, 0, 1]),
                        CountVector::from_slice(&[0, 2, 0, 0]),
                        CountVector::from_slice(&[1, 0, 2, 0]),
                    ]
                ]
                .into_iter()
            )
        );
        /*
        Valid Guided GS necklaces for 010201020:
        10 20 10 20 01 02 01 02 (closing 00) -> abababab -> GS(a, b)
        1020 1020 0102 0102 0010 2010 2001 0201 (closing 0200) -> aaaabaaac -> GS(a, a, a, a, b)
        2010 2001 0201 0200 1020 1020 0102 0102 (closing 0010) -> aaabaaaac -> GS(a, a, a, b, a)
         */
    }
    #[test]
    fn test_gets_len_2_guided_gs_for_diasem() {
        let diasem: [usize; 9] = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let gener_seqs = guided_gs_list_of_len(2, &diasem);
        assert_eq!(
            gener_seqs,
            vec![vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2])
            ]]
        );
    }

    #[test]
    fn test_gets_guided_gs_for_diamech() {
        let diamech_4sl: [usize; 11] = [1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2];
        let gener_seqs = guided_gs_list(&diamech_4sl);
        assert_eq!(
            BTreeSet::from_iter(gener_seqs.into_iter()),
            BTreeSet::from_iter(
                vec![
                    vec![
                        CountVector::from_slice(&[0, 2]),
                        CountVector::from_slice(&[0, 2]),
                        CountVector::from_slice(&[0, 1]),
                    ],
                    vec![
                        CountVector::from_slice(&[2, 0, 2]),
                        CountVector::from_slice(&[1, 0, 2]),
                        CountVector::from_slice(&[0, 2, 0]),
                        CountVector::from_slice(&[1, 0, 2]),
                        CountVector::from_slice(&[1, 0, 2]),
                        CountVector::from_slice(&[0, 2, 0]),
                        CountVector::from_slice(&[1, 0, 2]),
                        CountVector::from_slice(&[0, 2, 0]),
                        CountVector::from_slice(&[1, 0, 2]),
                    ],
                    vec![
                        CountVector::from_slice(&[0, 0, 1, 2, 2]),
                        CountVector::from_slice(&[0, 0, 1, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 1, 2]),
                        CountVector::from_slice(&[0, 0, 1, 2, 2]),
                    ],
                    vec![
                        CountVector::from_slice(&[0, 0, 0, 1, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 1, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 1, 2, 2]),
                        CountVector::from_slice(&[0, 0, 1, 2, 2, 2]),
                    ],
                    vec![
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 1, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 1, 1, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 1, 1, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 1, 1, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 2, 2, 2]),
                    ],
                    vec![
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 1, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 2, 2, 2, 2]),
                        CountVector::from_slice(&[0, 0, 0, 0, 1, 1, 2, 2, 2]),
                    ],
                ]
                .into_iter()
            )
        );
        /*
                    Valid Guided GS necklaces for 10202010202:
                    02 02 01 02 02 10 20 20 10 20 (21)

                    202 102 020 102 021 020 201 020 210 202 (010) -- abcbbcbcbad -> GS(a,b,c,b,b,c,b,c,b)

                    10202 10202 01020 21020 20102 02102 02010 20210 20201 02021 02020
                    00122 00122 00012 00122 00122 00122 00012 00122 00122 00122 00022-- aaba aaba aac -> GS(a,a,b,a)

                    020201 020210 202010 202102 020102 021020 201020 210202 010202 102020 102021
                    000122 000122 000122 001222 000122 000122 000122 001222 000122 000122 001122 -- aaab aaab aac -> GS(a,a,a,b)

                    20210202 01020210 20201020 21020201 02021020 20102021 02020102 02102020 10202102 02010202 10202010
                    00012222 00001122 00001222 00011222 00001222 00011222 00001222 00001222 00011222 00001222 00001122
                    b        a        c        d        c        d        c        c        d        c        a

                    102020102 021020201 020210202 010202102 020102021 020201020 210202010 202102020 102021020 201020210 202010202
                    000011222 000011222 000012222 000011222 000011222 000001222 000011222 000012222 000011222 000011222 000012222
                    aabaacabaab
        */
    }

    #[test]
    fn test_9l6m10s_unimodular_bug() {
        // LLsmsLmsLsLmsLsmLsLsmLsms - This scale has a unimodular basis:
        // equave: [9, 6, 10]
        // gener_1: [3, 2, 3] (stacked 4 times)
        // gener_2: [2, 1, 2]
        // The determinant of these vectors is ±1, making them a unimodular basis.
        let scale: [usize; 25] = [
            0, 0, 2, 1, 2, 0, 1, 2, 0, 2, 0, 1, 2, 0, 2, 1, 0, 2, 0, 2, 1, 0, 2, 1, 2,
        ];
        let gfs = guide_frames(&scale);
        assert!(!gfs.is_empty(), "Should find at least one guide frame");

        // Check that the unimodular basis is found
        use crate::word_to_profile;
        let profile = word_to_profile(&scale);
        assert!(
            profile.structure.is_some(),
            "Should find a unimodular basis for 9L6m10s scale LLsmsLmsLsLmsLsmLsLsmLsms"
        );

        // Verify the lattice basis contains the expected generators
        if let Some(lattice_basis) = profile.lattice_basis {
            // The basis should contain gener_1=[3,2,3] or gener_2=[2,1,2] or their multiples
            let has_expected_basis = lattice_basis.iter().any(|v| {
                (v[0] == 3 && v[1] == 2 && v[2] == 3) || (v[0] == 2 && v[1] == 1 && v[2] == 2)
            });
            assert!(
                has_expected_basis,
                "Lattice basis should contain one of the expected generators"
            );
        }
    }
}
