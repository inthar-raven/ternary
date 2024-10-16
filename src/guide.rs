use std::collections::BTreeSet;

use itertools::Itertools;
use serde::Serialize;
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
        dyad_on_degree, offset_vec, rotate, rotations, weak_period, word_on_degree, CountVector,
        Letter, Subtendable,
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
        .map(|chain| weak_period(&chain[..len - 1]))
        .collect::<Vec<_>>()
}
/*
    A guided generator sequence (Guided GS) is a generator sequence made of stacked `k`-steps,
    where `k` is fixed and `gcd(k, scale.len()) == 1`.
    The interval left over after stacking (which is a `k`-step) is different from all the others,
    and where the generators in the generator sequence are distinct from any non-k-step interval.
    There is no need to check the last condition for step vectors.
*/
fn k_step_guided_gs_list(k: usize, neck: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    guided_gs_chains(&stacked_k_steps(k, neck))
}

pub fn stacked_k_steps<T>(k: usize, neck: &[T]) -> Vec<T::Interval>
where
    T: Subtendable + std::fmt::Debug,
{
    // println!("scale: {:?}", neck);
    (0..neck.len())
        .map(|i| word_on_degree(neck, k * i, k))
        .map(|subword| <T as Subtendable>::interval_from_slice(&subword))
        .collect()
}

/// TODO: To get all Guided GSes, we actually need to iterate through all the step classes, not just ones <= len / 2.
/// All Guided GSes that generate a given abstract necklace.
pub fn guided_gs_list(neck: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    let len = neck.len();
    (2..=len - 2) // Don't include 1-step GSes
        .filter(|&k| gcd(k as u64, len as u64) == 1)
        .flat_map(|k| k_step_guided_gs_list(k, neck))
        .collect()
}

/// All Guided GSes of length `l` that generate a given abstract necklace.
pub fn guided_gs_list_of_len(l: usize, neck: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    let neck_len = neck.len();
    (2..=neck_len / 2) // Don't include 1-step GSes
        .filter(|&k| gcd(k as u64, neck_len as u64) == 1)
        .flat_map(|k| k_step_guided_gs_list(k, neck))
        .filter(|vs| l == vs.len())
        .collect()
}

#[allow(unused)]
fn k_step_guided_gs_list_for_subscale(
    k: usize,
    subscale: &[CountVector<usize>],
) -> Vec<Vec<CountVector<usize>>> {
    guided_gs_chains(&stacked_k_steps(k, subscale))
}

// Guided GS of a chain which is represented as `Vec<CountVector<usize>>` rather than `Vec<usize>`.
fn guided_gs_list_for_subscale(subscale: &[CountVector<usize>]) -> Vec<Vec<CountVector<usize>>> {
    if subscale.len() == 2 {
        vec![vec![subscale[0].clone()]]
    } else {
        let len = subscale.len();
        (2..=len / 2) // Don't include 1-step GSes
            .filter(|&k| gcd(k as u64, len as u64) == 1)
            .flat_map(|k| k_step_guided_gs_list_for_subscale(k, subscale))
            .collect()
    }
}

/// This is [Ploidacot](https://en.xen.wiki/w/Ploidacot) but with possible "contorsion".
#[derive(Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub struct Ploidacot {
    /// The number of parts the 2/1 is split into.
    pub ploid: usize,
    /// The number of parts the voicing of 3/2 is split into
    /// If the 3/2 is comprised of several ploids, this is 0 ("acot").
    pub cot: usize,
    /// Determines the voicing of 3/2 used.
    /// 0 indicates 3/2 is used; each increment elevates the voicing by a (2^(1/ploid)).
    pub shear: usize,
}

impl Ploidacot {
    pub fn try_get_ploidacot(scale: &[Letter]) -> Option<Self> {
        println!("gf: {:?}", guide_frames(scale));
        if let Some(gf) = guide_frames(scale).first() {
            let n = scale.len();
            println!("gf: {:?}", gf.clone());
            let po = &gf.polyoffset;
            let gs = &gf.gs;
            let multigen_class = gs.len() * gs[0].len();
            let patent_3_mapping =
                f64::round(n as f64 * f64::log(3.0, std::f64::consts::E) / std::f64::consts::LN_2)
                    as usize;
            let ploid = if po.len() == 1 || gcd(patent_3_mapping as u64, n as u64) > 1 {
                1 // do this instead of polyploid acot
            } else {
                n / po[1].len()
            };
            let patent_fifth_mapping = patent_3_mapping - n;
            let patent_fourth_mapping = 2 * n - patent_3_mapping;
            // Check if one generator in the GS can be interpreted as a fifth
            if gs[0].len() == patent_fifth_mapping || gs[0].len() == patent_fourth_mapping {
                Some(Self {
                    ploid,
                    cot: 1,
                    shear: 0,
                })
            } else {
                // Check if aggregate generator can be interpreted as a fifth
                for i in 0..gs.len() {
                    if patent_fifth_mapping + i * n / ploid == multigen_class {
                        let cot = gs.len();
                        let shear = i;
                        return Some(Self { ploid, cot, shear });
                    }
                }
                let patent_fourth_mapping = 2 * n - patent_3_mapping;
                for i in 0..gs.len() {
                    if patent_fourth_mapping + i * n / ploid == multigen_class {
                        let cot = gs.len();
                        let shear = (cot - ploid - i) % cot;
                        return Some(Self { ploid, cot, shear });
                    }
                }
                None
            }
        } else {
            None
        }
    }
}

/// A guide frame structure for a scale word, consisting of a generator sequence together with a set of offsets or a multiplicity.
/// Multiplicity greater than 1 is a generalization of diregular MV3s;
/// scales of this type always have the number of notes divisible by the multiplicity.
/// It consists of m copies of the same generator sequence offset by `scale.len()` / m steps,
/// where `m` is called the *multiplicity*.
///
/// (Actually, diregular scales are also interleaved; the strand is a detempered 2edo and the polyoffset has a generator.
/// Diachrome is currently considered to have multiplicity 2 rather than having a polyoffset of 6 notes.
/// As ad-hoc as it may seem, this complexity measure rests on the "thin" generator-offset structure of nice ternary scales
/// (i.e. one side is the generator which there is usually much more of than there are offsets)
/// and yields reasonable complexities for quasi-diatonic aberrismic scales. I might rename what I call these guide frame structures in light of this.)
#[derive(Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct GuideFrame {
    /// Either Guided GS or multiple interleaved GSes that are Guided GSes when considered individually.
    /// `gs` generates a guided generator sequence (detempered single-period MOS) subscale.
    pub gs: Vec<CountVector<usize>>,
    /// `polyoffset` is the set of intervals that each guided generator sequence chain is based on. Always includes `CountVector::ZERO`.
    pub polyoffset: Vec<CountVector<usize>>,
}

impl GuideFrame {
    pub fn new_simple(gs: Vec<CountVector<usize>>) -> Self {
        Self {
            gs,
            polyoffset: vec![CountVector::ZERO],
        }
    }
    pub fn new_multiple(gs: Vec<CountVector<usize>>, polyoffset: Vec<CountVector<usize>>) -> Self {
        Self { gs, polyoffset }
    }
    // The comoplexity of a guide frame.
    pub fn complexity(&self) -> usize {
        self.gs.len() * self.polyoffset.len()
    }
    // The multiplicity
    pub fn multiplicity(&self) -> usize {
        self.polyoffset.len()
    }
    // Try to get simple or interleaved guide frames with k-step generators.
    pub fn try_simple(scale: &[usize], k: usize) -> Vec<Self> {
        if scale.is_empty() || gcd(scale.len() as u64, k as u64) != 1 {
            vec![]
        } else {
            k_step_guided_gs_list(k, scale)
                .into_iter()
                .map(|gs| Self {
                    gs,
                    polyoffset: vec![CountVector::ZERO],
                })
                .sorted()
                .dedup()
                .collect::<Vec<_>>()
        }
    }
    pub fn try_multiple(scale: &[usize], multiplicity: usize, k: usize) -> Vec<Self> {
        // The scale cannot be empty and its size must be divisible by `multiplicity`.
        if multiplicity == 1 || scale.is_empty() || scale.len() % multiplicity != 0 {
            vec![]
        } else {
            let d = gcd(k as u64, scale.len() as u64) as usize;
            let co_d = scale.len() / d;
            if co_d % multiplicity != 0 {
                if d == multiplicity {
                    // println!("interleaved");
                    // It's an interleaved scale.
                    let subscales = (0..d)
                        .map(|degree| rotate(scale, degree))
                        .map(|rotation| {
                            stacked_k_steps(d, &rotation[..scale.len()])[..scale.len() / d].to_vec()
                        })
                        .collect::<Vec<_>>();
                    let subscales_cloned = subscales.clone();
                    let subscale_on_root = subscales_cloned
                        .split_first()
                        .expect("since we checked that `scale` is nonempty, this operation should be infallible")
                        .0;
                    // println!("subscale_on_root: {:?}", subscale_on_root);

                    // All of the subscales must be rotations of one another.
                    // `offset_vec()` returns a witness to rotational equivalence (an offset) if there is any;
                    // the offsets are combined to form the polyoffset.
                    // If it returns `None` for any subscale, the whole procedure fails.
                    let maybe_offsets = subscales
                        .into_iter()
                        .enumerate()
                        .map(|(i, subscale)| {
                            offset_vec(subscale_on_root, &subscale).map(|offset| {
                                // `.map()` returns `None` if the previous result is `None` and functorially applies the closure to `Some`s.
                                CountVector::from_slice(&word_on_degree(scale, 0, offset * d + i))
                            })
                        })
                        // collect returns `None` if there is any `None` returned by `map`.
                        .collect::<Option<Vec<CountVector<usize>>>>();
                    // println!("subscale_on_root: {:?}", subscale_on_root);
                    // println!("maybe_offsets: {:?}", maybe_offsets);
                    if let Some(offsets) = maybe_offsets {
                        // sort list of offsets by step class
                        let offsets: Vec<CountVector<usize>> =
                            offsets.into_iter().sorted_by_key(|v| v.len()).collect();
                        // If polyoffset is {0} use multiplicity 1
                        if offsets == [CountVector::ZERO] {
                            guided_gs_list(scale)
                                .into_iter()
                                .map(|gs| Self {
                                    gs,
                                    polyoffset: offsets.to_owned(),
                                })
                                .sorted()
                                .dedup()
                                .collect::<Vec<_>>()
                        } else {
                            guided_gs_list_for_subscale(subscale_on_root)
                                .into_iter()
                                .map(|gs| Self {
                                    gs,
                                    polyoffset: offsets.to_owned(),
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
                    // For every degree of `scale`, get stack of gs_length_limit-many k-steps on that degree.
                    let gen_chains_enumerated = (0..scale.len())
                        .map(|degree| {
                            let mode = rotate(scale, degree);
                            stacked_k_steps(k, &mode)[0..chain_length].to_vec()
                        })
                        .enumerate()
                        .filter(|(_, stack)| {
                            // Each stack is generated by a GS,
                            // but for the GS to be a guided GS, the last element must differ from all previous elements.
                            let mut init = stack.iter().take(chain_length - 1);
                            let last = stack
                                .last()
                                .expect("last exists because gs_length_limit >= 2");
                            // `init` will be nonempty, so the following check won't be vacuous.
                            init.all(|k_step| *k_step != *last)
                        });
                    let gses: Vec<Vec<CountVector<Letter>>> = gen_chains_enumerated
                        .clone()
                        // Take prefix of gs_length_limit - 1 elements and get what GS it is generated by
                        .map(|(_, chain)| weak_period(&chain[0..chain_length - 1]))
                        .sorted()
                        .dedup()
                        .collect();
                    // println!("{:?}", gses);
                    gses.iter()
                        .map(|gs| {
                            (
                                gs,
                                gen_chains_enumerated
                                    .clone()
                                    // Check only the prefix of gs_length_limit - 1 elements, because that's what the guided GS is based on.
                                    .filter(|(_, gen_chain)| {
                                        // println!("gen_chain: {:?}", gen_chain);
                                        //println!(
                                        //     "weak_period: {:?}",
                                        //     weak_period(&gen_chain[..chain_length - 1])
                                        // );
                                        weak_period(&gen_chain[..chain_length - 1]) == *gs.clone()
                                    })
                                    // Get all indices on which this particular GS occurs.
                                    .map(|(i, _)| i)
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .filter(|(_, polyoffset_indices)| {
                            // Filter by whether the number of gen chains each GS occurs on is equal to the multiplcity.
                            // We also need to check that all of the chains are disjoint.
                            let mut union_of_chains: Vec<_> = polyoffset_indices
                                .iter()
                                .flat_map(|first| {
                                    (0..scale.len() / multiplicity)
                                        .map(|i| (first + i * k) % scale.len())
                                        .collect::<Vec<_>>()
                                })
                                .collect();
                            union_of_chains.sort();
                            union_of_chains.dedup();
                            // println!("union_of_chains: {:?}", union_of_chains);
                            let chains_are_disjoint: bool = union_of_chains.len() == scale.len();
                            chains_are_disjoint && polyoffset_indices.len() == multiplicity
                        })
                        .map(|(gs, polyoffset_indices)| {
                            // println!("gs: {:?}", gs);
                            // println!("polyoffset_indices: {:?}", polyoffset_indices);
                            let first_deg = polyoffset_indices[0];
                            let polyoffset: Vec<CountVector<Letter>> = polyoffset_indices
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
                                gs: gs.clone(),
                                polyoffset,
                            }
                        })
                        .collect()
                }
            }
        }
    }
    // Get both Simple and Multiple guide frames.
    fn try_all_variants(scale: &[usize], k: usize) -> Vec<Self> {
        // Let's just do primes for now.
        let prime_factors: Vec<usize> = factorize(scale.len() as u64)
            .into_iter()
            .dedup()
            .map(|p| p as usize)
            .collect();
        let simple_guide_moses: Vec<GuideFrame> = Self::try_simple(scale, k);
        let multiple_guide_moses: Vec<GuideFrame> = if BTreeSet::from_iter(scale.iter()).len() > 1 {
            prime_factors
                .into_iter()
                .flat_map(|p| Self::try_multiple(scale, p, k))
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
        .flat_map(|k| GuideFrame::try_all_variants(scale, k))
        .sorted_by_key(GuideFrame::complexity)
        .collect()
}
#[cfg(test)]
mod tests {
    #[allow(unused)]
    use std::collections::{BTreeMap, BTreeSet, HashSet};

    use crate::words::{CountVector, Letter};

    use super::*;
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
            // Diamech is tricot
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
    #[test]
    fn test_lllmllms() {
        let bad_scale: [usize; 8] = [0, 0, 0, 1, 0, 0, 1, 2];
        let complexity_2_gses = GuideFrame::try_multiple(&bad_scale, 2, 2);
        println!("{:?}", complexity_2_gses);
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
    fn test_stacked_k_steps() {
        let diasem: [usize; 9] = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let one_steps = stacked_k_steps(1, &diasem);
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
        let two_steps = stacked_k_steps(2, &diasem);
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
        let three_steps = stacked_k_steps(3, &diasem);
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
        let two_steps = stacked_k_steps(2, &diasem);
        let chains = guided_gs_chains(two_steps.as_slice());
        assert_eq!(
            chains,
            vec![vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2]),
            ]]
        );
        let four_steps = stacked_k_steps(4, &diasem);
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
        let gen_seqs = guided_gs_list_of_len(2, &diasem);
        assert_eq!(
            gen_seqs,
            vec![vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2])
            ]]
        );
    }

    #[test]
    fn test_gets_guided_gs_for_diamech() {
        let diamech_4sl: [usize; 11] = [1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2];
        let gen_seqs = guided_gs_list(&diamech_4sl);
        assert_eq!(
            BTreeSet::from_iter(gen_seqs.into_iter()),
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
}
