use std::collections::BTreeSet;

use itertools::Itertools;

use crate::{
    helpers::gcd,
    primes::{factorize, is_prime},
    words::{offset_vec, rotate, rotations, weak_period, word_on_degree, CountVector, Subtendable},
};

// Given a necklace of stacked k-steps, where k is fixed,
// get all valid WFGSes using k-steps on any rotation.
fn wfgs_chains<T>(chain: &[T]) -> Vec<Vec<T>>
where
    T: core::fmt::Debug + PartialEq + Clone + Eq + Send,
{
    let len = chain.len();
    rotations(chain)
        .into_iter()
        .filter(|list| !(list[..len - 1].contains(&list[len - 1])))
        .map(|chain| weak_period(&chain[..len - 1]))
        .collect::<Vec<_>>()
}
/*
    A well-formed generator sequence (WFGS) is a generator sequence made of stacked `k`-steps,
    where `k` is fixed and `gcd(k, scale.len()) == 1`.
    The interval left over after stacking (which is a `k`-step) is different from all the others,
    and where the generators in the generator sequence are distinct from any non-k-step interval.
    There is no need to check the last condition for step vectors.
*/
fn k_step_wfgs_list(k: usize, neck: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    wfgs_chains(&stacked_k_steps(k, neck))
}

pub fn stacked_k_steps<T>(k: usize, neck: &[T]) -> Vec<T::Interval>
where
    T: Subtendable + std::fmt::Debug,
{
    (0..neck.len())
        .map(|i| word_on_degree(neck, k * i, k))
        .map(|subword| <T as Subtendable>::interval_from_slice(&subword))
        .collect()
}

/// TODO: To get all WFGSes, we actually need to iterate through all the step classes, not just ones <= len / 2.
/// All WFGSes that generate a given abstract necklace.
pub fn wfgs_list(neck: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    let len = neck.len();
    (2..=len - 2) // Don't include 1-step GSes
        .filter(|&k| gcd(k as u64, len as u64) == 1)
        .flat_map(|k| k_step_wfgs_list(k, neck))
        .collect()
}

/// All WFGSes of length `l` that generate a given abstract necklace.
pub fn wfgs_list_of_len(l: usize, neck: &[usize]) -> Vec<Vec<CountVector<usize>>> {
    let neck_len = neck.len();
    (2..=neck_len / 2) // Don't include 1-step GSes
        .filter(|&k| gcd(k as u64, neck_len as u64) == 1)
        .flat_map(|k| k_step_wfgs_list(k, neck))
        .filter(|vs| l == vs.len())
        .collect()
}

// WFGS of a chain which is represented as `Vec<CountVector<usize>>` rather than `Vec<usize>`.
fn wfgs_list_for_subscale(subscale: &[CountVector<usize>]) -> Vec<Vec<CountVector<usize>>> {
    if subscale.len() == 2 {
        vec![vec![subscale[0].clone()]]
    } else {
        let len = subscale.len();
        (2..=len / 2) // Don't include 1-step GSes
            .filter(|&k| gcd(k as u64, len as u64) == 1)
            .flat_map(|k| k_step_wfgs_list_for_subscale(k, subscale))
            .collect()
    }
}

fn k_step_wfgs_list_for_subscale(
    k: usize,
    subscale: &[CountVector<usize>],
) -> Vec<Vec<CountVector<usize>>> {
    wfgs_chains(&stacked_k_steps(k, subscale))
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
    /// Either WFGS or multiple interleaved GSes that are WFGSes when considered individually.
    /// `gs` generates a well-formed generator sequence (detempered single-period MOS) subscale.
    pub gs: Vec<CountVector<usize>>,
    /// `polyoffset` is the set of intervals that each well-formed generator sequence chain is based on. Always includes `CountVector::ZERO`.
    pub polyoffset: Vec<CountVector<usize>>,
    /// The base GS chains in a multiple GS structure don't form interleaved scales. Instead they form a detempered copy of m-edo.
    pub multiplicity: usize,
}

impl GuideFrame {
    pub fn new_multiple(gs: Vec<CountVector<usize>>, multiplicity: usize) -> Self {
        Self {
            gs,
            multiplicity,
            polyoffset: vec![CountVector::ZERO],
        }
    }
    pub fn new_simple(gs: Vec<CountVector<usize>>, polyoffset: Vec<CountVector<usize>>) -> Self {
        Self {
            gs,
            polyoffset,
            multiplicity: 1,
        }
    }
    // The comoplexity of a guide frame.
    pub fn complexity(&self) -> usize {
        self.gs.len() * self.polyoffset.len() * self.multiplicity
    }
    // Try to get simple or interleaved guide frames with k-step generators.
    pub fn try_simple_or_interleaved(scale: &[usize], k: usize) -> Vec<Self> {
        if scale.is_empty() {
            vec![]
        } else {
            let d = gcd(scale.len() as u64, k as u64) as usize;
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
                    wfgs_list(scale)
                        .into_iter()
                        .map(|gs| Self {
                            gs,
                            polyoffset: offsets.to_owned(),
                            multiplicity: 1,
                        })
                        .sorted()
                        .dedup()
                        .collect::<Vec<_>>()
                } else {
                    wfgs_list_for_subscale(subscale_on_root)
                        .into_iter()
                        .map(|gs| Self {
                            gs,
                            polyoffset: offsets.to_owned(),
                            multiplicity: 1,
                        })
                        .sorted()
                        .dedup()
                        .collect::<Vec<_>>()
                }
            } else {
                vec![]
            }
        }
    }
    pub fn try_multiple(scale: &[usize], multiplicity: usize, k: usize) -> Vec<Self> {
        // The scale cannot be empty and its size must be divisible by `multiplicity`.
        if is_prime(scale.len() as u64) || scale.is_empty() || scale.len() % multiplicity != 0 {
            vec![]
        } else {
            let d = gcd(k as u64, scale.len() as u64) as usize;
            if d == 1 {
                // One-strand multi-GS scales.
                // For each rotation (there are `scale.len() / multiplicity` of them we want to consider),
                // we will get a collection with `multiplicity` elements
                // all of which have to be equal.
                let valid_gses: Vec<Vec<CountVector<usize>>> = (0..scale.len() / multiplicity)
                    .map(|degree| {
                        // Stack k-steps and split the result into `multiplicity` vecs of equal length
                        let s = stacked_k_steps(k, &rotate(scale, degree));
                        (0..multiplicity)
                            .map(|i| {
                                s[scale.len() / multiplicity * i
                                    ..scale.len() / multiplicity * (i + 1)]
                                    .to_vec()
                            })
                            .collect::<Vec<_>>()
                    })
                    .filter(|gses| {
                        // To qualify as a valid GS, the last element cannot be in the GS.
                        !gses[0][0..scale.len() / multiplicity - 1]
                            .iter()
                            .contains(&gses[0][scale.len() / multiplicity - 1])
                        // Ignoring the last element, all of the vecs have to be equal.
                        && gses.iter()
                            .map(|gs|
                                gs[0..scale.len() / multiplicity - 1].to_vec()
                            )
                            .all_equal()
                    })
                    // Get the first generator chain, which should exist and be equal to all the other GSes in the list.
                    .map(|gses| gses[0][0..scale.len() / multiplicity - 1].to_vec())
                    // Turn the chain into a generator sequence recipe.
                    .map(|gs| weak_period(&gs))
                    .collect();
                // Convert each valid multi-GS into a `GuideFrame` struct.
                valid_gses
                    .into_iter()
                    .map(|gs| Self {
                        gs,
                        polyoffset: vec![CountVector::ZERO],
                        multiplicity,
                    })
                    .sorted()
                    .dedup()
                    .collect()
            } else {
                // Interleaved multi-GS scales are not handled yet.
                vec![]
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
        let simple_guide_moses: Vec<GuideFrame> = Self::try_simple_or_interleaved(scale, k);
        let multiple_guide_moses: Vec<GuideFrame> = if BTreeSet::from_iter(scale.iter()).len() > 1 {
            prime_factors
                .into_iter()
                .flat_map(|p| Self::try_multiple(scale, p, k))
                .collect()
        } else {
            vec![]
        };
        [simple_guide_moses, multiple_guide_moses].concat()
    }
}

/// Return the collection of guide frames for the given scale word, sorted by complexity.
pub fn guide_structures(scale: &[usize]) -> Vec<GuideFrame> {
    (2..=scale.len() - 2) // steps subtended by generator used for the guided generator sequence
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
    fn test_fix_bug_for_4sr() {
        let diamech_4sr: [Letter; 11] = [0, 2, 0, 1, 0, 2, 0, 2, 0, 1, 2];

        assert!(wfgs_list_of_len(3, &diamech_4sr).contains(&vec![
            CountVector::from_slice(&[0, 2]),
            CountVector::from_slice(&[0, 1]),
            CountVector::from_slice(&[0, 2]),
        ]));
        assert!(wfgs_list(&diamech_4sr).contains(&vec![
            CountVector::from_slice(&[0, 2]),
            CountVector::from_slice(&[0, 1]),
            CountVector::from_slice(&[0, 2]),
        ]));
        let guide_frames = GuideFrame::try_simple_or_interleaved(&diamech_4sr, 2);
        println!("{:?}", guide_frames);
        assert!(guide_frames.contains(&GuideFrame::new_simple(
            vec![
                CountVector::from_slice(&[0, 2]),
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2]),
            ],
            vec![CountVector::ZERO],
        )));
    }

    #[test]
    fn test_try_simple_or_interleaved() {
        let diachrome_5sc = [0, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 1];
        let should_be_nonempty = GuideFrame::try_simple_or_interleaved(&diachrome_5sc, 6);
        assert_ne!(should_be_nonempty, vec![]);
        println!("{:?}", should_be_nonempty);
        println!(
            "{:?}",
            should_be_nonempty
                .into_iter()
                .map(|gf| gf.complexity())
                .collect::<Vec<_>>()
        );
    }
    #[test]
    fn test_wfgs_based_guide_structure() {
        let pinedye = [0, 0, 1, 0, 1, 0, 0, 2];
        let guide_moses = guide_structures(&pinedye);
        println!("Pinedye has guide MOS structures: {:?}", guide_moses);
        assert!(guide_moses.contains(&GuideFrame::new_simple(
            vec![
                CountVector::from_slice(&[0, 0, 2]),
                CountVector::from_slice(&[0, 0, 1]),
                CountVector::from_slice(&[0, 0, 1]),
            ],
            vec![CountVector::ZERO],
        )));

        let diasem = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let guide_moses = guide_structures(&diasem);
        println!("Diasem has guide MOS structures: {:?}", guide_moses);
        assert!(guide_moses.contains(&GuideFrame::new_simple(
            vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2])
            ],
            vec![CountVector::ZERO],
        )));
        assert_eq!(
            GuideFrame::new_simple(
                vec![
                    CountVector::from_slice(&[0, 1]),
                    CountVector::from_slice(&[0, 2])
                ],
                vec![CountVector::ZERO],
            )
            .complexity(),
            2
        );

        let blackdye: [usize; 10] = [0, 1, 0, 2, 0, 1, 0, 2, 0, 2];
        let guide_moses = guide_structures(&blackdye);
        println!("Blackdye has guide MOS structures: {:?}", guide_moses);
        assert!(guide_moses.contains(&GuideFrame::new_simple(
            vec![CountVector::from_slice(&[0, 0, 1, 2])],
            vec![CountVector::ZERO, CountVector::from_slice(&[0])],
        )));
        assert_eq!(
            GuideFrame::new_simple(
                vec![CountVector::from_slice(&[0, 0, 1, 2])],
                vec![CountVector::ZERO, CountVector::from_slice(&[0])],
            )
            .complexity(),
            2
        );

        let diamech_4sl: [usize; 11] = [1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2];
        let guide_moses = guide_structures(&diamech_4sl);
        println!("Diamech has guide MOS structures: {:?}", guide_moses);
        assert!(guide_moses.contains(&GuideFrame::new_simple(
            vec![
                CountVector::from_slice(&[0, 2]),
                CountVector::from_slice(&[0, 2]),
                CountVector::from_slice(&[0, 1]),
            ],
            vec![CountVector::ZERO],
        )));
        assert_eq!(
            GuideFrame::new_simple(
                vec![
                    CountVector::from_slice(&[0, 2]),
                    CountVector::from_slice(&[0, 2]),
                    CountVector::from_slice(&[0, 1]),
                ],
                vec![CountVector::ZERO],
            )
            .complexity(),
            3
        );

        let diachrome_5sc = [0, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2, 1];
        let guide_moses = guide_structures(&diachrome_5sc);
        println!("Diachrome has guide MOS structures: {:?}", guide_moses);
        assert!(guide_moses.contains(&GuideFrame::new_multiple(
            vec![CountVector::from_slice(&[0, 0, 1, 2, 2]),],
            2,
        )));
        assert_eq!(
            GuideFrame::new_multiple(vec![CountVector::from_slice(&[0, 0, 1, 2, 2])], 2)
                .complexity(),
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
    fn test_wfgs_chains() {
        let diasem: [usize; 9] = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let two_steps = stacked_k_steps(2, &diasem);
        let chains = wfgs_chains(two_steps.as_slice());
        assert_eq!(
            chains,
            vec![vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2]),
            ]]
        );
        let four_steps = stacked_k_steps(4, &diasem);
        let chains = wfgs_chains(four_steps.as_slice());
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
        Valid WFGS necklaces for 010201020:
        10 20 10 20 01 02 01 02 (closing 00) -> abababab -> GS(a, b)
        1020 1020 0102 0102 0010 2010 2001 0201 (closing 0200) -> aaaabaaac -> GS(a, a, a, a, b)
        2010 2001 0201 0200 1020 1020 0102 0102 (closing 0010) -> aaabaaaac -> GS(a, a, a, b, a)
         */
    }
    #[test]
    fn test_gets_len_2_wfgs_for_diasem() {
        let diasem: [usize; 9] = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let gen_seqs = wfgs_list_of_len(2, &diasem);
        assert_eq!(
            gen_seqs,
            vec![vec![
                CountVector::from_slice(&[0, 1]),
                CountVector::from_slice(&[0, 2])
            ]]
        );
    }

    #[test]
    fn test_gets_wfgs_for_diamech() {
        let diamech_4sl: [usize; 11] = [1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 2];
        let gen_seqs = wfgs_list(&diamech_4sl);
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
                    Valid WFGS necklaces for 10202010202:
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
