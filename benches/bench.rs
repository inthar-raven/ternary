use criterion::{black_box, criterion_group, criterion_main, Criterion};

use itertools::Itertools;

use ternary::guide::{stacked_k_steps, GuideFrame};
use ternary::primes::is_prime;
use ternary::words::{rotate, rotations, weak_period, CountVector, Letter};
use ternary::{comb::necklaces_fixed_content, interval::Dyad, ji::odd_limit, ji_ratio::RawJiRatio};

use ternary::utils::gcd as my_gcd;

// Returns the octave-reduced intervals of a specified odd limit.
pub fn odd_limit_itertools(limit: u64) -> Vec<RawJiRatio> {
    let odds = (0..=(limit - 1) / 2).map(|i| 2 * i + 1);
    odds.clone()
        .cartesian_product(odds)
        .map(|(m, n)| {
            RawJiRatio::try_new(m, n)
                .expect("every element of `odds` should be positive")
                .rd(RawJiRatio::OCTAVE)
        })
        .sorted_unstable()
        .dedup()
        .collect()
}

// Compare ways to get Cartesian product iterators: `pairs()` vs. `itertools::cartesian_product()`
pub fn cartesian_bench(c: &mut Criterion) {
    c.bench_function("pairs()", |b| b.iter(|| odd_limit(black_box(81))));
    c.bench_function("itertools Cartesian product", |b| {
        b.iter(|| odd_limit_itertools(black_box(81)))
    });
}

pub fn bench_gcd(c: &mut Criterion) {
    c.bench_function("my gcd", |b| {
        b.iter(|| {
            for i in 2..100 {
                for j in 2..100 {
                    my_gcd(i, j);
                }
            }
        })
    });
}

pub fn bench_prime_sieve(c: &mut Criterion) {
    /*
    c.bench_function("my eratosthenes", |b| {
        b.iter(|| {
            for i in black_box(100..=1000) {
                ternary::primes::eratosthenes(black_box(i));
            }
        })
    });
    */
    c.bench_function("my atkin", |b| {
        b.iter(|| {
            for i in black_box(100..=1000) {
                ternary::primes::atkin(black_box(i));
            }
        })
    });
}

pub fn bench_factorize(c: &mut Criterion) {
    c.bench_function("my factorize", |b| {
        b.iter(|| {
            for i in black_box(1..=1000) {
                ternary::primes::factorize(black_box(i));
            }
        })
    });
}

/*
pub fn threads(c: &mut Criterion) {
    c.bench_function("single-threaded", |b| {
        b.iter(|| {
            ternary::comb::necklaces_fixed_content(&[6, 4, 3]);
            ternary::ji::solve_step_sig_81_odd_limit(
                &[5, 2, 3],
                ternary::monzo::Monzo::OCTAVE,
                false,
            );
        })
    });
    c.bench_function("multithreaded", |b| {
        b.iter(|| {
            let ji_tunings_thread = std::thread::spawn(move || {
                ternary::ji::solve_step_sig_81_odd_limit(
                    &[5, 2, 3],
                    ternary::monzo::Monzo::OCTAVE,
                    false,
                )
            });
            let scales_thread =
                std::thread::spawn(move || ternary::comb::necklaces_fixed_content(&[5, 2, 3]));

            ji_tunings_thread.join().unwrap();
            scales_thread.join().unwrap();
        })
    });
}
pub fn bench_gs(c: &mut Criterion) {
    c.bench_function("Finding WFGSes", |b| {
        b.iter(|| {
            for neck in necklaces_fixed_content(&[5, 2, 4]) {
                ternary::guidemos::wfgs_list(&neck);
            }
        })
    });
}
*/
pub fn cs_checking(c: &mut Criterion) {
    c.bench_function("exact", |b| {
        b.iter(|| ternary::ji::is_cs_ji_scale(black_box(&ternary::ji_ratio::RawJiRatio::EURYBIA)))
    });
    c.bench_function("approximate", |b| {
        b.iter(|| {
            ternary::ji::is_cs_ji_scale_fast(black_box(&ternary::ji_ratio::RawJiRatio::EURYBIA))
        })
    });
}

pub fn bench_guide(c: &mut Criterion) {
    fn try_multiple_1(scale: &[usize], multiplicity: usize, k: usize) -> Vec<GuideFrame> {
        // The scale cannot be empty and its size must be divisible by `multiplicity`.
        if is_prime(scale.len() as u64) || scale.len() == 0 || scale.len() % multiplicity != 0 {
            vec![]
        } else {
            let d = my_gcd(k as u64, scale.len() as u64) as usize;
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
                            .into_iter()
                            .map(|i| {
                                s[scale.len() / multiplicity * i
                                    ..scale.len() / multiplicity * (i + 1)]
                                    .to_vec()
                            })
                            .collect::<Vec<_>>()
                    })
                    .filter(|gses| {
                        // Ignoring the last element, all of the vecs have to be equal.
                        let mut iter = gses.iter();
                        iter.clone().map(|gs|
                                gs[0..scale.len() / multiplicity - 1].to_vec()
                            )
                            .all_equal()
                            // To qualify as a valid GS, the last element cannot be in the GS.
                        && iter.all(|gs| !gs[0..scale.len() / multiplicity - 1]
                                .into_iter()
                                .contains(&gs[scale.len() / multiplicity - 1]))
                    })
                    // Get the first generator chain, which should exist and be equal to all the other GSes in the list.
                    .map(|gses| gses[0][0..scale.len() / multiplicity - 1].to_vec())
                    // Turn the chain into a GS recipe.
                    .map(|gs| weak_period(&gs))
                    .collect();
                // Convert each valid multi-GS into a `GuideFrame` struct.
                valid_gses
                    .into_iter()
                    .map(|gs| GuideFrame::new_multiple(gs, multiplicity))
                    .sorted()
                    .dedup()
                    .collect()
            } else {
                // Interleaved multi-GS scales are not handled yet.
                vec![]
            }
        }
    }
    c.bench_function("iterator practice 1", |b| {
        b.iter(|| GuideFrame::try_multiple(black_box(&[2, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1]), 2, 5))
    });
    c.bench_function("iterator practice 2", |b| {
        b.iter(|| try_multiple_1(black_box(&[2, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1]), 2, 5))
    });
}

pub fn bench_lyndon(c: &mut Criterion) {
    pub fn least_mode_naive(scale: &[Letter]) -> Vec<Letter> {
        let sorted_modes = rotations(scale).into_iter().sorted().collect::<Vec<_>>();
        let result = sorted_modes[0].to_owned();
        result
    }
    c.bench_function("naive", |b| {
        b.iter(|| {
            for scale in black_box(ternary::words::mos_substitution_scales(&[5, 2, 10])) {
                ternary::words::least_mode_naive(&scale);
            }
        });
    });
    c.bench_function("Booth", |b| {
        b.iter(|| {
            for scale in black_box(ternary::words::mos_substitution_scales(&[5, 2, 10])) {
                ternary::words::least_mode_booth(&scale);
            }
        })
    });
}
criterion_group!(
    benches,
    // cartesian_bench,
    // bench_gcd,
    // cs_checking,
    // bench_prime_sieve,
    // bench_factorize,
    // bench_guide,
    bench_lyndon,
);
criterion_main!(benches);
