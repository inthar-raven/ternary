#![deny(warnings)]

const STEP_LETTERS: [&'static str; 12] = [
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
use std::collections::{BTreeMap, HashSet};

use actix_rt::task::yield_now;
use actix_web::http::header::ContentType;
use actix_web::middleware::Logger;
use actix_web::web::Query;
use actix_web::HttpResponse;
use actix_web::Result;
use actix_web::{get, App, HttpServer};

use actix_rt::time::timeout;
use itertools::Itertools;
use serde::Deserialize;

use guide::guide_structures;
use guide::GuideFrame;
use interval::JiRatio;
use serde::Serialize;
use serde_json::json;
use ternary::words::{maximum_variety, monotone_lm, monotone_ms, rotation_to_least_mode_booth};
use words::{monotone_s0, CountVector, Letter};

pub mod comb;
#[macro_use]
pub mod equal;
pub mod guide;
pub mod interval;
pub mod ji;
pub mod ji_ratio;
#[macro_use]
pub mod monzo;
pub mod primes;
pub mod utils;
pub mod words;

#[derive(Debug, Serialize)]
pub struct GuideResult {
    /// Either WFGS or multiple interleaved WFGSes
    /// `wfgs` generates a well-formed generator sequence (detempered single-period MOS) subscale.
    gs: Vec<CountVector<char>>,
    /// The aggregate generator
    aggregate: CountVector<char>,
    /// `polyoffset` is the set of intervals that each well-formed generator sequence chain is based on. Always includes `CountVector::ZERO`.
    polyoffset: Vec<CountVector<char>>,
    /// complexity result
    /// The base GS chains in a multiple GS structure don't form interleaved scales. Instead they form a detempered copy of m-edo.
    multiplicity: usize,
    complexity: usize,
}

#[derive(Clone, Debug, Deserialize)]
pub struct GuideMosRequest {
    // The key must match the key used in the URL query
    word: String,
}

#[derive(Clone, Debug, Deserialize)]
pub struct StepSigRequest {
    // The key must match the key used in the URL query
    step_sig: String,
    mos_subst: String,
    // The following 3 are either "on" or "off"
    monotone_lm: String,
    monotone_ms: String,
    monotone_s0: String,
    // The constraints are either "exactly" or "at-most"; 0 for usize is used for blank.
    ggs_len_constraint: String,
    ggs_len: usize,
    complexity_constraint: String,
    complexity: usize,
    mv_constraint: String,
    mv: usize,
}

const TIMEOUT_DURATION: std::time::Duration = std::time::Duration::from_secs(6);

#[actix_web::main]
async fn main() -> Result<(), std::io::Error> {
    // access logs are printed with the INFO level so ensure it is enabled by default
    env_logger::init_from_env(env_logger::Env::new().default_filter_or("info"));
    std::env::set_var("RUST_BACKTRACE", "full");
    println!("Listening on http://127.0.0.1:8080/static/index.html");

    HttpServer::new(|| {
        App::new()
            // Serve static files
            .service(actix_files::Files::new("/static", "static").show_files_listing())
            // Use handlers to serve endpoints
            .service(step_sig_result)
            .service(word_result)
            // Middleware logger
            .wrap(Logger::default())
    })
    .client_request_timeout(TIMEOUT_DURATION) // set timeout
    // Bind to URL
    .bind(("127.0.0.1", 8080))?
    .run()
    .await
}

fn guide_structure_to_json(structure: &GuideFrame, arity: usize) -> GuideResult {
    match structure {
        GuideFrame {
            gs,
            polyoffset,
            multiplicity,
        } => {
            if *multiplicity == 1 {
                GuideResult {
                    gs: gs
                        .into_iter()
                        .map(|cv| {
                            CountVector::from_iter(cv.into_inner().into_iter().map(|(k, v)| {
                                (STEP_LETTERS[min(arity, 12)].chars().nth(k).unwrap(), v)
                            }))
                        })
                        .collect(),
                    aggregate: gs
                        .into_iter()
                        .map(|cv| {
                            CountVector::from_iter(cv.into_inner().into_iter().map(|(k, v)| {
                                (STEP_LETTERS[min(arity, 12)].chars().nth(k).unwrap(), v)
                            }))
                        })
                        .fold(CountVector::<char>::ZERO, |acc, v| acc.add(&v)),
                    polyoffset: polyoffset
                        .into_iter()
                        .map(|cv| {
                            CountVector::from_iter(cv.into_inner().into_iter().map(|(k, v)| {
                                (STEP_LETTERS[min(arity, 12)].chars().nth(k).unwrap(), v)
                            }))
                        })
                        .collect(),
                    multiplicity: 1,
                    complexity: structure.complexity(),
                }
            } else {
                GuideResult {
                    gs: gs
                        .into_iter()
                        .map(|cv| {
                            CountVector::from_iter(cv.into_inner().into_iter().map(|(k, v)| {
                                (STEP_LETTERS[min(arity, 12)].chars().nth(k).unwrap(), v)
                            }))
                        })
                        .collect(),
                    aggregate: gs
                        .into_iter()
                        .map(|cv| {
                            CountVector::from_iter(cv.into_inner().into_iter().map(|(k, v)| {
                                (STEP_LETTERS[min(arity, 12)].chars().nth(k).unwrap(), v)
                            }))
                        })
                        .fold(CountVector::<char>::ZERO, |acc, v| acc.add(&v)),
                    polyoffset: vec![CountVector::<char>::ZERO],
                    multiplicity: *multiplicity,
                    complexity: structure.complexity(),
                }
            }
        }
    }
}

#[get("/word")]
async fn word_result(req: Query<GuideMosRequest>) -> HttpResponse {
    let start = std::time::Instant::now();
    let task = async {
        let string = req.word.to_owned();
        let arity = string.chars().collect::<HashSet<_>>().len();
        let mut word = vec![];
        for c in string.chars() {
            if let Some(letter) = STEP_LETTERS[min(arity, 12)].find(c) {
                word.push(letter);
            } else {
                return HttpResponse::InternalServerError().into();
            }
        }
        let k = rotation_to_least_mode_booth(&word);
        let brightest = [&string[k..], &string[0..k]].concat(); // canonicalize
        println!("{:?}", word);
        let step_sig = {
            let mut r = vec![0; 3];
            for letter in word.clone() {
                if letter < 3 {
                    r[letter] += 1;
                }
            }
            r
        };
        let ji_tunings =
            crate::ji::solve_step_sig_81_odd_limit(&step_sig, crate::monzo::Monzo::OCTAVE, false);
        let ji_tunings: Vec<BTreeMap<char, String>> = ji_tunings
            .into_iter()
            .map(|v| {
                v.into_iter()
                    .enumerate()
                    .map(|(i, monzo)| {
                        (
                            STEP_LETTERS[min(12, arity)].chars().nth(i).unwrap(),
                            format!("{}/{}", monzo.numer(), monzo.denom()),
                        )
                    })
                    .collect::<BTreeMap<_, _>>()
            })
            .collect();
        let structures = guide_structures(&word);
        yield_now().await;
        if let Some(structure) = structures.first() {
            let result = guide_structure_to_json(structure, arity);

            if step_sig.len() == 3 {
                let ed_tunings = crate::equal::ed_tunings_for_ternary(
                    &step_sig,
                    crate::ji_ratio::RawJiRatio::OCTAVE,
                    53,
                    20.0,
                    60.0,
                );
                yield_now().await;
                let profile: serde_json::Value = json!({
                        "structure": Some(result),
                        "lm": monotone_lm(&word),
                        "ms": monotone_ms(&word),
                        "s0": monotone_s0(&word),
                        "mv": maximum_variety(&word),
                });
                yield_now().await;
                HttpResponse::Ok().json(json!({
                    "profile": profile,
                    "brightest": brightest,
                    "ji_tunings": ji_tunings,
                    "ed_tunings": ed_tunings,
                }))
            } else {
                let profile: serde_json::Value = json!({
                        "structure": Some(result),
                        "mv": maximum_variety(&word),
                });
                yield_now().await;
                HttpResponse::Ok().json(json!({
                    "profile": profile,
                    "brightest": brightest,
                    "ji_tunings": ji_tunings,
                }))
            }
        } else {
            if step_sig.len() == 3 {
                let ed_tunings = crate::equal::ed_tunings_for_ternary(
                    &step_sig,
                    crate::ji_ratio::RawJiRatio::OCTAVE,
                    53,
                    20.0,
                    60.0,
                );
                yield_now().await;
                let profile = json!({
                    "structure": None::<GuideFrame>,
                    "lm": monotone_lm(&word),
                    "ms": monotone_ms(&word),
                    "s0": monotone_s0(&word),
                    "mv": maximum_variety(&word),
                });
                yield_now().await;
                HttpResponse::Ok().json(json!({
                    "profile": profile,
                    "brightest": brightest,
                    "ji_tunings": ji_tunings,
                    "ed_tunings": ed_tunings,
                }))
            } else {
                let profile = json!({
                    "structure": None::<GuideFrame>,
                    "mv": maximum_variety(&word),
                });
                yield_now().await;
                HttpResponse::Ok().json(json!({
                    "brightest": brightest,
                    "profile": profile,
                    "ji_tunings": ji_tunings,
                }))
            }
        }
    };
    timeout(TIMEOUT_DURATION, task).await.unwrap_or({
        let elapsed = start.elapsed().as_secs_f64();
        println!("{:?}", elapsed);
        HttpResponse::GatewayTimeout()
            .content_type(ContentType::plaintext())
            .body(format!("{:.2} seconds", elapsed))
    })
}

#[get("/stepsig")]
async fn step_sig_result(req: Query<StepSigRequest>) -> HttpResponse {
    let start = std::time::Instant::now();
    let task = async {
        // Use monotone MOS conditions, GGS length, guide frame complexity, and maximum variety.
        // The `mos_subst` query will be used to determine which function is used to generate the scales before they are filtered by this closure.
        let filtering_cond = |scale: &[Letter]| {
            ((req.monotone_lm != "on") || monotone_lm(scale))
                && ((req.monotone_ms != "on") || monotone_ms(scale))
                && ((req.monotone_s0 != "on") || monotone_s0(scale))
                && (match req.ggs_len {
                    0 => true,
                    l => {
                        let guide_frames = guide_structures(scale);
                        if req.ggs_len_constraint == "exactly" {
                            !guide_frames.is_empty() && guide_frames[0].gs.len() == l
                        } else {
                            !guide_frames.is_empty() && guide_frames[0].gs.len() <= l
                        }
                    }
                })
                && (match req.mv {
                    0 => true,
                    mv => {
                        if req.mv_constraint == "exactly" {
                            maximum_variety(scale) == mv
                        } else {
                            maximum_variety(scale) <= mv
                        }
                    }
                })
                && (match req.complexity {
                    0 => true,
                    c => {
                        let guide_frames = guide_structures(scale);
                        if req.complexity_constraint == "exactly" {
                            !guide_frames.is_empty() && guide_frames[0].complexity() == c
                        } else {
                            !guide_frames.is_empty() && guide_frames[0].complexity() <= c
                        }
                    }
                })
        };
        let query = req.step_sig.to_owned();
        let step_sig = query
            .split(" ")
            .map(|digit| digit.parse().unwrap())
            .collect::<Vec<_>>();
        let arity = step_sig.len();

        let ji_tunings =
            crate::ji::solve_step_sig_81_odd_limit(&step_sig, crate::monzo::Monzo::OCTAVE, false);
        let scales = if req.mos_subst == "on" {
            words::mos_substitution_scales(&step_sig)
        } else {
            crate::comb::necklaces_fixed_content(&step_sig).await
        };
        // Now filter
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
        // Convert scale words in usize to words in char. The client side shouldn't know about the server side representation.
        let scales_data: Vec<serde_json::Value> = {
            let mut result = vec![];
            for scale in scales {
                result.push(json!({
                    "word": String::from_iter(
                        scale.iter()
                            .map(|i| STEP_LETTERS[min(12, arity)].chars().nth(*i).unwrap()),
                    ),
                    "profile": json!({
                        "structure": guide_structures(&scale)
                            .first()
                            .map(|s| guide_structure_to_json(s, arity)),
                        "lm": words::monotone_lm(&scale),
                        "ms": words::monotone_ms(&scale),
                        "s0": words::monotone_s0(&scale),
                        "mv": words::maximum_variety(&scale),
                    }),
                }));
                yield_now().await;
            }
            result
        };
        let ji_tunings: Vec<BTreeMap<char, String>> = ji_tunings
            .into_iter()
            .map(|v| {
                v.into_iter()
                    .enumerate()
                    .map(|(i, monzo)| {
                        (
                            STEP_LETTERS[min(12, arity)].chars().nth(i).unwrap(),
                            format!("{}/{}", monzo.numer(), monzo.denom()),
                        )
                    })
                    .collect::<BTreeMap<_, _>>()
            })
            .collect();
        if step_sig.len() == 3 {
            let ed_tunings = crate::equal::ed_tunings_for_ternary(
                &step_sig,
                crate::ji_ratio::RawJiRatio::OCTAVE,
                53,
                20.0,
                60.0,
            );
            HttpResponse::Ok().json(json!({
                "scales": scales_data,
                "ji_tunings": ji_tunings,
                "ed_tunings": ed_tunings,
            }))
        } else {
            HttpResponse::Ok().json(json!({
                "scales": scales_data,
                "ji_tunings": ji_tunings,
                "ed_tunings": Vec::<Vec<i32>>::new(),
            }))
        }
    };
    timeout(TIMEOUT_DURATION, task).await.unwrap_or({
        let elapsed = start.elapsed().as_secs_f64();
        HttpResponse::GatewayTimeout()
            .content_type(ContentType::plaintext())
            .body(format!("{:.2} seconds", elapsed))
    })
}
