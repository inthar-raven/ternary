use crate::helpers::gcd;
use crate::plane_geometry::{ConvexPolygon, advance, project_and_partition};
use crate::words::{Letter, least_mode};

/// Billiard scales with the given signature up to rotation.
pub fn billiard_scales(a: usize, b: usize, c: usize) -> Vec<Vec<Letter>> {
    let (a, b, c) = (a as u32, b as u32, c as u32);
    let d = gcd(gcd(a, b), c);
    let mut result = Vec::<Vec<Letter>>::new();
    let regions: Vec<ConvexPolygon> = project_and_partition(a / d, b / d, c / d);
    let n = a / d + b / d + c / d;
    for region in regions {
        let mut word = Vec::<Letter>::new();
        let mut next_point = region.centroid();
        let first_point = region.centroid();
        for i in 0..n {
            let res = advance(a, b, c, next_point).unwrap();
            next_point = res.0;
            debug_assert!(i == n - 1 || first_point != next_point); // There should not be a subperiod in the orbit of the "billiard ball".
            word.push(res.1);
        }
        result.push({
            // Repeat each word d times
            let mut r = Vec::with_capacity(word.len() * d as usize); // Pre-allocate memory
            for _ in 0..d {
                r.extend_from_slice(&word); // Append a slice of the original vec
            }
            r
        });
        debug_assert_eq!(first_point, next_point); // The last point should come back to the first point. When i = n - 1, next_point gets set to point number n.
    }
    result = result
        .into_iter()
        .map(|scale| least_mode(&scale))
        .collect::<Vec<_>>();
    result.sort();
    result.dedup();
    result
}
