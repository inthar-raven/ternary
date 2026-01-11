//! Algorithm for generating ternary [billiard scales](https://en.xen.wiki/w/Hypercubic_billiard_word).
//!
//! This module implements a geometric approach to scale generation: a "billiard ball"
//! bounces inside a 3D cube, and the sequence of faces it hits determines
//! the scale's step pattern. A 2D projection of the cube is used to simplify the computation.
//!
//! # Algorithm
//!
//! 1. Project the unit cube in R³ via a map with `(a, b, c)` in its kernel
//! 2. Partition the projected hexagon by constraint planes
//! 3. For each region, trace a trajectory and record which coordinate changes
//! 4. The sequence of changes (0=L, 1=m, 2=s) forms the scale word
//!
//! # Examples
//!
//! ```
//! use ternary::billiard::billiard_scales;
//!
//! // Generate all billiard scales with signature 5L 2m 2s
//! let scales = billiard_scales(5, 2, 2);
//! assert!(!scales.is_empty());
//!
//! // Each scale has the correct number of steps
//! for scale in &scales {
//!     assert_eq!(scale.len(), 9);  // 5 + 2 + 2 = 9 notes
//! }
//! ```
//!
//! # Properties
//!
//! Billiard scales have special structural properties:
//! - They are always [deletion-MOS](https://en.xen.wiki/w/Deletion-MOS)
//!
//! Not all ternary scales are billiard scales — this is a strict subset.

use core::fmt;

use num_rational::Rational32 as r32;

use crate::helpers::gcd;
use crate::plane_geometry::{ConvexPolygon, Line, Point, PointLineConfiguration, Slope};
use crate::words::{Letter, least_mode};

// ERRORS

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
/// Error type for illegal billiard trajectory states.
pub enum BadBilliardState {
    // Projected cube is invalid.
    ProjectedCubeInvalid((u32, u32, u32)),
    // Point is not inside the projected cube,
    PointNotInCube(Point),
    // Point is inside the projected cube but falls on a projected constraint plane
    PointOnConstraintPlane(Point),
}

impl fmt::Display for BadBilliardState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ProjectedCubeInvalid((a, b, c)) => {
                write!(f, "Projected cube for signature {a}L{b}m{c}s is invalid")
            }
            Self::PointNotInCube(pt) => write!(f, "Point {} is not inside the projected cube", pt),
            Self::PointOnConstraintPlane(pt) => write!(
                f,
                "Point {} is inside the projected cube but falls on a projected constraint plane",
                pt
            ),
        }
    }
}

impl std::error::Error for BadBilliardState {}

fn projected_constraint_planes(a: u32, b: u32, c: u32) -> (Vec<Line>, Vec<Line>, Vec<Line>) {
    let first = (0..b + c + 1)
        .map(|i| {
            Line::from_slope_and_point(
                &Slope::raw(0, 1),
                &Point::new(r32::new(0, 1), r32::new(-(b as i32) + (i as i32), c as i32)),
            )
        })
        .collect(); // horizontal lines going up (Left)
    let second = (0..a + c + 1)
        .map(|j| {
            Line::from_slope_and_point(
                &Slope::raw(1, 0),
                &Point::new(r32::new((c as i32) - (j as i32), c as i32), r32::new(0, 1)),
            )
        })
        .collect(); // vertical lines going left (Left)
    let third = (0..a + b + 1)
        .map(|k| {
            Line::from_slope_and_point(
                &Slope::raw(b as i32, a as i32),
                &Point::new(r32::new(0, 1), r32::new(-(b as i32) + (k as i32), a as i32)),
            )
        })
        .collect(); // lines of pos. slope going up (Left)
    (first, second, third)
}

/// Project the unit cube in R^3 via a linear map that has (a,b,c) in its kernel.
fn projected_cube(a: u32, b: u32, c: u32) -> Option<ConvexPolygon> {
    if a > 0 && b > 0 && c > 0 {
        let projected_cube_vertices = vec![
            Point::new(r32::new(0, 1), r32::new(1, 1)),
            Point::new(r32::new(1, 1), r32::new(1, 1)),
            Point::new(r32::new(1, 1), r32::new(0, 1)),
            Point::new(
                r32::new(-(a as i32), c as i32),
                r32::new(-(b as i32), c as i32) + r32::new(1, 1),
            ),
            Point::new(
                r32::new(-(a as i32), c as i32),
                r32::new(-(b as i32), c as i32),
            ),
            Point::new(
                r32::new(-(a as i32), c as i32) + r32::new(1, 1),
                r32::new(-(b as i32), c as i32),
            ),
        ];
        ConvexPolygon::new(&projected_cube_vertices)
    } else {
        None
    }
}

/// Partition the projected unit cube with the projected constraint planes.
fn project_and_partition(a: u32, b: u32, c: u32) -> Option<Vec<ConvexPolygon>> {
    projected_cube(a, b, c).map(|projected_cube| {
        let lineses = projected_constraint_planes(a, b, c); // it's a collection of collections of `Line`s
        let mut first_shapes = Vec::<ConvexPolygon>::new();
        let mut second_shapes = Vec::<ConvexPolygon>::new();
        let mut third_shapes = Vec::<ConvexPolygon>::new();
        let mut i: usize = 0;
        let mut remaining: ConvexPolygon = projected_cube;
        loop {
            let (maybe_remaining, maybe_next) = remaining.subdivide(&(lineses.0)[i]);
            // The remaining portion is to the Left.
            match (maybe_remaining, maybe_next) {
                (Some(r), None) => {
                    remaining = r;
                }
                (Some(r), Some(next)) => {
                    remaining = r;
                    first_shapes.push(next);
                }
                (None, Some(next)) => {
                    // This is the last `next` you need to look at.
                    first_shapes.push(next);
                    break;
                }
                (_, _) => {
                    break;
                }
            }
            i += 1;
        }
        for polygon in first_shapes {
            let mut j: usize = 0;
            let mut remaining2: ConvexPolygon = polygon;
            loop {
                let (maybe_remaining, maybe_next) = remaining2.subdivide(&(lineses.1)[j]); // The remaining portion is to the Left.
                match (maybe_remaining, maybe_next) {
                    (Some(r), None) => {
                        remaining2 = r;
                    }
                    (Some(r), Some(next)) => {
                        remaining2 = r;
                        second_shapes.push(next);
                    }
                    (None, Some(next)) => {
                        second_shapes.push(next);
                        break;
                    }
                    (_, _) => {
                        break;
                    }
                }
                j += 1;
            }
        }
        for polygon in second_shapes {
            let mut j: usize = 0;
            let mut remaining3: ConvexPolygon = polygon;
            loop {
                let (maybe_remaining, maybe_next) = remaining3.subdivide(&(lineses.2)[j]);
                // The remaining portion is to the Left.
                match (maybe_remaining, maybe_next) {
                    (Some(r), None) => {
                        remaining3 = r;
                    }
                    (Some(r), Some(next)) => {
                        remaining3 = r;
                        third_shapes.push(next);
                    }
                    (None, Some(next)) => {
                        third_shapes.push(next);
                        break;
                    }
                    (_, _) => {
                        break;
                    }
                }
                j += 1;
            }
        }
        third_shapes
    })
}

/// Advance a projected point, `pt`, by one letter according to what integer coordinate it next hits as it traverses in the given velocity,
/// and return the next point and the corresponding letter, one of "0", "1", and "2".
fn advance(a: u32, b: u32, c: u32, pt: Point) -> Result<(Point, Letter), BadBilliardState> {
    if let Some(projected_cube) = projected_cube(a, b, c) {
        // Consider: _ /
        //            |
        if !projected_cube.has_point_inside(pt) {
            Err(BadBilliardState::PointNotInCube(pt))
        } else {
            let projected_x_y_equals_1_1 = Line::from_slope_and_point(
                &Slope::raw(b as i32, a as i32),
                &Point::new(r32::new(0, 1), r32::new(-(b as i32) + (a as i32), a as i32)),
            ); // the /
            let projected_x_z_equals_1_1 = Line::from_slope_and_point(
                &Slope::raw(0, 1),
                &Point::new(r32::new(0, 1), r32::new(-(b as i32) + (c as i32), c as i32)),
            ); // the _, oriented right
            let projected_y_z_equals_1_1 = Line::from_slope_and_point(
                &Slope::raw(1, 0),
                &Point::new(r32::new((c as i32) - (a as i32), c as i32), r32::new(0, 1)),
            ); // the |, oriented up
            if projected_x_y_equals_1_1.point_line_config(pt) == PointLineConfiguration::Right
                && projected_y_z_equals_1_1.point_line_config(pt) == PointLineConfiguration::Right
            {
                // Below the / and to the right of the | => x has changed, change x coordinate to 0 and project; move x coordinate to the left
                Ok((pt - Point::new(r32::new(1, 1), r32::new(0, 1)), 0)) // new point is one unit to the left
            } else if projected_x_z_equals_1_1.point_line_config(pt) == PointLineConfiguration::Left
                && projected_x_y_equals_1_1.point_line_config(pt) == PointLineConfiguration::Left
            {
                // Above the _ and to the left of the / => y has changed, change y coordinate to 0 and project; move y coordinate down
                Ok((pt - Point::new(r32::new(0, 1), r32::new(1, 1)), 1)) // new point is one unit below
            } else if projected_x_z_equals_1_1.point_line_config(pt)
                == PointLineConfiguration::Right
                && projected_y_z_equals_1_1.point_line_config(pt) == PointLineConfiguration::Left
            {
                // Below the _ and to the left of the | => z has changed, change y coordinate to 0 and project
                Ok((
                    pt + Point::new(r32::new(a as i32, c as i32), r32::new(b as i32, c as i32)),
                    2,
                )) // new point is one unit below in the z-direction, and e3 is projected to (-a/c, -b/c)
            } else {
                Err(BadBilliardState::PointOnConstraintPlane(pt))
            }
        }
    } else {
        Err(BadBilliardState::ProjectedCubeInvalid((a, b, c)))
    }
}

/// Generate all billiard scales with the given step signature.
///
/// Returns scales in canonical form (lexicographically first rotation),
/// sorted and deduplicated.
///
/// # Arguments
///
/// * `a` - Number of L (large) steps
/// * `b` - Number of m (medium) steps
/// * `c` - Number of s (small) steps
///
/// # Examples
///
/// ```
/// use ternary::billiard::billiard_scales;
///
/// // Generate billiard scales for 5L 2m 2s
/// let scales = billiard_scales(5, 2, 2);
/// assert!(!scales.is_empty());
///
/// // All scales are in canonical form (sorted)
/// for scale in &scales {
///     assert_eq!(scale.iter().filter(|&&x| x == 0).count(), 5);  // 5 L steps
/// }
/// ```
///
/// # Returns
///
/// Empty vector if the signature is invalid (any component is 0).
pub fn billiard_scales(a: usize, b: usize, c: usize) -> Vec<Vec<Letter>> {
    let (a, b, c) = (a as u32, b as u32, c as u32);
    let d = gcd(gcd(a, b), c);
    if d != 0 {
        let mut result = Vec::<Vec<Letter>>::new();
        let n = a / d + b / d + c / d;
        if let Some(regions) = project_and_partition(a / d, b / d, c / d) {
            for region in regions {
                let mut word = Vec::<Letter>::new();
                let mut next_point = region.centroid();
                let first_point = region.centroid();
                for i in 0..n {
                    let res = advance(a, b, c, next_point)
                    .expect("This should never fail as long as you begin from a point inside a polygonal region");
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
        } else {
            vec![]
        }
    } else {
        vec![]
    }
}
