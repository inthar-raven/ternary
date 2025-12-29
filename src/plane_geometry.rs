// Geometry of points and lines in the plane, implemented with rational numbers. The implementation is janky, though.
// Can be used to enumerate ternary billiard scales.
use num_traits::sign::Signed;
use std::fmt;
use std::iter::zip;
use std::ops::Add;
use std::ops::Neg;
use std::ops::Sub;

use num_rational::Rational32 as r32;

use crate::helpers::gcd;

#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
/// The possible configurations of a point and an oriented line on the plane.
pub enum PointLineConfiguration {
    Left,      // The point is on the left half-plane.
    Right,     // The point is on the right half-plane.
    OnTheLine, // The point is on the line.
}

fn is_between_r32(a: r32, b: r32, c: r32) -> bool {
    if b == c {
        a == b
    } else if b < c {
        b < a && a < c
    } else {
        c < a && a < b
    }
}

/// A rational number or unsigned infinity; represents a valid slope for a line in Q^2.
#[derive(Copy, Clone, Debug)]
pub struct Slope {
    numer: i32,
    denom: i32,
}

impl PartialEq for Slope {
    fn eq(&self, other: &Self) -> bool {
        self.numer * other.denom == other.numer * self.denom
    }
}

impl fmt::Display for Slope {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}", self.numer, self.denom)
    }
}

impl Slope {
    /// Converts a `Slope` to a rational; returns Err if `denom` is 0.
    pub fn as_rational(&self) -> Result<r32, String> {
        if self.denom != 0 {
            Ok(r32::new(self.numer, self.denom))
        } else {
            Err("Attempted to convert a Slope with `denom` == 0 into a Rational64".to_string())
        }
    }

    /// Converts the projective rational to a rational; panics if the denominator is zero.
    pub fn as_rational_raw(&self) -> r32 {
        r32::new(self.numer, self.denom)
    }

    /// Create a new reduced Slope; panics if both numer and denom are 0.
    pub fn new(n: i32, d: i32) -> Result<Self, String> {
        let (mut numer, mut denom) = (n, d);
        if numer == 0 && denom == 0 {
            Err("Attempted to create a Slope with both `numer` and `denom` equal to 0".to_string())
        } else {
            if denom < 0 {
                numer = -numer;
                denom = -denom;
            }
            let d = gcd(numer as u32, denom as u32) as i32;
            if d > 1 {
                numer /= d;
                denom /= d;
            }
            Ok(Self { numer, denom })
        }
    }

    /// Create a new Slope without checking that `numer` == `denom` == 0, or reducing the rational afterwards.
    pub fn raw(numer: i32, denom: i32) -> Self {
        Self { numer, denom }
    }

    /// Create a new `Slope` from a rational.
    pub fn rational(rat: r32) -> Self {
        Self {
            numer: *rat.numer(),
            denom: *rat.denom(),
        }
    }

    /// Create a new `Slope` from an integer.
    pub fn integer(n: i32) -> Self {
        // use 'raw' since an integer is always reduced
        Self::raw(n, 1)
    }

    /// Shorthand for slope 0.
    pub fn zero() -> Self {
        Self::raw(0, 1)
    }

    /// Shorthand for slope 1.
    pub fn one() -> Self {
        Self::raw(1, 1)
    }

    /// Shorthand for infinite slope.
    pub fn infinity() -> Self {
        Self::raw(1, 0)
    }
}

/// A finite point; a member of Q^2.
#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
pub struct Point {
    // The x-coordinate.
    x: r32,
    // The y-coordinate.
    y: r32,
}

impl Add for Point {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for Point {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Neg for Point {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl Point {
    pub fn new(x: r32, y: r32) -> Self {
        Self { x, y }
    }

    pub fn zero() -> Self {
        Self {
            x: r32::new(0, 1),
            y: r32::new(0, 1),
        }
    }

    pub fn scalar_mult(lambda: r32, v: Point) -> Self {
        Self {
            x: lambda * v.x,
            y: lambda * v.y,
        }
    }

    pub fn average(points: &Vec<Self>) -> Self {
        let mut sum = Self::zero();
        for p in points {
            sum = sum + *p;
        }
        Self::scalar_mult(r32::new(1, 1) / r32::new(points.len() as i32, 1), sum)
    }

    pub fn theta(v: Point) -> f64 {
        (*v.y.numer() as f64 / *v.y.denom() as f64).atan2(*v.x.numer() as f64 / *v.x.denom() as f64)
    }

    /// Return a copy of `points` sorted using the direction they make from the centroid, according to theta(v-centroid) in (-pi, pi].
    fn ordered_ccw(mut points: Vec<Self>) -> Vec<Self> {
        let centroid = Self::average(&points);
        points.sort_by(|a, b| {
            Self::theta(*a - centroid)
                .partial_cmp(&Self::theta(*b - centroid))
                .unwrap()
        });
        points
    }

    pub fn midpoint(p1: &Self, p2: &Self) -> Self {
        Self {
            x: (p1.x + p2.x) / r32::new(2, 1),
            y: (p1.y + p2.y) / r32::new(2, 1),
        }
    }
    // Use the slope-intercept formula including the case where the slope is infinite.
    pub fn are_collinear(p1: Self, p2: Self, p3: Self) -> bool {
        (p1.y - p3.y) * (p1.x - p2.x) == (p1.y - p2.y) * (p1.x - p3.x)
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

/// An oriented line in Q^2 of rational or infinite slope.
#[derive(Copy, Clone, Debug)]
pub struct Line {
    point1: Point,
    point2: Point,
}

impl PartialEq for Line {
    fn eq(&self, other: &Self) -> bool {
        self.slope() == other.slope() && self.has_point(other.point())
    }
}

impl Line {
    pub fn slope(self: Line) -> Slope {
        Slope::new(
            (*(self.point2.y - self.point1.y).numer()) * (*(self.point2.x - self.point1.x).denom()),
            (*(self.point2.x - self.point1.x).numer()) * (*(self.point2.y - self.point1.y).denom()),
        )
        .unwrap()
    }

    pub fn point(self: Line) -> Point {
        self.point1
    }

    pub fn from_slope_and_point(sl: Slope, point: Point) -> Self {
        Self {
            point1: point,
            point2: point + Point::new(r32::new_raw(sl.denom, 1), r32::new_raw(sl.numer, 1)),
        }
    }
    pub fn from_points(point1: Point, point2: Point) -> Result<Self, String> {
        if point1 == point2 {
            Err("The points {point1} and {point2} are equal".to_string())
        } else {
            Ok(Line { point1, point2 })
        }
    }

    /// Test whether a line is vertical.
    fn is_vertical(&self) -> bool {
        self.point1.x == self.point2.x
    }

    /// Test whether `p` satisfies the equation of the line.
    fn has_point(&self, p: Point) -> bool {
        (p.y - self.point().y) * r32::new_raw(self.slope().denom, 1)
            == r32::new_raw(self.slope().numer, 1) * (p.x - self.point().x)
    }

    /// Determine whether `p1` and `p2` are on the same side of the line.
    fn on_same_side(&self, p1: Point, p2: Point) -> bool {
        if self.is_vertical() {
            (p1.x < self.point1.x && p2.x < self.point1.x)
                || (p1.x > self.point1.x && p2.x > self.point1.x)
        } else {
            // both `p1` and `p2` are above `self`
            ((p1.y - self.point1.y)*r32::new_raw(self.slope().denom,1) > r32::new_raw(self.slope().numer,1)*(p1.x - self.point1.x) &&
                (p2.y - self.point1.y)*r32::new_raw(self.slope().denom,1) > r32::new_raw(self.slope().numer,1)*(p2.x - self.point1.x))
            || // ... or both are below `self`
            ((p1.y - self.point1.y)*r32::new_raw(self.slope().denom,1) < r32::new_raw(self.slope().numer,1)*(p1.x - self.point1.x) &&
                (p2.y - self.point1.y)*r32::new_raw(self.slope().denom,1) < r32::new_raw(self.slope().numer,1)*(p2.x - self.point1.x))
        }
    }

    /// Depending on the orientation of the line, determine whether the point `p` is to the left, to the right, or on the line.
    pub fn point_line_config(&self, p: Point) -> PointLineConfiguration {
        if self.point2.x != self.point1.x {
            // For non-vertical lines
            if ((p.y - self.point().y) * r32::new_raw(self.slope().denom, 1)
                - r32::new_raw(self.slope().numer, 1) * (p.x - self.point().x))
                * (self.point2.x - self.point1.x).signum()
                > r32::new_raw(0, 1)
            {
                PointLineConfiguration::Left
            } else if ((p.y - self.point().y) * r32::new_raw(self.slope().denom, 1)
                - r32::new_raw(self.slope().numer, 1) * (p.x - self.point().x))
                * (self.point2.x - self.point1.x).signum()
                == r32::new_raw(0, 1)
            {
                PointLineConfiguration::OnTheLine
            } else {
                PointLineConfiguration::Right
            }
        } else {
            // For vertical lines
            // If point2.y < point1.y, then the line is oriented downwards, so multiply by -1 to compensate.
            if (p.x - self.point().x) * (self.point2.y - self.point1.y) < r32::new_raw(0, 1) {
                // upwards vertical line and point to the left of it, or downwards vertical line and point to the right of it.
                PointLineConfiguration::Left
            } else if (p.x - self.point().x) * (self.point2.y - self.point1.y) == r32::new_raw(0, 1)
            {
                PointLineConfiguration::OnTheLine
            } else {
                PointLineConfiguration::Right
            }
        }
    }
}

impl fmt::Display for Line {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(slope: {}, point: {})", self.slope(), self.point1)
    }
}

/// For two non-parallel ones, find the intersection; for two parallel lines, return `None`.
pub fn intersection(l1: Line, l2: Line) -> Option<Point> {
    let (x1, y1, x2, y2) = (l1.point().x, l1.point().y, l2.point().x, l2.point().y);
    match (l1.slope(), l2.slope()) {
        (s1, s2) if s1 == s2 => None, // Same slope
        (s1, s2) if s1 == Slope::infinity() && s2 != Slope::infinity() => Some(Point {
            x: x1,
            y: s2.as_rational_raw() * (x1 - x2) + y2,
        }),
        (s1, s2) if s1 != Slope::infinity() && s2 == Slope::infinity() => Some(Point {
            x: x2,
            y: s1.as_rational_raw() * (x2 - x1) + y1,
        }),
        (s1, s2) => Some(Point {
            x: (y1 - s1.as_rational_raw() * x1 + s2.as_rational_raw() * x2 - y2)
                / (s2.as_rational_raw() - s1.as_rational_raw()),
            y: (s1.as_rational_raw() * s2.as_rational_raw() * (x2 - x1)
                + s2.as_rational_raw() * y1
                - s1.as_rational_raw() * y2)
                / (s2.as_rational_raw() - s1.as_rational_raw()),
        }),
    }
}

#[derive(Clone, Debug)]
pub struct ConvexPolygon {
    vertices: Vec<Point>, // Assumes that the vertices are ordered in clockwise or counterclockwise order.
}

impl PartialEq for ConvexPolygon {
    fn eq(&self, other: &Self) -> bool {
        // If all the points in both polygons are equal then they are equal; we ensure this by enforcing a standard order when a `ConvexPolygon` is created.
        if other.vertices.len() != self.vertices.len() {
            return false;
        }
        for (a, b) in zip(self.vertices.iter(), other.vertices.iter()) {
            if a != b {
                return false;
            }
        }
        true
    }
}

impl ConvexPolygon {
    // Assumes that the vertices form a convex polygon. Does not check if the `Vec` of vertices consists of distinct elements.
    pub fn new(vertices: Vec<Point>) -> Option<Self> {
        if vertices.len() >= 3 {
            Some(Self {
                vertices: Point::ordered_ccw(vertices), // This will work for vertices of a convex polygon.
            })
        } else {
            None
        }
    }

    pub fn centroid(&self) -> Point {
        Point::average(&self.vertices)
    }

    pub fn has_point_inside(&self, point: Point) -> bool {
        // Do `point` and edge_i+1 x edge_i+2 lie on the same side of edge_i?
        for i in 0..self.vertices.len() - 2 {
            if !Line::from_points(self.vertices[i], self.vertices[i + 1])
                .unwrap()
                .on_same_side(point, self.vertices[i + 2])
            {
                return false;
            }
        }
        // Check last two triples.
        if !Line::from_points(
            self.vertices[self.vertices.len() - 2],
            self.vertices[self.vertices.len() - 1],
        )
        .unwrap()
        .on_same_side(point, self.vertices[0])
        {
            false
        } else {
            Line::from_points(self.vertices[self.vertices.len() - 1], self.vertices[0])
                .unwrap()
                .on_same_side(point, self.vertices[1])
        }
    }

    /// A vector storing the edges of the polygon in order.
    fn lines_for_edges(&self) -> Vec<Line> {
        let mut result: Vec<Line> = (0..self.vertices.len() - 1)
            .map(|i| Line::from_points(self.vertices[i], self.vertices[i + 1]).unwrap())
            .collect();
        result.push(
            Line::from_points(self.vertices[self.vertices.len() - 1], self.vertices[0]).unwrap(),
        );
        result
    }

    /// Subdivide a polygon with the given line and return the two pieces.
    pub fn subdivide(&self, line: Line) -> (Option<ConvexPolygon>, Option<ConvexPolygon>) {
        // Collect vertices that fall to the left, on thhe line, and to the right, respectively. The left and right collections are missing two vertices.
        let mut vertices_left: Vec<Point> = self
            .vertices
            .iter()
            .filter(|&x| line.point_line_config(*x) == PointLineConfiguration::Left)
            .cloned()
            .collect();
        let mut vertices_right: Vec<Point> = self
            .vertices
            .iter()
            .filter(|&x| line.point_line_config(*x) == PointLineConfiguration::Right)
            .cloned()
            .collect();
        let vertices_on_line: Vec<Point> = self
            .vertices
            .iter()
            .filter(|&x| line.point_line_config(*x) == PointLineConfiguration::OnTheLine)
            .cloned()
            .collect();
        if !vertices_left.is_empty() || !vertices_right.is_empty() {
            // Add any vertices on the line both to the new left polygon and the new right polygon.
            for point in &vertices_on_line {
                vertices_left.push(*point);
                vertices_right.push(*point);
            }
            // Look for intersection points on edges with the line.
            for edge in &self.lines_for_edges() {
                if let Some(intsn) = intersection(*edge, line)
                    && is_between_r32(intsn.x, edge.point1.x, edge.point2.x)
                    && is_between_r32(intsn.y, edge.point1.y, edge.point2.y)
                {
                    vertices_left.push(intsn);
                    vertices_right.push(intsn);
                }
            }
            // At most two vertices are added, so if there was no vertex to one side originally, that side doesn't have a polygon.
        }
        (
            ConvexPolygon::new(vertices_left),
            ConvexPolygon::new(vertices_right),
        )
    }
}
