use std::collections::BTreeSet;
use std::iter::Iterator;
use {std::cmp::Ord, std::cmp::Ordering};

#[derive(Default, Debug, PartialEq, Hash)]
/// "Top-level errors" for scale construction and analysis.
pub enum ScaleError {
    /// When the required generator step class for an analysis is not coprime with scale length
    NonCoprimeGenError,
    /// When a subset's size is required to divide the size of the whole scale but does not
    NonDivisibleSubsetError,
    /// When an offset fails to meet the interleavability condition, resulting in non-interleaved scales
    NotInterleavable,
    /// Default error value for generic scale construction failures
    #[default]
    CannotMakeScale,
}

impl std::fmt::Display for ScaleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonCoprimeGenError => {
                write!(f, "generator class used is not coprime")
            }
            &Self::NonDivisibleSubsetError => {
                write!(f, "chose a subset size that does not divide the scale size")
            }
            &Self::NotInterleavable => {
                write!(
                    f,
                    "the chosen interval does not interleave the chosen strand scale"
                )
            }
            &Self::CannotMakeScale => {
                write!(f, "erorr making scale")
            }
        }
    }
}
impl std::error::Error for ScaleError {}

/// The vector of all pairs with the first component from `v1` and second from `v2`.
/// Produces the Cartesian product v1 × v2.
pub fn pairs<T, U>(v1: &[T], v2: &[U]) -> Vec<(T, U)>
where
    T: Clone,
    U: Clone,
{
    v1.iter()
        .flat_map(|t| v2.iter().cloned().map(move |u| (t.clone(), u)))
        .collect()
}

/// The powerset of the given collection.
/// Returns all 2^n possible subsets (including empty set and full set).
pub fn powerset<T>(s: &[T]) -> Vec<Vec<T>>
where
    T: Clone,
{
    (0..2usize.pow(s.len() as u32))
        .map(|i| {
            s.iter()
                .enumerate()
                .filter(|&(t, _)| (i >> t) % 2 == 1)
                .map(|(_, element)| element.clone())
                .collect()
        })
        .collect()
}

/// Check if given subset is a maximal set in the collection (not properly contained in any other set).
pub fn is_maximal_in<T>(set: &[T], sets: &[Vec<T>]) -> bool
where
    T: Clone + Eq + Ord,
{
    let set_as_btree: BTreeSet<&T> = set.iter().collect();
    for set2 in sets {
        let set2_as_btree: BTreeSet<&T> = set2.iter().collect();
        if set_as_btree.is_subset(&set2_as_btree) && set.len() < set2.len() {
            return false;
        }
    }
    true
}

/// Given a `&[Vec]`, convert to borrowed slice references `Vec<&[T]>`.
pub fn slicify_each<T>(vecs: &[Vec<T>]) -> Vec<&[T]> {
    vecs.iter().map(|x| x.as_slice()).collect()
}

/// Given a `&[&[T]]`, convert to owned vectors `Vec<Vec<T>>`.
pub fn vectorize_each<T>(vecs: &[&[T]]) -> Vec<Vec<T>>
where
    T: Clone,
{
    (vecs.iter().map(|x| (*x).to_vec())).collect::<Vec<_>>()
}

/// Given a descending sorted vector, get the first index i where v[i] == t.
/// Returns None if not found or if earlier elements are smaller (due to sorting).
pub fn first_index_desc<T>(v: &[T], t: T) -> Option<usize>
where
    T: PartialEq + Ord,
{
    for (i, value) in v.iter().enumerate() {
        match (*value).cmp(&t) {
            Ordering::Equal => return Some(i),
            Ordering::Less => return None,
            _ => {}
        }
    }
    None
}

/// Get the first index i where v[i] < t. Returns None if no such index exists.
pub fn first_index_smaller<T>(v: &[T], t: T) -> Option<usize>
where
    T: PartialEq + Ord,
{
    for (i, value) in v.iter().enumerate() {
        if *value < t {
            return Some(i);
        }
    }
    None
}

/// Type alias for a binary relation (predicate function on pairs).
pub type Relation<T> = fn(&T, &T) -> bool;

/// Assuming that `equiv` is an equivalence relation,
/// return one representative from each equivalence class in the set under `equiv`.
pub fn equivalence_class_representatives<T>(set: &[T], equiv: Relation<T>) -> Vec<T>
where
    T: Clone + PartialEq,
{
    let mut result: Vec<T> = Vec::new();
    'b: for elem in set {
        let mut found_class_equivalent = false;
        'a: for item in &result {
            // All vecs in result are nonempty.
            if equiv(item, elem) {
                found_class_equivalent = true;
                continue 'a;
            }
        }
        if !found_class_equivalent {
            result.push(elem.to_owned());
            break 'b;
        }
    }
    result
}

/// Whether `a` is between `b` and `c` (both ends inclusive; works regardless of whether `b > c`, `b == c`, or `b < c`).
pub fn is_between<T>(a: T, b: T, c: T) -> bool
where
    T: Ord,
{
    match b.cmp(&c) {
        Ordering::Equal => a == b,
        Ordering::Less => b < a && a < c,
        Ordering::Greater => c < a && a < b,
    }
}

/// Compute `n % abs(m)` (modulo with absolute value).
fn modulo(n: i32, m: i32) -> i32 {
    ((n % m) + m) % m
}

/// Return modular inverse: find x such that (x*a) mod b = 1.
/// Returns error if a and b are not coprime.
pub fn modinv(a: i32, b: i32) -> Result<i32, String> {
    let (gcd, x, _) = extended_gcd(a, b);
    if gcd == 1 {
        Ok(modulo(x, b))
    } else {
        Err("Non-coprime generator error".to_string())
    }
}

/// Return gcd(a, b) and the Bézout coefficients: (g, x, y) such that ax + by = g.
/// Uses the extended Euclidean algorithm.
pub fn extended_gcd(a: i32, b: i32) -> (i32, i32, i32) {
    let (mut r, mut old_r) = (a, b);
    let (mut s, mut old_s) = (1, 0);
    let (mut t, mut old_t) = (0, 1);
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
        (old_t, t) = (t, old_t - quotient * t);
    }
    (old_r, old_s, old_t)
}

/// Gives the greatest common divisor of the two inputs.
/// Uses the binary GCD algorithm (Stein's algorithm) for efficiency.
pub fn gcd(mut u: u32, mut v: u32) -> u32 {
    // Base cases: gcd(n, 0) = gcd(0, n) = n
    if u == 0 {
        return v;
    } else if v == 0 {
        return u;
    }

    // gcd(2ⁱ u, 2ʲ v) = 2ᵏ gcd(u, v) with u, v odd and k = min(i, j)
    // 2ᵏ is the greatest power of two that divides both 2ⁱ u and 2ʲ v
    let i = u.trailing_zeros();
    u >>= i;
    let j = v.trailing_zeros();
    v >>= j;
    let k = std::cmp::min(i, j);

    loop {
        // u and v are odd at the start of the loop

        // Swap if necessary so u ≤ v
        if u > v {
            std::mem::swap(&mut u, &mut v);
        }

        // Identity 4: gcd(u, v) = gcd(u, v-u) as u ≤ v and u, v are both odd
        v -= u;
        // v is now even

        if v == 0 {
            // Identity 1: gcd(u, 0) = u
            // The shift by k is necessary to add back the 2ᵏ factor that was removed before the loop
            return u << k;
        }

        // Identity 3: gcd(u, 2ʲ v) = gcd(u, v) as u is odd
        v >>= v.trailing_zeros();
    }
}

/// Find the lcm (least common multiple) of two numbers.
pub fn lcm(m: u32, n: u32) -> u32 {
    m * n / gcd(m, n)
}

/// Finds a Bézout solution for a linear Diophantine equation with multiple coefficients.
/// Returns (gcd, solutions) where gcd is gcd(coeffs) and solutions satisfies:
/// `coeffs[0]*solutions[0] + coeffs[1]*solutions[1] + ... = gcd(coeffs)`
pub fn bezout(coeffs: &[i32]) -> (i32, Vec<i32>) {
    debug_assert!(!coeffs.is_empty(), "`coeffs` should not be empty");
    if coeffs.len() == 1 {
        (coeffs[0], vec![1])
    } else if coeffs.len() == 2 {
        let (d, x, y) = extended_gcd(coeffs[0], coeffs[1]);
        (d, vec![x, y])
    } else {
        let (d_, xs) = bezout(&coeffs[..coeffs.len() - 1]);
        let (d, x, t) = extended_gcd(*coeffs.last().expect("ok because coeffs.len() > 3"), d_);
        let mut result_vec = xs.into_iter().map(|x| t * x).collect::<Vec<i32>>();
        result_vec.push(x);
        (d, result_vec)
    }
}

/// Whether a slice is in strictly descending order.
pub fn is_sorted_strictly_desc<T>(ts: &[T]) -> bool
where
    T: Ord + Copy,
{
    match ts.len() {
        0 | 1 => true,
        _ => {
            let mut prev = ts[0];
            for t in &ts[1..] {
                if *t >= prev {
                    return false;
                } else {
                    prev = *t;
                }
            }
            true
        }
    }
}

#[cfg(test)]
mod tests {
    #[allow(unused)]
    use super::*;
    #[test]
    fn test_modinv() {
        assert_eq!(modinv(1, 2), Ok(1));
    }
    #[test]
    fn test_bezout() {
        assert_eq!(1, bezout(&[5, 2, 3]).0);
        assert_eq!(3, bezout(&[9, 6, 15]).0);
        assert_eq!(1, bezout(&[7, 5]).0);
    }
}
