/// Calculate the determinant of a 3x3 matrix.
/// Takes three column vectors and computes the determinant using the rule of Sarrus.
///
/// # Arguments
/// * `col0` - First column of the matrix (3 elements)
/// * `col1` - Second column of the matrix (3 elements)
/// * `col2` - Third column of the matrix (3 elements)
///
/// # Returns
/// The determinant as an i32
///
/// # Panics
/// May panic if any column has fewer than 3 elements
pub fn det3(col0: &[u16], col1: &[u16], col2: &[u16]) -> i32 {
    // Cast once per element to use WASM-native i32 arithmetic
    let (a0, a1, a2) = (col0[0] as i32, col0[1] as i32, col0[2] as i32);
    let (b0, b1, b2) = (col1[0] as i32, col1[1] as i32, col1[2] as i32);
    let (c0, c1, c2) = (col2[0] as i32, col2[1] as i32, col2[2] as i32);

    a0 * b1 * c2 + a1 * b2 * c0 + a2 * b0 * c1 - a2 * b1 * c0 - a1 * b0 * c2 - a0 * b2 * c1
}

/// Calculate the inverse of a unimodular 3x3 matrix.
/// # Arguments
/// * `col0` - First column of the matrix (3 elements)
/// * `col1` - Second column of the matrix (3 elements)
/// * `col2` - Third column of the matrix (3 elements)
///
/// # Returns
/// The inverse as a `Vec<Vec<i32>>` (column-major order)
pub fn unimodular_inv(col0: &[i32], col1: &[i32], col2: &[i32]) -> Vec<Vec<i32>> {
    let (a, b, c, d, e, f, g, h, i) = (
        col0[0], col1[0], col2[0], col0[1], col1[1], col2[1], col0[2], col1[2], col2[2],
    );
    let inv0 = vec![e * i - f * h, f * g - d * i, d * h - e * g];
    let inv1 = vec![c * h - b * i, a * i - c * g, b * g - a * h];
    let inv2 = vec![b * f - c * e, c * d - a * f, a * e - b * d];
    vec![inv0, inv1, inv2]
}

/// Calculate the product of a 3x3 matrix M and a vector v.
///
/// # Arguments
/// * `col0` - First column of M (3 elements)
/// * `col1` - Second column of M (3 elements)
/// * `col2` - Third column of M (3 elements)
/// * `v` - The vector premultiplied by M
///
/// # Returns
/// The vector Mv
///
/// # Panics
/// May panic if any vector has fewer than 3 elements
pub fn matrix_times_vector(col0: &[i32], col1: &[i32], col2: &[i32], v: &[i32]) -> Vec<i32> {
    let mv0 = col0[0] * v[0] + col1[0] * v[1] + col2[0] * v[2];
    let mv1 = col0[1] * v[0] + col1[1] * v[1] + col2[1] * v[2];
    let mv2 = col0[2] * v[0] + col1[2] * v[1] + col2[2] * v[2];
    vec![mv0, mv1, mv2]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_det3_identity() {
        // Identity matrix has determinant 1
        let col0 = [1, 0, 0];
        let col1 = [0, 1, 0];
        let col2 = [0, 0, 1];
        assert_eq!(det3(&col0, &col1, &col2), 1);
    }

    #[test]
    fn test_det3_zero() {
        // Matrix with two identical rows has determinant 0
        let col0 = [1, 2, 3];
        let col1 = [1, 2, 3];
        let col2 = [4, 5, 6];
        assert_eq!(det3(&col0, &col1, &col2), 0);
    }

    #[test]
    fn test_det3_unimodular() {
        // Test the unimodular basis from the 9L6m10s scale
        let equave = [9, 6, 10];
        let gener_1 = [3, 2, 3];
        let gener_2 = [2, 1, 2];
        let det = det3(&equave, &gener_1, &gener_2);
        assert_eq!(
            det.abs(),
            1,
            "Determinant should be ±1 for unimodular basis"
        );
    }

    #[test]
    fn test_det3_negative() {
        // Test a case with negative determinant
        let col0 = [1, 0, 0];
        let col1 = [0, 0, 1];
        let col2 = [0, 1, 0];
        assert_eq!(det3(&col0, &col1, &col2), -1);
    }

    #[test]
    fn test_matrix_inv() {
        let col0 = [1, 0, 0];
        let col1 = [2, 1, 0];
        let col2 = [3, 4, 1];
        assert_eq!(
            unimodular_inv(&col0, &col1, &col2),
            vec![vec![1, 0, 0], vec![-2, 1, 0], vec![5, -4, 1],]
        );
    }
}
