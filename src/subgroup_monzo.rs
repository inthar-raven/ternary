use nalgebra::Const;
use nalgebra::Matrix;
use nalgebra::Vector3;

use crate::monzo::Monzo;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum SubgroupMonzoErr {
    BadDim,
    BadBasis,
}

#[derive(Debug, PartialEq)]
pub struct SubgroupMonzo(((Monzo, i32), (Monzo, i32), (Monzo, i32))); // Using just 3-member vectors for now

impl SubgroupMonzo {
    /// Tries to create a new `SubgroupMonzo` with length 3 enforced.
    #[allow(clippy::result_unit_err)]
    pub fn try_new(basis: &[Monzo], vector: &[i32]) -> Result<Self, SubgroupMonzoErr> {
        // placeholder Err type
        if basis.len() < 3 || vector.len() < 3 {
            Err(SubgroupMonzoErr::BadDim)
        } else {
            // convert entries to f64
            let f64_basis: Vec<Matrix<_, _, _, _>> = basis
                .iter()
                .map(|v| {
                    let real_v: Vec<f64> = v.into_inner().iter().map(|x| *x as f64).collect();
                    Matrix::<f64, Const<22>, Const<1>, _>::from_column_slice(&real_v)
                })
                .collect();
            // Check that `basis` is actually linearly independent
            let matrix = Matrix::<f64, Const<22>, Const<3>, _>::from_columns(&f64_basis);
            if matrix.rank(f64::EPSILON) < 3 {
                Err(SubgroupMonzoErr::BadBasis)
            } else {
                Ok(Self((
                    (basis[0], vector[0]),
                    (basis[1], vector[1]),
                    (basis[2], vector[2]),
                )))
            }
        }
    }
    /// Gets an `nalgebra::Vector` from a `SubgroupMonzo`.
    pub fn extract_vector(&self) -> Vector3<i32> {
        let (v1, v2, v3) = (self.0 .0 .1, self.0 .1 .1, self.0 .2 .1);
        Vector3::new(v1, v2, v3)
    }
    /// Convert the `SubgroupMonzo` to a 79-limit `Monzo` representation
    pub fn to_monzo(&self) -> Monzo {
        let (v1, v2, v3) = (self.0 .0 .0, self.0 .1 .0, self.0 .2 .0);
        let (a1, a2, a3) = (self.0 .0 .1, self.0 .1 .1, self.0 .2 .1);
        v1 * a1 + v2 * a2 + v3 * a3
    }
}
