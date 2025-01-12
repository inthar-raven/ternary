use nalgebra::Const;
use nalgebra::Vector3;
use nalgebra::Matrix;

use crate::monzo::Monzo;

pub struct SubgroupMonzo(((Monzo, i32), (Monzo, i32), (Monzo, i32))); // Using just 3-member vectors for now

impl SubgroupMonzo {
    /// Tries to create a new `SubgroupMonzo` with length 3 enforced.
    pub fn try_new(basis: &[Monzo], vector: &[i32]) -> Result<Self, ()> { // placeholder Err type
        if basis.len() < 3 || vector.len() < 3 {
            Err(())
        } else {
            let real_basis: Vec<Matrix<_, _, _, _>> = basis.into_iter().map(|v| {
                let real_v: Vec<f64> = v.into_inner().iter().map(|x| *x as f64).collect();
                Matrix::<f64, Const<22>, Const<1>, _>::from_column_slice(&real_v)
            })
                .collect();
            // Check that `basis` is actually linearly independent
            let matrix = Matrix::<f64, Const<22>, Const<3>, _>::from_columns(&real_basis);
            if matrix.rank(f64::EPSILON) < 3 {
                Err(())
            } else {
                Ok(Self(((basis[0], vector[0]), (basis[1], vector[1]), (basis[2], vector[2]))))
            }
        }
    }
    /// Gets an `nalgebra::Vector` from a `SubgroupMonzo`.
    pub fn extract_vector(&self) -> Vector3<i32> {
        let (v1, v2, v3) = (self.0.0.1, self.0.1.1, self.0.2.1);
        Vector3::new(v1, v2, v3)
    }
}