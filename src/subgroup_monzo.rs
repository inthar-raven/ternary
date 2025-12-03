use crate::monzo::Monzo;

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
    // Use `js_namespace` here to bind `console.log(..)` instead of just
    // `log(..)`
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}
#[allow(unused_macros)]
macro_rules! console_log {
    // Note that this is using the `log` function imported above during
    // `bare_bones`
    ($($t:tt)*) => (log(&format_args!($($t)*).to_string()))
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
/// Error types for invalid `SubgroupMonzo` construction
pub enum SubgroupMonzoErr {
    /// Dimension mismatch (need at least 3 elements)
    BadDim,
    /// Basis vectors are not linearly independent
    BadBasis,
}

#[derive(Debug, PartialEq)]
pub struct SubgroupMonzo(((Monzo, i32), (Monzo, i32), (Monzo, i32))); // Using just 3-member vectors for now

impl SubgroupMonzo {
    /// Tries to create a new `SubgroupMonzo` with length 3 enforced.
    /// Validates that basis vectors are linearly independent.
    pub fn try_new(basis: &[Monzo], vector: &[i32]) -> Result<Self, SubgroupMonzoErr> {
        if basis.len() < 3 || vector.len() < 3 {
            console_log!("SubgroupMonzoErr::BadDim");
            Err(SubgroupMonzoErr::BadDim)
        } else {
            // Simple linear independence check: check if the three vectors are not all parallel
            // (more sophisticated checks could use Gram-Schmidt or determinant, but this is sufficient)
            let v1 = basis[0].into_inner();
            let v2 = basis[1].into_inner();
            let v3 = basis[2].into_inner();

            // Check they're not collinear by computing rank using Gram-Schmidt
            let v2_proj = {
                let dot_v1_v2 = v1.dot(&v2);
                let dot_v1_v1 = v1.dot(&v1);
                if dot_v1_v1 == 0 {
                    v2
                } else {
                    let scalar = (dot_v1_v2 as f64) / (dot_v1_v1 as f64);
                    v2 - (v1 * scalar as i32)
                }
            };

            let dot_v2proj_v2proj = v2_proj.dot(&v2_proj);
            if dot_v2proj_v2proj == 0 {
                console_log!("SubgroupMonzoErr::BadBasis");
                Err(SubgroupMonzoErr::BadBasis)
            } else {
                let v3_proj = {
                    let dot_v1_v3 = v1.dot(&v3);
                    let dot_v1_v1 = v1.dot(&v1);
                    let dot_v2proj_v3 = v2_proj.dot(&v3);
                    let scalar1 = if dot_v1_v1 == 0 {
                        0.0
                    } else {
                        (dot_v1_v3 as f64) / (dot_v1_v1 as f64)
                    };
                    let scalar2 = if dot_v2proj_v2proj == 0 {
                        0.0
                    } else {
                        (dot_v2proj_v3 as f64) / (dot_v2proj_v2proj as f64)
                    };
                    v3 - (v1 * scalar1 as i32) - (v2_proj * scalar2 as i32)
                };

                if v3_proj.dot(&v3_proj) == 0 {
                    console_log!("SubgroupMonzoErr::BadBasis");
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
    }
    /// Gets a tuple of the three coefficients from a `SubgroupMonzo`.
    pub fn extract_vector(&self) -> (i32, i32, i32) {
        (self.0.0.1, self.0.1.1, self.0.2.1)
    }
    /// Convert the `SubgroupMonzo` to a 79-limit `Monzo` representation
    /// by combining the basis vectors with the subgroup vector.
    pub fn to_monzo(&self) -> Monzo {
        let (v1, v2, v3) = (self.0.0.0, self.0.1.0, self.0.2.0);
        let (a1, a2, a3) = (self.0.0.1, self.0.1.1, self.0.2.1);
        v1 * a1 + v2 * a2 + v3 * a3
    }
}
