use std::cmp::max;

use crate::GuideResult;
use crate::countvector_to_u16_vec;
use crate::guide::*;
use crate::guide_frame_to_result;
use crate::matrix;
use crate::matrix::unimodular_inv;
use crate::word_to_sig;
use crate::words::CountVector;

#[allow(dead_code)]
/// A struct representing a sub-traversal of a full row-by-row traversal of a lattice parallelogram.
#[derive(Clone, Debug)]
pub struct QuasiParallelogram {
    row_count: i32,
    full_row_len: i32,
    first_row_len: i32,
    last_row_len: i32,
}

impl QuasiParallelogram {
    fn new(row_count: i32, full_row_len: i32, first_row_len: i32, last_row_len: i32) -> Self {
        QuasiParallelogram {
            row_count,
            full_row_len,
            first_row_len,
            last_row_len,
        }
    }
}

pub fn get_unimodular_basis(
    structures: &[GuideFrame],
    step_sig: &[u16],
) -> Option<(Vec<Vec<u16>>, GuideResult)> {
    for structure in structures {
        if structure.multiplicity() == 1 {
            let result = guide_frame_to_result(structure);
            let gs = &result.gs;
            for i in 0..gs.len() {
                for j in (i + 1)..gs.len() {
                    if matrix::det3(step_sig, &gs[i], &gs[j]).abs() == 1 {
                        return Some((vec![gs[i].clone(), gs[j].clone()], result));
                    }
                }
            }
            for v in &result.offset_chord {
                for w in gs {
                    if matrix::det3(step_sig, v, w).abs() == 1 {
                        return Some((vec![v.clone(), w.clone()], result));
                    }
                }
            }
        } else {
            // this branch handles multiplicity > 1 scales
            let result = guide_frame_to_result(structure);
            let gs = &result.gs;
            // Check all pairs from gs and offset_chord
            for v in &result.offset_chord {
                for w in gs {
                    if matrix::det3(step_sig, w, v).abs() == 1 {
                        return Some((vec![w.clone(), v.clone()], result));
                    }
                }
            }
        }
    }
    None
}

// Try to get a lattice of pitch classes in the scale
// using the unimodular basis found with guide_frames.
pub fn try_pitch_class_lattice(query: &[usize]) -> Option<Vec<Vec<i32>>> {
    let gfs = guide_frames(query);
    let sig = word_to_sig(query)
        .iter()
        .map(|x| *x as u16)
        .collect::<Vec<u16>>();
    get_unimodular_basis(&gfs, &sig).map(|(basis_, _)| {
        // For every pitch in the scale expressed as a CountVector,
        // do a change of basis from scale steps basis
        // to (equave, generator1, generator2) basis.
        // basis_ doesn't have the equave, so we add it back first.
        let (equave, gener_1, gener_2) = (
            sig.iter().map(|x| *x as i32).collect::<Vec<_>>(),
            basis_[0].iter().map(|x| *x as i32).collect::<Vec<_>>(),
            basis_[1].iter().map(|x| *x as i32).collect::<Vec<_>>(),
        );
        // Invert [equave, gener_1, gener_2] to get basis change matrix
        let basis_change = unimodular_inv(&equave, &gener_1, &gener_2);
        // Now project all vectors to the (generator1, generator2)-plane
        // Remove the first coordinate
        let mut pitch_classes = vec![];
        let mut count_vector = CountVector::ZERO;
        for step in query {
            // Add the current step
            count_vector = count_vector.add(&CountVector::from_slice(&[*step]));
            let v_i = countvector_to_u16_vec(&count_vector)
                .iter()
                .map(|x| *x as i32)
                .collect::<Vec<_>>();
            let v_i_transformed = matrix::matrix_times_vector(
                &basis_change[0],
                &basis_change[1],
                &basis_change[2],
                &v_i,
            );
            let v_i_projected = v_i_transformed[1..3].to_vec();
            pitch_classes.push(v_i_projected);
        }
        pitch_classes
    })
}

/// Whether the result is `Some` or `None` depends on
/// whether the ternary scale `query` has the "quasi-parallelogram property",
/// i.e. its pitch classes forming a substring of a traversal
/// ```text
/// (0,0), (0,1), (0,2), ..., (0,n),
/// (1,0), (1,1), (1,2), ..., (1,n),
/// (2,0), (2,1), (2,2), ..., (2,n),
/// ...
/// (m,0), (m,1), (m,2), ..., (m,n)
/// ```
/// under some choice of coordinate vectors (v, w).
/// If `Some`, returns a `QuasiParallelogram` struct containing the following info:
/// row count, length of a full row, length of first row, length of last row.
pub fn quasi_parallelogram_info(query: &[usize]) -> Option<QuasiParallelogram> {
    // If no unimodular basis, return false
    try_pitch_class_lattice(query).and_then(|pitch_classes| {
        // Get all pairwise differences between distinct points
        let mut pairwise_differences: Vec<Vec<i32>> = vec![];
        for i in 0..query.len() {
            for j in i + 1..query.len() {
                let diff = vec![
                    pitch_classes[j][0] - pitch_classes[i][0],
                    pitch_classes[j][1] - pitch_classes[i][1],
                ];
                pairwise_differences.push(diff);
            }
        }
        // Look for a unimodular basis that witnesses the quasi-parallelogram condition
        for (i, v1) in pairwise_differences.iter().enumerate() {
            for v2 in pairwise_differences.iter().skip(i + 1) {
                if (v1[0] * v2[1] - v1[1] * v2[0]).abs() == 1 {
                    // Change coordinates to basis (v1, v2)
                    let basis_change: Vec<Vec<i32>> =
                        vec![vec![v2[1], -v1[1]], vec![-v2[0], v1[0]]];
                    let mut pitch_classes_transformed = pitch_classes
                        .iter()
                        .map(|v| {
                            vec![
                                basis_change[0][0] * v[0] + basis_change[1][0] * v[1],
                                basis_change[0][1] * v[0] + basis_change[1][1] * v[1],
                            ]
                        })
                        .collect::<Vec<_>>();
                    // Get window dimensions: x_min, x_max, y_min, y_max
                    let mut xs: Vec<_> = pitch_classes_transformed.iter().map(|v| v[0]).collect();
                    let mut ys: Vec<_> = pitch_classes_transformed.iter().map(|v| v[1]).collect();
                    xs.sort();
                    ys.sort();
                    let x_min = xs[0];
                    let x_max = xs[xs.len() - 1];
                    let y_min = ys[0];
                    let y_max = ys[ys.len() - 1];
                    // Check all 4 possible traversals:
                    // 1. each row LTR (increases in x), rows go BTT (increases in y) (equivalently each row RTL, rows go TTB)
                    // 2. each row RTL (decreases in x), rows go BTT (increases in y) (equivalently each row LTR, rows go TTB)
                    // 3. each row BTT (increases in y), rows go LTR (increases in x) (equivalently each row TTB, rows go RTL)
                    // 4. each row TTB (decreases in y), rows go LTR (increases in x) (equivalently each row BTT, rows go RTL)
                    'traversal12: {
                        // Sort pitch_classes_transformed in lex order for traversal 1
                        pitch_classes_transformed.sort_by(|v1, v2| {
                            // Sort by *ascending* y values, if y values are equal sort by *ascending* x values
                            v1[1].cmp(&v2[1]).then(v1[0].cmp(&v2[0]))
                        });
                        let mut index = 0; // index into pitch_classes_transformed
                        // Check if middle rows are fully occupied; if not break out of block early
                        let ys_middle = (y_min + 1)..=(y_max - 1);
                        let full_row_len = x_max - x_min + 1; // Required length of each middle row
                        for y in ys_middle {
                            let mut row_counter = 0; // Count pitches with this y value
                            while pitch_classes_transformed[index][1] < y {
                                index += 1;
                            }
                            while pitch_classes_transformed[index][1] == y {
                                row_counter += 1;
                                index += 1;
                            }
                            if row_counter != full_row_len {
                                break 'traversal12;
                            }
                        }
                        'traversal1: {
                            // Check outer rows for traversal 1
                            // Last row must be a prefix of a row traversal
                            let mut last_row = vec![];
                            while pitch_classes_transformed[index][1] < y_max {
                                index += 1;
                            }
                            while index < pitch_classes_transformed.len() {
                                last_row.push(pitch_classes_transformed[index].clone());
                                index += 1;
                            }
                            // First x value in last_row == x_min AND
                            // |last x value in row - first x value in row| + 1 == last_row.len()
                            // (last_row should be sorted by *ascending* x values)
                            let last_row_is_prefix = last_row[0][0] == x_min
                                && (last_row[last_row.len() - 1][0] - x_min + 1) as usize
                                    == last_row.len();
                            if !last_row_is_prefix {
                                break 'traversal1;
                            }
                            // First row must be a suffix of a row traversal
                            index = 0;
                            let mut first_row = vec![];
                            while pitch_classes_transformed[index][1] == y_min {
                                first_row.push(pitch_classes_transformed[index].clone());
                                index += 1;
                            }
                            // Last x value in first_row == x_max AND
                            // |last x value in row - first x value in row| + 1 == first_row.len()
                            // (first_row should be sorted by *ascending* x values)
                            let first_row_is_suffix = first_row[first_row.len() - 1][0] == x_max
                                && (x_max - first_row[0][0] + 1) as usize == first_row.len();
                            if first_row_is_suffix {
                                let row_count = y_max - y_min + 1;
                                let first_row_len = first_row.len() as i32;
                                let last_row_len = last_row.len() as i32;
                                return Some(QuasiParallelogram::new(
                                    row_count,
                                    if row_count == 2 {
                                        max(first_row_len, last_row_len)
                                    } else {
                                        full_row_len
                                    },
                                    first_row_len,
                                    last_row_len,
                                ));
                            }
                        }
                        // Sort pitch_classes_transformed in lex order for traversal 2
                        pitch_classes_transformed.sort_by(|v1, v2| {
                            // Sort by *ascending* y values, if y values are equal sort by *descending* x values
                            v1[1].cmp(&v2[1]).then(v1[0].cmp(&v2[0]).reverse())
                        });
                        'traversal2: {
                            // Check outer rows for traversal 2
                            // Last row must be a prefix of a row traversal
                            let mut last_row = vec![];
                            index = 0;
                            while pitch_classes_transformed[index][1] > y_min {
                                index += 1;
                            }
                            while index < pitch_classes_transformed.len() {
                                last_row.push(pitch_classes_transformed[index].clone());
                                index += 1;
                            }
                            // First x value in last_row == x_max AND
                            // |last x value in row - first x value in row| + 1 == last_row.len()
                            // (last_row should be sorted by *descending* x values)
                            let last_row_is_prefix = last_row[0][0] == x_max
                                && (x_max - last_row[last_row.len() - 1][0] + 1) as usize
                                    == last_row.len();
                            if !last_row_is_prefix {
                                break 'traversal2;
                            }
                            // First row must be a suffix of a row traversal
                            index = 0;
                            let mut first_row = vec![];
                            while pitch_classes_transformed[index][1] == y_max {
                                first_row.push(pitch_classes_transformed[index].clone());
                                index += 1;
                            }
                            // Last x value in first_row == x_min AND
                            // |last x value in row - first x value in row| + 1 == first_row.len()
                            // (first_row should be sorted by *descending* x values)
                            let first_row_is_suffix = first_row[first_row.len() - 1][0] == x_min
                                && (first_row[0][0] - x_min + 1) as usize == first_row.len();
                            if first_row_is_suffix {
                                let row_count = y_max - y_min + 1;
                                let first_row_len = first_row.len() as i32;
                                let last_row_len = last_row.len() as i32;
                                return Some(QuasiParallelogram::new(
                                    row_count,
                                    if row_count == 2 {
                                        max(first_row_len, last_row_len)
                                    } else {
                                        full_row_len
                                    },
                                    first_row_len,
                                    last_row_len,
                                ));
                            }
                        }
                    }
                    'traversal34: {
                        // Sort pitch_classes_transformed in lex order for traversal 3
                        pitch_classes_transformed.sort_by(|v1, v2| {
                            // Sort by *ascending* x values, if x values are equal sort by *ascending* y values
                            v1[0].cmp(&v2[0]).then(v1[1].cmp(&v2[1]))
                        });
                        let mut index = 0; // index into pitch_classes_transformed
                        // Check if middle rows are fully occupied; if not break out of block early
                        let xs_middle = (x_min + 1)..=(x_max - 1);
                        let full_row_len = y_max - y_min + 1; // Required length of each middle row
                        for x in xs_middle {
                            let mut row_counter = 0; // Count pitches with this y value
                            while pitch_classes_transformed[index][0] < x {
                                index += 1;
                            }
                            while pitch_classes_transformed[index][0] == x {
                                row_counter += 1;
                                index += 1;
                            }
                            if row_counter != full_row_len {
                                break 'traversal34;
                            }
                        }
                        'traversal3: {
                            // Check outer rows for traversal 3
                            // Last row must be a prefix of a row traversal
                            let mut last_row = vec![];
                            while pitch_classes_transformed[index][0] < x_max {
                                index += 1;
                            }
                            while index < pitch_classes_transformed.len() {
                                last_row.push(pitch_classes_transformed[index].clone());
                                index += 1;
                            }
                            // First y value in last_row == y_min AND
                            // |last y value in row - first y value in row| + 1 == last_row.len()
                            // (last_row should be sorted by *ascending* y values)
                            let last_row_is_prefix = last_row[0][1] == y_min
                                && (last_row[last_row.len() - 1][1] - y_min + 1) as usize
                                    == last_row.len();
                            if !last_row_is_prefix {
                                break 'traversal3;
                            }
                            // First row must be a suffix of a row traversal
                            index = 0;
                            let mut first_row = vec![];
                            while pitch_classes_transformed[index][0] == x_min {
                                first_row.push(pitch_classes_transformed[index].clone());
                                index += 1;
                            }
                            // Last y value in first_row == y_max AND
                            // |last y value in row - first y value in row| + 1 == first_row.len()
                            // (first_row should be sorted by *ascending* y values)
                            let first_row_is_suffix = first_row[first_row.len() - 1][1] == y_max
                                && (y_max - first_row[0][1] + 1) as usize == first_row.len();
                            if first_row_is_suffix {
                                let row_count = x_max - x_min + 1;
                                let first_row_len = first_row.len() as i32;
                                let last_row_len = last_row.len() as i32;
                                return Some(QuasiParallelogram::new(
                                    row_count,
                                    if row_count == 2 {
                                        max(first_row_len, last_row_len)
                                    } else {
                                        full_row_len
                                    },
                                    first_row_len,
                                    last_row_len,
                                ));
                            }
                        }
                        // Sort pitch_classes_transformed in lex order for traversal 4
                        pitch_classes_transformed.sort_by(|v1, v2| {
                            // Sort by *ascending* x values, if y values are equal sort by *descending* y values
                            v1[0].cmp(&v2[0]).then(v1[1].cmp(&v2[1]).reverse())
                        });
                        'traversal4: {
                            // Check outer rows for traversal 4
                            // Last row must be a prefix of a row traversal
                            let mut last_row = vec![];
                            index = 0;
                            while pitch_classes_transformed[index][0] > x_min {
                                index += 1;
                            }
                            while index < pitch_classes_transformed.len() {
                                last_row.push(pitch_classes_transformed[index].clone());
                                index += 1;
                            }
                            // First y value in last_row == y_max AND
                            // |last y value in row - first y value in row| + 1 == last_row.len()
                            // (last_row should be sorted by descending y values)
                            let last_row_is_prefix = last_row[0][1] == y_max
                                && (y_max - last_row[last_row.len() - 1][1] + 1) as usize
                                    == last_row.len();
                            if !last_row_is_prefix {
                                break 'traversal4;
                            }
                            // First row must be a suffix of a row traversal
                            index = 0;
                            let mut first_row = vec![];
                            while pitch_classes_transformed[index][1] == y_max {
                                first_row.push(pitch_classes_transformed[index].clone());
                                index += 1;
                            }
                            // Last y value in first_row == y_min AND
                            // |last y value in row - first y value in row| + 1 == first_row.len()
                            // (last_row should be sorted by descending x values)
                            let first_row_is_suffix = first_row[first_row.len() - 1][1] == y_min
                                && (first_row[0][1] - y_min + 1) as usize == first_row.len();
                            if first_row_is_suffix {
                                let row_count = x_max - x_min + 1;
                                let first_row_len = first_row.len() as i32;
                                let last_row_len = last_row.len() as i32;
                                return Some(QuasiParallelogram::new(
                                    row_count,
                                    if row_count == 2 {
                                        max(first_row_len, last_row_len)
                                    } else {
                                        full_row_len
                                    },
                                    first_row_len,
                                    last_row_len,
                                ));
                            }
                        }
                    }
                }
            }
        }
        None
    })
}

#[cfg(test)]
mod tests {
    use crate::{
        // helpers::gcd,
        lattice::quasi_parallelogram_info,
    };
    // use std::fs;
    #[test]
    fn test_quasi_parallelogram() {
        let diasem = [0, 1, 0, 2, 0, 1, 0, 2, 0];
        let blackdye = [0, 1, 0, 2, 0, 1, 0, 2, 0, 2];
        let diaslen_4sc = [2, 0, 1, 0, 2, 0, 2, 0, 1, 0, 2]; // sLmLsLsLmLs
        let fourth_example = [
            0, 0, 2, 1, 2, 0, 1, 2, 0, 2, 0, 1, 2, 0, 2, 1, 0, 2, 0, 2, 1, 0, 2, 1, 2,
        ]; // LLsmsLmsLsLmsLsmLsLsmLsms
        assert!(quasi_parallelogram_info(&diasem).is_some());
        assert!(quasi_parallelogram_info(&blackdye).is_some());
        assert!(quasi_parallelogram_info(&diaslen_4sc).is_some());
        assert!(quasi_parallelogram_info(&fourth_example).is_some());

        let nonexample = [0, 0, 0, 0, 2, 0, 1, 0, 0, 2, 0, 0, 0, 1, 2]; // LLLLsLmLLsLLLms
        let nonexample_2 = [0, 2, 0, 2, 0, 2, 1, 2, 0, 2, 0, 2, 1, 2, 2]; // LsLsLsmsLsLsmss
        assert!(quasi_parallelogram_info(&nonexample).is_none());
        assert!(quasi_parallelogram_info(&nonexample_2).is_none());
    }
    /*
    #[test]
    fn test_() {
        let file_path = "output.txt";
        let mut content = "Quasi-parallelogram MOS substitution scales\n".to_string();
        for scale_size in 7..=34 {
            content = format!("{content}=== {scale_size} notes ===\n");
            for a in 2..=(scale_size - 2) {
                for b in 1..=(scale_size - a) / 2 {
                    let c = scale_size - a - b;
                    if gcd(a as u64, gcd(b as u64, c as u64)) == 1 {
                        println!("Now iterating for [a, b, c] = [{a}, {b}, {c}]");
                        let mos_subst_scales =
                            crate::words::mos_substitution_scales_one_perm(a, b, c);
                        for scale in mos_subst_scales {
                            if let Some(QuasiParallelogram {
                                row_count,
                                full_row_len,
                                first_row_len,
                                last_row_len,
                            }) = quasi_parallelogram_info(&scale)
                            {
                                let scale_string = {
                                    let mut result = "".to_string();
                                    for step in scale {
                                        if step == 0 {
                                            result += "x";
                                        } else if step == 1 {
                                            result += "y";
                                        } else if step == 2 {
                                            result += "z";
                                        }
                                    }
                                    result
                                };
                                println!(
                                    "{scale_string} is QP\trow_count: {row_count}\tfull_row_len: {full_row_len}\tfirst_row_len: {first_row_len}\tlast_row_len: {last_row_len}"
                                );
                                content = format!(
                                    "{content}* {a}x({b}y{c}z)\t{scale_string}\trow_count: {row_count}\tfull_row_len: {full_row_len}\tfirst_row_len: {first_row_len}\tlast_row_len: {last_row_len}\n"
                                );
                            }
                        }
                    }
                }
            }
        }
        let _ = fs::write(file_path, content);
    }
    */
}
