// Replacement for nalgebra::SVector<i32, SMALL_PRIMES_COUNT> - a fixed-size vector wrapper
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, Index, Mul, Neg, Sub, SubAssign};

use crate::primes::SMALL_PRIMES_COUNT;

/// A fixed-size vector of `SMALL_PRIMES_COUNT` i32 elements, replacing nalgebra::SVector
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Vector([i32; SMALL_PRIMES_COUNT]);

impl Vector {
    /// Create from array
    pub const fn new(arr: [i32; SMALL_PRIMES_COUNT]) -> Self {
        Vector(arr)
    }

    /// Create from slice, padding with zeros
    pub fn from_slice(slice: &[i32]) -> Self {
        let mut arr = [0i32; SMALL_PRIMES_COUNT];
        for (i, &val) in slice.iter().take(SMALL_PRIMES_COUNT).enumerate() {
            arr[i] = val;
        }
        Vector(arr)
    }

    /// Get the underlying array
    pub fn as_array(&self) -> &[i32; SMALL_PRIMES_COUNT] {
        &self.0
    }

    /// Get mutable access to underlying array
    pub fn as_array_mut(&mut self) -> &mut [i32; SMALL_PRIMES_COUNT] {
        &mut self.0
    }

    /// Iterator over elements
    pub fn iter(&self) -> std::slice::Iter<'_, i32> {
        self.0.iter()
    }

    /// Mutable iterator
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, i32> {
        self.0.iter_mut()
    }

    /// Dot product with another vector
    pub fn dot(&self, other: &Vector) -> i32 {
        let mut sum = 0i32;
        for i in 0..SMALL_PRIMES_COUNT {
            sum += self.0[i] * other.0[i];
        }
        sum
    }

    /// Element-wise multiplication (Hadamard product)
    pub fn hadamard(&self, other: &Vector) -> Vector {
        let mut result = [0i32; SMALL_PRIMES_COUNT];
        for (i, item) in result.iter_mut().enumerate() {
            *item = self.0[i] * other.0[i];
        }
        Vector(result)
    }

    /// Sum of all elements
    pub fn sum(&self) -> i32 {
        self.0.iter().sum()
    }

    /// Check if all elements are zero
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|&x| x == 0)
    }
}

impl Index<usize> for Vector {
    type Output = i32;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Add for Vector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = self.0;
        for (i, item) in result.iter_mut().enumerate() {
            *item += rhs.0[i];
        }
        Vector(result)
    }
}

impl AddAssign for Vector {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..SMALL_PRIMES_COUNT {
            self.0[i] += rhs.0[i];
        }
    }
}

impl Sub for Vector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut result = self.0;
        for (i, item) in result.iter_mut().enumerate() {
            *item -= rhs.0[i];
        }
        Vector(result)
    }
}

impl SubAssign for Vector {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..SMALL_PRIMES_COUNT {
            self.0[i] -= rhs.0[i];
        }
    }
}

impl Neg for Vector {
    type Output = Self;
    fn neg(self) -> Self {
        let mut result = self.0;
        for item in &mut result {
            *item = -*item;
        }
        Vector(result)
    }
}

impl Mul<i32> for Vector {
    type Output = Self;
    fn mul(self, scalar: i32) -> Self {
        let mut result = self.0;
        for item in &mut result {
            *item *= scalar;
        }
        Vector(result)
    }
}

impl Mul<i32> for &Vector {
    type Output = Vector;
    fn mul(self, scalar: i32) -> Vector {
        (*self) * scalar
    }
}

impl Mul<Vector> for i32 {
    type Output = Vector;
    fn mul(self, vec: Vector) -> Vector {
        vec * self
    }
}

impl Mul<&Vector> for i32 {
    type Output = Vector;
    fn mul(self, vec: &Vector) -> Vector {
        (*vec) * self
    }
}

impl Div<i32> for Vector {
    type Output = Self;
    fn div(self, scalar: i32) -> Self {
        if scalar == 0 {
            panic!("Division by zero");
        }
        let mut result = self.0;
        for item in &mut result {
            *item /= scalar;
        }
        Vector(result)
    }
}

impl Sum for Vector {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Vector([0; SMALL_PRIMES_COUNT]), |a, b| a + b)
    }
}

impl IntoIterator for Vector {
    type Item = i32;
    type IntoIter = std::array::IntoIter<i32, SMALL_PRIMES_COUNT>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a Vector {
    type Item = &'a i32;
    type IntoIter = std::slice::Iter<'a, i32>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

/// A fixed-size vector of SMALL_PRIMES_COUNT f64 elements
#[derive(Debug, Clone, Copy)]
pub struct Vectorf64([f64; SMALL_PRIMES_COUNT]);

impl Vectorf64 {
    /// Create from array
    pub const fn new(arr: [f64; SMALL_PRIMES_COUNT]) -> Self {
        Vectorf64(arr)
    }

    /// Create from slice, padding with zeros
    pub fn from_slice(slice: &[f64]) -> Self {
        let mut arr = [0.0f64; SMALL_PRIMES_COUNT];
        for (i, &val) in slice.iter().take(SMALL_PRIMES_COUNT).enumerate() {
            arr[i] = val;
        }
        Vectorf64(arr)
    }

    /// Create from i32 iterator
    pub fn from_vec(vec: Vec<f64>) -> Self {
        let mut arr = [0.0f64; SMALL_PRIMES_COUNT];
        for (i, &val) in vec.iter().take(SMALL_PRIMES_COUNT).enumerate() {
            arr[i] = val;
        }
        Vectorf64(arr)
    }

    /// Get the underlying array
    pub fn as_array(&self) -> &[f64; SMALL_PRIMES_COUNT] {
        &self.0
    }

    /// Iterator over elements
    pub fn iter(&self) -> std::slice::Iter<'_, f64> {
        self.0.iter()
    }

    /// Dot product with another vector
    pub fn dot(&self, other: &Vectorf64) -> f64 {
        let mut sum = 0.0;
        for i in 0..SMALL_PRIMES_COUNT {
            sum += self.0[i] * other.0[i];
        }
        sum
    }

    /// Element-wise multiplication (Hadamard product)
    pub fn hadamard(&self, other: &Vectorf64) -> Vectorf64 {
        let mut result = [0.0f64; SMALL_PRIMES_COUNT];
        for (i, item) in result.iter_mut().enumerate() {
            *item = self.0[i] * other.0[i];
        }
        Vectorf64(result)
    }
}

impl Index<usize> for Vectorf64 {
    type Output = f64;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Add for Vectorf64 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = self.0;
        for (i, item) in result.iter_mut().enumerate() {
            *item += rhs.0[i];
        }
        Vectorf64(result)
    }
}

impl Mul<f64> for Vectorf64 {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        let mut result = self.0;
        for item in &mut result {
            *item *= scalar;
        }
        Vectorf64(result)
    }
}

impl IntoIterator for Vectorf64 {
    type Item = f64;
    type IntoIter = std::array::IntoIter<f64, SMALL_PRIMES_COUNT>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a Vectorf64 {
    type Item = &'a f64;
    type IntoIter = std::slice::Iter<'a, f64>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

/// A fixed-size row vector of SMALL_PRIMES_COUNT i32 elements, replacing nalgebra::RowSVector
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RowVector([i32; SMALL_PRIMES_COUNT]);

impl RowVector {
    /// Create from array
    pub const fn new(arr: [i32; SMALL_PRIMES_COUNT]) -> Self {
        RowVector(arr)
    }

    /// Create from slice, padding with zeros
    pub fn from_slice(slice: &[i32]) -> Self {
        let mut arr = [0i32; SMALL_PRIMES_COUNT];
        for (i, &val) in slice.iter().take(SMALL_PRIMES_COUNT).enumerate() {
            arr[i] = val;
        }
        RowVector(arr)
    }

    /// Get the underlying array
    pub fn as_array(&self) -> &[i32; SMALL_PRIMES_COUNT] {
        &self.0
    }

    /// Iterator over elements
    pub fn iter(&self) -> std::slice::Iter<'_, i32> {
        self.0.iter()
    }

    /// Dot product with Vector (row * column)
    pub fn dot(&self, other: &Vector) -> i32 {
        let mut sum = 0i32;
        for i in 0..SMALL_PRIMES_COUNT {
            sum += self.0[i] * other.0[i];
        }
        sum
    }
}

impl Index<usize> for RowVector {
    type Output = i32;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Add for RowVector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = self.0;
        for (i, item) in result.iter_mut().enumerate() {
            *item += rhs.0[i];
        }
        RowVector(result)
    }
}

impl AddAssign for RowVector {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..SMALL_PRIMES_COUNT {
            self.0[i] += rhs.0[i];
        }
    }
}

impl Sub for RowVector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut result = self.0;
        for (i, item) in result.iter_mut().enumerate() {
            *item -= rhs.0[i];
        }
        RowVector(result)
    }
}

impl Neg for RowVector {
    type Output = Self;
    fn neg(self) -> Self {
        let mut result = self.0;
        for item in &mut result {
            *item = -*item;
        }
        RowVector(result)
    }
}

impl Mul<i32> for RowVector {
    type Output = Self;
    fn mul(self, scalar: i32) -> Self {
        let mut result = self.0;
        for item in &mut result {
            *item *= scalar;
        }
        RowVector(result)
    }
}

impl Mul<i32> for &RowVector {
    type Output = RowVector;
    fn mul(self, scalar: i32) -> RowVector {
        (*self) * scalar
    }
}

impl Sum for RowVector {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(RowVector([0; SMALL_PRIMES_COUNT]), |a, b| a + b)
    }
}
