//! The `objects` modules contains types that are used for interacting
//! with the FHE-based OPRF that we construct
use tfhe::shortint::{Ciphertext, ClientKey};
use rand_core::{OsRng, RngCore};
use std::ops::Mul;

/// Trait for defining shared behaviour for `Matrix` types.
pub trait Matrix {
  fn width(&self) -> usize;
  fn height(&self) -> usize;
  fn cols(&self) -> Vec<Vec<u64>>;
}

/// Functions for sampling and operating over binary matrices.
#[derive(Clone, Debug)]
pub struct BinMatrix {
  cols: Vec<Vec<u64>>,
  height: usize,
  width: usize,
}
impl BinMatrix {
  pub fn new(cols: Vec<Vec<u64>>, height: usize, width: usize) -> Self {
    Self {
      cols,
      height,
      width,
    }
  }

  pub fn sample(height: usize, width: usize) -> Self {
    let cols = (0..width)
      .into_iter()
      .map(|_| sample_bin_vec_u64(height))
      .collect();
    Self {
      cols,
      height,
      width,
    }
  }
}
impl Matrix for BinMatrix {
  fn width(&self) -> usize {
    self.width
  }
  fn height(&self) -> usize {
    self.height
  }
  fn cols(&self) -> Vec<Vec<u64>> {
    self.cols.clone()
  }
}

/// Functions for sampling and operating over ternary matrices, with
/// entries in `{0,1,2}`.
#[derive(Debug, Clone)]
pub struct TernaryMatrix {
  cols: Vec<Vec<u64>>,
  height: usize,
  width: usize,
}
impl TernaryMatrix {
  pub fn sample(height: usize, width: usize) -> Self {
    let cols = (0..width)
      .into_iter()
      .map(|_| {
        let row = vec![0u64; height];
        row
          .iter()
          .map(|_| {
            let mut r = OsRng.next_u32();
            loop {
              // we need to perform rejection sampling as division rounds down
              let too_big = r > 3 * (u32::MAX / 3);
              if !too_big {
                break
              } 
              r = OsRng.next_u32();
            }
            (r % 3) as u64
          })
          .collect()
      })
      .collect();
    Self {
      cols,
      height,
      width,
    }
  }
}
impl Matrix for TernaryMatrix {
  fn width(&self) -> usize {
    self.width
  }
  fn height(&self) -> usize {
    self.height
  }
  fn cols(&self) -> Vec<Vec<u64>> {
    self.cols.clone()
  }
}

/// Functions for building and operating over vectors that can be
/// multiplied with matrices.
#[derive(Debug, Clone, PartialEq)]
pub struct Vector {
  value: Vec<u64>,
  width: usize,
}
impl Vector {
  pub fn encrypt(&self, client_key: &ClientKey) -> Vec<Ciphertext> {
    self
      .as_vec()
      .iter()
      .map(|&x| client_key.encrypt(x))
      .collect()
  }

  pub fn to_modulo(&self, modulus: u32) -> Vector {
    let mod2_vals: Vec<u64> = self
      .value
      .iter()
      .map(|v| modulo(*v, modulus) as u64)
      .collect();
    Vector::from(mod2_vals)
  }

  pub fn append_one(&mut self) {
    self.value.push(1);
    self.width += 1;
  }

  pub fn as_vec(&self) -> &Vec<u64> {
    &self.value
  }

  pub fn get(&self, index: usize) -> u64 {
    self.value[index]
  }

  pub fn set(&mut self, index: usize, value: u64) {
    self.value[index] = value;
  }
}
impl From<Vec<u64>> for Vector {
  fn from(v: Vec<u64>) -> Self {
    Self {
      value: v.clone(),
      width: v.len(),
    }
  }
}
impl From<Vec<u8>> for Vector {
  fn from(v: Vec<u8>) -> Self {
    let v64 = v.iter().map(|&x| x as u64).collect();
    Self {
      value: v64,
      width: v.len(),
    }
  }
}

impl<T> Mul<T> for Vector
where
  T: Matrix,
{
  type Output = Self;
  fn mul(self, m: T) -> Self {
    let h = m.height();
    let w = m.width();
    if h != self.width {
      panic!(
        "Width of vector ({}) does not equal height of matrix ({})",
        self.width, h
      );
    }
    Vector::from(
      (0..w)
        .into_iter()
        .map(|i| {
          let mut acc = 0;
          let row = &m.cols()[i];
          for j in 0..h {
            acc += self.get(j) * row[j];
          }
          acc
        })
        .collect::<Vec<u64>>(),
    )
  }
}

impl<T> Mul<&T> for &Vector
where
  T: Matrix,
{
  type Output = Vector;
  fn mul(self, m: &T) -> Vector {
    let h = m.height();
    let w = m.width();
    if h != self.width {
      panic!(
        "Width of vector ({}) does not equal height of matrix ({})",
        self.width, h
      );
    }
    Vector::from(
      (0..w)
        .into_iter()
        .map(|i| {
          let mut acc = 0;
          let row = &m.cols()[i];
          for j in 0..h {
            acc += self.get(j) * row[j];
          }
          acc
        })
        .collect::<Vec<u64>>(),
    )
  }
}

/// Sample a random vector of dimension `dim` where each element has a
/// 50% chance of being `1u64` or `0u64`.
fn sample_bin_vec_u64(dim: usize) -> Vec<u64> {
  let v = vec![0u64; dim];
  v.iter().map(|_| (OsRng.next_u32() % 2) as u64).collect()
}

/// Sample a random vector of dimension `dim` where each element has a
/// 50% chance of being `1u8` or `0u8`.
pub fn sample_bin_vec(dim: usize) -> Vec<u8> {
  let v = vec![0u8; dim];
  v.iter().map(|_| (OsRng.next_u32() % 2) as u8).collect()
}

// transforms a zero-centered `value` into a positive value modulo the `modulus`
pub fn modulo(val: u64, modulus: u32) -> u32 {
  let m = modulus as u64;
  (((val % m) + m) % m) as u32
}
