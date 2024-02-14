use num_traits::{NumCast, ToPrimitive};
use std::fmt::{self, Display};
use std::ops::{Index, IndexMut};

/// Describes data layout of a single row in `Band2d`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Stripe {
  pub begin: usize,
  pub end: usize,
}

impl Stripe {
  #[inline]
  pub fn new<T, U>(begin: T, end: U) -> Stripe
  where
    T: NumCast + Copy,
    U: NumCast + Copy,
  {
    Stripe {
      begin: begin.to_usize().unwrap(),
      end: end.to_usize().unwrap(),
    }
  }

  #[inline]
  pub const fn len(&self) -> usize {
    self.end - self.begin
  }

  #[inline]
  pub const fn is_empty(&self) -> bool {
    self.end == self.begin
  }
}

pub fn simple_stripes(mean_shift: i32, band_width: usize, ref_len: usize, qry_len: usize) -> Vec<Stripe> {
  // Begin runs diagonally, with max(0, mean_shift - band_width + i)
  // End runs diagonally, with min(qry_len, mean_shift + band_width + i)

  // TODO: Increase by start/end bands
  let mut stripes = Vec::<Stripe>::with_capacity(ref_len + 1);
  let band_width_i32 = band_width.to_i32().unwrap();
  let ref_len_i32 = ref_len.to_i32().unwrap();
  let qry_len_i32 = qry_len.to_i32().unwrap();
  for i in 0..=ref_len_i32 {
    let begin = num::clamp(-mean_shift - band_width_i32 + i, 0, qry_len_i32);
    let end = num::clamp(-mean_shift + band_width_i32 + i + 1, 1, qry_len_i32 + 1);
    stripes.push(Stripe::new(begin, end));
  }
  // Make sure first and last stripe can reach origin/end
  stripes[0].begin = 0;
  stripes[ref_len].end = qry_len + 1;
  stripes
}

pub fn full_matrix(ref_len: usize, qry_len: usize) -> Vec<Stripe> {
  let mut stripes = Vec::<Stripe>::with_capacity(ref_len + 1);
  let ref_len_i32 = ref_len.to_i32().unwrap();
  let qry_len_i32 = qry_len.to_i32().unwrap();
  for i in 0..=ref_len_i32 {
    stripes.push(Stripe::new(0, qry_len + 1));
  }
  stripes
}

/// Represents a diagonal band in a matrix.
///
/// The underlying storage is sparse - the row storage consists of `Stripe`s, each of a given size (`stripe.length`)
/// and shifted by a given amount (`stripe.begin`) relative to the left boundary of the matrix. In each row, the cells
/// which are outside of the corresponding stripe are not allocated and accessing them is illegal.
///
/// Stripe begins must increase monotonically
#[derive(Clone, PartialEq, Eq)]
pub struct Band2d<T>
where
  T: Default + Clone,
{
  data: Vec<T>,
  stripes: Vec<Stripe>,
  row_start_points: Vec<usize>,
  n_rows: usize,
  n_cols: usize,
}

impl<T> Band2d<T>
where
  T: Default + Clone,
{
  pub fn new(stripes: &[Stripe]) -> Self {
    let (n_rows, _, row_start_points) = calculate_dimensions(stripes);
    let data: Vec<T> = vec![T::default(); row_start_points[n_rows]];
    Self::with_data(stripes, &data)
  }

  pub fn with_data(stripes: &[Stripe], data: &[T]) -> Self {
    let (n_rows, n_cols, row_start_points) = calculate_dimensions(stripes);
    Self {
      data: data.to_vec(),
      stripes: stripes.to_vec(),
      row_start_points,
      n_rows,
      n_cols,
    }
  }

  #[inline]
  pub const fn num_rows(&self) -> usize {
    self.n_rows
  }

  #[inline]
  pub const fn num_cols(&self) -> usize {
    self.n_cols
  }

  #[inline]
  pub fn data_len(&self) -> usize {
    self.data.len()
  }

  #[inline]
  fn get_index<I: NumCast + Copy, J: NumCast + Copy>(&self, index2d: (I, J)) -> usize {
    let row = index2d.0.to_usize().unwrap();
    let col = index2d.1.to_usize().unwrap();
    let stripe = &self.stripes[row];
    assert!(
      stripe.begin <= col && col < stripe.end,
      "Stripe col out of bounds: stripe.begin = {}, stripe.end = {}, col = {}",
      stripe.begin,
      stripe.end,
      col
    );
    self.row_start_points[row] + (col - stripe.begin)
  }
}

fn calculate_dimensions(stripes: &[Stripe]) -> (usize, usize, Vec<usize>) {
  let n_rows = stripes.len();
  let mut row_start_points = vec![0_usize; n_rows + 1];
  let mut n_cols = 0_usize;
  row_start_points[0] = 0;
  for (i, stripe) in stripes.iter().enumerate() {
    row_start_points[i + 1] = row_start_points[i] + stripe.len();
    n_cols = n_cols.max(stripe.end);
  }
  (n_rows, n_cols, row_start_points)
}

/// Allows 2-dimensional indexing using a tuple
impl<T: Default + Clone, I: NumCast + Copy, J: NumCast + Copy> Index<(I, J)> for Band2d<T> {
  type Output = T;

  #[inline]
  fn index(&self, index2d: (I, J)) -> &Self::Output {
    self.data.index(self.get_index(index2d))
  }
}

/// Allows 2-dimensional mutable indexing using a tuple
impl<T: Default + Clone, I: NumCast + Copy, J: NumCast + Copy> IndexMut<(I, J)> for Band2d<T> {
  #[inline]
  fn index_mut(&mut self, index2d: (I, J)) -> &mut Self::Output {
    self.data.index_mut(self.get_index(index2d))
  }
}

impl<T: Default + Clone + Display> fmt::Debug for Band2d<T> {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    for i in 0..self.n_rows {
      for j in 0..self.n_cols {
        if self.stripes[i].begin <= j && j < self.stripes[i].end {
          write!(f, "{:2} ", self[(i, j)])?;
        } else {
          write!(f, " - ")?;
        }
      }
      writeln!(f)?;
    }
    Ok(())
  }
}
