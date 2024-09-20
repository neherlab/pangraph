use crate::make_internal_report;
use eyre::Report;
use ndarray::{ArrayBase, ArrayView, Data, Dimension, IntoDimension};
use std::fmt::Debug;

/// Broadcasts array to given shape or fails with an error.
///
/// This is a thin wrapper around ndarray:ArrayBase::broadcast(), which returns an error instead of an option.
#[allow(clippy::needless_pass_by_value)]
#[inline]
pub fn broadcast<A, S, D1, D2, E2>(arr: &ArrayBase<S, D1>, dim: E2) -> Result<ArrayView<'_, A, E2::Dim>, Report>
where
  S: Data<Elem = A>,
  D1: Dimension,
  D2: Dimension,
  E2: IntoDimension<Dim = D2> + Debug + Clone,
{
  arr.broadcast(dim.clone()).ok_or_else(|| {
    make_internal_report!(
      "Unable to broadcast array of shape '{:#?}' to shape '{:#?}'",
      arr.shape(),
      &dim
    )
  })
}
