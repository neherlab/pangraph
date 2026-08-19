use crate::pangraph::strand::Strand;
use crate::utils::interval::Interval;
use serde::{Deserialize, Serialize};

/// A genome annotation feature, normalized from an annotation file (currently GFF)
/// into a single, format-agnostic representation.
///
/// The model is deliberately format-neutral: additional readers (e.g. GenBank) can be
/// added later without touching anything downstream. Coordinates live in `interval` as
/// 0-based, half-open `[start, end)` over the genome/contig identified by `seqid`.
/// Everything downstream (matching, the annotation lift, writers) consumes `Feature`
/// and never the original file format.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Feature {
  /// Identifier of the genome/contig the feature belongs to (GFF `seqid` column).
  /// Matched against pangraph path names.
  pub seqid: String,

  /// Free-text provenance of the annotation (GFF `source` column), if any.
  pub source: Option<String>,

  /// Feature type, e.g. `"CDS"` or `"gene"` (GFF `type` column).
  pub feature_type: String,

  /// Location on `seqid`, as a 0-based half-open interval `[start, end)`.
  pub interval: Interval,

  /// Strand of the feature, or `None` when unstranded (e.g. GFF `.` or `?`).
  pub strand: Option<Strand>,

  /// Stable feature identifier (GFF `ID` attribute), if present.
  pub id: Option<String>,

  /// Human-readable name (GFF `Name` attribute), if present.
  pub name: Option<String>,

  /// All remaining key→value metadata, preserving order and duplicate keys.
  pub attributes: Vec<(String, String)>,
}

/// Build a 0-based half-open [`Interval`] from 1-based, fully-closed coordinates,
/// the convention used by GFF.
///
/// For example `(1, 3)` (three bases, 1-based inclusive) becomes `[0, 3)`.
pub fn interval_from_one_based_inclusive(start: usize, end: usize) -> Interval {
  debug_assert!(start >= 1, "1-based coordinates must be >= 1, got start={start}");
  Interval::new(start - 1, end)
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_interval_from_one_based_inclusive_multibase() {
    assert_eq!(interval_from_one_based_inclusive(1, 3), Interval::new(0, 3));
  }

  #[test]
  fn test_interval_from_one_based_inclusive_single_base() {
    assert_eq!(interval_from_one_based_inclusive(5, 5), Interval::new(4, 5));
  }

  #[test]
  fn test_interval_from_one_based_inclusive_offset() {
    assert_eq!(interval_from_one_based_inclusive(10, 20), Interval::new(9, 20));
  }
}
