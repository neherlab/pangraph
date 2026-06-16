use crate::pangraph::strand::Strand;
use crate::utils::interval::Interval;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

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

/// Keep only the features whose `feature_type` passes the type filters.
///
/// When `only` is non-empty, retain only features whose type is in it (whitelist); independently,
/// drop any feature whose type is in `exclude` (blacklist). Empty slices are no-ops, so the default
/// (no flags) passes everything through. Matching is exact (case-sensitive). The two filters are
/// mutually exclusive at the CLI, but applying both here is order-independent and safe.
pub fn filter_features_by_type(features: Vec<Feature>, only: &[String], exclude: &[String]) -> Vec<Feature> {
  let only: BTreeSet<&str> = only.iter().map(String::as_str).collect();
  let exclude: BTreeSet<&str> = exclude.iter().map(String::as_str).collect();
  features
    .into_iter()
    .filter(|f| {
      (only.is_empty() || only.contains(f.feature_type.as_str())) && !exclude.contains(f.feature_type.as_str())
    })
    .collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  /// Build a minimal feature carrying only the `feature_type` the type-filter tests care about.
  fn typed(feature_type: &str) -> Feature {
    Feature {
      seqid: "chr1".to_owned(),
      source: None,
      feature_type: feature_type.to_owned(),
      interval: Interval::new(0, 1),
      strand: Some(Strand::Forward),
      id: None,
      name: None,
      attributes: vec![],
    }
  }

  /// Collect just the `feature_type`s after filtering, for compact assertions.
  fn types_after(only: &[&str], exclude: &[&str]) -> Vec<String> {
    let features = vec![typed("gene"), typed("CDS"), typed("region")];
    let only: Vec<String> = only.iter().map(|s| (*s).to_owned()).collect();
    let exclude: Vec<String> = exclude.iter().map(|s| (*s).to_owned()).collect();
    filter_features_by_type(features, &only, &exclude)
      .into_iter()
      .map(|f| f.feature_type)
      .collect()
  }

  #[test]
  fn filter_no_filters_passes_everything() {
    assert_eq!(types_after(&[], &[]), vec!["gene", "CDS", "region"]);
  }

  #[test]
  fn filter_only_keeps_whitelisted_types() {
    assert_eq!(types_after(&["gene"], &[]), vec!["gene"]);
    assert_eq!(types_after(&["gene", "CDS"], &[]), vec!["gene", "CDS"]);
  }

  #[test]
  fn filter_exclude_drops_blacklisted_types() {
    assert_eq!(types_after(&[], &["region"]), vec!["gene", "CDS"]);
  }

  #[test]
  fn filter_is_case_sensitive() {
    // `cds` does not match `CDS`, so the whitelist keeps nothing.
    assert!(types_after(&["cds"], &[]).is_empty());
  }

  #[test]
  fn filter_unknown_type_is_a_silent_no_op() {
    // Whitelisting a type that no feature has yields an empty result rather than an error.
    assert!(types_after(&["mRNA"], &[]).is_empty());
    // Excluding a type that no feature has leaves everything in place.
    assert_eq!(types_after(&[], &["mRNA"]), vec!["gene", "CDS", "region"]);
  }

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
