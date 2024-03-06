use crate::align::alignment::Alignment;
use clap::{Args, ValueHint};
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Clone, Debug, SmartDefault, Args, Serialize, Deserialize)]
pub struct SplitMatchesArgs {
  /// Minimum block size for alignment graph (in nucleotides)
  #[default = 100]
  #[clap(long = "len", short = 'l', default_value_t = SplitMatchesArgs::default().indel_len_threshold)]
  #[clap(value_hint = ValueHint::Other)]
  pub indel_len_threshold: usize, // Pass `PangraphBuildArgs::len` here
}

/// Split the alignments whenever an alignment contains an in/del longer than the threshold length.
pub fn split_matches(aln: &Alignment, args: &SplitMatchesArgs) -> Result<Vec<Alignment>, Report> {
  Ok(vec![])
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::align::bam::cigar::parse_cigar_str;
  use crate::o;
  use crate::pangraph::strand::Strand;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_split_matches_simple_case_forward() {
    let aln = Alignment {
      qry: Hit {
        name: o!("qry"),
        length: 500,
        start: 200,
        stop: 255,
      },
      reff: Hit {
        name: o!("ref"),
        length: 500,
        start: 100,
        stop: 140,
      },
      matches: 0,
      length: 0,
      quality: 10,
      cigar: parse_cigar_str("3I 6M 3I 3M 4D 5M 14I 7M 3D 4I 5M 5D 3M 3I".replace(' ', "")).unwrap(),
      orientation: Strand::Forward,
      divergence: Some(0.1),
      align: None,
    };

    let actual = split_matches(
      &aln,
      &SplitMatchesArgs {
        indel_len_threshold: 10,
      },
    )
    .unwrap();

    let expected = vec![
      Alignment {
        qry: Hit {
          name: o!("qry"),
          length: 500,
          start: 203,
          stop: 220,
        },
        reff: Hit {
          name: o!("ref"),
          length: 500,
          start: 100,
          stop: 118,
        },
        matches: 14,
        length: 21,
        quality: 10,
        cigar: parse_cigar_str("6M 3I 3M 4D 5MM".replace(' ', "")).unwrap(),
        orientation: Strand::Forward,
        divergence: Some(0.1),
        align: None,
      },
      Alignment {
        qry: Hit {
          name: o!("qry"),
          length: 500,
          start: 234,
          stop: 253,
        },
        reff: Hit {
          name: o!("ref"),
          length: 500,
          start: 118,
          stop: 141,
        },
        matches: 15,
        length: 27,
        quality: 10,
        cigar: parse_cigar_str("7M 3D 4I 5M 5D 3M".replace(' ', "")).unwrap(),
        orientation: Strand::Forward,
        divergence: Some(0.1),
        align: None,
      },
    ];

    assert_eq!(expected, actual);
  }

  #[rstest]
  fn test_split_matches_simple_case_reverse() {
    let aln = Alignment {
      qry: Hit {
        name: o!("qry"),
        length: 500,
        start: 200,
        stop: 256,
      },
      reff: Hit {
        name: o!("ref"),
        length: 500,
        start: 100,
        stop: 141,
      },
      matches: 0,
      length: 0,
      quality: 10,
      cigar: parse_cigar_str("3I 6M 3I 3M 4D 5M 14I 7M 3D 4I 5M 5D 3M 3I".replace(' ', "")).unwrap(),
      orientation: Strand::Reverse,
      divergence: Some(0.1),
      align: None,
    };

    let actual = split_matches(
      &aln,
      &SplitMatchesArgs {
        indel_len_threshold: 10,
      },
    )
    .unwrap();

    let expected = vec![
      Alignment {
        qry: Hit {
          name: o!("qry"),
          length: 500,
          start: 236,
          stop: 253,
        },
        reff: Hit {
          name: o!("ref"),
          length: 500,
          start: 100,
          stop: 118,
        },
        matches: 14,
        length: 21,
        quality: 10,
        cigar: parse_cigar_str("6M 3I 3M 4D 5M".replace(' ', "")).unwrap(),
        orientation: Strand::Reverse,
        divergence: Some(0.1),
        align: None,
      },
      Alignment {
        qry: Hit {
          name: o!("qry"),
          length: 500,
          start: 203,
          stop: 222,
        },
        reff: Hit {
          name: o!("ref"),
          length: 500,
          start: 118,
          stop: 141,
        },
        matches: 15,
        length: 27,
        quality: 10,
        cigar: parse_cigar_str("7M 3D 4I 5M 5D 3M".replace(' ', "")).unwrap(),
        orientation: Strand::Reverse,
        divergence: Some(0.1),
        align: None,
      },
    ];

    assert_eq!(expected, actual);
  }
}
