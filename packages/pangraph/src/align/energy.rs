use crate::align::alignment::{Alignment, Hit};
use crate::align::alignment_args::AlignmentArgs;
use clap::Args;
use serde::{Deserialize, Serialize};

/// Calculate the energy of the alignment.
///
/// This is a function of alignment length, identity, and whether it will create a block split.
pub fn alignment_energy(aln: &Alignment, args: &AlignmentArgs) -> f64 {
  let minblock = args.indel_len_threshold;
  if aln.length < minblock {
    return f64::INFINITY;
  }

  let ncuts = (cuts(&aln.qry, args) + cuts(&aln.reff, args)) as f64;
  let nmuts = aln.divergence.unwrap_or_default() * aln.length as f64;

  -(aln.length as f64) + args.alpha * ncuts + args.beta * nmuts
}

fn cuts(hit: &Hit, args: &AlignmentArgs) -> usize {
  let minblock = args.indel_len_threshold;
  // FIXME: Addition of booleans?
  (hit.start > minblock) as usize + ((hit.length - hit.stop) > minblock) as usize
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::bam::cigar::parse_cigar_str;
  use crate::o;
  use crate::pangraph::strand::Strand;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_alignment_energy_simple_case() {
    let aln = Alignment {
      qry: Hit {
        name: o!("qry_3"),
        length: 997,
        start: 0,
        stop: 980,
      },
      reff: Hit {
        name: o!("ref_3"),
        length: 1000,
        start: 18,
        stop: 1000,
      },
      matches: 965,
      length: 982,
      quality: 60,
      orientation: Strand::Reverse,
      cigar: parse_cigar_str("124M1D416M1D440M").unwrap(),
      divergence: Some(0.0173),
      align: Some(889.0),
    };

    let args = AlignmentArgs::default();

    assert_ulps_eq!(alignment_energy(&aln, &args), -812.114);
  }
}
