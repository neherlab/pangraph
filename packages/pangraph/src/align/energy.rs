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
  (hit.interval.start > minblock) as usize + ((hit.length - hit.interval.end) > minblock) as usize
}

/// Calculate the energy of the alignment
///
/// The formula for the energy is:
///
/// $ E = -L + \alpha C + \beta M $
///
/// where:
/// - `L` is the alignment length. We can measure it as number of matches, to avoid counting indels.
/// - `C` is the number of *cuts* we need to make if we want to merge the block on this alignment. This depend on whether the alignment extends to the beginning/end of the query/reference sequences. If not we need to introduce a cut. This number is between 0 and 4.
/// - `M` is the number of mismatches in the alignment. It can be estimated from the divergence.
/// - `alpha` and `beta` are the CLI parameters of pangraph. They act as weight, and control how much penalty is given for each mismatch and each cut.
#[allow(non_snake_case)]
pub fn alignment_energy2(aln: &Alignment, args: &AlignmentArgs) -> f64 {
  let L = aln.matches;
  let M = aln.divergence.unwrap_or_default() * L as f64;
  let mut C = 4;
  if aln.qry.interval.start == 0 {
    C -= 1;
  }
  if aln.qry.interval.end == aln.qry.length {
    C -= 1;
  }
  if aln.reff.interval.start == 0 {
    C -= 1;
  }
  if aln.reff.interval.end == aln.reff.length {
    C -= 1;
  }
  -(L as f64) + (C as f64) * args.alpha + M * args.beta
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::bam::cigar::parse_cigar_str;
  use crate::o;
  use crate::pangraph::strand::Strand;
  use crate::utils::interval::Interval;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_alignment_energy_simple_case() {
    let aln = Alignment {
      qry: Hit::new("qry_3", 997, (0, 980)),
      reff: Hit::new("ref_3", 1000, (18, 1000)),
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

  #[rstest]
  fn test_alignment_energy2_simple_case() {
    let aln = Alignment {
      qry: Hit::new("bl_1", 100, (0, 50)),
      reff: Hit::new("bl_2", 200, (120, 200)),
      matches: 40,
      length: 60,
      quality: 100,
      orientation: Strand::Forward,
      cigar: parse_cigar_str("10I40M10D").unwrap(),
      divergence: Some(0.02),
      align: Some(0.1),
    };

    let args = AlignmentArgs {
      alpha: 10.0,
      beta: 10.0,
      ..Default::default()
    };

    assert_ulps_eq!(alignment_energy2(&aln, &args), -12.0);
  }
}
