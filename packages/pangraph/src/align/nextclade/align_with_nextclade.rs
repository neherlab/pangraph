use crate::align::nextclade::align::align::align_nuc;
use crate::align::nextclade::align::backtrace::AlignmentOutput;
use crate::align::nextclade::align::gap_open::get_gap_open_close_scores_flat;
use crate::align::nextclade::align::insertions_strip::{insertions_strip, Insertion};
pub use crate::align::nextclade::align::params::NextalignParams;
use crate::align::nextclade::align::seed_match::CodonSpacedIndex;
use crate::align::nextclade::alphabet::nuc::{from_nuc_seq, to_nuc_seq, Nuc};
use crate::align::nextclade::analyze::nuc_changes::{find_nuc_changes, FindNucChangesOutput};
use crate::align::nextclade::analyze::nuc_del::NucDelRange;
use crate::align::nextclade::analyze::nuc_sub::NucSub;
use crate::io::fasta::FastaRecord;
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Serialize, Deserialize)]
pub struct AlignWithNextcladeOutput {
  pub qry_aln: String,
  pub ref_aln: String,
  pub substitutions: Vec<NucSub>,
  pub deletions: Vec<NucDelRange>,
  pub insertions: Vec<Insertion<Nuc>>,
  pub is_reverse_complement: bool,
}

pub fn align_with_nextclade(
  reff: impl AsRef<str>,
  qry: impl AsRef<str>,
  params: &NextalignParams,
) -> Result<AlignWithNextcladeOutput, Report> {
  let ref_seq = to_nuc_seq(reff.as_ref()).wrap_err("When converting reference sequence")?;
  let qry_seq = to_nuc_seq(qry.as_ref()).wrap_err("When converting query sequence")?;
  let seed_index = CodonSpacedIndex::from_sequence(&ref_seq);
  let gap_open_close = get_gap_open_close_scores_flat(&ref_seq, params);

  let alignment =
    align_nuc(0, "", &qry_seq, &ref_seq, &seed_index, &gap_open_close, params).wrap_err("When aligning sequences")?;

  let stripped = insertions_strip(&alignment.qry_seq, &alignment.ref_seq);

  let FindNucChangesOutput {
    substitutions,
    deletions,
    ..
  } = find_nuc_changes(&stripped.qry_seq, &ref_seq);

  Ok(AlignWithNextcladeOutput {
    qry_aln: from_nuc_seq(&stripped.qry_seq),
    ref_aln: from_nuc_seq(&alignment.ref_seq),
    substitutions,
    deletions,
    insertions: stripped.insertions,
    is_reverse_complement: alignment.is_reverse_complement,
  })
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::nextclade::coord::position::NucRefGlobalPosition;
  use crate::align::nextclade::coord::range::NucRefGlobalRange;
  use crate::o;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::{fixture, rstest};

  #[rstest]
  fn test_align_with_nextclade_general_case() -> Result<(), Report> {
    //                0         1         2         3         4         5         6         7         8         9
    //                0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    let ref_seq = o!("CTTGGAGGTTCCGTGGCTAGATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTGAGCCTTTGT");
    let qry_seq = o!("CTTGGAGGTTCCGTGGCTATAAAGATAACAGAACATTCTTGGAATGCTGATCAAGCTCATGGGACANNNNNCATGGTGGACAGCCTTTGT");
    let qry_aln = o!("CTTGGAGGTTCCGTGGCTAGATAACAGAACATTCTTGGAATGCTGATC-----AAGCTCATGGGACANNNNNCATGGTGAGCCTTTGT");
    let ref_aln = o!("CTTGGAGGTTCCGTGGCTA----GATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTG---AGCCTTTGT");

    let params = NextalignParams {
      min_length: 3,
      ..NextalignParams::default()
    };

    let actual = align_with_nextclade(ref_seq, qry_seq, &params)?;

    let expected = AlignWithNextcladeOutput {
      qry_aln,
      ref_aln,
      substitutions: vec![NucSub {
        pos: NucRefGlobalPosition::new(62),
        ref_nuc: Nuc::C,
        qry_nuc: Nuc::G,
      }],
      deletions: vec![NucDelRange::from_usize(48, 53)],
      insertions: vec![
        Insertion {
          pos: 18,
          ins: to_nuc_seq("TAAA")?,
        },
        Insertion {
          pos: 78,
          ins: to_nuc_seq("GAC")?,
        },
      ],

      is_reverse_complement: false,
    };

    assert_eq!(expected, actual);

    Ok(())
  }
}
