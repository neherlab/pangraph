use crate::align::nextclade::alphabet::letter::Letter;
use crate::align::nextclade::alphabet::nuc::Nuc;
use crate::align::nextclade::analyze::nuc_del::NucDelRange;
use crate::align::nextclade::analyze::nuc_sub::NucSub;
use crate::align::nextclade::coord::range::NucRefGlobalRange;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FindNucChangesOutput {
  pub substitutions: Vec<NucSub>,
  pub deletions: Vec<NucDelRange>,
  pub alignment_range: Option<NucRefGlobalRange>,
}

/// Finds nucleotide changes (nucleotide substitutions and deletions) as well
/// as the beginning and end of the alignment range.
///
/// @pre Precondition: sequences are expected to be aligned and stripped from insertions.
pub fn find_nuc_changes(qry_aln: &[Nuc], ref_aln: &[Nuc]) -> FindNucChangesOutput {
  assert_eq!(ref_aln.len(), qry_aln.len());

  let mut n_del: i64 = 0;
  let mut del_pos: i64 = -1;
  let mut before_alignment = true;

  let mut substitutions = Vec::<NucSub>::new();
  let mut deletions = Vec::<NucDelRange>::new();
  let mut alignment_start: i64 = -1;
  let mut alignment_end: i64 = -1;

  for i in 0..qry_aln.len() {
    let d = qry_aln[i];

    if !d.is_gap() {
      if before_alignment {
        alignment_start = i as i64;
        before_alignment = false;
      } else if n_del > 0 {
        deletions.push(NucDelRange::from_usize(del_pos as usize, (del_pos + n_del) as usize));
        n_del = 0;
      }
      alignment_end = (i + 1) as i64;
    }

    let ref_nuc = ref_aln[i];
    if !d.is_gap() && (d != ref_nuc) {
      substitutions.push(NucSub {
        ref_nuc,
        pos: i.into(),
        qry_nuc: d,
      });
    } else if d.is_gap() && !before_alignment {
      if n_del == 0 {
        del_pos = i as i64;
      }
      n_del += 1;
    }
  }

  substitutions.sort();
  deletions.sort();

  let alignment_range = (alignment_start >= 0 && alignment_end >= 0)
    .then(|| NucRefGlobalRange::from_usize(alignment_start as usize, alignment_end as usize));

  FindNucChangesOutput {
    substitutions,
    deletions,
    alignment_range,
  }
}
#[cfg(test)]
mod tests {

  use super::*;
  use crate::align::nextclade::alphabet::nuc::to_nuc_seq;

  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  // Simple test with no differences between query and reference.
  #[rstest]
  fn test_find_nuc_changes_no_differences() -> Result<(), Report> {
    // Assuming Nuc can be created from characters and that a gap is represented by Nuc::gap().
    // The implementation of Nuc::from(char) and Nuc::gap() is expected.
    let ref_aln = to_nuc_seq("GGGGGGGGGGGGGGGGGG")?;
    let qry_aln = to_nuc_seq("GGGGGGGGGGGGGGGGGG")?;

    let result = find_nuc_changes(&qry_aln, &ref_aln);
    let expected_result = FindNucChangesOutput {
      substitutions: vec![],
      deletions: vec![],
      alignment_range: Some(NucRefGlobalRange::from_usize(0, 18)),
    };
    assert_eq!(result, expected_result);

    Ok(())
  }

  // Test with a single substitution.
  #[rstest]
  fn test_find_nuc_changes_single_substitution() -> Result<(), Report> {
    let ref_aln = to_nuc_seq("GGGGGGGGGGGGGGGGGG")?;
    let qry_aln = to_nuc_seq("GGGAGGGGGGGGGGGGGG")?;

    let result = find_nuc_changes(&qry_aln, &ref_aln);
    let expected_result = FindNucChangesOutput {
      substitutions: vec![NucSub {
        ref_nuc: Nuc::G,
        pos: 3.into(),
        qry_nuc: Nuc::A,
      }],
      deletions: vec![],
      alignment_range: Some(NucRefGlobalRange::from_usize(0, 18)),
    };
    assert_eq!(result, expected_result);

    Ok(())
  }

  // test with deletion in query
  #[rstest]
  fn test_find_nuc_changes_single_deletion() -> Result<(), Report> {
    let ref_aln = to_nuc_seq("GGGGGGGGGGGGGGGGGG")?;
    let qry_aln = to_nuc_seq("GGG--GGGGGGGGGGGGG")?;

    let result = find_nuc_changes(&qry_aln, &ref_aln);
    let expected_result = FindNucChangesOutput {
      substitutions: vec![],
      deletions: vec![NucDelRange::from_usize(3, 5)],
      alignment_range: Some(NucRefGlobalRange::from_usize(0, 18)),
    };
    assert_eq!(result, expected_result);

    Ok(())
  }

  // test for deletion at the beginning / end of the query
  #[rstest]
  fn test_find_nuc_changes_deletion_at_edges() -> Result<(), Report> {
    let ref_aln = to_nuc_seq("GGGGGGGGGGGGGGGGGGG")?;
    let qry_aln = to_nuc_seq("--GGGGGGGGGGGGGGG--")?;

    let result = find_nuc_changes(&qry_aln, &ref_aln);
    let expected_result = FindNucChangesOutput {
      substitutions: vec![],
      deletions: vec![],
      alignment_range: Some(NucRefGlobalRange::from_usize(2, 17)),
    };
    assert_eq!(result, expected_result);

    Ok(())
  }

  // test with fully deleted query
  #[rstest]
  fn test_find_nuc_changes_full_deletion() -> Result<(), Report> {
    let ref_aln = to_nuc_seq("GGGGGGGGGGGGGGGGGG")?;
    let qry_aln = to_nuc_seq("------------------")?;

    let result = find_nuc_changes(&qry_aln, &ref_aln);
    let expected_result = FindNucChangesOutput {
      substitutions: vec![],
      deletions: vec![],
      alignment_range: None,
    };
    assert_eq!(result, expected_result);

    Ok(())
  }
}
