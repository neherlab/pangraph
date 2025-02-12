use crate::align::nextclade::align::band_2d::Band2d;
use crate::align::nextclade::align::score_matrix::{
  BOUNDARY, MATCH, QRY_GAP_EXTEND, QRY_GAP_MATRIX, REF_GAP_EXTEND, REF_GAP_MATRIX,
};
use crate::align::nextclade::alphabet::letter::Letter;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct AlignmentOutput<T> {
  pub qry_seq: Vec<T>,
  pub ref_seq: Vec<T>,
  pub alignment_score: i32,
  pub is_reverse_complement: bool,
  pub hit_boundary: bool,
}

pub fn backtrace<T: Letter<T>>(
  qry_seq: &[T],
  ref_seq: &[T],
  scores: &Band2d<i32>,
  paths: &Band2d<i8>,
) -> AlignmentOutput<T> {
  let num_cols = scores.num_cols();
  let num_rows = scores.num_rows();

  // max length of the alignment is the sum of query and reference length
  let aln_capacity = scores.num_cols() + scores.num_rows();
  let mut aln_ref = Vec::<T>::with_capacity(aln_capacity);
  let mut aln_qry = Vec::<T>::with_capacity(aln_capacity);

  // Add right overhang, i.e. unaligned parts of the query or reference
  let mut r_pos = num_rows - 1;
  let mut q_pos = num_cols - 1;

  let mut origin: i8;
  let mut current_matrix = 0;
  let mut hit_boundary = false;
  // Do backtrace in the aligned region
  while r_pos > 0 || q_pos > 0 {
    origin = paths[(r_pos, q_pos)];
    if (origin & BOUNDARY) > 0 {
      hit_boundary = true;
    }

    if (origin & MATCH) != 0 && (current_matrix == 0) {
      // Match -- decrement both strands and add match to alignment
      q_pos -= 1;
      r_pos -= 1;
      aln_qry.push(qry_seq[q_pos]);
      aln_ref.push(ref_seq[r_pos]);
    } else if ((origin & REF_GAP_MATRIX) != 0 && current_matrix == 0) || current_matrix == REF_GAP_MATRIX {
      // Insertion in ref -- decrement query, increase shift
      q_pos -= 1;
      aln_qry.push(qry_seq[q_pos]);
      aln_ref.push(T::GAP);
      current_matrix = if (origin & REF_GAP_EXTEND) != 0 {
        // Remain in gap-extension mode and ignore best-overall score
        REF_GAP_MATRIX
      } else {
        // Close gap, return to best-overall score
        0
      }
    } else if ((origin & QRY_GAP_MATRIX) != 0 && current_matrix == 0) || current_matrix == QRY_GAP_MATRIX {
      // Deletion in query -- decrement reference, reduce shift
      aln_qry.push(T::GAP);
      r_pos -= 1;
      aln_ref.push(ref_seq[r_pos]);
      current_matrix = if (origin & QRY_GAP_EXTEND) != 0 {
        // Remain in gap-extension mode and ignore best-overall score
        QRY_GAP_MATRIX
      } else {
        // Close gap, return to best-overall score
        0
      }
    } else {
      // This should never be reached
      // origin = 0 and current_matrix = 0
      // Why would this ever happen?
      // Mistake in score_matrix?
      // TODO: This actually does seem to be reachable, at least when band is width 0, i.e. a line
      unreachable!("Problem in backtrace: origin = 0 and current_matrix = 0 before (0,0) reached. Please share the sequence with the developers.\nr_pos = {}, q_pos = {}, origin = {}, current_matrix = {}", r_pos, q_pos, origin, current_matrix);
    }
  }

  aln_qry.reverse();
  aln_ref.reverse();

  AlignmentOutput {
    qry_seq: aln_qry,
    ref_seq: aln_ref,
    alignment_score: scores[(num_rows - 1, num_cols - 1)],
    is_reverse_complement: false,
    hit_boundary,
  }
}
