use crate::align::nextclade::align::backtrace::{AlignmentOutput, backtrace};
use crate::align::nextclade::align::band_2d::Stripe;
use crate::align::nextclade::align::band_2d::{full_matrix, simple_stripes};
use crate::align::nextclade::align::params::NextalignParams;
use crate::align::nextclade::align::score_matrix::{ScoreMatrixResult, score_matrix};
use crate::align::nextclade::align::seed_alignment::create_alignment_band;
use crate::align::nextclade::align::seed_match::{
  CodonSpacedIndex, SeedMatchesResult, get_seed_matches_maybe_reverse_complement,
};
use crate::align::nextclade::alphabet::letter::Letter;
use crate::align::nextclade::alphabet::nuc::Nuc;
use crate::make_error;
use eyre::{Report, WrapErr};
use log::{debug, trace, warn};
use std::cmp::max;

fn align_pairwise<T: Letter<T>>(
  qry_seq: &[T],
  ref_seq: &[T],
  gap_open_close: &[i32],
  params: &NextalignParams,
  stripes: &[Stripe],
) -> AlignmentOutput<T> {
  trace!("Align pairwise: started. Params: {params:?}");

  let ScoreMatrixResult { scores, paths } = score_matrix(qry_seq, ref_seq, gap_open_close, stripes, params);

  backtrace(qry_seq, ref_seq, &scores, &paths)
}

pub fn align_nuc_simplestripe(
  qry_seq: &[Nuc],
  ref_seq: &[Nuc],
  gap_open_close: &[i32],
  mean_shift: i32,
  params: &NextalignParams,
) -> Result<AlignmentOutput<Nuc>, Report> {
  let qry_len = qry_seq.len();
  let ref_len = ref_seq.len();
  let min_len = params.min_length;
  if qry_len < min_len {
    return make_error!(
      "Unable to align: sequence is too short. Details: sequence length: {qry_len}, min length allowed: {min_len}. This is likely due to a low quality of the provided sequence, or due to using incorrect reference sequence."
    );
  }

  // mean shift equal to half difference between ref and query lengths
  // let mean_shift = (ref_len as i32 - qry_len as i32) / 2;
  // band width equal to 2 times absolute difference + 10
  let mut band_width = (mean_shift.abs() * 4 + 10) as usize;
  let mut stripes = simple_stripes(mean_shift, band_width, ref_len, qry_len);

  let mut attempt = 0;
  let mut alignment = align_pairwise(qry_seq, ref_seq, gap_open_close, params, &stripes);

  while alignment.hit_boundary && attempt < params.max_alignment_attempts {
    debug!(
      "In nucleotide alignment: Band boundary is hit on attempt {}. Retrying with relaxed parameters. Alignment score was: {}",
      attempt + 1,
      alignment.alignment_score
    );
    // double bandwidth parameters or increase to one if 0
    band_width = max(2 * band_width, 1);
    stripes = simple_stripes(mean_shift, band_width, ref_len, qry_len);
    // realign
    attempt += 1;
    alignment = align_pairwise(qry_seq, ref_seq, gap_open_close, params, &stripes);
  }

  // report success/failure of broadening of band width
  if alignment.hit_boundary {
    warn!(
      "In nucleotide alignment: Attempted to relax band parameters {attempt} times, but still hitting the band boundary. Returning last attempt with score: {}",
      alignment.alignment_score
    );
  } else if attempt > 0 {
    debug!(
      "In nucleotide alignment: Succeeded without hitting band boundary on attempt {}. Alignment score was: {}",
      attempt + 1,
      alignment.alignment_score
    );
  }
  Ok(alignment)
}

/// align nucleotide sequences via seed alignment and banded smith watermann without penalizing terminal gaps
pub fn align_nuc(
  index: usize,
  seq_name: &str,
  qry_seq: &[Nuc],
  ref_seq: &[Nuc],
  seed_index: &CodonSpacedIndex,
  gap_open_close: &[i32],
  params: &NextalignParams,
) -> Result<Option<AlignmentOutput<Nuc>>, Report> {
  let qry_len = qry_seq.len();
  let ref_len = ref_seq.len();
  let min_len = params.min_length;
  if qry_len < min_len {
    return make_error!(
      "Unable to align: sequence is too short. Details: sequence length: {qry_len}, min length allowed: {min_len}. This is likely due to a low quality of the provided sequence, or due to using incorrect reference sequence."
    );
  }

  if ref_len + qry_len < (20 * params.kmer_length) {
    // for very short sequences, use full square
    let stripes = full_matrix(ref_len, qry_len);
    trace!(
      "When processing sequence #{index} '{seq_name}': In nucleotide alignment: Band construction: short sequences, using full matrix"
    );
    return Ok(Some(align_pairwise(qry_seq, ref_seq, gap_open_close, params, &stripes)));
  }

  // otherwise, determine seed matches roughly regularly spaced along the query sequence
  let SeedMatchesResult {
    qry_seq,
    seed_matches,
    is_reverse_complement,
  } = get_seed_matches_maybe_reverse_complement(qry_seq, ref_seq, seed_index, params)
    .wrap_err("When calculating seed matches")?;

  if seed_matches.is_empty() {
    return Ok(None);
  }

  let mut terminal_bandwidth = params.terminal_bandwidth as isize;
  let mut excess_bandwidth = params.excess_bandwidth as isize;
  let mut minimal_bandwidth = max(1, params.allowed_mismatches as isize);
  let max_band_area = params.max_band_area;
  let mut attempt = 0;

  let (mut stripes, mut band_area) = create_alignment_band(
    &seed_matches,
    qry_len as isize,
    ref_len as isize,
    terminal_bandwidth,
    excess_bandwidth,
    minimal_bandwidth,
  );
  if band_area > max_band_area {
    return make_error!(
      "Alignment matrix size {band_area} exceeds maximum value {max_band_area}. The threshold can be adjusted using CLI flag '--max-band-area' or using 'maxBandArea' field in the dataset's pathogen.json"
    );
  }

  let mut alignment = align_pairwise(&qry_seq, ref_seq, gap_open_close, params, &stripes);

  while alignment.hit_boundary && attempt < params.max_alignment_attempts {
    debug!(
      "When processing sequence #{index} '{seq_name}': In nucleotide alignment: Band boundary is hit on attempt {}. Retrying with relaxed parameters. Alignment score was: {}",
      attempt + 1,
      alignment.alignment_score
    );
    // double bandwidth parameters or increase to one if 0
    terminal_bandwidth = max(2 * terminal_bandwidth, 1);
    excess_bandwidth = max(2 * excess_bandwidth, 1);
    minimal_bandwidth = max(2 * minimal_bandwidth, 1);
    attempt += 1;
    // make new band
    (stripes, band_area) = create_alignment_band(
      &seed_matches,
      qry_len as isize,
      ref_len as isize,
      terminal_bandwidth,
      excess_bandwidth,
      minimal_bandwidth,
    );
    // discard stripes and break to return previous alignment
    if band_area > max_band_area {
      break;
    }
    // realign
    alignment = align_pairwise(&qry_seq, ref_seq, gap_open_close, params, &stripes);
  }
  // report success/failure of broadening of band width
  if alignment.hit_boundary {
    debug!(
      "When processing sequence #{index} '{seq_name}': In nucleotide alignment: Attempted to relax band parameters {attempt} times, but still hitting the band boundary. Returning last attempt with score: {}",
      alignment.alignment_score
    );
    if band_area > max_band_area {
      debug!(
        "When processing sequence #{index} '{seq_name}': final band area {band_area} exceeded the cutoff {max_band_area}"
      );
    }
  } else if attempt > 0 {
    debug!(
      "When processing sequence #{index} '{seq_name}': In nucleotide alignment: Succeeded without hitting band boundary on attempt {}. Alignment score was: {}",
      attempt + 1,
      alignment.alignment_score
    );
  }
  alignment.is_reverse_complement = is_reverse_complement;
  Ok(Some(alignment))
}
