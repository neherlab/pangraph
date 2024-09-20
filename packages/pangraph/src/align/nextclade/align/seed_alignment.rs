use crate::align::nextclade::align::band_2d::Stripe;
use crate::align::nextclade::align::seed_match::SeedMatch2;
use crate::align::nextclade::alphabet::letter::Letter;

use num_traits::abs;
use num_traits::clamp;
use std::cmp::max;
use std::cmp::min;

/// generate a vector of query sequence positions that are followed by at least `seed_length`
/// valid characters. Positions in this vector are thus "good" positions to start a query k-mer.
fn get_map_to_good_positions<L: Letter<L>>(qry_seq: &[L], seed_length: usize) -> Vec<usize> {
  let qry_len = qry_seq.len();

  let mut map_to_good_positions = Vec::<usize>::with_capacity(qry_len);
  let mut distance_to_last_bad_pos: i32 = 0;

  for (i, letter) in qry_seq.iter().enumerate() {
    // TODO: Exclude ambiguous letters
    if letter.is_unknown() {
      distance_to_last_bad_pos = 0;
    } else if distance_to_last_bad_pos >= seed_length as i32 {
      map_to_good_positions.push(i - seed_length);
    }
    distance_to_last_bad_pos += 1;
  }

  map_to_good_positions
}

#[derive(Debug, Clone, Copy)]
pub struct SeedMatch {
  pub qry_pos: usize,
  pub ref_pos: usize,
  pub score: usize,
}

fn abs_shift(seed1: &SeedMatch2, seed2: &SeedMatch2) -> isize {
  abs(seed2.offset - seed1.offset)
}

// Function that processes a seed, adds a band, rewind the band vector as necessary
// to accommodate any extended band width implied by look-back-length
struct RewindResult {
  look_back_length: isize,
  current_ref_end: isize,
}
fn extend_and_rewind(
  current_seed: &SeedMatch2,
  current_band: TrapezoidDirectParams,
  bands: &mut Vec<TrapezoidDirectParams>,
  mean_offset: isize,
  look_forward_length: isize,
  look_back_length: isize,
  minimal_bandwidth: isize,
  ref_len: isize,
) -> RewindResult {
  let mut current_band = current_band;
  let mut look_back_length = look_back_length;
  // generate new current trapezoid for the body of the current seed
  let current_seed_end = (current_seed.ref_pos + current_seed.length) as isize;
  if current_seed_end > current_band.ref_end {
    bands.push(current_band);
    current_band = TrapezoidDirectParams {
      ref_start: current_band.ref_end,
      ref_end: current_seed_end,
      min_offset: current_seed.offset - minimal_bandwidth,
      max_offset: current_seed.offset + minimal_bandwidth,
    };
  }
  look_back_length = max(
    max(look_back_length, mean_offset - current_band.min_offset),
    current_band.max_offset - mean_offset,
  );

  // rewind the bands until the ref_start of the last one preceeds the one to add
  while current_band.ref_start > max(0, current_seed_end - look_back_length - 1) {
    current_band = if let Some(tmp_band) = bands.pop() {
      tmp_band
    } else {
      // we rewound all the way to the beginning, add a new terminal
      TrapezoidDirectParams {
        ref_start: 0,
        ref_end: min(current_seed.ref_pos as isize + look_forward_length, ref_len + 1),
        min_offset: current_seed.offset - look_back_length,
        max_offset: current_seed.offset + look_back_length,
      }
    };
    // increase look-back-length to cover absorbed previous bands
    look_back_length = max(
      max(look_back_length, mean_offset - current_band.min_offset),
      current_band.max_offset - mean_offset,
    );
  }
  // terminate previous trapezoid where the new one will start and push unless band is empty
  current_band.ref_end = max(0, current_seed_end - look_back_length);
  if current_band.ref_end > current_band.ref_start {
    bands.push(current_band);
  }
  // return the updated look-back-length and the end of the current band
  RewindResult {
    look_back_length,
    current_ref_end: current_band.ref_end,
  }
}

/// Takes in seed matches and returns a vector of stripes
/// Stripes define the query sequence range for each reference position
pub fn create_alignment_band(
  chain: &[SeedMatch2],
  qry_len: isize,
  ref_len: isize,
  terminal_bandwidth: isize,
  excess_bandwidth: isize,
  minimal_bandwidth: isize,
) -> (Vec<Stripe>, usize) {
  // This function steps through the chained seeds and determines and appropriate band
  // defined via stripes in query coordinates. These bands will later be chopped to reachable ranges

  // the broad idea is the following:
  // pre: deal with a special case the beginning and allow for terminal bandwidth
  // within: for each pair of chained seed (current and next), add a Trapezoid centered at the junction
  //         and a body Trapezoid for next. Extend the gap Trapezoid into the current and next
  // post: deal with the terminal trapezoid and allow of terminal bandwidth

  let mut bands = Vec::<TrapezoidDirectParams>::with_capacity(2 * chain.len() + 2);
  // make initial trapezoid starting at 0 and extending into match by terminal_bandwidth
  let mut current_seed = &chain[0];
  let mut look_back_length = terminal_bandwidth;
  let mut look_forward_length = terminal_bandwidth;
  let mut current_ref_end = min(current_seed.ref_pos as isize + look_forward_length, ref_len + 1);
  let mut current_band = TrapezoidDirectParams {
    ref_start: 0,
    ref_end: current_ref_end,
    min_offset: current_seed.offset - terminal_bandwidth,
    max_offset: current_seed.offset + terminal_bandwidth,
  };

  // loop over remaining seeds in chain
  for next_seed in chain.iter().skip(1) {
    let mean_offset = (next_seed.offset + current_seed.offset) / 2; // offset of gap seed
    let shift = abs_shift(current_seed, next_seed) / 2; // distance from mean offset
    look_forward_length = shift + excess_bandwidth;

    // attempt to add new a band for current_seed, then rewind as necessary to accommodate shift
    RewindResult {
      look_back_length,
      current_ref_end,
    } = extend_and_rewind(
      current_seed,
      current_band,
      &mut bands,
      mean_offset,
      look_forward_length,
      shift + excess_bandwidth,
      minimal_bandwidth,
      ref_len,
    );

    // generate trapezoid for the gap between seeds that goes look-forward-length into the next seed
    current_band = TrapezoidDirectParams {
      ref_start: current_ref_end,
      ref_end: next_seed.ref_pos as isize + look_forward_length,
      min_offset: mean_offset - max(look_back_length, excess_bandwidth),
      max_offset: mean_offset + max(look_back_length, excess_bandwidth),
    };
    current_seed = next_seed;
  }

  // process the final seed (different offset)
  RewindResult {
    look_back_length,
    current_ref_end,
  } = extend_and_rewind(
    current_seed,
    current_band,
    &mut bands,
    current_seed.offset,
    look_forward_length,
    max(terminal_bandwidth, look_back_length),
    minimal_bandwidth,
    ref_len,
  );
  // add band that extends all the way to the end
  current_band = TrapezoidDirectParams {
    ref_start: current_ref_end,
    ref_end: ref_len + 1,
    min_offset: current_seed.offset - look_back_length,
    max_offset: current_seed.offset + look_back_length,
  };
  bands.push(current_band);

  // construct strips from the ordered bands
  let mut stripes = Vec::<Stripe>::with_capacity(ref_len as usize + 1);
  for band in bands {
    for ref_pos in band.ref_start..band.ref_end {
      stripes.push(Stripe {
        begin: (ref_pos + band.min_offset).clamp(0, qry_len - minimal_bandwidth) as usize,
        end: (ref_pos + band.max_offset).clamp(0, qry_len) as usize + 1,
      });
    }
  }

  // trim stripes to reachable regions
  regularize_stripes(stripes, qry_len as usize)
}

#[derive(Clone, Copy, Debug)]
struct TrapezoidDirectParams {
  ref_start: isize,
  ref_end: isize,
  min_offset: isize,
  max_offset: isize,
}

/// Chop off unreachable parts of the stripes.
/// Overhanging parts are pruned
fn regularize_stripes(mut stripes: Vec<Stripe>, qry_len: usize) -> (Vec<Stripe>, usize) {
  // assure stripe begin are non-decreasing -- such states would be unreachable in the alignment
  let stripes_len = stripes.len();
  stripes[0].begin = 0;
  for i in 1..stripes_len {
    stripes[i].begin = clamp(stripes[i].begin, stripes[i - 1].begin, qry_len);
  }

  // analogously, assure that strip ends are non-decreasing. this needs to be done in reverse.
  stripes[stripes_len - 1].end = qry_len + 1;
  let mut band_area = stripes[stripes_len - 1].end - stripes[stripes_len - 1].begin;
  for i in (0..(stripes_len - 1)).rev() {
    stripes[i].end = clamp(stripes[i].end, stripes[i].begin + 1, stripes[i + 1].end);
    band_area += stripes[i].end - stripes[i].begin;
  }

  (stripes, band_area)
}
