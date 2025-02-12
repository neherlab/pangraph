use crate::align::nextclade::align::params::NextalignParams;
use crate::align::nextclade::alphabet::nuc::Nuc;

pub type GapScoreMap = Vec<i32>;

pub fn get_gap_open_close_scores_flat(ref_seq: &[Nuc], params: &NextalignParams) -> GapScoreMap {
  let value = params.penalty_gap_open;
  let len = ref_seq.len() + 2;
  vec![value; len]
}
