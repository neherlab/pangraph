#![allow(unsafe_code)]
use crate::ptr::is_null;
use eyre::{eyre, Report};
use minimap2_sys::{mm_check_opt, mm_idxopt_t, mm_mapopt_t, mm_set_opt, MM_F_CIGAR, MM_F_NO_DIAG};
use std::ffi::{c_char, CString};
use std::mem::MaybeUninit;
use std::ptr::null;

#[derive(Copy, Clone, Debug)]
pub enum Minimap2Preset {
  LrHqae,
  LrHq,
  Splice,
  SpliceHq,
  Asm,
  Asm5,
  Asm10,
  Asm20,
  Sr,
  MapPb,
  MapHifi,
  MapOnt,
  AvaPb,
  AvaOnt,
  Short,
  Map10k,
  Cdna,
}

impl Minimap2Preset {
  pub fn as_str(&self) -> &str {
    match self {
      Minimap2Preset::LrHqae => "lr:hqae",
      Minimap2Preset::LrHq => "lr:hq",
      Minimap2Preset::Splice => "splice",
      Minimap2Preset::SpliceHq => "splice:hq",
      Minimap2Preset::Asm => "asm",
      Minimap2Preset::Asm5 => "asm5",
      Minimap2Preset::Asm10 => "asm10",
      Minimap2Preset::Asm20 => "asm20",
      Minimap2Preset::Sr => "sr",
      Minimap2Preset::MapPb => "map-pb",
      Minimap2Preset::MapHifi => "map-hifi",
      Minimap2Preset::MapOnt => "map-ont",
      Minimap2Preset::AvaPb => "ava-pb",
      Minimap2Preset::AvaOnt => "ava-ont",
      Minimap2Preset::Short => "short",
      Minimap2Preset::Map10k => "map10k",
      Minimap2Preset::Cdna => "cdna",
    }
  }

  pub fn to_c_string(&self) -> CString {
    CString::new(self.as_str()).expect("CString::new failed")
  }
}

#[must_use]
#[derive(Clone, Debug)]
pub struct Minimap2Options {
  idx_opt: mm_idxopt_t,
  map_opt: mm_mapopt_t,
}

impl Minimap2Options {
  pub fn new() -> Result<Self, Report> {
    init_options(null())
  }

  pub fn with_preset(preset: Minimap2Preset) -> Result<Self, Report> {
    init_options(preset.to_c_string().as_ptr())
  }

  pub(crate) fn idx_opt(&self) -> &mm_idxopt_t {
    &self.idx_opt
  }

  pub(crate) fn idx_opt_mut(&mut self) -> *mut mm_idxopt_t {
    &mut self.idx_opt
  }

  pub(crate) fn map_opt(&self) -> &mm_mapopt_t {
    &self.map_opt
  }

  pub(crate) fn map_opt_mut(&mut self) -> *mut mm_mapopt_t {
    &mut self.map_opt
  }
}

fn init_options(preset: *const c_char) -> Result<Minimap2Options, Report> {
  let mut idx_opt = MaybeUninit::<mm_idxopt_t>::uninit();
  let mut map_opt = MaybeUninit::<mm_mapopt_t>::uninit();

  // SAFETY: calling unsafe function
  unsafe { mm_set_opt(null(), idx_opt.as_mut_ptr(), map_opt.as_mut_ptr()) };

  // SAFETY: calling unsafe function
  if unsafe { mm_set_opt(preset, idx_opt.as_mut_ptr(), map_opt.as_mut_ptr()) } != 0 {
    return Err(eyre!("minimap2: mm_set_opt(): failed to set options: incorrect preset"));
  }

  {
    // SAFETY: calling unsafe function
    let idx_opt = unsafe { idx_opt.assume_init_mut() };

    // SAFETY: calling unsafe function
    let map_opt = unsafe { map_opt.assume_init_mut() };

    apply_additional_options(idx_opt, map_opt)?;
  }

  if is_null(&idx_opt) || is_null(&map_opt) {
    return Err(eyre!("minimap2: mm_set_opt(): null pointer returned"));
  }

  // SAFETY: calling unsafe function
  if unsafe { mm_check_opt(idx_opt.as_ptr(), map_opt.as_ptr()) } != 0 {
    return Err(eyre!("minimap2: mm_check_opt(): options are invalid"));
  }

  if is_null(&idx_opt) || is_null(&map_opt) {
    return Err(eyre!("minimap2: mm_check_opt(): null pointer returned"));
  }

  // SAFETY: calling unsafe function
  let idx_opt = unsafe { idx_opt.assume_init() };

  // SAFETY: calling unsafe function
  let map_opt = unsafe { map_opt.assume_init() };

  Ok(Minimap2Options { idx_opt, map_opt })
}

fn apply_additional_options(idx_opt: &mut mm_idxopt_t, map_opt: &mut mm_mapopt_t) -> Result<(), Report> {
  let minblock = 10;

  map_opt.flag |= MM_F_CIGAR as i64; // Output CIGAR (in minimap2 CLI: `-c`)

  // FIXME: make these configurable
  idx_opt.k = 15; // k-mer size (no larger than 28) (in minimap2 CLI: `-k`, default: 15)

  // This is `-X`, but with `MM_F_NO_DUAL` it returns empty results when e.g. ref name is "0" and query name is "1"
  // map_opt.flag |= (MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN) as i64; // (in minimap2 CLI: `-X`)
  // In julia version only this flag is set.
  map_opt.flag |= MM_F_NO_DIAG as i64;

  // FIXME: Overrides from Julia version
  map_opt.min_dp_max = minblock - 10; // minimal peak DP alignment score (in minimap2 CLI: `-s`, default: 80)
  idx_opt.bucket_bits = 14;

  // FIXME: from minimap2/python/mappy.pyx
  // idx_opt.batch_size = 0x7fffffffffffffff; // always build a uni-part index
  // map_opt.mid_occ = 1000; // don't filter high-occ seeds

  Ok(())
}

//   Indexing:
//     -H           use homopolymer-compressed k-mer (preferrable for PacBio)
//     -k INT       k-mer size (no larger than 28) [15]
//     -w INT       minimizer window size [10]
//     -I NUM       split index for every ~NUM input bases [8G]
//     -d FILE      dump index to FILE []
//   Mapping:
//     -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [0.0002]
//     -g NUM       stop chain enlongation if there are no minimizers in INT-bp [5000]
//     -G NUM       max intron length (effective with -xsplice; changing -r) [200k]
//     -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]
//     -r NUM[,NUM] chaining/alignment bandwidth and long-join bandwidth [500,20000]
//     -n INT       minimal number of minimizers on a chain [3]
//     -m INT       minimal chaining score (matching bases minus log gap penalty) [40]
//     -X           skip self and dual mappings (for the all-vs-all mode)
//     -p FLOAT     min secondary-to-primary score ratio [0.8]
//     -N INT       retain at most INT secondary alignments [5]

// fn main() {
//   let mut idxopt = mm_idxopt_t {
//     k: 15, // default k-mer size
//     w: 10, // default minimizer window size
//     flag: 0,
//     batch_size: 8_000_000_000, // default batch size
//     bucket_bits: 0,
//     mini_batch_size: 0,
//   };
//
//   let mut mapopt = mm_mapopt_t {
//     flag: 0,
//     min_chain_score: 40, // default minimal chaining score
//     max_gap: 5000,       // default stop chain elongation gap
//     max_frag_len: 800,   // default max fragment length
//     pri_ratio: 0.8,      // default min secondary-to-primary score ratio
//     best_n: 5,           // retain at most INT secondary alignments
//     min_cnt: 3,          // minimal number of minimizers on each chain
//     bw: 500,
//     bw_long: 20000,
//     max_gap_ref: 0,
//     ..Default::default()
//   };
//
//   let args = std::env::args().collect::<Vec<String>>();
//   let mut i = 1;
//   while i < args.len() {
//     match args[i].as_str() {
//       "-k" => idxopt.k = args[i + 1].parse().unwrap(),
//       "-w" => idxopt.w = args[i + 1].parse().unwrap(),
//       "-H" => idxopt.flag |= MM_I_HPC as i16,
//       "-I" => idxopt.batch_size = mm_parse_num(&args[i + 1]),
//       "-d" => {} // Dump index file handling not included
//       "-m" => mapopt.min_chain_score = args[i + 1].parse().unwrap(),
//       "-g" => mapopt.max_gap = args[i + 1].parse().unwrap(),
//       "-G" => {
//         let max_intron_len = args[i + 1].parse::<i32>().unwrap();
//         mm_mapopt_max_intron_len(&mut mapopt, max_intron_len);
//       }
//       "-F" => mapopt.max_frag_len = args[i + 1].parse().unwrap(),
//       "-r" => {
//         let (bw, bw_long) = parse_bandwidths(&args[i + 1]);
//         mapopt.bw = bw;
//         mapopt.bw_long = bw_long;
//       }
//       "-n" => mapopt.min_cnt = args[i + 1].parse().unwrap(),
//       "-p" => mapopt.pri_ratio = args[i + 1].parse().unwrap(),
//       "-N" => mapopt.best_n = args[i + 1].parse().unwrap(),
//       "-X" => mapopt.flag |= (MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN) as i64,
//       _ => {}
//     }
//     i += 2; // move to the next flag
//   }
// }
//
// fn mm_parse_num(s: &str) -> u64 {
//   let mut multiplier = 1.0;
//   let value_str = s.trim_end_matches(|c: char| !c.is_digit(10));
//   let unit = s[value_str.len()..].to_lowercase();
//   match unit.as_str() {
//     "g" => multiplier = 1e9,
//     "m" => multiplier = 1e6,
//     "k" => multiplier = 1e3,
//     _ => {}
//   }
//   (value_str.parse::<f64>().unwrap() * multiplier).round() as u64
// }
//
// fn parse_bandwidths(s: &str) -> (i32, i32) {
//   let parts: Vec<&str> = s.split(',').collect();
//   let bw = parts[0].parse().unwrap();
//   let bw_long = if parts.len() > 1 { parts[1].parse().unwrap() } else { 0 };
//   (bw, bw_long)
// }
//
// fn mm_mapopt_max_intron_len(opt: &mut mm_mapopt_t, max_intron_len: i32) {
//   if (opt.flag & MM_F_SPLICE as i64) != 0 && max_intron_len > 0 {
//     opt.max_gap_ref = max_intron_len;
//     opt.bw = max_intron_len;
//     opt.bw_long = max_intron_len;
//   }
// }
