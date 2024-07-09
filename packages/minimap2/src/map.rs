#![allow(unsafe_code)]

use crate::buf::Minimap2Buffer;
use crate::index::Minimap2Index;
use eyre::Report;
use itertools::Itertools;
use minimap2_sys::{free, mm_extra_t, mm_idx_t, mm_map, mm_reg1_t, MM_CIGAR_CHARS};
use ordered_float::OrderedFloat;
use std::ffi::CString;
use std::os::raw::c_int;
use std::slice::from_raw_parts;

#[must_use]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Minimap2Result {
  pub regs: Vec<Minimap2Reg>,
}

#[must_use]
#[derive(Debug)]
pub struct Minimap2Mapper<'i> {
  buf: Minimap2Buffer,
  idx: &'i Minimap2Index,
}

impl<'i> Minimap2Mapper<'i> {
  pub fn new(idx: &'i Minimap2Index) -> Result<Self, Report> {
    Ok(Self {
      buf: Minimap2Buffer::new()?,
      idx,
    })
  }

  pub fn run_map(&mut self, seq: &str, name: &str) -> Result<Minimap2Result, Report> {
    Minimap2Result::new(seq, name, self.idx, &mut self.buf)
  }
}

impl Minimap2Result {
  pub fn new(seq: &str, name: &str, idx: &Minimap2Index, buf: &mut Minimap2Buffer) -> Result<Self, Report> {
    let regs = Minimap2RawRegs::new(seq, name, idx, buf)?;

    dbg!(&regs.as_slice());

    let regs = regs
      .as_slice()
      .iter()
      .map(Minimap2Reg::from_raw)
      .collect::<Result<Vec<Minimap2Reg>, Report>>()?;

    Ok(Minimap2Result { regs })
  }
}

#[must_use]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Minimap2Reg {
  pub id: i32,     // ID for internal uses (see also parent below)
  pub cnt: i32,    // number of minimizers; if on the reverse strand
  pub rid: i32,    // reference index; if this is an alignment from inversion rescue
  pub score: i32,  // DP alignment score
  pub qs: i32,     // query start
  pub qe: i32,     // query end
  pub rs: i32,     // reference start
  pub re: i32,     // reference end
  pub parent: i32, // parent==id if primary
  pub subsc: i32,  // best alternate mapping score
  pub as_: i32,    // offset in the a[] array (for internal uses only)
  pub mlen: i32,   // seeded exact match length
  pub blen: i32,   // seeded alignment block length
  pub n_sub: i32,  // number of suboptimal mappings
  pub score0: i32, // initial chaining score (before chain merging/spliting)
  pub hash: u32,
  pub div: OrderedFloat<f32>,
  pub p: Option<MinimapExtra>,
  pub mapq: u32,            // map quality
  pub split: u32,           // split
  pub rev: u32,             // reverse
  pub inv: u32,             // inversion
  pub sam_pri: u32,         // sam primary
  pub proper_frag: u32,     // proper fragment
  pub pe_thru: u32,         // paired-end through
  pub seg_split: u32,       // segment split
  pub seg_id: u32,          // segment ID
  pub split_inv: u32,       // split inversion
  pub is_alt: u32,          // is alternate
  pub strand_retained: u32, // strand retained
  pub dummy: u32,           // dummy
}

#[must_use]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MinimapExtra {
  pub capacity: u32, // the capacity of cigar[]
  pub dp_score: i32, // DP score
  pub dp_max: i32,   // score of the max-scoring segment
  pub dp_max2: i32,  // score of the best alternate mappings
  pub n_cigar: u32,  // number of cigar operations in cigar[]
  pub cigar: String,
  pub n_ambi: u32,       // number of ambiguous bases
  pub trans_strand: u32, // transcript strand: 0 for unknown, 1 for +, 2 for -
}

impl Minimap2Reg {
  pub fn from_raw(reg: &mm_reg1_t) -> Result<Self, Report> {
    // SAFETY: dereferencing a raw pointer
    let p = unsafe { reg.p.as_ref() };
    let p = p.map(MinimapExtra::from_raw).transpose()?;
    Ok(Self {
      id: reg.id,
      cnt: reg.cnt,
      rid: reg.rid,
      score: reg.score,
      qs: reg.qs,
      qe: reg.qe,
      rs: reg.rs,
      re: reg.re,
      parent: reg.parent,
      subsc: reg.subsc,
      as_: reg.as_,
      mlen: reg.mlen,
      blen: reg.blen,
      n_sub: reg.n_sub,
      score0: reg.score0,
      hash: reg.hash,
      div: OrderedFloat(reg.div),
      p,
      mapq: reg.mapq(),
      split: reg.split(),
      rev: reg.rev(),
      inv: reg.inv(),
      sam_pri: reg.sam_pri(),
      proper_frag: reg.proper_frag(),
      pe_thru: reg.pe_thru(),
      seg_split: reg.seg_split(),
      seg_id: reg.seg_id(),
      split_inv: reg.split_inv(),
      is_alt: reg.is_alt(),
      strand_retained: reg.strand_retained(),
      dummy: reg.dummy(),
    })
  }
}

impl MinimapExtra {
  pub fn from_raw(p: &mm_extra_t) -> Result<Self, Report> {
    dbg!(&p);
    let n_cigar = p.n_cigar as usize;

    // SAFETY: dereferencing raw pointer and calling unsafe function
    let cigar = unsafe { p.cigar.as_slice(n_cigar) };

    let cigar = cigar
      .iter()
      .map(|&op| {
        let len = op >> 4;
        let op = (op & 0xf) as usize;
        let ch = char::from(MM_CIGAR_CHARS[op]);
        format!("{len}{ch}")
      })
      .join("");

    Ok(Self {
      capacity: p.capacity,
      dp_score: p.dp_score,
      dp_max: p.dp_max,
      dp_max2: p.dp_max2,
      n_cigar: p.n_cigar,
      cigar,
      n_ambi: p.n_ambi(),
      trans_strand: p.trans_strand(),
    })
  }
}

// impl Minimap2Reg {
//   pub fn from(reg: &mm_reg1_t) -> Self {
//     // SAFETY: dereferencing a raw pointer
//     let p = unsafe { reg.p.as_ref() };
//     let p = p.map(MinimapExtra::from);
//     Self {
//       id: reg.id,
//       cnt: reg.cnt,
//       rid: reg.rid,
//       score: reg.score,
//       qs: reg.qs,
//       qe: reg.qe,
//       rs: reg.rs,
//       re: reg.re,
//       parent: reg.parent,
//       subsc: reg.subsc,
//       as_: reg.as_,
//       mlen: reg.mlen,
//       blen: reg.blen,
//       n_sub: reg.n_sub,
//       score0: reg.score0,
//       hash: reg.hash,
//       div: OrderedFloat(reg.div),
//       p,
//       mapq: reg.mapq(),
//       split: reg.split(),
//       rev: reg.rev(),
//       inv: reg.inv(),
//       sam_pri: reg.sam_pri(),
//       proper_frag: reg.proper_frag(),
//       pe_thru: reg.pe_thru(),
//       seg_split: reg.seg_split(),
//       seg_id: reg.seg_id(),
//       split_inv: reg.split_inv(),
//       is_alt: reg.is_alt(),
//       strand_retained: reg.strand_retained(),
//       dummy: reg.dummy(),
//     }
//   }
// }
//
// impl MinimapExtra {
//   pub fn from(p: &mm_extra_t) -> Self {
//     dbg!(&p);
//     let n_cigar = p.n_cigar as usize;
//
//     // SAFETY: dereferencing raw pointer and calling unsafe function
//     let cigar = unsafe { p.cigar.as_slice(n_cigar) };
//
//     let cigar = cigar
//       .iter()
//       .map(|&op| {
//         let len = op >> 4;
//         let op = (op & 0xf) as usize;
//         let ch = char::from(MM_CIGAR_CHARS[op]);
//         format!("{len}{ch}")
//       })
//       .join("");
//
//     Self {
//       capacity: p.capacity,
//       dp_score: p.dp_score,
//       dp_max: p.dp_max,
//       dp_max2: p.dp_max2,
//       n_cigar: p.n_cigar,
//       cigar,
//       n_ambi: p.n_ambi(),
//       trans_strand: p.trans_strand(),
//     }
//   }
// }

// impl From<&mm_reg1_t> for Minimap2Reg {
//   fn from(reg: &mm_reg1_t) -> Self {
//     // SAFETY: dereferencing a raw pointer
//     let p = unsafe { reg.p.as_ref() };
//     let p = p.map(MinimapExtra::from);
//     Self {
//       id: reg.id,
//       cnt: reg.cnt,
//       rid: reg.rid,
//       score: reg.score,
//       qs: reg.qs,
//       qe: reg.qe,
//       rs: reg.rs,
//       re: reg.re,
//       parent: reg.parent,
//       subsc: reg.subsc,
//       as_: reg.as_,
//       mlen: reg.mlen,
//       blen: reg.blen,
//       n_sub: reg.n_sub,
//       score0: reg.score0,
//       hash: reg.hash,
//       div: OrderedFloat(reg.div),
//       p,
//       mapq: reg.mapq(),
//       split: reg.split(),
//       rev: reg.rev(),
//       inv: reg.inv(),
//       sam_pri: reg.sam_pri(),
//       proper_frag: reg.proper_frag(),
//       pe_thru: reg.pe_thru(),
//       seg_split: reg.seg_split(),
//       seg_id: reg.seg_id(),
//       split_inv: reg.split_inv(),
//       is_alt: reg.is_alt(),
//       strand_retained: reg.strand_retained(),
//       dummy: reg.dummy(),
//     }
//   }
// }
//
// impl From<&mm_extra_t> for MinimapExtra {
//   fn from(p: &mm_extra_t) -> Self {
//     dbg!(&p);
//     let n_cigar = p.n_cigar as usize;
//
//     // SAFETY: dereferencing raw pointer and calling unsafe function
//     let cigar = unsafe { p.cigar.as_slice(n_cigar) };
//
//     let cigar = cigar
//       .iter()
//       .map(|&op| {
//         let len = op >> 4;
//         let op = (op & 0xf) as usize;
//         let ch = char::from(MM_CIGAR_CHARS[op]);
//         format!("{len}{ch}")
//       })
//       .join("");
//
//     Self {
//       capacity: p.capacity,
//       dp_score: p.dp_score,
//       dp_max: p.dp_max,
//       dp_max2: p.dp_max2,
//       n_cigar: p.n_cigar,
//       cigar,
//       n_ambi: p.n_ambi(),
//       trans_strand: p.trans_strand(),
//     }
//   }
// }

// #[must_use]
// #[derive(Clone, Debug, Eq, PartialEq)]
// pub struct Minimap2Reg {
//   pub id: i32,
//   pub cnt: i32,
//   pub rid: i32,
//   pub score: i32,
//   pub qs: i32,
//   pub qe: i32,
//   pub rs: i32,
//   pub re: i32,
//   pub parent: i32,
//   pub subsc: i32,
//   pub as_: i32,
//   pub mlen: i32,
//   pub blen: i32,
//   pub n_sub: i32,
//   pub score0: i32,
//   pub hash: u32,
//   pub div: OrderedFloat<f32>,
//   pub cigar: String,
// }
//
// impl Minimap2Reg {
//   fn from_raw(r: &mm_reg1_t) -> Result<Self, Report> {
//     let mm_reg1_t {
//       id,
//       cnt,
//       rid,
//       score,
//       qs,
//       qe,
//       rs,
//       re,
//       parent,
//       subsc,
//       as_,
//       mlen,
//       blen,
//       n_sub,
//       score0,
//       hash,
//       div,
//       p,
//       ..
//     } = r;
//
//     // let mm_extra_t {
//     //   capacity,
//     //   dp_score,
//     //   dp_max,
//     //   dp_max2,
//     //   dp_max0,
//     //   n_cigar,
//     //   cigar,
//     //   ..
//     // } = p;
//
//     // #[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
//     // pub struct ExtractedHit {
//     //   pub hit: Hit,
//     //   pub new_block_id: BlockId,
//     //   pub is_anchor: bool,
//     //   pub orientation: Strand,
//     // }
//     //
//     // /// Pairwise homologous alignment between two sequences.
//     // #[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
//     // pub struct Alignment {
//     //   pub qry: Hit,
//     //   pub reff: Hit,
//     //   pub matches: usize,
//     //   pub length: usize,
//     //   pub quality: usize,
//     //   pub orientation: Strand,
//     //
//     //   pub new_block_id: Option<BlockId>, // FIXME: it looks like this does not belong here and is a "partially-initialized object" anti-pattern
//     //   pub anchor_block: Option<AnchorBlock>, // FIXME: it looks like this does not belong here and is a "partially-initialized object" anti-pattern
//     //
//     //   #[serde(serialize_with = "serde_serialize_cigar")]
//     //   #[serde(deserialize_with = "serde_deserialize_cigar")]
//     //   pub cigar: Cigar, // TODO: We probably want Edits here instead?
//     //
//     //   pub divergence: Option<f64>,
//     //   pub align: Option<f64>,
//     // }
//
//     // SAFETY: dereferencing raw pointer
//     let p = unsafe { r.p.as_ref() };
//     let cigar = if let Some(p) = p {
//       let n_cigar = p.n_cigar as usize;
//
//       // SAFETY: dereferencing raw pointer and calling unsafe function
//       let cigar = unsafe { p.cigar.as_slice(n_cigar) };
//
//       cigar
//         .iter()
//         .map(|&op| {
//           let len = op >> 4;
//           let op = (op & 0xf) as usize;
//           let ch = char::from(MM_CIGAR_CHARS[op]);
//           format!("{len}{ch}")
//         })
//         .join("")
//     } else {
//       return Err(eyre!("minimap2: cigar is null"));
//     };
//
//     Ok(Self {
//       id: *id,
//       cnt: *cnt,
//       rid: *rid,
//       score: *score,
//       qs: *qs,
//       qe: *qe,
//       rs: *rs,
//       re: *re,
//       parent: *parent,
//       subsc: *subsc,
//       as_: *as_,
//       mlen: *mlen,
//       blen: *blen,
//       n_sub: *n_sub,
//       score0: *score0,
//       hash: *hash,
//       div: OrderedFloat(*div),
//       cigar,
//     })
//   }
// }

// RAII wrapper around array of `mm_reg1_t`
#[must_use]
#[derive(Clone, Debug)]
pub struct Minimap2RawRegs {
  regs: *mut mm_reg1_t,
  n_regs: usize,
}

impl Minimap2RawRegs {
  pub fn new(seq: &str, name: &str, idx: &Minimap2Index, buf: &mut Minimap2Buffer) -> Result<Self, Report> {
    let mi: *const mm_idx_t = idx.get();
    let map_opt = idx.options().map_opt();

    let l_seq: c_int = seq.len() as c_int;

    let seq = CString::new(seq)?;
    let seq = seq.as_ptr();

    let name = CString::new(name)?;
    let name = name.as_ptr();

    // let seq: *const c_char = to_c_char_ptr(seq)?;
    // let name: *const c_char = to_c_char_ptr(name)?;

    let mut n_regs: c_int = 0;

    // SAFETY: calling unsafe function. The `mm_map()` allocates memory for results array and the `.p` field each of
    // its elements' using `malloc()`. This memory needs to be deallocated using `free()`.
    let regs: *mut mm_reg1_t = unsafe { mm_map(mi, l_seq, seq, &mut n_regs, buf.get_mut(), map_opt, name) };

    Ok(Self {
      regs,
      n_regs: n_regs as usize,
    })
  }

  pub fn as_slice(&self) -> &[mm_reg1_t] {
    if self.regs.is_null() || self.n_regs == 0 {
      return &[];
    }
    // SAFETY: converting pointer to slice
    unsafe { from_raw_parts(self.regs, self.n_regs) }
  }
}

impl Drop for Minimap2RawRegs {
  fn drop(&mut self) {
    if self.regs.is_null() || self.n_regs == 0 {
      return;
    }
    for reg in self.as_slice() {
      if !reg.p.is_null() {
        // SAFETY: deallocating memory
        unsafe { free(reg.p.cast()) };
      }
    }
    // SAFETY: deallocating memory
    unsafe { free(self.regs.cast()) };
  }
}
