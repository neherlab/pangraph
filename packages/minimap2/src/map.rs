#![allow(unsafe_code)]

use crate::buf::Minimap2Buffer;
use crate::index::Minimap2Index;
use eyre::{eyre, Report};
use itertools::Itertools;
use minimap2_sys::{free, mm_idx_t, mm_map, mm_reg1_t, MM_CIGAR_CHARS};
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
  pub id: i32,
  pub cnt: i32,
  pub rid: i32,
  pub score: i32,
  pub qs: i32,
  pub qe: i32,
  pub rs: i32,
  pub re: i32,
  pub parent: i32,
  pub subsc: i32,
  pub as_: i32,
  pub mlen: i32,
  pub blen: i32,
  pub n_sub: i32,
  pub score0: i32,
  pub hash: u32,
  pub div: OrderedFloat<f32>,
  pub cigar: String,
}

// pub capacity: u32,
// pub dp_score: i32,
// pub dp_max: i32,
// pub dp_max2: i32,
// pub dp_max0: i32,
// pub _bitfield_align_1: [u32; 0],
// pub _bitfield_1: __BindgenBitfieldUnit<[u8; 4usize]>,
// pub n_cigar: u32,
// pub cigar: __IncompleteArrayField<u32>,

impl Minimap2Reg {
  fn from_raw(r: &mm_reg1_t) -> Result<Self, Report> {
    let mm_reg1_t {
      id,
      cnt,
      rid,
      score,
      qs,
      qe,
      rs,
      re,
      parent,
      subsc,
      as_,
      mlen,
      blen,
      n_sub,
      score0,
      hash,
      div,
      p,
      ..
    } = r;

    // let mm_extra_t {
    //   capacity,
    //   dp_score,
    //   dp_max,
    //   dp_max2,
    //   dp_max0,
    //   n_cigar,
    //   cigar,
    //   ..
    // } = p;

    // #[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
    // pub struct ExtractedHit {
    //   pub hit: Hit,
    //   pub new_block_id: BlockId,
    //   pub is_anchor: bool,
    //   pub orientation: Strand,
    // }
    //
    // /// Pairwise homologous alignment between two sequences.
    // #[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
    // pub struct Alignment {
    //   pub qry: Hit,
    //   pub reff: Hit,
    //   pub matches: usize,
    //   pub length: usize,
    //   pub quality: usize,
    //   pub orientation: Strand,
    //
    //   pub new_block_id: Option<BlockId>, // FIXME: it looks like this does not belong here and is a "partially-initialized object" anti-pattern
    //   pub anchor_block: Option<AnchorBlock>, // FIXME: it looks like this does not belong here and is a "partially-initialized object" anti-pattern
    //
    //   #[serde(serialize_with = "serde_serialize_cigar")]
    //   #[serde(deserialize_with = "serde_deserialize_cigar")]
    //   pub cigar: Cigar, // TODO: We probably want Edits here instead?
    //
    //   pub divergence: Option<f64>,
    //   pub align: Option<f64>,
    // }

    // SAFETY: dereferencing raw pointer
    let p = unsafe { r.p.as_ref() };
    let cigar = if let Some(p) = p {
      let n_cigar = p.n_cigar as usize;

      // SAFETY: dereferencing raw pointer and calling unsafe function
      let cigar = unsafe { p.cigar.as_slice(n_cigar) };

      cigar
        .iter()
        .map(|&op| {
          let len = op >> 4;
          let op = (op & 0xf) as usize;
          let ch = char::from(MM_CIGAR_CHARS[op]);
          format!("{len}{ch}")
        })
        .join("")
    } else {
      return Err(eyre!("minimap2: cigar is null"));
    };

    Ok(Self {
      id: *id,
      cnt: *cnt,
      rid: *rid,
      score: *score,
      qs: *qs,
      qe: *qe,
      rs: *rs,
      re: *re,
      parent: *parent,
      subsc: *subsc,
      as_: *as_,
      mlen: *mlen,
      blen: *blen,
      n_sub: *n_sub,
      score0: *score0,
      hash: *hash,
      div: OrderedFloat(*div),
      cigar,
    })
  }
}

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
