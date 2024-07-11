#![allow(non_snake_case)]
#![allow(unsafe_code)]

use crate::buf::Minimap2Buffer;
use crate::index::Minimap2Index;
use eyre::{eyre, Report};
use itertools::Itertools;
use minimap2_sys::{free, mm_event_identity, mm_extra_t, mm_idx_t, mm_map, mm_reg1_t, MM_CIGAR_CHARS};
use ordered_float::OrderedFloat;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};
use std::slice::from_raw_parts;

#[must_use]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Minimap2Result {
  // This contains data in a format which is close to minimap2 internal format (`struct mm_reg1_t`)
  pub regs: Vec<Minimap2Reg>,

  /// This contains data in a digested format that is close to the 'PAF' output of minimap2 CLI
  pub pafs: Vec<Minimap2PafRow>,
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

  pub fn run_map(&mut self, seq: impl AsRef<str>, name: impl AsRef<str>) -> Result<Minimap2Result, Report> {
    Minimap2Result::new(seq, name, self.idx, &mut self.buf)
  }
}

impl Minimap2Result {
  pub fn new(
    seq: impl AsRef<str>,
    name: impl AsRef<str>,
    idx: &Minimap2Index,
    buf: &mut Minimap2Buffer,
  ) -> Result<Self, Report> {
    let seq = seq.as_ref();
    let name = name.as_ref();

    let regs = Minimap2RawRegs::new(seq, name, idx, buf)?;

    dbg!(&regs.as_slice());

    let pafs = regs
      .as_slice()
      .iter()
      .map(|reg| Minimap2PafRow::from_raw(seq, name, idx.get_ref()?, reg))
      .collect::<Result<Vec<Minimap2PafRow>, Report>>()?;

    let regs = regs
      .as_slice()
      .iter()
      .map(Minimap2Reg::from_raw)
      .collect::<Result<Vec<Minimap2Reg>, Report>>()?;

    Ok(Minimap2Result { regs, pafs })
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
  pub p: Option<Minimap2RegExtra>,
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
pub struct Minimap2RegExtra {
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
    let p = p.map(Minimap2RegExtra::from_raw).transpose()?;
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

impl Minimap2RegExtra {
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

#[derive(Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Minimap2PafRowSeq {
  /// (1) Query sequence name
  pub name: String,
  /// (2) Query sequence length
  pub len: usize,
  /// (3) Query start coordinate (0-based)
  pub start: i32,
  /// (4) Query end coordinate (0-based)
  pub end: i32,
}

/// Represents a row in the minimap2 flavor of PAF (Minimap2 Pairwise mApping) format.
/// PAF is a TAB-delimited text format with each line consisting of at least 12 fields.
/// For more details, see the full spec in the
/// ['Output' section of the man pages](https://lh3.github.io/minimap2/minimap2.html)
#[derive(Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Minimap2PafRow {
  /// (1-4) Query sequence
  pub q: Minimap2PafRowSeq,
  /// (6-9) Target sequence
  pub t: Minimap2PafRowSeq,
  /// (5) ‘+’ if query/target on the same strand; ‘-’ if opposite
  pub strand: char,
  /// (10) Number of matching bases in the mapping
  pub mlen: i32,
  /// (11) Number bases, including gaps, in the mapping
  pub blen: i32,
  /// (12) Mapping quality (0-255 with 255 for missing)
  pub mapq: u8,
  /// (13) tp:A Type of aln: P/primary, S/secondary and I,i/inversion
  pub tp: Option<char>,
  /// (14) cm:i Number of minimizers on the chain
  pub cm: Option<i32>,
  /// (15) s1:i Chaining score
  pub s1: Option<i32>,
  /// (16) s2:i Chaining score of the best secondary chain
  pub s2: Option<i32>,
  /// (17) NM:i Total number of mismatches and gaps in the alignment
  pub NM: Option<i32>,
  /// (18) MD:Z To generate the ref sequence in the alignment
  pub MD: Option<String>,
  /// (19) AS:i DP alignment score
  pub AS: Option<i32>,
  /// (20) SA:Z List of other supplementary alignments
  pub SA: Option<String>,
  /// (21) ms:i DP score of the max scoring segment in the alignment
  pub ms: Option<i32>,
  /// (22) nn:i Number of ambiguous bases in the alignment
  pub nn: Option<u32>,
  /// (23) ts:A Transcript strand (splice mode only)
  pub ts: Option<char>,
  /// (24) cg:Z CIGAR string (only in PAF)
  pub cg: Option<String>,
  /// (25) cs:Z Difference string
  pub cs: Option<String>,
  /// (26) dv:f Approximate per-base sequence divergence
  pub dv: Option<OrderedFloat<f32>>,
  /// (27) de:f Gap-compressed per-base sequence divergence
  pub de: Option<OrderedFloat<f64>>,
  /// (28) rl:i Length of query regions harboring repetitive seeds
  pub rl: Option<i32>,
}

pub fn c_char_to_str<'a>(ptr: *mut c_char) -> Result<&'a str, std::str::Utf8Error> {
  // SAFETY: Converting raw pointer to string
  unsafe { CStr::from_ptr(ptr) }.to_str()
}

impl Minimap2PafRow {
  fn from_raw(
    seq: impl AsRef<str>,
    name: impl AsRef<str>,
    idx: &mm_idx_t,
    res: &mm_reg1_t,
  ) -> Result<Minimap2PafRow, Report> {
    let query_seq = Minimap2PafRowSeq {
      name: name.as_ref().to_owned(),
      len: seq.as_ref().len(),
      start: res.qs,
      end: res.qe,
    };
    let strand = if res.rev() == 0 { '+' } else { '-' };

    if idx.n_seq > 0 && idx.seq.is_null() {
      return Err(eyre!("minimap2: unable to access query sequence in the result"));
    }

    // SAFETY: Converting raw pointer to slice
    let mi_seq = unsafe { from_raw_parts(idx.seq, idx.n_seq as usize) };

    let target_seq = Minimap2PafRowSeq {
      name: c_char_to_str(mi_seq[res.rid as usize].name)?.to_owned(),
      len: mi_seq[res.rid as usize].len as usize,
      start: res.rs,
      end: res.re,
    };

    let num_matching_bases = res.mlen;
    let num_bases_with_gaps = res.blen;
    let mapping_quality = res.mapq() as u8;

    let tp = Some(match (res.id == res.parent, res.inv()) {
      (true, 0) => 'P',
      (true, _) => 'I',
      (false, 0) => 'S',
      (false, _) => 'i',
    });

    let cm = Some(res.cnt);
    let s1 = Some(res.score);
    let s2 = (res.parent == res.id).then_some(res.subsc);

    // SAFETY: dereferencing raw pointer
    let p = unsafe { res.p.as_ref() };
    let NM = p.map(|p| res.blen - res.mlen + p.n_ambi() as i32);
    let AS = p.map(|p| p.dp_score);
    let ms = p.map(|p| p.dp_max);
    let nn = p.map(|p| p.n_ambi());

    let ts = p.and_then(|p| match p.trans_strand() {
      1 => Some('+'),
      2 => Some('-'),
      _ => None,
    });

    let dv = (0.0..=1.0).contains(&res.div).then_some(OrderedFloat(res.div));
    let de = p.map(|_| {
      // SAFETY: calling unsafe function
      let ei = unsafe { mm_event_identity(res) };
      OrderedFloat(1.0 - ei)
    });

    let extra = p.map(Minimap2RegExtra::from_raw).transpose()?;
    let cg = extra.map(|extra| extra.cigar);

    Ok(Minimap2PafRow {
      q: query_seq,
      t: target_seq,
      strand,
      mlen: num_matching_bases,
      blen: num_bases_with_gaps,
      mapq: mapping_quality,
      tp,
      cm,
      s1,
      s2,
      NM,
      MD: None, // TODO
      AS,
      SA: None, // TODO
      ms,
      nn,
      ts,
      cg,
      cs: None, // TODO
      dv,
      de,
      rl: None, // TODO
    })
  }
}

/// RAII wrapper around array of `mm_reg1_t`
#[must_use]
#[derive(Clone, Debug)]
pub struct Minimap2RawRegs {
  regs: *mut mm_reg1_t,
  n_regs: usize,
}

impl Minimap2RawRegs {
  pub fn new(
    seq: impl AsRef<str>,
    name: impl AsRef<str>,
    idx: &Minimap2Index,
    buf: &mut Minimap2Buffer,
  ) -> Result<Self, Report> {
    let mi: *const mm_idx_t = idx.get();
    let map_opt = idx.options().map_opt();

    let l_seq: c_int = seq.as_ref().len() as c_int;

    let seq = CString::new(seq.as_ref())?;
    let seq = seq.as_ptr();

    let name = CString::new(name.as_ref())?;
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
