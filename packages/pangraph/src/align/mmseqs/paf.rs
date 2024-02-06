use crate::align::alignment::{Alignment, Hit};
use crate::align::bam::cigar::parse_cigar_str;
use crate::pangraph::pangraph::Strand;
use crate::utils::collections::remove_exactly_one;
use csv::ReaderBuilder as CsvReaderBuilder;
use eyre::Report;
use serde::Deserialize;
use serde_aux::serde_introspection::serde_introspect;
use std::io::Cursor;

/// Represents one row in the PAF file emitted by mmseqs
#[derive(Clone, Debug, Deserialize)]
pub struct PafTsvRecord {
  /* 01 */ query: String,
  /* 02 */ qlen: usize,
  /* 03 */ qstart: usize,
  /* 04 */ qend: usize,
  /* 05 */ empty: String,
  /* 06 */ target: String,
  /* 07 */ tlen: usize,
  /* 08 */ tstart: usize,
  /* 09 */ tend: usize,
  /* 10 */ nident: usize,
  /* 11 */ alnlen: usize,
  /* 12 */ bits: usize,
  /* 13 */ cigar: String,
  /* 14 */ fident: f64,
  /* 15 */ raw: f64,
}

impl PafTsvRecord {
  pub fn fields_names() -> &'static [&'static str] {
    serde_introspect::<PafTsvRecord>()
  }
}

impl Alignment {
  pub fn from_paf_str(paf_str: impl AsRef<str>) -> Result<Self, Report> {
    let paf = {
      let mut rdr = CsvReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(Cursor::new(paf_str.as_ref()));
      let records: Result<Vec<PafTsvRecord>, _> = rdr.deserialize().collect();
      remove_exactly_one(records?)?
    };

    let (qstart, qend, strand) = {
      let q1 = paf.qstart;
      let q2 = paf.qend;
      if q1 < q2 {
        (q1, q2, Strand::Forward)
      } else {
        (q2, q1, Strand::Reverse)
      }
    };

    Ok(Alignment {
      qry: Hit {
        name: paf.query,
        length: paf.qlen,
        start: qstart,
        stop: qend,
        seq: None,
      },
      reff: Hit {
        name: paf.target,
        length: paf.tlen,
        start: paf.tstart,
        stop: paf.tend,
        seq: None,
      },
      matches: paf.nident,
      length: paf.alnlen,
      quality: paf.bits,
      orientation: strand,
      cigar: parse_cigar_str(paf.cigar)?,
      divergence: Some(1.0 - paf.fident),
      align: Some(paf.raw),
    })
  }
}
