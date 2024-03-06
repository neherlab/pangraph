use crate::align::alignment::{Alignment, Hit};
use crate::align::bam::cigar::parse_cigar_str;
use crate::pangraph::strand::Strand;
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

    let (qstart, qend, strand) = order_range(paf.qstart, paf.qend);

    let (tstart, tend, _) = order_range(paf.tstart, paf.tend);

    Ok(Alignment {
      qry: Hit {
        name: paf.query,
        length: paf.qlen,
        start: qstart,
        stop: qend,
      },
      reff: Hit {
        name: paf.target,
        length: paf.tlen,
        start: tstart,
        stop: tend,
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

fn order_range(start: usize, end: usize) -> (usize, usize, Strand) {
  if start < end {
    (start, end, Strand::Forward)
  } else {
    (end, start, Strand::Reverse)
  }
}

#[cfg(test)]
mod tests {

  use super::*;
  use crate::o;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_paf_parse_forward() {
    // forward alignment
    let paf_content = "qry	507	1	497	-	ref	500	500	24	440	508	622	67M10D18M20I235M10I22M1I5M1D119M	0.866	693";
    let aln = Alignment {
      qry: Hit {
        name: o!("qry"),
        length: 507,
        start: 1,
        stop: 497,
      },
      reff: Hit {
        name: o!("ref"),
        length: 500,
        start: 24,
        stop: 500,
      },
      matches: 440,
      length: 508,
      quality: 622,
      orientation: Strand::Forward,
      cigar: parse_cigar_str("67M10D18M20I235M10I22M1I5M1D119M").unwrap(),
      divergence: Some(0.134),
      align: Some(693.0),
    };
    assert_eq!(Alignment::from_paf_str(paf_content).unwrap(), aln);
  }

  #[rstest]
  fn test_paf_parse_reverse() {
    // reverse alignment
    let paf_content = "rev_qry	507	507	11	-	ref	500	500	24	440	508	622	67M10D18M20I235M10I22M1I5M1D119M	0.866	693";
    let aln = Alignment {
      qry: Hit {
        name: o!("rev_qry"),
        length: 507,
        start: 11,
        stop: 507,
      },
      reff: Hit {
        name: o!("ref"),
        length: 500,
        start: 24,
        stop: 500,
      },
      matches: 440,
      length: 508,
      quality: 622,
      orientation: Strand::Reverse,
      cigar: parse_cigar_str("67M10D18M20I235M10I22M1I5M1D119M").unwrap(),
      divergence: Some(0.134),
      align: Some(693.0),
    };
    assert_eq!(Alignment::from_paf_str(paf_content).unwrap(), aln);
  }
}
