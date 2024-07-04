use crate::align::alignment::{Alignment, Hit};
use crate::align::bam::cigar::parse_cigar_str;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::strand::Strand;
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

#[allow(clippy::multiple_inherent_impl)]
impl Alignment {
  pub fn from_paf_str(paf_str: impl AsRef<str>) -> Result<Vec<Self>, Report> {
    let mut rdr = CsvReaderBuilder::new()
      .delimiter(b'\t')
      .has_headers(false)
      .from_reader(Cursor::new(paf_str.as_ref()));

    rdr
      .deserialize()
      .map(|paf| {
        let paf: PafTsvRecord = paf?;

        let (qstart, qend, strand) = order_range(paf.qstart, paf.qend);

        let (tstart, tend, _) = order_range(paf.tstart, paf.tend);

        let block_id_qry = BlockId::from_str(&paf.query)?;
        let block_id_ref = BlockId::from_str(&paf.target)?;

        Ok(Alignment {
          qry: Hit::new(block_id_qry, paf.qlen, (qstart, qend)),
          reff: Hit::new(block_id_ref, paf.tlen, (tstart, tend)),
          matches: paf.nident,
          length: paf.alnlen,
          quality: paf.bits,
          orientation: strand,
          new_block_id: None, // FIXME: initialize?
          anchor_block: None, // FIXME: initialize?
          cigar: parse_cigar_str(paf.cigar)?,
          divergence: Some(1.0 - paf.fident),
          align: Some(paf.raw),
        })
      })
      .collect()
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
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_paf_parse_forward() {
    // forward alignment
    let paf_content = "qry	507	1	497	-	ref	500	500	24	440	508	622	67M10D18M20I235M10I22M1I5M1D119M	0.866	693";
    let aln = vec![Alignment {
      qry: Hit::new(BlockId(1), 507, (1, 497)),
      reff: Hit::new(BlockId(2), 500, (24, 500)),
      matches: 440,
      length: 508,
      quality: 622,
      orientation: Strand::Forward,
      new_block_id: None,
      anchor_block: None,
      cigar: parse_cigar_str("67M10D18M20I235M10I22M1I5M1D119M").unwrap(),
      divergence: Some(0.134),
      align: Some(693.0),
    }];
    assert_eq!(Alignment::from_paf_str(paf_content).unwrap(), aln);
  }

  #[rstest]
  fn test_paf_parse_reverse() {
    // reverse alignment
    let paf_content = "rev_qry	507	507	11	-	ref	500	500	24	440	508	622	67M10D18M20I235M10I22M1I5M1D119M	0.866	693";
    let aln = vec![Alignment {
      qry: Hit::new(BlockId(3), 507, (11, 507)),
      reff: Hit::new(BlockId(4), 500, (24, 500)),
      matches: 440,
      length: 508,
      quality: 622,
      orientation: Strand::Reverse,
      new_block_id: None,
      anchor_block: None,
      cigar: parse_cigar_str("67M10D18M20I235M10I22M1I5M1D119M").unwrap(),
      divergence: Some(0.134),
      align: Some(693.0),
    }];
    assert_eq!(Alignment::from_paf_str(paf_content).unwrap(), aln);
  }
}
