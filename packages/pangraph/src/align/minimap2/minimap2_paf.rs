use crate::align::alignment::{Alignment, Hit};
use crate::align::bam::cigar::parse_cigar_str;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::strand::Strand;
use csv::ReaderBuilder as CsvReaderBuilder;
use eyre::Report;
use serde::Deserialize;
use std::io::Cursor;

/// Represents one row in the PAF file emitted by mmseqs
#[derive(Clone, Debug, Deserialize)]
pub struct MinimapPafTsvRecord {
  /* 01 */ query: BlockId,
  /* 02 */ qlen: usize,
  /* 03 */ qstart: usize,
  /* 04 */ qend: usize,
  /* 05 */ strand: Strand,
  /* 06 */ target: BlockId,
  /* 07 */ tlen: usize,
  /* 08 */ tstart: usize,
  /* 09 */ tend: usize,
  /* 10 */ nident: usize,
  /* 11 */ alnlen: usize,
  /* 12 */ mapq: usize,
  /* 13 */ nm: String,
  /* 14 */ ms: String,
  /* 15 */ ascore: String,
  /* 16 */ nn: String,
  /* 17 */ tp: String,
  /* 18 */ cm: String,
  /* 19 */ s1: String,
  /* 20 */ s2: String,
  /* 21 */ de: String,
  /* 22 */ rl: String,
  /* 23 */ cg: String,
}

#[allow(clippy::multiple_inherent_impl)]
impl Alignment {
  pub fn from_minimap_paf_str(paf_str: impl AsRef<str>) -> Result<Vec<Self>, Report> {
    let mut rdr = CsvReaderBuilder::new()
      .delimiter(b'\t')
      .has_headers(false)
      .from_reader(Cursor::new(paf_str.as_ref()));

    rdr
      .deserialize()
      .map(|paf| {
        let paf: MinimapPafTsvRecord = paf?;

        // strip divergence prefix and then parse as float
        let div = paf.de.strip_prefix("de:f:").unwrap_or_default().parse::<f64>()?;

        // strip prefix from alignment score and then parse
        let aln_score = paf.ascore.strip_prefix("AS:i:").unwrap_or_default().parse::<f64>()?;

        // strip prefix from cigar string and then parse
        let cigar = paf.cg.strip_prefix("cg:Z:").unwrap_or_default();
        let cigar = parse_cigar_str(cigar)?;

        Ok(Alignment {
          qry: Hit::new(paf.query, paf.qlen, (paf.qstart, paf.qend)),
          reff: Hit::new(paf.target, paf.tlen, (paf.tstart, paf.tend)),
          matches: paf.nident,
          length: paf.alnlen,
          quality: paf.mapq,
          orientation: paf.strand,
          new_block_id: None, // FIXME: initialize?
          anchor_block: None, // FIXME: initialize?
          cigar,
          divergence: Some(div),
          align: Some(aln_score),
        })
      })
      .collect()
  }
}

#[cfg(test)]
mod tests {

  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  // TODO: add test cases for parsing multiple entries at once

  #[rstest]
  fn test_minimap2_paf_parse_one_forward() {
    // forward alignment
    let paf_content = "qry_0	998	0	996	+	ref_0	1000	0	998	969	998	60	NM:i:29	ms:i:845	AS:i:845	nn:i:0	tp:A:P	cm:i:145	s1:i:808	s2:i:0	de:f:0.0291	rl:i:0	cg:Z:545M1D225M1D226M";
    let aln = vec![Alignment {
      qry: Hit::new(BlockId(1), 998, (0, 996)),
      reff: Hit::new(BlockId(2), 1000, (0, 998)),
      matches: 969,
      length: 998,
      quality: 60,
      orientation: Strand::Forward,
      new_block_id: None,
      anchor_block: None,
      cigar: parse_cigar_str("545M1D225M1D226M").unwrap(),
      divergence: Some(0.0291),
      align: Some(845.),
    }];
    assert_eq!(Alignment::from_minimap_paf_str(paf_content).unwrap(), aln);
  }

  #[rstest]
  fn test_minimap2_paf_parse_one_reverse() {
    // reverse alignment
    let paf_content = "qry_3	997	0	980	-	ref_3	1000	18	1000	965	982	60	NM:i:17	ms:i:889	AS:i:889	nn:i:0	tp:A:P	cm:i:151	s1:i:815	s2:i:0	de:f:0.0173	rl:i:0	cg:Z:124M1D416M1D440M";
    let aln = vec![Alignment {
      qry: Hit::new(BlockId(3), 997, (0, 980)),
      reff: Hit::new(BlockId(4), 1000, (18, 1000)),
      matches: 965,
      length: 982,
      quality: 60,
      orientation: Strand::Reverse,
      new_block_id: None,
      anchor_block: None,
      cigar: parse_cigar_str("124M1D416M1D440M").unwrap(),
      divergence: Some(0.0173),
      align: Some(889.0),
    }];
    assert_eq!(Alignment::from_minimap_paf_str(paf_content).unwrap(), aln);
  }
}
