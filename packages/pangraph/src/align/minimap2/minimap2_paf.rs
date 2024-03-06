use crate::align::alignment::{Alignment, Hit};
use crate::align::bam::cigar::parse_cigar_str;
use crate::pangraph::strand::Strand;
use csv::ReaderBuilder as CsvReaderBuilder;
use eyre::Report;
use serde::Deserialize;
use serde_aux::serde_introspection::serde_introspect;
use std::io::Cursor;

/// Represents one row in the PAF file emitted by mmseqs
#[derive(Clone, Debug, Deserialize)]
pub struct MinimapPafTsvRecord {
  /* 01 */ query: String,
  /* 02 */ qlen: usize,
  /* 03 */ qstart: usize,
  /* 04 */ qend: usize,
  /* 05 */ strand: String,
  /* 06 */ target: String,
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

impl MinimapPafTsvRecord {
  pub fn fields_names() -> &'static [&'static str] {
    serde_introspect::<MinimapPafTsvRecord>()
  }
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

        let strand = match paf.strand.as_str() {
          "+" => Strand::Forward,
          "-" => Strand::Reverse,
          _ => return Err(eyre::eyre!("Invalid strand")),
        };

        // strip divergence prefix and then parse as float
        let div = paf.de.strip_prefix("de:f:").unwrap_or_default().parse::<f64>()?;

        // strip prefix from alignment score and then parse
        let aln_score = paf.ascore.strip_prefix("AS:i:").unwrap_or_default().parse::<f64>()?;

        // strip prefix from cigar string and then parse
        let cigar = paf.cg.strip_prefix("cg:Z:").unwrap_or_default();
        let cigar = parse_cigar_str(cigar)?;

        Ok(Alignment {
          qry: Hit {
            name: paf.query,
            length: paf.qlen,
            start: paf.qstart,
            stop: paf.qend,
          },
          reff: Hit {
            name: paf.target,
            length: paf.tlen,
            start: paf.tstart,
            stop: paf.tend,
          },
          matches: paf.nident,
          length: paf.alnlen,
          quality: paf.mapq,
          orientation: strand,
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
  use crate::o;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  // TODO: add test cases for parsing multiple entries at once

  #[rstest]
  fn test_minimap2_paf_parse_one_forward() {
    // forward alignment
    let paf_content = "qry_0	998	0	996	+	ref_0	1000	0	998	969	998	60	NM:i:29	ms:i:845	AS:i:845	nn:i:0	tp:A:P	cm:i:145	s1:i:808	s2:i:0	de:f:0.0291	rl:i:0	cg:Z:545M1D225M1D226M";
    let aln = vec![Alignment {
      qry: Hit {
        name: o!("qry_0"),
        length: 998,
        start: 0,
        stop: 996,
      },
      reff: Hit {
        name: o!("ref_0"),
        length: 1000,
        start: 0,
        stop: 998,
      },
      matches: 969,
      length: 998,
      quality: 60,
      orientation: Strand::Forward,
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
      qry: Hit {
        name: o!("qry_3"),
        length: 997,
        start: 0,
        stop: 980,
      },
      reff: Hit {
        name: o!("ref_3"),
        length: 1000,
        start: 18,
        stop: 1000,
      },
      matches: 965,
      length: 982,
      quality: 60,
      orientation: Strand::Reverse,
      cigar: parse_cigar_str("124M1D416M1D440M").unwrap(),
      divergence: Some(0.0173),
      align: Some(889.0),
    }];
    assert_eq!(Alignment::from_minimap_paf_str(paf_content).unwrap(), aln);
  }
}
