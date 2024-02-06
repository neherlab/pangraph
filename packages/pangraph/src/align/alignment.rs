use crate::pangraph::pangraph::Strand;
use noodles::sam::record::Cigar;
use serde::{Deserialize, Serialize};

/// Hit is one side of a pairwise alignment between homologous sequences.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Hit {
  pub name: String,
  pub length: usize,
  pub start: usize,
  pub stop: usize,
  pub seq: Option<Vec<u8>>,
}

/// Alignment is a pairwise homologous alignment between two sequences.
#[derive(Clone, Debug, Default)]
pub struct Alignment {
  pub qry: Hit,
  pub reff: Hit,
  pub matches: usize,
  pub length: usize,
  pub quality: usize,
  pub orientation: Strand,
  pub cigar: Cigar,
  pub divergence: Option<f64>,
  pub align: Option<f64>,
}
