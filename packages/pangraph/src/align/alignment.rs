use crate::align::bam::cigar::parse_cigar_str;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::strand::Strand;
use crate::utils::interval::Interval;
use color_eyre::{Section, SectionExt};
use eyre::WrapErr;
use noodles::sam::record::Cigar;
use serde::de::Error;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::str::FromStr;

/// Hit is one side of a pairwise alignment between homologous sequences.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Hit {
  pub name: String,
  pub length: usize,
  pub interval: Interval,
}

impl Hit {
  pub fn new(name: impl Into<String>, length: usize, (start, end): (usize, usize)) -> Self {
    Self {
      name: name.into(),
      length,
      interval: Interval::new(start, end),
    }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Hit2 {
  pub name: BlockId,
  pub length: usize,
  pub interval: Interval,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ExtractedHit {
  pub hit: Hit2,
  pub new_block_id: BlockId,
  pub is_anchor: bool,
  pub orientation: bool,
}

/// Alignment is a pairwise homologous alignment between two sequences.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Alignment {
  pub qry: Hit,
  pub reff: Hit,
  pub matches: usize,
  pub length: usize,
  pub quality: usize,
  pub orientation: Strand,

  #[serde(serialize_with = "serde_serialize_cigar")]
  #[serde(deserialize_with = "serde_deserialize_cigar")]
  pub cigar: Cigar, // TODO: We probably want Edits here instead?

  pub divergence: Option<f64>,
  pub align: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Alignment2 {
  pub qry: Hit2,
  pub reff: Hit2,
  pub matches: usize,
  pub length: usize,
  pub quality: usize,
  pub orientation: bool,

  pub new_block_id: Option<BlockId>, // FIXME: it looks like this does not belong here and is a "partially-initialized object" anti-pattern
  pub anchor_block: Option<AnchorBlock>, // FIXME: it looks like this does not belong here and is a "partially-initialized object" anti-pattern

  #[serde(serialize_with = "serde_serialize_cigar")]
  #[serde(deserialize_with = "serde_deserialize_cigar")]
  pub cigar: Cigar, // TODO: We probably want Edits here instead?

  pub divergence: Option<f64>,
  pub align: Option<f64>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum AnchorBlock {
  Ref,
  Qry,
}

pub fn serde_serialize_cigar<S: Serializer>(cigar: &Cigar, s: S) -> Result<S::Ok, S::Error> {
  s.serialize_str(&cigar.to_string())
}

pub fn serde_deserialize_cigar<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Cigar, D::Error> {
  let s = String::deserialize(deserializer)?;
  let cigar = parse_cigar_str(&s)
    .wrap_err("When parsing Cigar string")
    .with_section(|| s.header("Cigar string:"))
    .map_err(Error::custom)?;
  Ok(cigar)
}
