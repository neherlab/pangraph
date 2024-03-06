use crate::align::bam::cigar::parse_cigar_str;
use crate::pangraph::strand::Strand;
use color_eyre::{Section, SectionExt};
use eyre::WrapErr;
use noodles::sam::record::Cigar;
use serde::de::Error;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::str::FromStr;

/// Hit is one side of a pairwise alignment between homologous sequences.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct Hit {
  pub name: String,
  pub length: usize, // TODO: having all of {start, stop, length, seq.length} seems redundant. Some of these can be functions.
  pub start: usize,
  pub stop: usize,
}

/// Alignment is a pairwise homologous alignment between two sequences.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
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
