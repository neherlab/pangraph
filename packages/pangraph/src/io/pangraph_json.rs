use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_reader_to_string;
use crate::io::json::json_parse;
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;
use std::str::FromStr;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphJson {
  pub paths: Vec<PangraphPath>,
  pub blocks: Vec<PangraphBlock>,
}

impl FromStr for PangraphJson {
  type Err = Report;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    json_parse(s).wrap_err("When parsing Pangraph JSON contents")
  }
}

impl PangraphJson {
  pub fn from_path<P: AsRef<Path>>(filepath: &Option<P>) -> Result<Self, Report> {
    let reader = open_file_or_stdin(filepath)?;
    let data = read_reader_to_string(reader).wrap_err("When reading Pangraph JSON")?;
    Self::from_str(&data).wrap_err("When parsing Pangraph JSON")
  }

  pub fn to_string_pretty(&self) -> Result<String, Report> {
    let mut tree_str = serde_json::to_string_pretty(self)?;
    tree_str += "\n";
    Ok(tree_str)
  }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphPath {
  pub name: String,
  pub offset: isize,
  pub circular: bool,
  pub position: Vec<usize>,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphBlock {
  pub id: String,
  pub sequence: String,
  pub gaps: BTreeMap<String, usize>,
  pub mutate: Vec<PangraphMutate>,
  pub insert: Vec<PangraphInsert>,
  pub delete: Vec<PangraphDelete>,
  pub positions: Vec<PangraphPositions>,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphBlockEntry {
  pub name: String,
  pub number: usize,
  pub strand: bool,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphMutate {
  #[serde(flatten)]
  pub entry: PangraphBlockEntry,
  pub muts: Vec<Mutation>,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Mutation {
  pub pos: usize,
  pub nuc: char,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphInsert {
  #[serde(flatten)]
  pub entry: PangraphBlockEntry,
  pub inss: Vec<Insertion>,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Insertion;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphDelete {
  #[serde(flatten)]
  pub entry: PangraphBlockEntry,
  pub dels: Vec<Insertion>,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Deletion;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphPositions {
  #[serde(flatten)]
  pub entry: PangraphBlockEntry,
  pub poss: Vec<usize>,
}
