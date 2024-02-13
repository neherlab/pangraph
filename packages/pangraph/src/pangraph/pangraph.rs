use crate::io::fasta::FastaRecord;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_reader_to_string;
use crate::io::json::json_parse;
use crate::utils::id::random_id;
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::Path;
use std::str::FromStr;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Pangraph {
  pub paths: Vec<PangraphPath>,
  pub blocks: Vec<PangraphBlock>,
}

impl Pangraph {
  pub fn singleton(fasta: FastaRecord, circular: bool) -> Self {
    let block = PangraphBlock::new(fasta.seq);
    let path = PangraphPath::new(fasta.seq_name, block.id, circular);
    Self {
      blocks: vec![block],
      paths: vec![path],
    }
  }

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

  pub fn sequences(&self) {}

  pub fn names(&self) {}

  pub fn sequence_map(&self) {}
}

impl FromStr for Pangraph {
  type Err = Report;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    json_parse(s).wrap_err("When parsing Pangraph JSON contents")
  }
}

#[derive(Clone, Serialize, Deserialize, Debug, SmartDefault, PartialEq, Eq)]
pub enum Strand {
  #[default]
  Forward,
  Reverse,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphNode {
  block_id: usize,
  strand: Strand,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphPath {
  pub name: String,
  pub nodes: Vec<PangraphNode>,
  pub offset: Option<isize>,
  pub circular: bool,
  pub position: Vec<usize>,
}

impl PangraphPath {
  pub fn new(name: String, block_id: usize, circular: bool) -> Self {
    Self {
      name,
      nodes: vec![PangraphNode {
        block_id,
        strand: Strand::default(), // FIXME: should we assume forward strand here?
      }],
      offset: circular.then_some(0),
      circular,
      position: vec![],
    }
  }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphBlock {
  pub id: usize,
  pub sequence: String,
  pub gaps: BTreeMap<String, usize>,
  pub mutate: Vec<PangraphMutate>,
  pub insert: Vec<PangraphInsert>,
  pub delete: Vec<PangraphDelete>,
  pub positions: Vec<PangraphPositions>,
}

impl PangraphBlock {
  pub fn new(sequence: String) -> Self {
    Self {
      id: random_id(),
      sequence,
      gaps: BTreeMap::default(),
      mutate: vec![],
      insert: vec![],
      delete: vec![],
      positions: vec![],
    }
  }
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
