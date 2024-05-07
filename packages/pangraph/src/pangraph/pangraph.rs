use crate::io::fasta::FastaRecord;
use crate::io::file::open_file_or_stdin;
use crate::io::fs::read_reader_to_string;
use crate::io::json::json_parse;
use crate::pangraph::pangraph_block::PangraphBlock;
use crate::pangraph::pangraph_node::PangraphNode;
use crate::pangraph::strand::Strand;
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
    let block = PangraphBlock::from_consensus(fasta.seq);
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

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PangraphPath {
  pub name: String,
  pub nodes: Vec<PangraphNode>,
  // pub offset: Option<isize>,
  pub circular: bool,
  // pub position: Vec<usize>,
}

impl PangraphPath {
  pub fn new(name: String, block_id: usize, circular: bool) -> Self {
    Self {
      name,
      nodes: vec![PangraphNode {
        id: 0,
        block_id,
        path_id: 0,
        strand: Strand::default(), // FIXME: should we assume forward strand here?
        position: (0, 0),
      }],
      // offset: circular.then_some(0),
      circular,
      // position: vec![],
    }
  }
}
