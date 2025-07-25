use crate::io::fasta::FastaRecord;
use crate::io::json::{JsonPretty, json_write_str};
use crate::io::seq::reverse_complement;
use crate::pangraph::edits::{Del, Edit, Ins, Sub};
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_node::NodeId;
use crate::pangraph::pangraph_path::PathId;
use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;
use crate::utils::collections::has_duplicates;
use crate::utils::interval::positions_to_intervals;
use derive_more::{Display, From};
use eyre::{Report, WrapErr};
use getset::{CopyGetters, Getters};
use itertools::Itertools;
use maplit::btreemap;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::hash::Hash;

#[derive(
  Copy, Clone, Debug, Display, From, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize, JsonSchema,
)]
pub struct BlockId(pub usize);

impl BlockId {
  pub fn from_str(block_id_str: &String) -> Result<Self, Report> {
    let block_id = block_id_str
      .parse::<usize>()
      .wrap_err_with(|| format!("When parsing Block ID: expected unsigned integer, but got '{block_id_str}'"))?;
    Ok(Self(block_id))
  }
}

#[must_use]
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash, Getters, CopyGetters, JsonSchema)]
pub struct PangraphBlock {
  #[get_copy = "pub"]
  id: BlockId,
  consensus: Seq,
  alignments: BTreeMap<NodeId, Edit>,
}

impl PangraphBlock {
  pub fn from_consensus(consensus: impl Into<Seq>, block_id: BlockId, nid: NodeId) -> Self {
    PangraphBlock::new(block_id, consensus, btreemap! {nid => Edit::empty()})
  }

  pub fn new(block_id: BlockId, consensus: impl Into<Seq>, alignments: BTreeMap<NodeId, Edit>) -> Self {
    let consensus = consensus.into();
    let id = block_id;
    Self {
      id,
      consensus,
      alignments,
    }
  }

  pub fn reverse_complement(&self) -> Result<Self, Report> {
    let rev_cons = reverse_complement(&self.consensus)?;

    let len = self.consensus_len();

    let rev_aln = self
      .alignments()
      .iter()
      .map(|(&nid, e)| Ok((nid, e.reverse_complement(len)?)))
      .collect::<Result<_, Report>>()?;

    Ok(Self::new(self.id(), rev_cons, rev_aln))
  }

  pub fn depth(&self) -> usize {
    self.alignments.len()
  }

  pub fn consensus(&self) -> &Seq {
    &self.consensus
  }

  pub fn consensus_mut(&mut self) -> &mut Seq {
    &mut self.consensus
  }

  pub fn set_consensus(&mut self, consensus: Seq) {
    self.consensus = consensus;
  }

  pub fn set_consensus_char(&mut self, pos: usize, c: AsciiChar) {
    self.consensus[pos] = c;
  }

  pub fn consensus_len(&self) -> usize {
    self.consensus.len()
  }

  pub fn unaligned_len_for_edit(&self, edits: &Edit) -> usize {
    let total_dels: usize = edits.dels.iter().map(|del| del.len).sum();
    let total_inss: usize = edits.inss.iter().map(|ins| ins.seq.len()).sum();
    self.consensus_len() + total_inss - total_dels
  }

  pub fn unaligned_len_for_node(&self, node_id: &NodeId) -> usize {
    self.unaligned_len_for_edit(&self.alignments[node_id])
  }

  pub fn alignment(&self, nid: NodeId) -> &Edit {
    &self.alignments[&nid]
  }

  pub fn alignment_mut(&mut self, nid: NodeId) -> &mut Edit {
    self.alignments.get_mut(&nid).unwrap()
  }

  pub fn alignment_keys(&self) -> BTreeSet<NodeId> {
    self.alignments.keys().copied().collect()
  }

  pub fn alignments(&self) -> &BTreeMap<NodeId, Edit> {
    &self.alignments
  }

  pub fn alignments_mut(&mut self) -> &mut BTreeMap<NodeId, Edit> {
    &mut self.alignments
  }

  pub fn alignment_insert(&mut self, node_id: NodeId, edit: Edit) -> Option<Edit> {
    self.alignments.insert(node_id, edit)
  }

  pub fn alignment_remove(&mut self, node_id: NodeId) -> Edit {
    self.alignments.remove(&node_id).unwrap()
  }

  pub fn isolates<'a>(&'a self, graph: &'a Pangraph) -> impl Iterator<Item = PathId> + 'a {
    self.alignments.keys().map(|&node_id| graph.nodes[&node_id].path_id())
  }

  pub fn is_duplicated(&self, graph: &Pangraph) -> bool {
    has_duplicates(self.isolates(graph))
  }

  pub fn sequences<'a>(
    &'a self,
    graph: &'a Pangraph,
    aligned: bool,
    record_naming: RecordNaming,
  ) -> impl Iterator<Item = Result<FastaRecord, Report>> + 'a {
    self.alignments().iter().map(move |(node_id, edits)| {
      let (id, desc) = match record_naming {
        RecordNaming::Node => {
          let node = &graph.nodes[node_id];
          let block_id = node.block_id();
          let path_name = &graph.paths[&node.path_id()].name().as_ref().unwrap();
          let (start, end) = node.position();
          let strand = node.strand();

          let meta = json_write_str(
            &json! ({
              "path_name": path_name,
              "block_id": block_id,
              "start": start,
              "end": end,
              "strand": strand,
            }),
            JsonPretty(false),
          )
          .unwrap();

          let id = node_id.to_string();
          let descr = Some(meta);

          (id, descr)
        },
        RecordNaming::Path => {
          let path_id = graph.nodes[node_id].path_id();
          let path = &graph.paths[&path_id];
          let id = path.name.as_ref().map_or_else(|| path_id.to_string(), |p| p.to_owned());
          let desc = path.desc.clone();
          (id, desc)
        },
      };

      let seq = if aligned {
        edits.apply_aligned(self.consensus())
      } else {
        edits.apply(self.consensus())
      }?;

      Ok(FastaRecord {
        seq_name: id,
        desc,
        seq,
        index: 0,
      })
    })
  }

  /// Finds all majority edits (substitutions, deletions, insertions) in this block
  pub fn find_majority_edits(&self) -> Edit {
    Edit::new(
      self.find_majority_insertions(),
      self.find_majority_deletions(),
      self.find_majority_substitutions(),
    )
  }

  /// Helper method to check if a count represents a majority
  fn is_majority(&self, count: usize) -> bool {
    count > self.depth() / 2
  }

  /// Finds majority substitutions in this block
  pub fn find_majority_substitutions(&self) -> Vec<Sub> {
    let mut substitutions: Vec<_> = self
      .alignments()
      .values()
      .flat_map(|edit| &edit.subs)
      .map(|sub| (sub.pos, sub.alt))
      .into_group_map()
      .into_iter()
      .filter_map(|(pos, alts)| {
        let (alt, count) = alts.into_iter().counts().into_iter().max_by_key(|(_, count)| *count)?;
        self.is_majority(count).then_some(Sub::new(pos, alt))
      })
      .collect();

    substitutions.sort_by_key(|sub| sub.pos);
    substitutions
  }

  /// Finds majority deletions in this block
  pub fn find_majority_deletions(&self) -> Vec<Del> {
    let majority_positions: Vec<usize> = self
      .alignments()
      .values()
      .flat_map(|edit| edit.dels.iter().flat_map(|del| del.range()))
      .counts()
      .into_iter()
      .filter_map(|(pos, count)| self.is_majority(count).then_some(pos))
      .collect();

    positions_to_intervals(&majority_positions)
      .into_iter()
      .map(|interval| Del::new(interval.start, interval.len()))
      .collect()
  }

  /// Finds majority insertions in this block
  pub fn find_majority_insertions(&self) -> Vec<Ins> {
    let mut insertions: Vec<_> = self
      .alignments()
      .values()
      .flat_map(|edit| &edit.inss)
      .map(|insertion| (insertion.pos, insertion.seq.clone()))
      .counts()
      .into_iter()
      .filter_map(|((pos, seq), count)| self.is_majority(count).then_some(Ins::new(pos, seq)))
      .collect();

    insertions.sort_by_key(|ins| ins.pos);
    insertions
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum RecordNaming {
  Node,
  Path,
}
