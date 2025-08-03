use crate::align::map_variations::{BandParameters, map_variations};
use crate::commands::build::build_args::PangraphBuildArgs;
use crate::io::fasta::FastaRecord;
use crate::io::json::{JsonPretty, json_write_str};
use crate::io::seq::reverse_complement;
use crate::make_internal_error;
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
  #[inline]
  pub fn is_majority(&self, count: usize) -> bool {
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

  /// Change a nucleotide in the consensus sequence at a specific position
  /// and update the alignments accordingly, without changing the block sequences.
  pub fn change_consensus_nucleotide_at_pos(&mut self, pos: usize, c: AsciiChar) -> Result<(), Report> {
    if pos >= self.consensus_len() {
      return make_internal_error!(
        "Position {pos} is out of bounds for consensus of length {}",
        self.consensus_len()
      );
    }

    // get the original character
    let original_char = self.consensus[pos];
    // check: the two must be different
    if original_char == c {
      return make_internal_error!(
        "Cannot change consensus character at position {pos} to '{c}' because it is already '{original_char}'"
      );
    }

    // update the consensus
    self.consensus[pos] = c;

    // Update the alignments
    self.alignments_mut().values_mut().try_for_each(|edit| {
      edit
        .reconcile_substitution_with_consensus(&Sub::new(pos, c), original_char)
        .wrap_err_with(|| format!("When reconciling substitution at position {pos} with character '{c}'"))
    })?;

    Ok(())
  }

  /// Applies a set of edits to the block's consensus sequence and re-aligns the sequences
  /// to the new consensus. Returns a new `PangraphBlock` object with the same BlockId.
  pub fn edit_consensus_and_realign(self, edits: &Edit, args: &PangraphBuildArgs) -> Result<Self, Report> {
    // apply the edits to the consensus
    let new_consensus = edits.apply(&self.consensus)?;
    debug_assert!(!new_consensus.is_empty(), "Consensus cannot be empty");

    // calculate alignment band parameters
    let band_params = BandParameters::from_edits(edits, self.consensus_len())?;

    // realign the block sequences to the new consensus
    let new_alignments = self
      .alignments()
      .iter()
      .map(|(&nid, edit)| {
        // reconstruct the alignment sequence
        let seq = edit.apply(&self.consensus)?;
        debug_assert!(!seq.is_empty(), "Aligned sequence cannot be empty");

        // calculate the alignment band parameters from the orgiginal alignment plus the displacement
        // given by the edits applied to the original consensus
        let old_band_params = BandParameters::from_edits(edit, self.consensus_len())?;
        let updated_band_params = BandParameters::new(
          old_band_params.mean_shift() - band_params.mean_shift(),
          old_band_params.band_width() + band_params.band_width(),
        );

        // re-align the sequence and returns the new set of edits
        let new_edits = map_variations(&new_consensus, &seq, updated_band_params, args)?;
        Ok((nid, new_edits))
      })
      .collect::<Result<BTreeMap<NodeId, Edit>, Report>>()?;

    // create a new block with the updated alignments
    Ok(PangraphBlock {
      id: self.id(),
      consensus: new_consensus,
      alignments: new_alignments,
    })
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum RecordNaming {
  Node,
  Path,
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_node::NodeId;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  fn s(pos: usize, alt: char) -> Sub {
    Sub::new(pos, alt)
  }

  fn e_subs(subs: Vec<Sub>) -> Edit {
    Edit::new(vec![], vec![], subs)
  }

  fn d(pos: usize, len: usize) -> Del {
    Del::new(pos, len)
  }

  fn e_dels(dels: Vec<Del>) -> Edit {
    Edit::new(vec![], dels, vec![])
  }

  fn i(pos: usize, seq: &str) -> Ins {
    Ins::new(pos, seq)
  }

  fn e_inss(inss: Vec<Ins>) -> Edit {
    Edit::new(inss, vec![], vec![])
  }

  #[test]
  fn test_find_majority_substitutions_single_node() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_subs(vec![s(0, 'G'), s(2, 'A')])
      },
    );
    let result = block.find_majority_substitutions();
    // Single node is always majority (1 > 0)
    assert_eq!(result, vec![s(0, 'G'), s(2, 'A')]);
  }

  #[test]
  fn test_find_majority_substitutions_no_majority() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_subs(vec![s(0, 'G')]),
        NodeId(2) => e_subs(vec![s(0, 'C')]),
        NodeId(3) => e_subs(vec![s(0, 'T')]),
      },
    );
    let result = block.find_majority_substitutions();
    // No substitution has majority (1 is not > 3/2 = 1)
    assert!(result.is_empty());
  }

  #[test]
  fn test_find_majority_substitutions_clear_majority() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_subs(vec![s(0, 'G'), s(2, 'A')]),
        NodeId(2) => e_subs(vec![s(0, 'G'), s(3, 'A')]),
        NodeId(3) => e_subs(vec![s(0, 'C'), s(2, 'A')]),
      },
    );
    let result = block.find_majority_substitutions();
    assert_eq!(result, vec![s(0, 'G'), s(2, 'A')]);
  }

  #[test]
  fn test_find_majority_substitutions_tie_no_majority() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_subs(vec![]),
        NodeId(2) => e_subs(vec![]),
        NodeId(3) => e_subs(vec![s(0, 'C')]),
        NodeId(4) => e_subs(vec![s(0, 'C')]),
      },
    );
    let result = block.find_majority_substitutions();
    assert!(result.is_empty());
  }

  #[test]
  fn test_find_majority_deletions_single_node() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCGAA",
      btreemap! {
        NodeId(1) => e_dels(vec![d(1, 2), d(4, 1)])
      },
    );
    let result = block.find_majority_deletions();
    // Single node is always majority (1 > 0)
    assert_eq!(result, vec![d(1, 2), d(4, 1)]);
  }

  #[test]
  fn test_find_majority_deletions_no_majority() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCGAA",
      btreemap! {
        NodeId(1) => e_dels(vec![d(0, 1)]),
        NodeId(2) => e_dels(vec![d(1, 1)]),
        NodeId(3) => e_dels(vec![d(2, 1)]),
      },
    );
    let result = block.find_majority_deletions();
    assert!(result.is_empty());
  }

  #[test]
  fn test_find_majority_deletions_clear_majority() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCGAA",
      btreemap! {
        NodeId(1) => e_dels(vec![d(1, 2), d(4, 1)]),
        NodeId(2) => e_dels(vec![d(1, 2), d(5, 1)]),
        NodeId(3) => e_dels(vec![d(0, 1), d(4, 1)]),
      },
    );
    let result = block.find_majority_deletions();
    assert_eq!(result, vec![d(1, 2), d(4, 1)]);
  }

  #[test]
  fn test_find_majority_deletions_overlapping_intervals() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCGAATT",
      btreemap! {
        NodeId(1) => e_dels(vec![d(1, 3)]),    // deletes positions 1,2,3
        NodeId(2) => e_dels(vec![d(2, 3)]),    // deletes positions 2,3,4
        NodeId(3) => e_dels(vec![d(3, 2)]),    // deletes positions 3,4
        NodeId(4) => e_dels(vec![d(6, 1)]),    // deletes position 6
        NodeId(5) => e_dels(vec![d(6, 2)]),    // deletes positions 6,7
      },
    );
    let result = block.find_majority_deletions();
    assert_eq!(result, vec![d(3, 1)]);
  }

  #[test]
  fn test_find_majority_deletions_contiguous_intervals() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCGAATT",
      btreemap! {
        NodeId(1) => e_dels(vec![d(1, 1), d(2, 1), d(3, 1)]),  // deletes 1,2,3 separately
        NodeId(2) => e_dels(vec![d(1, 3)]),                     // deletes 1,2,3 as one interval
        NodeId(3) => e_dels(vec![d(1, 1), d(2, 2)]),           // deletes 1, then 2,3
        NodeId(4) => e_dels(vec![d(5, 1)]),                     // deletes 5
        NodeId(5) => e_dels(vec![d(5, 1), d(6, 1)]),           // deletes 5,6 separately
      },
    );
    let result = block.find_majority_deletions();
    assert_eq!(result, vec![d(1, 3)]);
  }

  #[test]
  fn test_find_majority_insertions_empty_block() {
    let block = PangraphBlock::new(BlockId(1), "ATCG", btreemap! {});
    let result = block.find_majority_insertions();
    assert!(result.is_empty());
  }

  #[test]
  fn test_find_majority_insertions_single_node() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_inss(vec![i(1, "GG"), i(3, "AA")])
      },
    );
    let result = block.find_majority_insertions();
    // Single node is always majority (1 > 0)
    assert_eq!(result, vec![i(1, "GG"), i(3, "AA")]);
  }

  #[test]
  fn test_find_majority_insertions_no_majority() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_inss(vec![i(1, "A")]),
        NodeId(2) => e_inss(vec![i(1, "T")]),
        NodeId(3) => e_inss(vec![i(1, "G")]),
      },
    );
    let result = block.find_majority_insertions();
    // No insertion has majority (1 is not > 3/2 = 1)
    assert!(result.is_empty());
  }

  #[test]
  fn test_find_majority_insertions_clear_majority() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_inss(vec![i(1, "GGG"), i(3, "A")]),
        NodeId(2) => e_inss(vec![i(1, "GGG"), i(2, "TT")]),
        NodeId(3) => e_inss(vec![i(1, "CC"), i(3, "A")]),
      },
    );
    let result = block.find_majority_insertions();
    assert_eq!(result, vec![i(1, "GGG"), i(3, "A")]);
  }

  #[test]
  fn test_find_majority_insertions_exact_sequence_match() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_inss(vec![i(1, "ATG")]),
        NodeId(2) => e_inss(vec![i(1, "ATG")]),
        NodeId(3) => e_inss(vec![i(1, "ATG")]),
        NodeId(4) => e_inss(vec![i(1, "GTA")]),
        NodeId(5) => e_inss(vec![i(1, "GTA")]),
      },
    );
    let result = block.find_majority_insertions();
    assert_eq!(result, vec![i(1, "ATG")]);
  }

  #[test]
  fn test_find_majority_insertions_different_positions() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCGAA",
      btreemap! {
        NodeId(1) => e_inss(vec![i(0, "G"), i(2, "T"), i(4, "C")]),
        NodeId(2) => e_inss(vec![i(0, "G"), i(3, "A"), i(5, "T")]),
        NodeId(3) => e_inss(vec![i(1, "A"), i(2, "T"), i(4, "C")]),
        NodeId(4) => e_inss(vec![i(0, "C"), i(2, "T"), i(6, "G")]),
        NodeId(5) => e_inss(vec![i(0, "G"), i(3, "A"), i(4, "C")]),
      },
    );
    let result = block.find_majority_insertions();
    assert_eq!(result, vec![i(0, "G"), i(2, "T"), i(4, "C")]);
  }

  #[test]
  fn test_find_majority_insertions_tie_no_majority() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_inss(vec![]),
        NodeId(2) => e_inss(vec![]),
        NodeId(3) => e_inss(vec![i(1, "AA")]),
        NodeId(4) => e_inss(vec![i(1, "AA")]),
      },
    );
    let result = block.find_majority_insertions();
    assert!(result.is_empty());
  }

  #[test]
  fn test_find_majority_edits() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => Edit::empty(),
        NodeId(2) => Edit::empty(),
        NodeId(3) => Edit::empty(),
      },
    );
    let result = block.find_majority_edits();
    assert!(result.is_empty());
  }

  #[test]
  fn test_find_majority_edits_comprehensive() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCGAATT",
      btreemap! {
        NodeId(1) => Edit::new(vec![i(1, "GG"), i(4, "C")], vec![d(2, 1), d(6, 1)], vec![s(0, 'G'), s(5, 'C')]),
        NodeId(2) => Edit::new(vec![i(1, "GG"), i(3, "A")], vec![d(2, 1), d(7, 1)], vec![s(0, 'G'), s(5, 'T')]),
        NodeId(3) => Edit::new(vec![i(1, "AA"), i(4, "C")], vec![d(2, 1), d(6, 1)], vec![s(0, 'C'), s(5, 'C')]),
        NodeId(4) => Edit::new(vec![i(1, "GG"), i(4, "C")], vec![d(1, 1), d(6, 1)], vec![s(0, 'G'), s(4, 'A')]),
        NodeId(5) => Edit::new(vec![i(1, "GG"), i(4, "C")], vec![d(2, 1), d(5, 1)], vec![s(0, 'G'), s(5, 'C')]),
      },
    );
    let result = block.find_majority_edits();

    // Depth = 5, majority threshold = 5/2 = 2, so need > 2 = at least 3

    // Insertions:
    // Position 1: GG appears 4 times (majority), AA appears 1 time
    // Position 4: C appears 4 times (majority)
    // Position 3: A appears 1 time (not majority)

    // Deletions:
    // Position 2: deleted by 4 nodes (majority)
    // Position 6: deleted by 4 nodes (majority)
    // Position 1,5,7: deleted by 1 node each (not majority)

    // Substitutions:
    // Position 0: G appears 4 times (majority), C appears 1 time
    // Position 5: C appears 3 times (majority), T appears 1 time
    // Position 4: A appears 1 time (not majority)

    assert_eq!(result.inss, vec![i(1, "GG"), i(4, "C")]);
    assert_eq!(result.dels, vec![d(2, 1), d(6, 1)]);
    assert_eq!(result.subs, vec![s(0, 'G'), s(5, 'C')]);
  }

  #[test]
  fn test_change_consensus_nucleotide_at_pos_basic() {
    let mut block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => Edit::empty(),
        NodeId(2) => e_subs(vec![s(1, 'G'), s(2, 'C')]),
        NodeId(3) => e_subs(vec![s(1, 'A')]),
      },
    );

    let expected_block = PangraphBlock::new(
      BlockId(1),
      "AGCG",
      btreemap! {
        NodeId(1) => e_subs(vec![s(1, 'T')]),
        NodeId(2) => e_subs(vec![s(2, 'C')]),
        NodeId(3) => e_subs(vec![s(1, 'A')]),
      },
    );

    // Change position 1 from T to G
    block.change_consensus_nucleotide_at_pos(1, 'G'.into()).unwrap();
    assert_eq!(block, expected_block);
  }

  #[test]
  fn test_change_consensus_nucleotide_at_pos_with_deletion() {
    let mut block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => e_dels(vec![d(1, 2)]), // Node 1 deletes positions 1-2 (TC)
        NodeId(2) => Edit::empty(),
        NodeId(3) => e_subs(vec![s(1, 'A')]),
        NodeId(4) => e_subs(vec![s(1, 'G')]),
      },
    );

    let expected_block = PangraphBlock::new(
      BlockId(1),
      "AGCG",
      btreemap! {
        NodeId(1) => e_dels(vec![d(1, 2)]),
        NodeId(2) => e_subs(vec![s(1, 'T')]), // Node 2 should get a reversion substitution T at position 1
        NodeId(3) => e_subs(vec![s(1, 'A')]), // Node 3's substitution remains unchanged
        NodeId(4) => Edit::empty(), // Node 4's substitution is removed (now matches consensus)
      },
    );

    // Change position 1 from T to G
    block.change_consensus_nucleotide_at_pos(1, 'G'.into()).unwrap();
    assert_eq!(block, expected_block);
  }

  #[test]
  fn test_change_consensus_nucleotide_at_pos_out_of_bounds() {
    let mut block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => Edit::empty(),
      },
    );

    // Try to change position 4 (out of bounds for length 4)
    let result = block.change_consensus_nucleotide_at_pos(4, 'A'.into());
    assert!(result.is_err());
  }

  #[test]
  fn test_change_consensus_nucleotide_at_pos_same_character() {
    let mut block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => Edit::empty(),
      },
    );

    // Try to change position 1 to the same character (T)
    let result = block.change_consensus_nucleotide_at_pos(1, 'T'.into());
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("already"));
  }

  #[test]
  fn test_reverse_complement() {
    let block = PangraphBlock::new(
      BlockId(1),
      "ATCG",
      btreemap! {
        NodeId(1) => Edit::new(vec![i(1, "AA")], vec![d(2, 1)], vec![s(0, 'G')]),
        NodeId(2) => Edit::new(vec![], vec![], vec![s(1, 'G'), s(3, 'A')]),
        NodeId(3) => Edit::empty(),
      },
    );

    let expected_block = PangraphBlock::new(
      BlockId(1),
      "CGAT",
      btreemap! {
        NodeId(1) => Edit::new(vec![i(3, "TT")], vec![d(1, 1)], vec![s(3, 'C')]),
        NodeId(2) => Edit::new(vec![], vec![], vec![s(0, 'T'), s(2, 'C')]),
        NodeId(3) => Edit::empty(),
      },
    );

    let rev_block = block.reverse_complement().unwrap();
    assert_eq!(rev_block, expected_block);
  }
}
