use crate::make_internal_report;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use eyre::Report;
use std::collections::{BTreeMap, HashMap};
use std::hash::{Hash, Hasher};
use twox_hash::XxHash64;

/// Number of hex digits used to render the consensus hash in a canonical name.
const HASH_WIDTH: usize = 16;

/// Number of decimal digits used to render the block id in a canonical name.
/// `usize::MAX` is 20 digits, so this never truncates.
const ID_WIDTH: usize = 20;

/// Content-derived ordering key for a block: hash of its consensus sequence, tie-broken by block id.
///
/// The tie-break only engages for blocks whose consensus is byte-identical, in which case the
/// alignment between them is symmetric and the choice cannot bias the result.
pub type BlockKey = (u64, BlockId);

/// Content-derived naming and ordering for the blocks handed to an alignment backend.
///
/// `BlockId`s are assigned from the position of a record in the input (see `Pangraph::singleton`),
/// so using them to name sequences leaks the input order into the aligner. `minimap2` run with
/// `-X` keeps only one direction of each pair, chosen by `strcmp` of the sequence names, which
/// makes the query/reference roles - and therefore the merged consensus - depend on the order the
/// FASTA files were listed on the command line.
///
/// This type replaces those names with a key derived from the block consensus, so that both the
/// order blocks are fed to the aligner and the names they carry are a pure function of the block
/// content. Names are zero-padded to a fixed width, so byte-wise (`strcmp`) comparison agrees with
/// the numeric `(consensus_hash, block_id)` order that [`BlockNames::canonical_order`] uses.
pub struct BlockNames {
  names: BTreeMap<BlockId, String>,
  ids: HashMap<String, BlockId>,
  keys: BTreeMap<BlockId, BlockKey>,
  order: Vec<BlockId>,
}

/// Hashes a block's consensus sequence.
///
/// Content-derived, and therefore independent of how `BlockId`s were assigned.
pub fn consensus_hash(block: &PangraphBlock) -> u64 {
  let mut hasher = XxHash64::with_seed(0);
  block.consensus().hash(&mut hasher);
  hasher.finish()
}

/// Computes the content-derived ordering key of a block.
fn consensus_key(id: BlockId, block: &PangraphBlock) -> BlockKey {
  (consensus_hash(block), id)
}

/// Renders a canonical, fixed-width name for a block.
fn canonical_name((hash, id): BlockKey) -> String {
  format!("{hash:0HASH_WIDTH$x}_{:0ID_WIDTH$}", id.0)
}

impl BlockNames {
  /// Builds the canonical naming and ordering for a set of blocks.
  pub fn from_blocks(blocks: &BTreeMap<BlockId, PangraphBlock>) -> Self {
    let keys: BTreeMap<BlockId, BlockKey> = blocks
      .iter()
      .map(|(&id, block)| (id, consensus_key(id, block)))
      .collect();

    let mut order: Vec<BlockId> = keys.keys().copied().collect();
    order.sort_unstable_by_key(|id| keys[id]);

    let names: BTreeMap<BlockId, String> = keys.iter().map(|(&id, &key)| (id, canonical_name(key))).collect();
    let ids: HashMap<String, BlockId> = names.iter().map(|(&id, name)| (name.clone(), id)).collect();

    Self {
      names,
      ids,
      keys,
      order,
    }
  }

  /// Blocks in canonical (content-derived) order.
  pub fn canonical_order(&self) -> impl Iterator<Item = BlockId> + '_ {
    self.order.iter().copied()
  }

  /// Canonical name of a block, as handed to the alignment backend.
  pub fn name(&self, id: BlockId) -> &str {
    &self.names[&id]
  }

  /// Recovers the block a canonical name refers to.
  ///
  /// Fails if the aligner returned a name we never fed it, which would otherwise surface as a
  /// silently wrong block id.
  pub fn id_of(&self, name: &str) -> Result<BlockId, Report> {
    self
      .ids
      .get(name)
      .copied()
      .ok_or_else(|| make_internal_report!("Aligner returned unknown sequence name '{name}'"))
  }

  /// Content-derived sort key of a block, for order-independent tie-breaking.
  pub fn sort_key(&self, id: BlockId) -> BlockKey {
    self.keys[&id]
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::pangraph_node::NodeId;
  use crate::representation::seq::Seq;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  /// Builds a block map from `(id, consensus)` pairs.
  fn blocks_of(spec: &[(usize, &str)]) -> BTreeMap<BlockId, PangraphBlock> {
    spec
      .iter()
      .map(|&(id, seq)| {
        let bid = BlockId(id);
        (bid, PangraphBlock::from_consensus(Seq::from_str(seq), bid, NodeId(id)))
      })
      .collect()
  }

  #[rstest]
  fn canonical_names_are_fixed_width_and_unique() {
    let blocks = blocks_of(&[(0, "ACGT"), (1, "TTTT"), (12345, "GGGG")]);
    let names = BlockNames::from_blocks(&blocks);

    let rendered = blocks.keys().map(|&id| names.name(id).to_owned()).collect_vec();
    assert!(rendered.iter().all(|n| n.len() == HASH_WIDTH + 1 + ID_WIDTH));
    assert_eq!(rendered.iter().unique().count(), rendered.len());
  }

  #[rstest]
  fn canonical_order_matches_lexicographic_name_order() {
    let blocks = blocks_of(&[(0, "ACGT"), (1, "TTTT"), (2, "GGGG"), (3, "CCCC")]);
    let names = BlockNames::from_blocks(&blocks);

    let by_order = names.canonical_order().map(|id| names.name(id)).collect_vec();
    let sorted = by_order.iter().copied().sorted().collect_vec();
    assert_eq!(by_order, sorted);
  }

  /// The property the fix exists for: permuting block ids over the same consensus set must not
  /// change the order or the relative naming the aligner sees.
  #[rstest]
  fn canonical_order_is_invariant_under_id_permutation() {
    let seqs = ["ACGTACGT", "TTTTGGGG", "GGGGCCCC", "CACACACA"];

    let forward = blocks_of(&seqs.iter().enumerate().map(|(i, s)| (i, *s)).collect_vec());
    let reversed = blocks_of(
      &seqs
        .iter()
        .enumerate()
        .map(|(i, s)| (seqs.len() - 1 - i, *s))
        .collect_vec(),
    );

    let n_fwd = BlockNames::from_blocks(&forward);
    let n_rev = BlockNames::from_blocks(&reversed);

    // The sequences appear in the same order regardless of how ids were assigned.
    let seq_of = |blocks: &BTreeMap<BlockId, PangraphBlock>, id: BlockId| blocks[&id].consensus().as_str().to_owned();
    let order_fwd = n_fwd.canonical_order().map(|id| seq_of(&forward, id)).collect_vec();
    let order_rev = n_rev.canonical_order().map(|id| seq_of(&reversed, id)).collect_vec();
    assert_eq!(order_fwd, order_rev);
  }

  #[rstest]
  fn names_round_trip_to_block_ids() {
    let blocks = blocks_of(&[(0, "ACGT"), (7, "TTTT")]);
    let names = BlockNames::from_blocks(&blocks);

    for &id in blocks.keys() {
      assert_eq!(names.id_of(names.name(id)).unwrap(), id);
    }
  }

  #[rstest]
  fn unknown_name_is_rejected() {
    let blocks = blocks_of(&[(0, "ACGT")]);
    let names = BlockNames::from_blocks(&blocks);
    let err = names.id_of("not-a-name").unwrap_err();
    assert!(err.to_string().contains("unknown sequence name"));
  }

  /// Blocks sharing a consensus must still receive distinct names: identical names would make
  /// `minimap2` treat the pair as a self-comparison and discard the diagonal anchors that carry
  /// their (perfect) alignment, so the two blocks could never merge.
  #[rstest]
  fn identical_consensus_blocks_get_distinct_names() {
    let blocks = blocks_of(&[(0, "ACGTACGT"), (1, "ACGTACGT")]);
    let names = BlockNames::from_blocks(&blocks);

    assert_ne!(names.name(BlockId(0)), names.name(BlockId(1)));
    assert_eq!(names.sort_key(BlockId(0)).0, names.sort_key(BlockId(1)).0);
  }
}
