use crate::align::map_variations::{BandParameters, map_variations};
use crate::commands::build::build_args::PangraphBuildArgs;
use crate::make_error;
use crate::pangraph::detach_unaligned::detach_unaligned_nodes;
use crate::pangraph::edits::Edit;
use crate::pangraph::edits::Sub;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_node::NodeId;
// use crate::reconsensus::remove_nodes::remove_emtpy_nodes;
use crate::reconsensus::remove_nodes::find_empty_nodes;
use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;
use eyre::Report;
use itertools::{Either, Itertools};
use rayon::prelude::*;
use std::collections::BTreeMap;

/// Applies the reconsensus operation to each updated block in the graph:
/// - updates the block consensus following a merge
/// - removes potentially empty nodes
///
/// Reconsensus function:
///   - for blocks that:
///     - originate from a new merger
///     - have depth > 2
///   - goes through each block and re-defines the consensus
///     - for mutations, changes any mutated position. Also update the alignment (this is straightforward)
///     - for any in/del present in > M/2 sites, appends it to the consensus
///   - if the consensus has updated indels, then re-aligns all the sequences to the new consensus
pub fn reconsensus_graph(
  graph: &mut Pangraph,
  ids_updated_blocks: &[BlockId],
  args: &PangraphBuildArgs,
) -> Result<(), Report> {
  // there should be no empty nodes in the graph
  debug_assert!(
    find_empty_nodes(graph, ids_updated_blocks).is_empty(),
    "Empty nodes found in the graph"
  );

  // Analyze blocks and determine which need realignment
  let (blocks_with_mutations_only, blocks_needing_realignment) =
    analyze_blocks_for_reconsensus(graph, ids_updated_blocks);

  // Apply mutation-only reconsensus (no realignment needed)
  blocks_with_mutations_only.into_iter().try_for_each(|block_id| {
    let block = graph.blocks.get_mut(&block_id).unwrap();
    let majority_edits = block.find_majority_edits();
    apply_mutation_reconsensus(block, &majority_edits.subs)
  })?;

  // Handle blocks requiring realignment
  if !blocks_needing_realignment.is_empty() {
    // Pop the realigned blocks from graph.blocks
    let mut realigned_blocks: Vec<_> = blocks_needing_realignment
      .iter()
      .filter_map(|block_id| graph.blocks.remove(block_id))
      .collect();

    // Apply full reconsensus with realignment
    realigned_blocks.iter_mut().try_for_each(|block| {
      let majority_edits = block.find_majority_edits();
      apply_full_reconsensus(block, &majority_edits, args)
    })?;

    // Apply detach_unaligned_nodes. This removes unaligned nodes and re-adds them to the list as new blocks.
    detach_unaligned_nodes(&mut realigned_blocks, &mut graph.nodes)?;

    // Re-add all the blocks (including potentially new singleton blocks) to graph.blocks
    realigned_blocks.into_iter().for_each(|block| {
      graph.blocks.insert(block.id(), block);
    });
  }

  Ok(())
}

/// Analyzes blocks to determine which need realignment vs. mutation-only reconsensus
fn analyze_blocks_for_reconsensus(graph: &Pangraph, block_ids: &[BlockId]) -> (Vec<BlockId>, Vec<BlockId>) {
  let (mutations_only, need_realignment): (Vec<_>, Vec<_>) = block_ids
    .iter()
    .filter_map(|&block_id| {
      let block = &graph.blocks[&block_id];
      let majority_edits = block.find_majority_edits();

      if majority_edits.has_indels() {
        Some(Either::Right(block_id))
      } else if majority_edits.has_subs() {
        Some(Either::Left(block_id))
      } else {
        None // Blocks with no variants are skipped
      }
    })
    .partition_map(|either| either);

  (mutations_only, need_realignment)
}

/// Updates alignment for a single mutation
fn update_alignment_for_mutation(
  edit: &mut Edit,
  pos: usize,
  alt: AsciiChar,
  original: AsciiChar,
) -> Result<(), Report> {
  let subs_at_pos: Vec<_> = edit.subs.iter().filter(|s| s.pos == pos).cloned().collect();

  match subs_at_pos.len() {
    0 => {
      // if the position is not in a deletion, append the mutation
      if !edit.dels.iter().any(|d| d.contains(pos)) {
        edit.subs.push(Sub::new(pos, original));
        edit.subs.sort_by_key(|s| s.pos);
      }
    },
    1 => {
      let s = &subs_at_pos[0];
      if s.alt == alt {
        edit.subs.retain(|sub| !(sub.pos == pos && sub.alt == alt));
      }
    },
    _ => {
      return make_error!(
        "At position {pos}: sequence states disagree: {:}",
        subs_at_pos
          .iter()
          .map(|sub| sub.alt.to_string())
          .collect_vec()
          .join(", ")
      );
    },
  }
  Ok(())
}

/// Applies only mutation reconsensus without realignment
fn apply_mutation_reconsensus(block: &mut PangraphBlock, subs: &[Sub]) -> Result<(), Report> {
  subs.iter().try_for_each(|sub| {
    let original = block.consensus()[sub.pos];

    // Update consensus
    block.set_consensus_char(sub.pos, sub.alt);

    // Update alignments
    block
      .alignments_mut()
      .values_mut()
      .try_for_each(|edit| update_alignment_for_mutation(edit, sub.pos, sub.alt, original))
  })
}

/// Applies full reconsensus including indels and realignment
fn apply_full_reconsensus(
  block: &mut PangraphBlock,
  majority_edits: &Edit,
  args: &PangraphBuildArgs,
) -> Result<(), Report> {
  // First apply mutations
  apply_mutation_reconsensus(block, &majority_edits.subs)?;

  // Then apply indels and realign if present
  if majority_edits.has_indels() {
    let consensus = block.consensus();
    let new_consensus = majority_edits.apply(consensus)?;

    let band_params = BandParameters::from_edits(majority_edits, block.consensus_len())?;

    debug_assert!(!new_consensus.is_empty(), "Consensus is empty after indels");

    update_block_consensus(block, &new_consensus, band_params, args)?;
  }

  Ok(())
}

/// Updates the consensus sequence of the block and re-aligns the sequences to the new consensus.
fn update_block_consensus(
  block: &mut PangraphBlock,
  consensus: &Seq,
  new_consensus_band_params: BandParameters,
  args: &PangraphBuildArgs,
) -> Result<(), Report> {
  // Reconstruct block sequences
  let seqs = block
    .alignments()
    .iter()
    .map(|(&nid, edit)| {
      let seq = edit.apply(block.consensus())?;
      let old_band_params = BandParameters::from_edits(edit, block.consensus_len())?;
      let updated_band_params = BandParameters::new(
        old_band_params.mean_shift() - new_consensus_band_params.mean_shift(),
        old_band_params.band_width() + new_consensus_band_params.band_width(),
      );
      Ok((nid, (seq, updated_band_params)))
    })
    .collect::<Result<BTreeMap<NodeId, (Seq, BandParameters)>, Report>>()?;

  // debug assets: all sequences are non-empty
  #[cfg(any(debug_assertions, test))]
  {
    if let Some((nid, (_seq, _))) = seqs.iter().find(|(_, (seq, _))| seq.is_empty()) {
      return make_error!(
        "node is empty!\nblock: {}\nnode: {}\nedits: {:?}\nconsensus: {}",
        block.id(),
        nid,
        block.alignments().get(nid),
        block.consensus()
      );
    }
  }

  // Re-align sequences
  let alignments = seqs
    .into_par_iter()
    .map(|(nid, (seq, band_params))| Ok((nid, map_variations(consensus, &seq, band_params, args)?)))
    .collect::<Result<_, Report>>()?;

  *block = PangraphBlock::new(block.id(), consensus, alignments);

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
  use crate::pangraph::pangraph_node::NodeId;
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::pangraph_path::{PangraphPath, PathId};
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use crate::utils::id::id;

  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  fn block_1() -> PangraphBlock {
    let consensus = "AGGACTTCGATCTATTCGGAGAA";
    //               0         1         2
    //               01234567890123456789012
    //    node 1)    .T...--..........A.....
    //    node 2)    .T...--...C......|.....
    //    node 3)    .T...--...C.....--.....
    //    node 4)    .C.......---.....A.....
    //    node 5)    .....|...........A.....
    //     L = 23, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![],                      vec![Del::new(5, 2)],     vec![Sub::new(1, 'T'), Sub::new(17, 'A')]),
      NodeId(2) => Edit::new(vec![Ins::new(17, "CAAT")],  vec![Del::new(5, 2)],     vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(3) => Edit::new(vec![],                      vec![Del::new(5, 2), Del::new(16,2)],     vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(4) => Edit::new(vec![],                      vec![Del::new(9, 3)],     vec![Sub::new(1, 'C'), Sub::new(17, 'A')]),
      NodeId(5) => Edit::new(vec![],                      vec![Del::new(5, 2)],     vec![Sub::new(17, 'A')]),
    };
    PangraphBlock::new(BlockId(0), consensus, aln)
  }

  fn block_1_mut_reconsensus() -> PangraphBlock {
    let consensus = "ATGACTTCGATCTATTCAGAGAA";
    //               0         1         2
    //               01234567890123456789012
    //    node 1)    .....--................
    //    node 2)    .....--...C......G|....
    //    node 3)    .....--...C.....--.....
    //    node 4)    .C.......---...........
    //    node 5)    .G...|.................
    //     L = 23, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![],                      vec![Del::new(5, 2)],     vec![]),
      NodeId(2) => Edit::new(vec![Ins::new(17, "CAAT")],  vec![Del::new(5, 2)],     vec![Sub::new(10, 'C'), Sub::new(17, 'G')]),
      NodeId(3) => Edit::new(vec![],                      vec![Del::new(5, 2), Del::new(16,2)],     vec![Sub::new(10, 'C')]),
      NodeId(4) => Edit::new(vec![],                      vec![Del::new(9, 3)],     vec![Sub::new(1, 'C')]),
      NodeId(5) => Edit::new(vec![],                      vec![Del::new(5, 2)],     vec![Sub::new(1, 'G')]),
    };
    PangraphBlock::new(BlockId(0), consensus, aln)
  }

  fn block_2() -> PangraphBlock {
    let consensus = "AGGACTTCGATCTATTCGGAGAA";
    //               0         1         2
    //               01234567890123456789012
    //    node 1)   |.T.|.--......|...A..-..
    //    node 2)   |.T...--...C..|......--.|
    //    node 3)    .T..----..C............|
    //    node 4)    .C.|.....---.....A.....|
    //    node 5)   |...|.........|...A.--..
    //     L = 23, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![Ins::new(0, "G"), Ins::new(3, "AA"), Ins::new(13, "AA")],   vec![Del::new(5, 2), Del::new(20, 1)],  vec![Sub::new(1, 'T'), Sub::new(17, 'A')]),
      NodeId(2) => Edit::new(vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")],  vec![Del::new(5, 2), Del::new(20, 2)],  vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(3) => Edit::new(vec![Ins::new(23, "TT")],                                        vec![Del::new(4, 4)],                   vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(4) => Edit::new(vec![Ins::new(3, "C"), Ins::new(23, "TT")],                      vec![Del::new(9, 3)],                   vec![Sub::new(1, 'C'), Sub::new(17, 'A')]),
      NodeId(5) => Edit::new(vec![Ins::new(0, "G"), Ins::new(3, "C"), Ins::new(13, "AA")],    vec![Del::new(19, 2)],                  vec![Sub::new(17, 'A')])
    };
    PangraphBlock::new(BlockId(0), consensus, aln)
  }

  fn block_3() -> PangraphBlock {
    let consensus = "GCCTCTTCCCGACCACGCGTTACAACATGGGACAGGCCTGCGCTTGAGGC";
    let aln = btreemap! {
        NodeId(1) => Edit::new(vec![],                                        vec![Del::new(19, 4)],  vec![Sub::new(5, 'A')]),
        NodeId(2) => Edit::new(vec![Ins::new(35, "AA"), Ins::new(50, "TT")],  vec![Del::new(20, 3)],  vec![Sub::new(5, 'A')]),
        NodeId(3) => Edit::new(vec![],                                        vec![],                 vec![Sub::new(14, 'G'), Sub::new(27, 'G')]),
        NodeId(4) => Edit::new(vec![Ins::new(50, "TT")],                      vec![Del::new(20, 3)],  vec![Sub::new(5, 'A')]),
        NodeId(5) => Edit::new(vec![Ins::new(50, "TT")],                      vec![],                 vec![])
    };
    PangraphBlock::new(BlockId(3), consensus, aln)
  }

  fn block_3_reconsensus() -> PangraphBlock {
    let consensus = "GCCTCATCCCGACCACGCGTAACATGGGACAGGCCTGCGCTTGAGGCTT";
    let aln = btreemap! {
        NodeId(1) => Edit::new(vec![],                    vec![Del::new(19, 1), Del::new(47, 2)],   vec![]),
        NodeId(2) => Edit::new(vec![Ins::new(32, "AA")],  vec![],                                   vec![]),
        NodeId(3) => Edit::new(vec![Ins::new(20, "TAC")], vec![Del::new(47, 2)],                    vec![Sub::new(5, 'T'), Sub::new(14, 'G'), Sub::new(24, 'G')]),
        NodeId(4) => Edit::new(vec![],                    vec![],                                   vec![]),
        NodeId(5) => Edit::new(vec![Ins::new(20, "TAC")], vec![],                                   vec![Sub::new(5, 'T')])
    };
    PangraphBlock::new(BlockId(3), consensus, aln)
  }

  fn block_mutations_only() -> PangraphBlock {
    let consensus = "ATGCGATCGATCGA";
    //               01234567890123
    //    node 1)    .C............  (mutation at pos 1: T->C)
    //    node 2)    .C............  (mutation at pos 1: T->C)
    //    node 3)    .C............  (mutation at pos 1: T->C)
    //    node 4)    ..........G...  (mutation at pos 10: C->G)
    //    node 5)    ..........G...  (mutation at pos 10: C->G)
    //     L = 14, N = 5
    // Position 1: T->C appears in 3/5 sequences (majority)
    // Position 10: C->G appears in 2/5 sequences (not majority)
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![], vec![], vec![Sub::new(1, 'C')]),
      NodeId(2) => Edit::new(vec![], vec![], vec![Sub::new(1, 'C')]),
      NodeId(3) => Edit::new(vec![], vec![], vec![Sub::new(1, 'C')]),
      NodeId(4) => Edit::new(vec![], vec![], vec![Sub::new(10, 'G')]),
      NodeId(5) => Edit::new(vec![], vec![], vec![Sub::new(10, 'G')]),
    };
    PangraphBlock::new(BlockId(10), consensus, aln)
  }

  fn block_mutations_only_after_reconsensus() -> PangraphBlock {
    let consensus = "ACGCGATCGATCGA";
    //               01234567890123
    //    node 1)    ..............  (no mutations after consensus update)
    //    node 2)    ..............  (no mutations after consensus update)
    //    node 3)    ..............  (no mutations after consensus update)
    //    node 4)    .T........G...  (gets T at pos 1, keeps G at pos 10)
    //    node 5)    .T........G...  (gets T at pos 1, keeps G at pos 10)
    //     L = 14, N = 5
    // After reconsensus: consensus[1] = 'C' (majority), original T becomes mutation in nodes 4,5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![], vec![], vec![]),
      NodeId(2) => Edit::new(vec![], vec![], vec![]),
      NodeId(3) => Edit::new(vec![], vec![], vec![]),
      NodeId(4) => Edit::new(vec![], vec![], vec![Sub::new(1, 'T'), Sub::new(10, 'G')]),
      NodeId(5) => Edit::new(vec![], vec![], vec![Sub::new(1, 'T'), Sub::new(10, 'G')]),
    };
    PangraphBlock::new(BlockId(10), consensus, aln)
  }

  #[test]
  fn test_reconsensus_mutations() {
    let mut block = block_1();
    let expected_block = block_1_mut_reconsensus();
    let subs = block.find_majority_substitutions();
    apply_mutation_reconsensus(&mut block, &subs).unwrap();
    assert_eq!(block, expected_block);
  }

  #[test]
  fn test_majority_deletions() {
    let dels: Vec<usize> = block_2()
      .find_majority_deletions()
      .iter()
      .flat_map(|d| d.range())
      .collect();
    assert_eq!(dels, vec![5, 6, 20]);
  }

  #[test]
  fn test_majority_insertions() {
    let ins = block_2().find_majority_insertions();
    assert_eq!(ins, vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")]);
  }

  #[test]
  fn test_apply_edits() {
    let consensus = "AGGACTTCGATCTATTCGGAGAA";
    let dels = vec![Del::new(5, 2), Del::new(20, 1)];
    let ins = vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")];
    let edit = Edit::new(ins, dels, vec![]);
    let cons = edit.apply(consensus).unwrap();
    assert_eq!(cons, "GAGGACCGATCTAAATTCGGAAATT");
  }

  #[test]
  fn test_reconsensus() {
    let mut block = block_3();
    let expected_block = block_3_reconsensus();

    let majority_edits = block.find_majority_edits();

    // Apply mutations
    apply_mutation_reconsensus(&mut block, &majority_edits.subs).unwrap();

    // Apply indels if present
    let has_indels = majority_edits.has_indels();
    if has_indels {
      let consensus = block.consensus();
      let new_consensus = majority_edits.apply(consensus).unwrap();
      let band_params = BandParameters::from_edits(&majority_edits, block.consensus_len()).unwrap();
      update_block_consensus(&mut block, &new_consensus, band_params, &PangraphBuildArgs::default()).unwrap();
    }

    assert!(has_indels); // This test expects realignment to have occurred
    assert_eq!(block.consensus(), expected_block.consensus());
    assert_eq!(block.alignments(), expected_block.alignments());
  }

  #[test]
  fn test_reconsensus_mutations_only_no_realignment() {
    let mut block = block_mutations_only();
    let expected_block = block_mutations_only_after_reconsensus();

    let majority_edits = block.find_majority_edits();

    // Apply mutations
    apply_mutation_reconsensus(&mut block, &majority_edits.subs).unwrap();

    // Check that no indels require re-alignment
    let has_indels = majority_edits.has_indels();
    assert!(!has_indels); // Should return false because no indels require re-alignment

    // But consensus should be updated and mutations healed
    assert_eq!(block.consensus(), expected_block.consensus());
    assert_eq!(block.alignments(), expected_block.alignments());
  }

  fn block_for_graph_test() -> PangraphBlock {
    let consensus = "GCCTCTTCCCGACCACGCGTTACAACATGGGACAGGCCTGCGCTTGAGGC";
    //               0         1         2         3         4
    //               01234567890123456789012345678901234567890123456789
    let aln = btreemap! {
        NodeId(1) => Edit::new(vec![],  vec![Del::new(0, 40)],   vec![]),  // deletion from position 0 to 40
        NodeId(2) => Edit::new(vec![],  vec![Del::new(35, 15)],  vec![]),  // deletion from position 35 to end (50-35=15)
        NodeId(3) => Edit::new(vec![],  vec![Del::new(35, 15)],  vec![]),  // deletion from position 35 to end
        NodeId(4) => Edit::new(vec![],  vec![Del::new(35, 15)],  vec![]),  // deletion from position 35 to end
        NodeId(5) => Edit::new(vec![],  vec![],                  vec![])   // no deletions
    };
    PangraphBlock::new(BlockId(20), consensus, aln)
  }

  fn block_for_graph_test_expected() -> PangraphBlock {
    // After reconsensus, the majority deletions (positions 35-49) should be applied
    // Original consensus: "GCCTCTTCCCGACCACGCGTTACAACATGGGACAGGCCTGCGCTTGAGGC" (50 chars)
    // Expected consensus: "GCCTCTTCCCGACCACGCGTTACAACATGGGACAG" (35 chars)
    let consensus = "GCCTCTTCCCGACCACGCGTTACAACATGGGACAG";
    let aln = btreemap! {
        NodeId(2) => Edit::new(vec![],                           vec![],                  vec![]),       // no edits (was majority deletion)
        NodeId(3) => Edit::new(vec![],                           vec![],                  vec![]),       // no edits (was majority deletion)
        NodeId(4) => Edit::new(vec![],                           vec![],                  vec![]),       // no edits (was majority deletion)
        NodeId(5) => Edit::new(vec![Ins::new(35, "GCCTGCGCTTGAGGC")], vec![],                  vec![])        // insertion to represent chars that were deleted from consensus
    };
    // NodeId(1) was detached because it was empty after deletions
    PangraphBlock::new(BlockId(20), consensus, aln)
  }

  fn singleton_block_expected() -> PangraphBlock {
    // This is the singleton block that should be created for NodeId(1) after reconsensus
    // this should be the reverse-complement of CGCTTGAGGC, because NodeId(1) was in Reverse strand
    let consensus = Seq::from("GCCTCAAGCG");
    PangraphBlock::from_consensus(consensus.clone(), BlockId(id((NodeId(1), &consensus))), NodeId(1))
  }

  #[test]
  fn test_reconsensus_graph() {
    // Create a block that will be modified by reconsensus
    let initial_block = block_for_graph_test();
    let expected_block = block_for_graph_test_expected();
    let singleton_block_exp = singleton_block_expected();

    // Create nodes for the block with lengths reflecting actual sequence lengths
    let nodes = btreemap! {
      NodeId(1) => PangraphNode::new(Some(NodeId(1)), initial_block.id(), PathId(1), Reverse, (0, 10)),   // 50 - 40 = 9 (deletes positions 0-39)
      NodeId(2) => PangraphNode::new(Some(NodeId(2)), initial_block.id(), PathId(2), Forward, (0, 35)),  // 50 - 15 = 35 (deletes positions 35-49)
      NodeId(3) => PangraphNode::new(Some(NodeId(3)), initial_block.id(), PathId(3), Forward, (0, 35)),  // 50 - 15 = 35 (deletes positions 35-49)
      NodeId(4) => PangraphNode::new(Some(NodeId(4)), initial_block.id(), PathId(4), Forward, (0, 35)),  // 50 - 15 = 35 (deletes positions 35-49)
      NodeId(5) => PangraphNode::new(Some(NodeId(5)), initial_block.id(), PathId(5), Forward, (0, 49)),  // no deletions
    };

    // Create paths
    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), [NodeId(1)], 49, false, None, None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), [NodeId(2)], 49, false, None, None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), [NodeId(3)], 49, false, None, None),
      PathId(4) => PangraphPath::new(Some(PathId(4)), [NodeId(4)], 49, false, None, None),
      PathId(5) => PangraphPath::new(Some(PathId(5)), [NodeId(5)], 49, false, None, None),
    };

    // Create blocks map
    let blocks = btreemap! {
      initial_block.id() => initial_block.clone()
    };

    // Create the graph
    let mut graph = Pangraph { paths, blocks, nodes };

    // Apply reconsensus_graph
    let result = reconsensus_graph(&mut graph, &[initial_block.id()], &PangraphBuildArgs::default());

    // Check that the operation succeeded
    result.unwrap();

    // Get the two created blocks, the modified and the singleton
    let final_block = &graph.blocks[&initial_block.id()];
    let singleton_block = &graph.blocks[&singleton_block_exp.id()];

    // Direct comparison with expected result
    assert_eq!(final_block.consensus(), expected_block.consensus());
    assert_eq!(final_block.alignments(), expected_block.alignments());

    // Check that the empty node was detached and a new singleton block was created
    assert_eq!(singleton_block.consensus(), singleton_block_exp.consensus());
    assert_eq!(singleton_block.alignments(), singleton_block_exp.alignments());

    // check that the node was updated correctly, flipping the strandedness
    let new_node1 = &graph.nodes[&NodeId(1)];
    let expected_node1 = PangraphNode::new(Some(NodeId(1)), singleton_block_exp.id(), PathId(1), Forward, (0, 10));
    assert_eq!(new_node1, &expected_node1);
  }
}
