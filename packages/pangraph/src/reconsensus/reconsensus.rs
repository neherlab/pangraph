use crate::commands::build::build_args::PangraphBuildArgs;
use crate::make_report;
use crate::pangraph::detach_unaligned::detach_unaligned_nodes;
use crate::pangraph::edits::Edit;
use crate::pangraph::edits::Sub;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::reconsensus::remove_nodes::find_empty_nodes;
use eyre::{Context, Report};
use itertools::{Either, Itertools};
use serde::{Deserialize, Serialize};

/// Result of analyzing blocks for reconsensus operation
#[derive(Clone, Debug, Serialize, Deserialize)]
struct BlockAnalysis {
  /// Blocks that need only mutation-based reconsensus (no realignment)
  mutations_only: Vec<(BlockId, Edit)>,
  /// Blocks that need full reconsensus with realignment due to indels
  need_realignment: Vec<(BlockId, Edit)>,
}

/// Applies the reconsensus operation to each updated block in the graph:
/// - updates the block consensus following a merge
/// - removes potentially empty nodes
///
/// Reconsensus function:
/// for blocks that originate from a new merger:
///   - check whether blocks have majority substitutions, deletions, or insertions
///   - for blocks with only majority substitutions, updates the consensus
///     and transfers the substitutions to the alignment
///   - for blocks with majority indels, updates the consensus and realigns the sequences.
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

  // Analyze blocks, find majority edits and determine which need realignment
  let analysis = analyze_blocks_for_reconsensus(graph, ids_updated_blocks);

  // Apply mutation-only reconsensus (no realignment needed)
  analysis.mutations_only.into_iter().try_for_each(|(block_id, edits)| {
    // get mutable block
    let block = graph
      .blocks
      .get_mut(&block_id)
      .ok_or_else(|| make_report!("Block {} not found in graph", block_id))?;
    // apply the substitions to the block consensus and alignment
    apply_substitutions_to_block(block, &edits.subs)
      .wrap_err_with(|| format!("When processing block {} for reconsensus via substitutions", block.id()))
  })?;

  // Handle blocks requiring realignment
  if !analysis.need_realignment.is_empty() {
    let mut realigned_blocks: Vec<PangraphBlock> = analysis
      .need_realignment
      .into_iter()
      .map(|(block_id, edits)| {
        // Pop the realigned blocks from graph.blocks
        let block = graph
          .blocks
          .remove(&block_id)
          .ok_or_else(|| make_report!("Block {} not found in graph", block_id))?;

        // apply edits to the block consensus and re-align sequences
        // returns a new block
        let new_block = block
          .edit_consensus_and_realign(&edits, args)
          .wrap_err_with(|| "When processing block for reconsensus via realignment")?;
        // Return the realigned block if successful
        Ok::<PangraphBlock, Report>(new_block)
      })
      .collect::<Result<Vec<_>, _>>()?;

    // Apply detach_unaligned_nodes. This removes unaligned nodes and re-adds them to
    // the `realigned_blocks` list as new blocks.
    // it also modifies the nodes dictionary accordingly.
    // Nb: This is done to handle the edge-case when re-alignment could generate
    // nodes that are not empty, but have zero aligned nucleotides to the consensus
    detach_unaligned_nodes(&mut realigned_blocks, &mut graph.nodes)?;

    // Re-add all the blocks (including potentially new singleton blocks) to graph.blocks
    realigned_blocks.into_iter().for_each(|block| {
      graph.blocks.insert(block.id(), block);
    });
  }

  Ok(())
}

/// Analyzes blocks to determine which need realignment vs. mutation-only reconsensus
fn analyze_blocks_for_reconsensus(graph: &Pangraph, block_ids: &[BlockId]) -> BlockAnalysis {
  let (mutations_only, need_realignment): (Vec<_>, Vec<_>) = block_ids
    .iter()
    .filter_map(|&block_id| {
      let block = &graph.blocks[&block_id];
      let majority_edits = block.find_majority_edits();

      if majority_edits.has_indels() {
        // the block needs full realignment
        Some(Either::Right((block_id, majority_edits)))
      } else if majority_edits.has_subs() {
        // no re-alignment needed, substitutions are sufficient
        Some(Either::Left((block_id, majority_edits)))
      } else {
        None // Blocks with no variants are skipped
      }
    })
    .partition_map(|either| either);

  BlockAnalysis {
    mutations_only,
    need_realignment,
  }
}

/// Applies a set of substitutions to the block's consensus sequence and alignment.
/// Does not re-align the sequences.
fn apply_substitutions_to_block(block: &mut PangraphBlock, subs: &[Sub]) -> Result<(), Report> {
  subs.iter().try_for_each(|sub| {
    block
      .change_consensus_nucleotide_at_pos(sub.pos, sub.alt)
      .wrap_err_with(|| format!("Failed to apply mutation at pos {}: {}", sub.pos, sub.alt))
  })
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
  use crate::representation::seq::Seq;
  use crate::utils::id::id;

  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  fn block_0() -> PangraphBlock {
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
    PangraphBlock::new(BlockId(0), consensus, aln)
  }

  fn block_0_reconsensus() -> PangraphBlock {
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
    PangraphBlock::new(BlockId(0), consensus, aln)
  }

  fn block_1() -> PangraphBlock {
    let consensus = "AGGACTTCGATCTATTCGGAGAA";
    //               0         1         2
    //               01234567890123456789012
    //    node 1)    .T...--.........|A.....
    //    node 2)    .T...--...C............
    //    node 3)    .T...--...C.....--.....
    //    node 4)    .C.......---.....A.....
    //    node 5)    .....|-..........A.....
    //     L = 23, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![Ins::new(17, "TTTT")],  vec![Del::new(5, 2)],                   vec![Sub::new(1, 'T'), Sub::new(17, 'A')]),
      NodeId(2) => Edit::new(vec![],                      vec![Del::new(5, 2)],                   vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(3) => Edit::new(vec![],                      vec![Del::new(5, 2), Del::new(16,2)],   vec![Sub::new(1, 'T'), Sub::new(10, 'C')]),
      NodeId(4) => Edit::new(vec![],                      vec![Del::new(9, 3)],                   vec![Sub::new(1, 'C'), Sub::new(17, 'A')]),
      NodeId(5) => Edit::new(vec![Ins::new(5, "AA")],     vec![Del::new(5, 2)],                   vec![Sub::new(17, 'A')]),
    };
    PangraphBlock::new(BlockId(1), consensus, aln)
  }

  fn block_1_mut_reconsensus() -> PangraphBlock {
    let consensus = "ATGACTTCGATCTATTCAGAGAA";
    //               0         1         2
    //               01234567890123456789012
    //    node 1)    .....--..........|.....
    //    node 2)    .....--...C......G.....
    //    node 3)    .....--...C.....--.....
    //    node 4)    .C.......---...........
    //    node 5)    .G...|-................
    //     L = 23, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![Ins::new(17, "TTTT")],  vec![Del::new(5, 2)],                  vec![]),
      NodeId(2) => Edit::new(vec![],                      vec![Del::new(5, 2)],                  vec![Sub::new(10, 'C'), Sub::new(17, 'G')]),
      NodeId(3) => Edit::new(vec![],                      vec![Del::new(5, 2), Del::new(16,2)],  vec![Sub::new(10, 'C')]),
      NodeId(4) => Edit::new(vec![],                      vec![Del::new(9, 3)],                  vec![Sub::new(1, 'C')]),
      NodeId(5) => Edit::new(vec![Ins::new(5, "AA")],     vec![Del::new(5, 2)],                  vec![Sub::new(1, 'G')]),
    };
    PangraphBlock::new(BlockId(1), consensus, aln)
  }

  fn block_1_reconsensus() -> PangraphBlock {
    let consensus = "ATGACCGATCTATTCAGAGAA";
    //               0         1         2
    //               012345678901234567890
    //    node 1)    .....|.........|.....
    //    node 2)    .....|..C......G.....
    //    node 3)    .....|..C.....--.....
    //    node 4)    .C...|.---...........
    //    node 5)    .G...|...............
    //     L = 21, N = 5
    let aln = btreemap! {
      NodeId(1) => Edit::new(vec![Ins::new(15, "TTTT")],  vec![],                   vec![]),
      NodeId(2) => Edit::new(vec![],                      vec![],                   vec![Sub::new(8, 'C'), Sub::new(15, 'G')]),
      NodeId(3) => Edit::new(vec![],                      vec![Del::new(14,2)],     vec![Sub::new(8, 'C')]),
      NodeId(4) => Edit::new(vec![Ins::new(5, "TT")],     vec![Del::new(7, 3)],     vec![Sub::new(1, 'C')]),
      NodeId(5) => Edit::new(vec![Ins::new(5, "AA")],     vec![],                   vec![Sub::new(1, 'G')]),
    };
    PangraphBlock::new(BlockId(1), consensus, aln)
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
    PangraphBlock::new(BlockId(2), consensus, aln)
  }

  // fn block_2_reconsensus() -> PangraphBlock {
  //   let consensus = "GATGACTTCGATCTATTCGGAGAA";
  //   //               0         1         2
  //   //               01234567890123456789012
  //   //    node 1)    ....|.--......|...A..-..
  //   //    node 2)    ......--...C..|......--.|
  //   //    node 3)    -....----..C............|
  //   //    node 4)    -.C.|.....---.....A.....|
  //   //    node 5)    ..G.|.........|...A.--..
  //   //     L = 23, N = 5
  //   let aln = btreemap! {
  //     NodeId(1) => Edit::new(vec![Ins::new(4, "AA"), Ins::new(13, "AA")],   vec![Del::new(5, 2), Del::new(20, 1)],  vec![Sub::new(17, 'A')]),
  //     NodeId(2) => Edit::new(vec![Ins::new(13, "AA"), Ins::new(23, "TT")],  vec![Del::new(5, 2), Del::new(20, 2)],  vec![Sub::new(10, 'C')]),
  //     NodeId(3) => Edit::new(vec![Ins::new(23, "TT")],                                        vec![Del::new(0,1), Del::new(4, 4)],                   vec![Sub::new(10, 'C')]),
  //     NodeId(4) => Edit::new(vec![Ins::new(4, "C"), Ins::new(23, "TT")],                      vec![Del::new(0,1), Del::new(9, 3)],                   vec![Sub::new(2, 'C'), Sub::new(17, 'A')]),
  //     NodeId(5) => Edit::new(vec![Ins::new(4, "C"), Ins::new(13, "AA")],    vec![Del::new(19, 2)],                  vec![Sub::new(2, 'G'), Sub::new(17, 'A')]),
  //   };
  //   PangraphBlock::new(BlockId(2), consensus, aln)
  // }

  fn block_3() -> PangraphBlock {
    let consensus = "GCCTCTTCCCGACCACGCGTTACAACATGGGACAGGCCTGCGCTTGAGGC";
    //               0         1         2         3         4
    //               01234567890123456789012345678901234567890123456789
    //    node 1)    .....A.............----...........................
    //    node 2)    .....A..............---.....................|....|
    //    node 3)    ..............G............G......................
    //    node 4)    .....A..............---..........................|
    //    node 5)    .................................................|
    //    L = 50, N = 5
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
    //               0         1         2         3         4
    //               0123456789012345678901234567890123456789012345678
    //    node 1)    ...................-...........................--
    //    node 2)    ................................|................
    //    node 3)    .....T........G.....|..G.......................--
    //    node 4)    .................................................
    //    node 5)    .....T..............|............................
    //     L = 49, N = 5
    let aln = btreemap! {
        NodeId(1) => Edit::new(vec![],                    vec![Del::new(19, 1), Del::new(47, 2)],   vec![]),
        NodeId(2) => Edit::new(vec![Ins::new(32, "AA")],  vec![],                                   vec![]),
        NodeId(3) => Edit::new(vec![Ins::new(20, "TAC")], vec![Del::new(47, 2)],                    vec![Sub::new(5, 'T'), Sub::new(14, 'G'), Sub::new(24, 'G')]),
        NodeId(4) => Edit::new(vec![],                    vec![],                                   vec![]),
        NodeId(5) => Edit::new(vec![Ins::new(20, "TAC")], vec![],                                   vec![Sub::new(5, 'T')])
    };
    PangraphBlock::new(BlockId(3), consensus, aln)
  }

  // TODO:
  // - test majority substitutions, insertions, deletions on the examples
  // x test analyze_blocks_for_reconsensus
  // - test only-substitutions reconsensus
  // - test full reconsensus

  #[test]
  fn test_analyze_block_reconsensus() {
    let graph = Pangraph {
      blocks: btreemap! (
        BlockId(0) => block_0(),
        BlockId(1) => block_1(),
        BlockId(2) => block_2(),
        BlockId(3) => block_3(),
      ),
      paths: btreemap! {},
      nodes: btreemap! {},
    };

    let block_ids = vec![BlockId(0), BlockId(1), BlockId(2), BlockId(3)];
    let results = analyze_blocks_for_reconsensus(&graph, &block_ids);

    let subs_blockids: Vec<BlockId> = results.mutations_only.iter().map(|(bid, _)| *bid).collect();
    let realign_blockids: Vec<BlockId> = results.need_realignment.iter().map(|(bid, _)| *bid).collect();

    assert_eq!(subs_blockids, vec![BlockId(0)]);
    assert_eq!(realign_blockids, vec![BlockId(1), BlockId(2), BlockId(3)]);
  }

  #[test]
  fn test_find_majority_edits_block0() {
    let edits = block_0().find_majority_edits();
    let expected_edits = Edit::new(vec![], vec![], vec![Sub::new(1, 'C')]);
    assert_eq!(edits, expected_edits);
  }

  #[test]
  fn test_find_majority_edits_block1() {
    let edits = block_1().find_majority_edits();
    let expected_edits = Edit::new(vec![], vec![Del::new(5, 2)], vec![Sub::new(1, 'T'), Sub::new(17, 'A')]);
    assert_eq!(edits, expected_edits);
  }

  #[test]
  fn test_find_majority_edits_block2() {
    let edits = block_2().find_majority_edits();
    let expected_edits = Edit::new(
      vec![Ins::new(0, "G"), Ins::new(13, "AA"), Ins::new(23, "TT")],
      vec![Del::new(5, 2), Del::new(20, 1)],
      vec![Sub::new(1, 'T'), Sub::new(17, 'A')],
    );
    assert_eq!(edits, expected_edits);
  }

  #[test]
  fn test_find_majority_edits_block3() {
    let edits = block_3().find_majority_edits();
    let expected_edits = Edit::new(vec![Ins::new(50, "TT")], vec![Del::new(20, 3)], vec![Sub::new(5, 'A')]);
    assert_eq!(edits, expected_edits);
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
  fn test_mutations_only_reconsensus_block0() {
    let mut block = block_0();
    let expected_block = block_0_reconsensus();

    // Apply mutations
    let majority_edits = block.find_majority_edits();
    assert!(!majority_edits.has_indels()); // This block has no indels requiring re-alignment
    apply_substitutions_to_block(&mut block, &majority_edits.subs).unwrap();

    // But consensus should be updated and mutations healed
    assert_eq!(block, expected_block);
  }

  #[test]
  fn test_mutations_only_reconsensus_block1() {
    let mut block = block_1();
    let expected_block = block_1_mut_reconsensus();

    // Apply mutations
    let majority_edits = block.find_majority_edits();
    apply_substitutions_to_block(&mut block, &majority_edits.subs).unwrap();

    // But consensus should be updated and mutations healed
    assert_eq!(block, expected_block);
  }

  #[test]
  fn test_realign_reconsensus_block1() {
    let block = block_1();
    let expected_block = block_1_reconsensus();

    // Apply majority edits
    let majority_edits = block.find_majority_edits();
    assert!(majority_edits.has_indels()); // This block has indels requiring re-alignment
    let block = block
      .edit_consensus_and_realign(&majority_edits, &PangraphBuildArgs::default())
      .unwrap();

    // Check that the re-alignment produced the expected result
    assert_eq!(block, expected_block);
  }

  // #[test]
  // fn test_realign_reconsensus_block2() {
  //   let block = block_2();
  //   let expected_block = block_2_reconsensus();

  //   // Apply majority edits
  //   let majority_edits = block.find_majority_edits();
  //   assert!(majority_edits.has_indels()); // This block has indels requiring re-alignment
  //   let block = block
  //     .edit_consensus_and_realign(&majority_edits, &PangraphBuildArgs::default())
  //     .unwrap();

  //   // Check that the re-alignment produced the expected result
  //   assert_eq!(block, expected_block);
  // }

  #[test]
  fn test_realign_reconsensus_block3() {
    let block = block_3();
    let expected_block = block_3_reconsensus();

    // Apply majority edits
    let majority_edits = block.find_majority_edits();
    assert!(majority_edits.has_indels()); // This block has indels requiring re-alignment
    let block = block
      .edit_consensus_and_realign(&majority_edits, &PangraphBuildArgs::default())
      .unwrap();

    // Check that the re-alignment produced the expected result
    assert_eq!(block, expected_block);
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

  fn signleton_block_expected() -> PangraphBlock {
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
    let singleton_block_exp = signleton_block_expected();

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
