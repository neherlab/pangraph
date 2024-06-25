use crate::pangraph::pangraph_block::PangraphBlock;
use eyre::Report;

/// Merge two blocks into a single block. This function modifies left block.
pub fn merge_blocks_inplace<'l>(
  left: &'l mut PangraphBlock,
  right: &PangraphBlock,
) -> Result<&'l PangraphBlock, Report> {
  // Append new sequences to left block
  for (node_id, edits) in &right.alignments {
    let seq = edits.apply(&right.consensus)?;
    left.append_sequence(&seq, *node_id)?;
  }
  Ok(left)
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_node::NodeId;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  #[rstest]
  fn test_merge_blocks_simple_case() -> Result<(), Report> {
    let cons_A = o!("GACTAAACCTGTCCGCTGAAACTGATCGGGGTACTGCAGC");
    let aln_A = BTreeMap::from([
      (
        NodeId(0),
        Edit {
          inss: vec![Ins::new(29, "GAT")],
          dels: vec![Del::new(7, 3)],
          subs: vec![Sub::new(33, 'T')],
        },
      ),
      (
        NodeId(1),
        Edit {
          inss: vec![Ins::new(15, "ACA")],
          dels: vec![Del::new(30, 3)],
          subs: vec![Sub::new(21, 'T')],
        },
      ),
      (
        NodeId(2),
        Edit {
          inss: vec![Ins::new(14, "AAT")],
          dels: vec![Del::new(0, 3)],
          subs: vec![Sub::new(0, 'C')],
        },
      ),
    ]);

    let mut block_A = PangraphBlock::new(cons_A, aln_A);

    let cons_B = o!("GACCAAACCTGTCCGCTGAAACTGCGGGGTACTGCAGC");
    let aln_B = BTreeMap::from([
      (
        NodeId(3),
        Edit {
          inss: vec![],
          dels: vec![Del::new(13, 3)],
          subs: vec![Sub::new(4, 'T')],
        },
      ),
      (
        NodeId(4),
        Edit {
          inss: vec![Ins::new(26, "AAA")],
          dels: vec![],
          subs: vec![Sub::new(3, 'T')],
        },
      ),
    ]);
    let block_B = PangraphBlock::new(cons_B, aln_B);
    let new_block = merge_blocks_inplace(&mut block_A, &block_B)?;

    let new_cons = o!("GACTAAACCTGTCCGCTGAAACTGATCGGGGTACTGCAGC");
    let new_aln = BTreeMap::from([
      (
        NodeId(0),
        Edit {
          inss: vec![Ins::new(29, "GAT")],
          dels: vec![Del::new(7, 3)],
          subs: vec![Sub::new(33, 'T')],
        },
      ),
      (
        NodeId(1),
        Edit {
          inss: vec![Ins::new(15, "ACA")],
          dels: vec![Del::new(30, 3)],
          subs: vec![Sub::new(21, 'T')],
        },
      ),
      (
        NodeId(2),
        Edit {
          inss: vec![Ins::new(14, "AAT")],
          dels: vec![Del::new(0, 3)],
          subs: vec![Sub::new(0, 'C')],
        },
      ),
      (
        NodeId(3),
        Edit {
          inss: vec![],
          dels: vec![Del::new(13, 3), Del::new(24, 2)],
          subs: vec![Sub::new(3, 'C'), Sub::new(4, 'T')],
        },
      ),
      (
        NodeId(4),
        Edit {
          inss: vec![Ins::new(28, "AAA")],
          dels: vec![Del::new(24, 2)],
          subs: vec![],
        },
      ),
    ]);

    let expected_new_block = PangraphBlock {
      consensus: new_cons,
      alignments: new_aln,
    };

    assert_eq!(&expected_new_block, new_block);
    Ok(())
  }
}
