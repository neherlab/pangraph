use crate::circularize::circularize::remove_transitive_edges;
use crate::commands::simplify::simplify_args::PangraphSimplifyArgs;
use crate::io::file::create_file_or_stdout;
use crate::io::json::{json_write, JsonPretty};
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_path::PathId;
use eyre::Report;
use std::collections::BTreeSet;

#[allow(unused_must_use)]
pub fn simplify_run(args: PangraphSimplifyArgs) -> Result<(), Report> {
  let PangraphSimplifyArgs { input, output, strains } = args;

  let mut graph = Pangraph::from_path(&input)?;

  let strain_set = strains.into_iter().collect();
  simplify(&mut graph, &strain_set);

  let f = create_file_or_stdout(&output)?;
  json_write(f, &graph, JsonPretty(true))
}

fn simplify(graph: &mut Pangraph, focal_paths: &BTreeSet<String>) -> Result<(), Report> {
  let path_ids_to_remove: Vec<PathId> = graph
    .paths
    .iter()
    .filter(|(_, path)| !focal_paths.contains(path.name().as_ref().unwrap()))
    .map(|(id, _)| *id)
    .collect();

  for pid in path_ids_to_remove {
    graph.remove_path(pid);
  }

  remove_transitive_edges(graph)?;

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
  use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
  use crate::pangraph::pangraph_path::PangraphPath;
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;

  const NID11: NodeId = NodeId(13172209629052542373);
  const NID12: NodeId = NodeId(16864511183100055928);

  fn block_a() -> PangraphBlock {
    //          0         1         2         3
    //          01234567890123456789012345678901
    // cons:    ACTATATTACGGCGATCGATCGATTACTCGCT
    //   n1:    ...G............................  l = 32
    //   n2:    .......|.....xxx................  l = 31
    //   n3:    ................................| l = 35
    PangraphBlock::new(
      BlockId(1),
      "ACTATATTACGGCGATCGATCGATTACTCGCT",
      btreemap! {
        NodeId(1) => Edit::new(vec![],                    vec![],                vec![Sub::new(3, 'G')]),
        NodeId(2) => Edit::new(vec![Ins::new(7, "AA")],   vec![Del::new(13, 3)], vec![]),
        NodeId(3) => Edit::new(vec![Ins::new(32, "CCC")], vec![],                vec![]),
      },
    )
  }

  fn block_b() -> PangraphBlock {
    //          0         1         2         3
    //          01234567890123456789012345678901
    // cons:    CATGCTACGCTACGCATTATCGATCGCATCGA
    //   n4:    ..........G.....................  l = 32
    //   n5:    .............xxx................  l = 29
    //   n6:    ................................| l = 35
    PangraphBlock::new(
      BlockId(2),
      "CATGCTACGCTACGCATTATCGATCGCATCGA",
      btreemap! {
        NodeId(4) => Edit::new(vec![],                    vec![],                vec![Sub::new(10, 'G')]),
        NodeId(5) => Edit::new(vec![],                    vec![Del::new(13, 3)], vec![]),
        NodeId(6) => Edit::new(vec![Ins::new(32, "AAA")], vec![],                vec![]),
      },
    )
  }

  fn block_c() -> PangraphBlock {
    //          0         1
    //          01234567890123456
    // cons:    ACGTGTACTAGTACTGC
    //   n7:    ................. l = 17
    //   n8:    ............C.... l = 17
    PangraphBlock::new(
      BlockId(3),
      "ACGTGTACTAGTACTGC",
      btreemap! {
        NodeId(7) => Edit::new(vec![], vec![], vec![]),
        NodeId(8) => Edit::new(vec![], vec![], vec![Sub::new(12, 'C')]),
      },
    )
  }

  fn block_ab() -> PangraphBlock {
    //          0         1         2         3         4         5         6
    //          0123456789012345678901234567890123456789012345678901234567890123
    // cons:    ACTATATTACGGCGATCGATCGATTACTCGCTCATGCTACGCTACGCATTATCGATCGCATCGA
    //  n10:    ...G......................................G.....................
    //  n11:    .......|.....xxx.............................xxx................
    PangraphBlock::new(
      BlockId(1),
      "ACTATATTACGGCGATCGATCGATTACTCGCTCATGCTACGCTACGCATTATCGATCGCATCGA",
      btreemap! {
        NID11 => Edit::new(vec![],                  vec![],                                 vec![Sub::new(3, 'G'), Sub::new(42, 'G')]),
        NID12 => Edit::new(vec![Ins::new(7, "AA")], vec![Del::new(13, 3), Del::new(45, 3)], vec![]),
      },
    )
  }

  fn graph() -> Pangraph {
    // n1+ -> n4+ -> n7+
    // n2+ -> n5+ -> n8-
    // n3+ -> n6-
    let nodes = btreemap! {
      NodeId(1) => PangraphNode::new(Some(NodeId(1)), BlockId(1), PathId(1), Forward, (0, 32)),
      NodeId(2) => PangraphNode::new(Some(NodeId(2)), BlockId(1), PathId(2), Forward, (0, 31)),
      NodeId(3) => PangraphNode::new(Some(NodeId(3)), BlockId(1), PathId(3), Forward, (0, 35)),
      NodeId(4) => PangraphNode::new(Some(NodeId(4)), BlockId(2), PathId(1), Forward, (32, 64)),
      NodeId(5) => PangraphNode::new(Some(NodeId(5)), BlockId(2), PathId(2), Forward, (31, 60)),
      NodeId(6) => PangraphNode::new(Some(NodeId(6)), BlockId(2), PathId(3), Forward, (35, 0)),
      NodeId(7) => PangraphNode::new(Some(NodeId(7)), BlockId(3), PathId(1), Forward, (64, 0)),
      NodeId(8) => PangraphNode::new(Some(NodeId(8)), BlockId(3), PathId(2), Reverse, (60, 0)),
    };
    let blocks = btreemap! {
      BlockId(1) => block_a(),
      BlockId(2) => block_b(),
      BlockId(3) => block_c(),
    };
    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), vec![NodeId(1), NodeId(4), NodeId(7)], 81, true, Some(o!("pathA")), None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), vec![NodeId(2), NodeId(5), NodeId(8)], 77, true, Some(o!("pathB")), None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), vec![NodeId(3), NodeId(6)],            70, true, Some(o!("pathC")), None),
    };
    Pangraph { paths, blocks, nodes }
  }

  fn expected_graph() -> Pangraph {
    let nodes = btreemap! {
      NID11 =>     PangraphNode::new(Some(NodeId(11)), BlockId(1), PathId(1), Forward, (0, 64)),
      NID12 =>     PangraphNode::new(Some(NodeId(12)), BlockId(1), PathId(2), Forward, (0, 60)),
      NodeId(7) => PangraphNode::new(Some(NodeId(7)),  BlockId(3), PathId(1), Forward, (64, 0)),
      NodeId(8) => PangraphNode::new(Some(NodeId(8)),  BlockId(3), PathId(2), Reverse, (60, 0)),
    };
    let blocks = btreemap! {
        BlockId(1) => block_ab(),
        BlockId(3) => block_c(),
    };
    let paths = btreemap! {
      PathId(1) => PangraphPath::new(Some(PathId(1)), vec![NID11, NodeId(7)], 81, true, Some(o!("pathA")), None),
      PathId(2) => PangraphPath::new(Some(PathId(2)), vec![NID12, NodeId(8)], 77, true, Some(o!("pathB")), None),
    };

    Pangraph { paths, blocks, nodes }
  }

  #[test]
  fn test_remove_path() {
    let mut graph = graph();

    graph.remove_path(PathId(1));

    let expected_paths = btreemap! {
      PathId(2) => PangraphPath::new(Some(PathId(2)), vec![NodeId(2), NodeId(5), NodeId(8)], 77, true, Some(o!("pathB")), None),
      PathId(3) => PangraphPath::new(Some(PathId(3)), vec![NodeId(3), NodeId(6)], 70, true, Some(o!("pathC")), None),
    };

    let expected_nodes = btreemap! {
      NodeId(2) => PangraphNode::new(Some(NodeId(2)), BlockId(1), PathId(2), Forward, (0, 31)),
      NodeId(3) => PangraphNode::new(Some(NodeId(3)), BlockId(1), PathId(3), Forward, (0, 35)),
      NodeId(5) => PangraphNode::new(Some(NodeId(5)), BlockId(2), PathId(2), Forward, (31, 60)),
      NodeId(6) => PangraphNode::new(Some(NodeId(6)), BlockId(2), PathId(3), Forward, (35, 0)),
      NodeId(8) => PangraphNode::new(Some(NodeId(8)), BlockId(3), PathId(2), Reverse, (60, 0)),
    };

    let expected_blocks = btreemap! {
      BlockId(1) => PangraphBlock::new(BlockId(1), "ACTATATTACGGCGATCGATCGATTACTCGCT", btreemap!{
        NodeId(2) => Edit::new(vec![Ins::new(7, "AA")],   vec![Del::new(13, 3)], vec![]),
        NodeId(3) => Edit::new(vec![Ins::new(32, "CCC")], vec![],                vec![]),
      }),
      BlockId(2) => PangraphBlock::new(BlockId(2), "CATGCTACGCTACGCATTATCGATCGCATCGA", btreemap!{
        NodeId(5) => Edit::new(vec![],                    vec![Del::new(13, 3)], vec![]),
        NodeId(6) => Edit::new(vec![Ins::new(32, "AAA")], vec![],                vec![]),
      }),
      BlockId(3) => PangraphBlock::new(BlockId(3), "ACGTGTACTAGTACTGC", btreemap!{
        NodeId(8) => Edit::new(vec![], vec![], vec![Sub::new(12, 'C')]),
      }),
    };

    assert_eq!(graph.paths, expected_paths);
    assert_eq!(graph.nodes, expected_nodes);
    assert_eq!(graph.blocks, expected_blocks);
  }

  #[test]
  fn test_simplify() {
    let mut graph = graph();
    let expected_graph = expected_graph();
    simplify(&mut graph, &btreeset! {o!("pathA"), o!("pathB")}).unwrap();
    assert_eq!(graph.paths, expected_graph.paths);
    assert_eq!(graph.blocks, expected_graph.blocks);
  }
}
