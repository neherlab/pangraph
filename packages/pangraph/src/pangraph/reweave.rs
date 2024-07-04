use crate::align::alignment::{Alignment2, AnchorBlock, ExtractedHit};
use crate::io::seq::reverse_complement;
use crate::pangraph::pangraph::{GraphUpdate, Pangraph};
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::pangraph_interval::extract_intervals;
use crate::pangraph::pangraph_node::NodeId;
use crate::pangraph::slice::block_slice;
use crate::utils::id::id;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::hash::Hash;

#[derive(Debug, Eq, PartialEq)]
struct MergePromise {
  anchor_block: PangraphBlock,
  append_block: PangraphBlock,
  orientation: bool,
}

impl MergePromise {
  fn new(anchor_block: PangraphBlock, append_block: PangraphBlock, orientation: bool) -> Self {
    Self {
      anchor_block,
      append_block,
      orientation,
    }
  }

  fn solve_promise(&mut self) -> Result<PangraphBlock, Report> {
    for (node_id, edits) in self.append_block.alignments() {
      let mut seq = edits.apply(self.append_block.consensus())?;
      if self.orientation {
        seq = reverse_complement(&seq)?;
      };
      self.anchor_block.append_sequence(&seq, *node_id)?;
    }
    Ok(self.anchor_block.clone())
  }
}

#[derive(Debug)]
struct ToMerge {
  pub block: PangraphBlock,
  pub is_anchor: bool,
  pub orientation: bool,
}

impl ToMerge {
  pub fn new(block: PangraphBlock, is_anchor: bool, orientation: bool) -> Self {
    Self {
      block,
      is_anchor,
      orientation,
    }
  }

  pub fn block_id(&self) -> BlockId {
    self.block.id()
  }
}

fn assign_new_block_ids(mergers: &mut [Alignment2]) {
  for a in mergers.iter_mut() {
    a.new_block_id = Some(BlockId(
      // FIXME: looks like this is trying to calculate its own hash id? It should probably not be done in random places like this.
      id((&a.qry.name, &a.qry.interval, &a.reff.name, &a.reff.interval)),
    ));
  }
}

fn assign_anchor_block(mergers: &mut [Alignment2], graph: &Pangraph) {
  for m in mergers.iter_mut() {
    let ref_block = &graph.blocks[&m.reff.name];
    let qry_block = &graph.blocks[&m.qry.name];
    if ref_block.depth() >= qry_block.depth() {
      m.anchor_block = Some(AnchorBlock::Ref);
    } else {
      m.anchor_block = Some(AnchorBlock::Qry);
    }
  }
}

fn target_blocks(mergers: &[Alignment2]) -> BTreeMap<BlockId, Vec<Alignment2>> {
  let mut target_blocks = BTreeMap::new();

  for merger in mergers {
    target_blocks
      .entry(merger.qry.name)
      .or_insert_with(Vec::new)
      .push(merger.clone());

    target_blocks
      .entry(merger.reff.name)
      .or_insert_with(Vec::new)
      .push(merger.clone());
  }

  target_blocks
}

fn extract_hits(bid: BlockId, mergers: &[Alignment2]) -> Vec<ExtractedHit> {
  let mut hits = Vec::with_capacity(mergers.len() * 2);

  for m in mergers {
    if m.reff.name == bid {
      hits.push(ExtractedHit {
        hit: m.reff.clone(),                                    // TODO: avoid copy
        new_block_id: m.new_block_id.unwrap(), // FIXME: "partially initialized object" anti-pattern shoots back
        is_anchor: m.anchor_block.unwrap() == AnchorBlock::Ref, // FIXME: "partially initialized object" anti-pattern shoots back
        orientation: m.orientation,
      });
    }

    if m.qry.name == bid {
      hits.push(ExtractedHit {
        hit: m.qry.clone(),
        new_block_id: m.new_block_id.unwrap(),
        is_anchor: m.anchor_block.unwrap() == AnchorBlock::Qry,
        orientation: m.orientation,
      });
    }
  }

  hits.shrink_to_fit();

  hits
}

fn group_promises(h: &[ToMerge]) -> Vec<MergePromise> {
  let mut promises = Vec::new();

  let groups = h
    .iter()
    .sorted_by_key(|x| x.block_id())
    .group_by(|x| x.block_id())
    .into_iter()
    .map(|(a, b)| (a, b.collect_vec()))
    .collect::<BTreeMap<_, _>>();

  for (_, bs) in groups {
    debug_assert_eq!(bs.len(), 2, "Only two blocks can be merged");

    let (b1, b2) = (&bs[0], &bs[1]);
    debug_assert!(b1.is_anchor ^ b2.is_anchor, "Only one block is the anchor");
    debug_assert_eq!(b1.orientation, b2.orientation, "Orientation must be the same");

    let anchor_block = if b1.is_anchor { &b1.block } else { &b2.block };
    let append_block = if b1.is_anchor { &b2.block } else { &b1.block };
    promises.push(MergePromise {
      anchor_block: anchor_block.clone(),
      append_block: append_block.clone(),
      orientation: b1.orientation,
    });
  }
  promises
}

fn split_block(bid: BlockId, mergers: &[Alignment2], graph: &Pangraph, thr_len: usize) -> (GraphUpdate, Vec<ToMerge>) {
  let extracted_hits = extract_hits(bid, mergers);
  let consensus_len = graph.blocks[&bid].consensus_len();
  let intervals = extract_intervals(&extracted_hits, consensus_len, thr_len);

  let mut u = GraphUpdate {
    b_old_id: bid,
    b_new: Vec::new(),
    n_new: graph.blocks[&bid]
      .alignment_keys()
      .iter()
      .map(|nid| (*nid, vec![]))
      .collect(),
  };

  let mut h = Vec::new();
  let b = &graph.blocks[&bid];
  for interval in intervals {
    let (b_slice, n_dict) = block_slice(b, &interval, graph);
    for (old_nid, new_node) in n_dict {
      u.n_new.entry(old_nid).or_default().push(new_node);
    }
    if interval.aligned {
      h.push(ToMerge {
        block: b_slice,
        is_anchor: interval.is_anchor.unwrap(), // FIXME: this field should never be uninitialized
        orientation: interval.orientation.unwrap(), // FIXME: this field should never be uninitialized
      });
    } else {
      u.b_new.push(b_slice);
    }
  }

  for (old_node_id, nodes) in &mut u.n_new {
    let strand = graph.nodes[old_node_id].strand();
    if !strand {
      *nodes = nodes.iter().cloned().rev().collect();
    }
  }
  (u, h)
}

fn reweave(mergers: &mut [Alignment2], graph: &mut Pangraph, thr_len: usize) -> (Pangraph, Vec<MergePromise>) {
  assign_new_block_ids(mergers);
  assign_anchor_block(mergers, graph);

  let tb = target_blocks(mergers);
  let mut u = vec![];
  let mut h = vec![];

  for (bid, m) in tb {
    let (graph_update, to_merge) = split_block(bid, &m, graph, thr_len);
    u.push(graph_update);
    h.extend(to_merge);
  }

  let merge_promises = group_promises(&h);
  for graph_update in u {
    graph.update(&graph_update);
  }
  (graph.clone(), merge_promises)
}

#[cfg(test)]
mod tests {
  #![allow(
    non_snake_case,
    clippy::redundant_clone,
    clippy::string_slice,
    clippy::many_single_char_names
  )]

  use super::*;
  use crate::align::alignment::Hit2;
  use crate::io::seq::generate_random_nuc_sequence;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::pangraph_path::{PangraphPath, PathId};
  use crate::utils::interval::Interval;
  use crate::utils::random::get_random_number_generator;
  use maplit::{btreemap, btreeset};
  use noodles::sam::record::Cigar;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_extract_hits() {
    fn new_hit(name: BlockId, start: usize) -> Hit2 {
      Hit2 {
        name,
        length: 0,                             // FIXME: `length` set to `None` in Python version
        interval: Interval::new(start, start), // FIXME: `end` was `None` in Python version
      }
    }

    fn create_hit(new_bid: BlockId, is_anchor: bool, strand: bool, hit: Hit2) -> ExtractedHit {
      ExtractedHit {
        new_block_id: new_bid,
        is_anchor,
        orientation: strand,
        hit,
      }
    }

    fn new_aln(new_block_id: BlockId, reff: &Hit2, qry: &Hit2, strand: bool, anchor_block: AnchorBlock) -> Alignment2 {
      Alignment2 {
        reff: reff.clone(),
        qry: qry.clone(),
        new_block_id: Some(new_block_id),
        orientation: strand,
        anchor_block: Some(anchor_block),

        // FIXME: these were all unset in Python version
        matches: 0,
        length: 0,
        quality: 0,
        cigar: Cigar::default(),
        divergence: None,
        align: None,
      }
    }

    let h1_a = new_hit(BlockId(1), 10);
    let h1_b = new_hit(BlockId(1), 20);
    let h1_c = new_hit(BlockId(1), 30);
    let h1_d = new_hit(BlockId(1), 40);
    let h2_e = new_hit(BlockId(2), 50);
    let h2_f = new_hit(BlockId(2), 60);
    let h2_g = new_hit(BlockId(2), 70);
    let h2_h = new_hit(BlockId(2), 80);

    let a1 = new_aln(BlockId(3), &h1_a, &h1_b, true, AnchorBlock::Ref);
    let a2 = new_aln(BlockId(4), &h1_c, &h2_e, true, AnchorBlock::Qry);

    let a3 = new_aln(BlockId(5), &h2_f, &h1_d, false, AnchorBlock::Ref);
    let a4 = new_aln(BlockId(6), &h2_g, &h2_h, false, AnchorBlock::Qry);

    let hits = extract_hits(BlockId(1), &[a1, a2, a3, a4]);

    assert_eq!(
      hits,
      vec![
        create_hit(BlockId(3), true, true, h1_a),
        create_hit(BlockId(3), false, true, h1_b),
        create_hit(BlockId(4), false, true, h1_c),
        create_hit(BlockId(5), false, false, h1_d),
      ]
    );
  }

  #[rustfmt::skip]
  #[test]
  fn test_group_promises() {
    let b1_anchor = PangraphBlock::new(Some(BlockId(1)), "A", btreemap! {} /* {1: None, 2: None, 3: None} */);
    let b1_append = PangraphBlock::new(Some(BlockId(1)), "B", btreemap! {} /* {4: None, 5: None} */);
    let b2_anchor = PangraphBlock::new(Some(BlockId(2)), "C", btreemap! {} /* {6: None, 7: None, 8: None} */);
    let b2_append = PangraphBlock::new(Some(BlockId(2)), "D", btreemap! {} /* {7: None, 8: None} */);
    let b3_anchor = PangraphBlock::new(Some(BlockId(3)), "E", btreemap! {} /* {11: None, 12: None} */);
    let b3_append = PangraphBlock::new(Some(BlockId(3)), "F", btreemap! {} /* {13: None} */);

    let h = &[
      ToMerge::new(b1_anchor.clone(), true, true),
      ToMerge::new(b1_append.clone(), false, true),
      ToMerge::new(b3_anchor.clone(), true, false),
      ToMerge::new(b2_append.clone(), false, true),
      ToMerge::new(b2_anchor.clone(), true, true),
      ToMerge::new(b3_append.clone(), false, false),
    ];

    let promises = group_promises(h);
    assert_eq!(
      promises,
      vec![
        MergePromise::new(b1_anchor, b1_append, true),
        MergePromise::new(b2_anchor, b2_append, true),
        MergePromise::new(b3_anchor, b3_append, false),
      ]
    );
  }

  #[test]
  fn test_assign_anchor_block() {
    fn new_hit(block_id: usize) -> Hit2 {
      Hit2 {
        name: BlockId(block_id),
        length: 0,                     // FIXME
        interval: Interval::default(), // FIXME
      }
    }

    fn new_aln(q: usize, r: usize) -> Alignment2 {
      Alignment2 {
        qry: new_hit(q),
        reff: new_hit(r),

        // FIXME: these were all unset in Python version
        matches: 0,
        length: 0,
        quality: 0,
        orientation: false,
        new_block_id: None,
        anchor_block: None,
        cigar: Cigar::default(),
        divergence: None,
        align: None,
      }
    }

    fn e(nids: &[usize]) -> BTreeMap<NodeId, Edit> {
      nids.iter().map(|nid| (NodeId(*nid), Edit::empty())).collect()
    }

    let b1 = PangraphBlock::new(Some(BlockId(1)), "A", e(&[1, 2, 3]));
    let b2 = PangraphBlock::new(Some(BlockId(2)), "B", e(&[4, 5]));
    let b3 = PangraphBlock::new(Some(BlockId(3)), "C", e(&[6]));
    let b4 = PangraphBlock::new(Some(BlockId(4)), "D", e(&[7, 8, 9, 10]));

    let pangraph = Pangraph {
      blocks: btreemap! {
        BlockId(1) => b1,
        BlockId(2) => b2,
        BlockId(3) => b3,
        BlockId(4) => b4,
      },
      paths: btreemap! {},
      nodes: btreemap! {},
    };

    let mut mergers = vec![new_aln(1, 2), new_aln(3, 4), new_aln(4, 1)];
    assign_anchor_block(&mut mergers, &pangraph);

    assert_eq!(mergers[0].anchor_block, Some(AnchorBlock::Qry));
    assert_eq!(mergers[1].anchor_block, Some(AnchorBlock::Ref));
    assert_eq!(mergers[2].anchor_block, Some(AnchorBlock::Qry));
  }

  #[test]
  fn test_target_blocks() {
    fn new_hit(block_id: usize) -> Hit2 {
      Hit2 {
        name: BlockId(block_id),
        length: 0,                     // FIXME
        interval: Interval::default(), // FIXME
      }
    }

    fn new_aln(qry: Hit2, reff: Hit2) -> Alignment2 {
      Alignment2 {
        qry,
        reff,

        // FIXME: these were all unset in Python version
        matches: 0,
        length: 0,
        quality: 0,
        orientation: false,
        new_block_id: None,
        anchor_block: None,
        cigar: Cigar::default(),
        divergence: None,
        align: None,
      }
    }

    let h1 = new_hit(1);
    let h2 = new_hit(2);
    let h3 = new_hit(3);
    let h4 = new_hit(4);

    let h5 = new_hit(1);
    let h6 = new_hit(2);
    let h7 = new_hit(3);
    let h8 = new_hit(4);

    let a1 = new_aln(h1, h2); // a1 : 1 -- 2
    let a2 = new_aln(h3, h4); // a2 : 3 -- 4
    let a3 = new_aln(h5, h8); // a3 : 1 -- 4
    let a4 = new_aln(h6, h7); // a4 : 2 -- 3

    let tb = target_blocks(&[a1.clone(), a2.clone(), a3.clone(), a4.clone()]);

    let expected = btreemap! {
        BlockId(1) => vec![a1.clone(), a3.clone()],
        BlockId(2) => vec![a1.clone(), a4.clone()],
        BlockId(3) => vec![a2.clone(), a4.clone()],
        BlockId(4) => vec![a2.clone(), a3.clone()],
    };

    assert_eq!(tb, expected);
  }

  #[test]
  fn test_split_block() {
    fn generate_example() -> (Pangraph, Vec<Alignment2>, BlockId) {
      let mut rng = get_random_number_generator(&Some(0));
      let consensus = generate_random_nuc_sequence(130, &mut rng);

      let bid = BlockId(1);
      let nid1 = NodeId(1000);
      let nid2 = NodeId(2000);
      let nid3 = NodeId(3000);

      let n1 = PangraphNode::new(Some(nid1), bid, PathId(100), true, (100, 230));
      let n2 = PangraphNode::new(Some(nid2), bid, PathId(200), false, (1000, 1130));
      let n3 = PangraphNode::new(Some(nid3), bid, PathId(300), false, (180, 110));

      let b1 = PangraphBlock::new(
        Some(bid),
        consensus,
        btreemap! {
          nid1 => Edit::empty(),
          nid2 => Edit::empty(),
          nid3 => Edit::empty(),
        },
      );

      let p1 = PangraphPath::new(Some(PathId(100)), [nid1], 2000, false);
      let p2 = PangraphPath::new(Some(PathId(200)), [nid2], 2000, false);
      let p3 = PangraphPath::new(Some(PathId(300)), [nid3], 200, false);

      let G = Pangraph {
        paths: btreemap! {
          p1.id() => p1,
          p2.id() => p2,
          p3.id() => p3,
        },
        blocks: btreemap! {
          bid => b1,
        },
        nodes: btreemap! {
          nid1 => n1,
          nid2 => n2,
          nid3 => n3,
        },
      };

      #[allow(clippy::items_after_statements)]
      fn h(name: usize, start: usize, stop: usize) -> Hit2 {
        Hit2 {
          name: BlockId(name),
          length: 0, // FIXME
          interval: Interval::new(start, stop),
        }
      }

      #[allow(clippy::items_after_statements)]
      fn a(qry: Hit2, reff: Hit2, strand: bool, new_block_id: usize, anchor_block: AnchorBlock) -> Alignment2 {
        Alignment2 {
          qry,
          reff,
          orientation: strand,
          new_block_id: Some(BlockId(new_block_id)),
          anchor_block: Some(anchor_block),

          // FIXME: these were all unset in Python version
          matches: 0,
          length: 0,
          quality: 0,
          cigar: Cigar::default(),
          divergence: None,
          align: None,
        }
      }

      let m = vec![
        a(h(1, 10, 50), h(2, 100, 150), true, 42, AnchorBlock::Qry),
        a(h(3, 1000, 1050), h(1, 80, 130), false, 43, AnchorBlock::Qry),
      ];

      (G, m, bid)
    }

    let thr_len = 20;
    let (G, m, bid) = generate_example();
    let (u, h) = split_block(bid, &m, &G, thr_len);

    assert_eq!(u.b_old_id, bid);

    let mut node_keys_0 = btreeset! {};
    let mut node_keys_1 = btreeset! {};
    let mut node_keys_2 = btreeset! {};

    let nodes = &u.n_new[&NodeId(1000)];
    node_keys_0.insert(nodes[0].id());
    node_keys_1.insert(nodes[1].id());
    node_keys_2.insert(nodes[2].id());
    assert_eq!(nodes[0].block_id(), BlockId(42));
    assert_eq!(nodes[0].strand(), true);
    assert_eq!(nodes[0].position(), (100, 150));
    assert_eq!(nodes[1].strand(), true);
    assert_eq!(nodes[1].position(), (150, 180));
    assert_eq!(nodes[2].block_id(), BlockId(43));
    assert_eq!(nodes[2].strand(), false);
    assert_eq!(nodes[2].position(), (180, 230));

    let nodes = &u.n_new[&NodeId(2000)];
    node_keys_0.insert(nodes[2].id());
    node_keys_1.insert(nodes[1].id());
    node_keys_2.insert(nodes[0].id());
    assert_eq!(nodes[0].block_id(), BlockId(43));
    assert_eq!(nodes[0].strand(), true);
    assert_eq!(nodes[0].position(), (1000, 1050));
    assert_eq!(nodes[1].strand(), false);
    assert_eq!(nodes[1].position(), (1050, 1080));
    assert_eq!(nodes[2].block_id(), BlockId(42));
    assert_eq!(nodes[2].strand(), false);
    assert_eq!(nodes[2].position(), (1080, 1130));

    let nodes = &u.n_new[&NodeId(3000)];
    node_keys_0.insert(nodes[2].id());
    node_keys_1.insert(nodes[1].id());
    node_keys_2.insert(nodes[0].id());
    assert_eq!(nodes[0].block_id(), BlockId(43));
    assert_eq!(nodes[0].strand(), true);
    assert_eq!(nodes[0].position(), (180, 30));
    assert_eq!(nodes[1].strand(), false);
    assert_eq!(nodes[1].position(), (30, 60));
    assert_eq!(nodes[2].block_id(), BlockId(42));
    assert_eq!(nodes[2].strand(), false);
    assert_eq!(nodes[2].position(), (60, 110));

    assert_eq!(u.b_new.len(), 1);
    let b = &u.b_new[0];
    assert_eq!(b.consensus(), &G.blocks[&bid].consensus()[50..80]);
    assert_eq!(b.alignment_keys(), node_keys_1);

    assert_eq!(h.len(), 2);

    assert_eq!(h[0].block.id(), BlockId(42));
    assert_eq!(h[0].is_anchor, true);
    assert_eq!(h[0].orientation, true);
    assert_eq!(h[1].block.id(), BlockId(43));
    assert_eq!(h[1].is_anchor, false);
    assert_eq!(h[1].orientation, false);

    assert_eq!(h[0].block.consensus(), &G.blocks[&bid].consensus()[0..50]);
    assert_eq!(h[1].block.consensus(), &G.blocks[&bid].consensus()[80..130]);

    assert_eq!(h[0].block.alignment_keys(), node_keys_0);
    assert_eq!(h[1].block.alignment_keys(), node_keys_2);
  }

  #[test]
  fn test_reweave() {
    fn i(pos: usize, len: usize, seq: impl AsRef<str>) -> Ins {
      Ins::new(pos, seq.as_ref().repeat(len))
    }

    fn d(pos: usize, len: usize) -> Del {
      Del::new(pos, len)
    }

    fn s(pos: usize, alt: char) -> Sub {
      Sub::new(pos, alt)
    }

    fn generate_example() -> (Pangraph, Vec<Alignment2>) {
      let nodes = btreemap! {
        NodeId(1) => PangraphNode::new(Some(NodeId(1)), BlockId(10), PathId(100), true, (700, 885)),
        NodeId(2) => PangraphNode::new(Some(NodeId(2)), BlockId(30), PathId(100), true, (885, 988)),
        NodeId(3) => PangraphNode::new(Some(NodeId(3)), BlockId(30), PathId(200), false, (100, 180)),
        NodeId(4) => PangraphNode::new(Some(NodeId(4)), BlockId(20), PathId(200), false, (180, 555)),
        NodeId(5) => PangraphNode::new(Some(NodeId(5)), BlockId(10), PathId(200), false, (555, 735)),
        NodeId(6) => PangraphNode::new(Some(NodeId(6)), BlockId(40), PathId(300), true, (600, 100)),
        NodeId(7) => PangraphNode::new(Some(NodeId(7)), BlockId(50), PathId(300), true, (100, 325)),
        NodeId(8) => PangraphNode::new(Some(NodeId(8)), BlockId(50), PathId(300), false, (325, 580)),
      };

      let paths = btreemap! {
        PathId(100) => PangraphPath::new(Some(PathId(100)), [NodeId(1), NodeId(2)], 1000, false),
        PathId(200) => PangraphPath::new(Some(PathId(200)), [NodeId(3), NodeId(4), NodeId(5)], 1000, false),
        PathId(300) => PangraphPath::new(Some(PathId(300)), [NodeId(6), NodeId(7), NodeId(8)], 1000, false),
      };

      #[rustfmt::skip]
      let ed = btreemap! {
        NodeId(1) => Edit::new([i(150, 10, "T")],                [d(50, 25)],              [s(125, 'G')]            ),
        NodeId(2) => Edit::new([i(50, 3, "G")],                  [],                       []                       ),
        NodeId(3) => Edit::new([i(25, 5, "G")],                  [d(50, 25)],              []                       ),
        NodeId(4) => Edit::new([i(250, 5, "A"), i(300, 5, "A")], [d(100, 25), d(350, 10)], [s(50, 'G'), s(225, 'T')]),
        NodeId(5) => Edit::new([i(200, 5, "A")],                 [d(100, 25)],             [s(25, 'T')]             ),
        NodeId(6) => Edit::new([i(200, 10, "T")],                [d(350, 10)],             [s(100, 'T')]            ),
        NodeId(7) => Edit::new([],                               [d(100, 25)],             [s(50, 'G')]             ),
        NodeId(8) => Edit::new([i(150, 5, "T")],                 [],                       []                       ),
      };

      let mut rng = get_random_number_generator(&Some(0));

      let bseq = btreemap! {
        BlockId(10) => generate_random_nuc_sequence(200, &mut rng),
        BlockId(20) => generate_random_nuc_sequence(400, &mut rng),
        BlockId(30) => generate_random_nuc_sequence(100, &mut rng),
        BlockId(40) => generate_random_nuc_sequence(500, &mut rng),
        BlockId(50) => generate_random_nuc_sequence(250, &mut rng),
      };

      let b = |bid: usize, nids: &[usize]| -> PangraphBlock {
        PangraphBlock::new(
          Some(BlockId(bid)),
          &bseq[&BlockId(bid)],
          nids
            .iter()
            .map(|&nid| (NodeId(nid), ed[&NodeId(nid)].clone()))
            .collect(),
        )
      };

      let blocks = btreemap! {
        BlockId(10) => b(10, &[1, 5]),
        BlockId(20) => b(20, &[4]),
        BlockId(30) => b(30, &[2, 3]),
        BlockId(40) => b(40, &[6]),
        BlockId(50) => b(50, &[7, 8]),
      };

      let G = Pangraph { paths, blocks, nodes };

      #[allow(clippy::items_after_statements)]
      fn h(name: usize, length: usize, start: usize, stop: usize) -> Hit2 {
        Hit2 {
          name: BlockId(name),
          length,
          interval: Interval::new(start, stop),
        }
      }

      #[allow(clippy::items_after_statements)]
      fn a(qry: Hit2, reff: Hit2, strand: bool) -> Alignment2 {
        Alignment2 {
          qry,
          reff,
          orientation: strand,

          // FIXME: these were all unset in Python version
          new_block_id: None,
          anchor_block: None,
          matches: 0,
          length: 0,
          quality: 0,
          cigar: Cigar::default(),
          divergence: None,
          align: None,
        }
      }

      let M = vec![
        a(h(10, 200, 10, 200), h(40, 500, 10, 200), true),
        a(h(20, 400, 0, 200), h(40, 500, 300, 500), false),
        a(h(20, 400, 300, 400), h(50, 250, 0, 100), true),
        a(h(30, 100, 0, 100), h(50, 250, 150, 250), true),
      ];

      (G, M)
    }

    let (mut G, mut M) = generate_example();
    let O = G.clone();
    let thr_len = 90;
    let (G, P) = reweave(&mut M, &mut G, thr_len);

    // new paths
    let p1 = &G.paths[&PathId(100)];
    let p2 = &G.paths[&PathId(200)];
    let p3 = &G.paths[&PathId(300)];

    // new nodes
    assert_eq!(p1.nodes.len(), 2);
    let (nid_100_1, nid_100_2) = (p1.nodes[0], p1.nodes[1]);
    assert_eq!(p2.nodes.len(), 5);
    let (nid_200_1, nid_200_2, nid_200_3, nid_200_4, nid_200_5) =
      (p2.nodes[0], p2.nodes[1], p2.nodes[2], p2.nodes[3], p2.nodes[4]);
    assert_eq!(p3.nodes.len(), 7);
    let (nid_300_1, nid_300_2, nid_300_3, nid_300_4, nid_300_5, nid_300_6, nid_300_7) = (
      p3.nodes[0],
      p3.nodes[1],
      p3.nodes[2],
      p3.nodes[3],
      p3.nodes[4],
      p3.nodes[5],
      p3.nodes[6],
    );

    // check node positions
    assert_eq!(G.nodes[&nid_100_1].position(), O.nodes[&NodeId(1)].position());
    assert_eq!(G.nodes[&nid_100_2].position(), O.nodes[&NodeId(2)].position());
    assert_eq!(G.nodes[&nid_200_1].position(), O.nodes[&NodeId(3)].position());
    assert_eq!(G.nodes[&nid_200_2].position(), (180, 275));
    assert_eq!(G.nodes[&nid_200_3].position(), (275, 380));
    assert_eq!(G.nodes[&nid_200_4].position(), (380, 555));
    assert_eq!(G.nodes[&nid_200_5].position(), O.nodes[&NodeId(5)].position());
    assert_eq!(G.nodes[&nid_300_1].position(), (600, 800));
    assert_eq!(G.nodes[&nid_300_2].position(), (800, 910));
    assert_eq!(G.nodes[&nid_300_3].position(), (910, 100));
    assert_eq!(G.nodes[&nid_300_4].position(), (100, 225));
    assert_eq!(G.nodes[&nid_300_5].position(), (225, 325));
    assert_eq!(G.nodes[&nid_300_6].position(), (325, 430));
    assert_eq!(G.nodes[&nid_300_7].position(), (430, 580));

    // check node orientation
    assert_eq!(G.nodes[&nid_100_1].strand(), true);
    assert_eq!(G.nodes[&nid_100_2].strand(), true);
    assert_eq!(G.nodes[&nid_200_1].strand(), false);
    assert_eq!(G.nodes[&nid_200_2].strand(), false);
    assert_eq!(G.nodes[&nid_200_3].strand(), false);
    assert_eq!(G.nodes[&nid_200_4].strand(), true);
    assert_eq!(G.nodes[&nid_200_5].strand(), false);
    assert_eq!(G.nodes[&nid_300_1].strand(), true);
    assert_eq!(G.nodes[&nid_300_2].strand(), true);
    assert_eq!(G.nodes[&nid_300_3].strand(), true);
    assert_eq!(G.nodes[&nid_300_4].strand(), true);
    assert_eq!(G.nodes[&nid_300_5].strand(), true);
    assert_eq!(G.nodes[&nid_300_6].strand(), false);
    assert_eq!(G.nodes[&nid_300_7].strand(), false);

    // block identity
    let bid10_1 = G.nodes[&nid_100_1].block_id();
    assert_eq!(G.nodes[&nid_200_5].block_id(), bid10_1);
    assert!(!G.blocks.contains_key(&bid10_1)); // missing block key, it is in merge promise
    assert!(P.iter().any(|p| p.anchor_block.id() == bid10_1));

    let bid20_2 = G.nodes[&nid_200_3].block_id();
    assert!(G.blocks.contains_key(&bid20_2)); // block is already in the graph
    assert!(!P
      .iter()
      .any(|p| p.anchor_block.id() == bid20_2 || p.append_block.id() == bid20_2));
    let ed20_2 = &G.blocks[&bid20_2].alignments()[&nid_200_3];
    assert_eq!(ed20_2, &Edit::new([i(50, 5, "A")], [], [s(25, 'T')]));

    let bid20_1 = G.nodes[&nid_200_1].block_id();
    for n in &[nid_100_2, nid_300_5, nid_300_6] {
      assert_eq!(G.nodes[&n].block_id(), bid20_1);
    }

    // merge promises
    assert_eq!(P.len(), 4);
    let p_dict: BTreeMap<_, _> = P.into_iter().map(|p| (p.anchor_block.id(), p)).collect();

    let p1 = &p_dict[&bid10_1];
    assert_eq!(p1.orientation, true);
    assert_eq!(p1.anchor_block.consensus(), O.blocks[&BlockId(10)].consensus());
    assert_eq!(p1.append_block.consensus(), &O.blocks[&BlockId(40)].consensus()[0..200]);
    assert_eq!(p1.append_block.id(), G.nodes[&nid_300_1].block_id());

    let bid40_3 = G.nodes[&nid_300_3].block_id();
    let p2 = &p_dict[&bid40_3];
    assert_eq!(p2.orientation, false);
    assert_eq!(
      p2.anchor_block.consensus(),
      &O.blocks[&BlockId(40)].consensus()[300..500]
    );
    assert_eq!(p2.append_block.id(), G.nodes[&nid_200_4].block_id());
    assert_eq!(p2.append_block.consensus(), &O.blocks[&BlockId(20)].consensus()[0..200]);
  }
}