use crate::annotation::feature::Feature;
use crate::make_internal_report;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_node::NodeId;
use crate::pangraph::pangraph_path::{PangraphPath, PathId};
use crate::pangraph::slice::consensus_coords_from_node_flagged;
use crate::pangraph::strand::Strand;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// A single annotation feature lifted onto **one graph node** (its block consensus).
///
/// A feature overlapping several nodes splits into several `LiftedAnnotation`s — its
/// *segments* — sharing a `parent_feature_id` and ordered by `segment_idx`. This node-level
/// view is the lossless source of truth: every feature instance is placed on the node(s) it
/// overlaps, with block-consensus coordinates obtained by inverting that node's edits.
///
/// Coordinates are 0-based half-open. `node_start`/`node_end` are in the node's
/// **consensus-oriented** sequence (after any reverse-strand flip), i.e. exactly the input to
/// the node→consensus inverse map; `cons_start`/`cons_end` are block-consensus coordinates.
///
/// The per-endpoint boolean flags refer to the **`cons_start` and `cons_end` endpoints of this
/// row** (not the genome 5'/3' ends). On a reverse-strand node the feature's genome-start maps
/// to `cons_end` and its genome-end to `cons_start`; `strand_on_consensus` lets a consumer map
/// back to genome orientation.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct LiftedAnnotation {
  /// Per-segment identifier (unique per row). Equal to `parent_feature_id` for a single-segment
  /// feature, suffixed with `.seg{idx}` when the feature spans multiple nodes.
  pub feature_id: String,

  /// Identifier of the original (un-split) feature, shared across all its segments. Taken from
  /// the GFF `ID` attribute; `None` when the source had no `ID`.
  pub parent_feature_id: Option<String>,

  /// 0-based index of this segment within the feature, in genome (low→high coordinate) order.
  pub segment_idx: usize,

  /// Total number of segments the feature was split into.
  pub n_segments: usize,

  /// Name of the genome (pangraph path) the feature belongs to.
  pub genome: String,

  /// Block whose consensus this segment is placed on.
  pub block_id: BlockId,

  /// Node (block instance on the path) this segment is placed on.
  pub node_id: NodeId,

  /// Strand of the feature relative to the block consensus (the source strand flipped iff the
  /// node is on the reverse strand). `None` when the source feature was unstranded.
  pub strand_on_consensus: Option<Strand>,

  /// Consensus-oriented node-local coordinates `[node_start, node_end)`.
  pub node_start: usize,
  pub node_end: usize,

  /// Block-consensus coordinates `[cons_start, cons_end)`.
  pub cons_start: usize,
  pub cons_end: usize,

  /// Whether the `cons_start` / `cons_end` endpoint is a real feature terminus (vs a fragment
  /// boundary introduced by a node/block split).
  pub start_is_terminus: bool,
  pub end_is_terminus: bool,

  /// Whether the `cons_start` / `cons_end` endpoint fell strictly inside an insertion (node
  /// bases absent from the consensus) and was snapped to the insertion anchor. A feature wholly
  /// inside an insertion has `cons_start == cons_end` with both flags set.
  pub start_in_insertion: bool,
  pub end_in_insertion: bool,

  /// Fraction of the original feature's length retained in this segment.
  pub frac_covered: f64,

  /// Feature type, e.g. `"CDS"` or `"gene"`.
  pub feature_type: String,

  /// Human-readable name (GFF `Name` attribute), if any.
  pub name: Option<String>,

  /// Flexible key→value metadata carried from the source feature.
  pub attributes: Vec<(String, String)>,
}

/// A feature's overlap with one node, in genome-oriented coordinates, before strand handling.
struct RawSegment {
  node_id: NodeId,
  /// Genome coordinates of the overlap `[genome_start, genome_end)`.
  genome_start: usize,
  genome_end: usize,
  /// Genome-oriented node-local coordinates of the overlap `[node_a, node_b)`.
  node_a: usize,
  node_b: usize,
}

/// Genome coverage piece(s) of a node on its path, each as `(genome_start, genome_end,
/// node_offset)` where `node_offset` is the genome-oriented node-local coordinate at
/// `genome_start`.
///
/// A node normally covers one interval `[p0, p1)`. On a circular path it may wrap the origin
/// (`p0 > p1`, covering `[p0, tot_len)` then `[0, p1)`) or span the whole circle (`p0 == p1`).
fn node_coverage_pieces(p0: usize, p1: usize, tot_len: usize, circular: bool) -> Vec<(usize, usize, usize)> {
  if p0 < p1 {
    return vec![(p0, p1, 0)];
  }
  if !circular {
    // A non-circular node always has p0 < p1; anything else is degenerate/empty.
    return vec![];
  }
  // Circular: `p0 == p1` marks a whole-circle node (its low piece runs up to p0); a genuine
  // origin wrap has `p0 > p1` (low piece runs up to p1).
  let second_end = if p0 == p1 { p0 } else { p1 };
  let mut pieces = vec![(p0, tot_len, 0)];
  if second_end > 0 {
    pieces.push((0, second_end, tot_len - p0));
  }
  pieces
}

/// Build the per-segment `feature_id` from the parent id (or a coordinate-based fallback when
/// the source had no `ID`), suffixing multi-segment features so each row is unique.
fn make_feature_id(parent: Option<&str>, feature: &Feature, segment_idx: usize, n_segments: usize) -> String {
  let base = parent.map_or_else(
    || format!("{}:{}-{}", feature.seqid, feature.interval.start, feature.interval.end),
    ToOwned::to_owned,
  );
  if n_segments > 1 {
    format!("{base}.seg{segment_idx}")
  } else {
    base
  }
}

/// Lift a single annotation [`Feature`] onto the nodes of the path it belongs to.
///
/// Returns one [`LiftedAnnotation`] per node the feature overlaps, in genome (low→high
/// coordinate) order. The feature interval is assumed to be a non-wrapping `[start, end)`
/// (origin-spanning features are encoded as separate GFF records); nodes, however, may wrap the
/// circular origin, which is handled here.
pub fn lift_feature(feature: &Feature, path: &PangraphPath, graph: &Pangraph) -> Result<Vec<LiftedAnnotation>, Report> {
  let tot_len = path.tot_len();
  let circular = path.circular();
  let genome = path.name().clone().unwrap_or_else(|| path.id().to_string());

  let f_s = feature.interval.start;
  let f_e = feature.interval.end;
  if f_e <= f_s {
    return Ok(vec![]); // empty feature: nothing to lift
  }

  // 1+2. Collect the feature's overlap with each node, then order by genome coordinate.
  let mut raw_segments: Vec<RawSegment> = Vec::new();
  for &node_id in path.nodes() {
    let node = graph
      .nodes
      .get(&node_id)
      .ok_or_else(|| make_internal_report!("When lifting feature: node {node_id} not found in graph"))?;
    let (p0, p1) = node.position();
    for (g_start, g_end, node_off) in node_coverage_pieces(p0, p1, tot_len, circular) {
      let ov_s = f_s.max(g_start);
      let ov_e = f_e.min(g_end);
      if ov_s < ov_e {
        raw_segments.push(RawSegment {
          node_id,
          genome_start: ov_s,
          genome_end: ov_e,
          node_a: node_off + (ov_s - g_start),
          node_b: node_off + (ov_e - g_start),
        });
      }
    }
  }
  raw_segments.sort_by_key(|s| s.genome_start);

  let n_segments = raw_segments.len();
  let feature_len = (f_e - f_s) as f64;
  let parent = feature.id.as_deref();

  let mut out = Vec::with_capacity(n_segments);
  for (segment_idx, seg) in raw_segments.into_iter().enumerate() {
    let node = &graph.nodes[&seg.node_id];
    let block_id = node.block_id();
    let block = &graph.blocks[&block_id];
    let block_l = block.consensus_len();
    let edits = block.alignment(seg.node_id);
    let len_node = block.unaligned_len_for_node(&seg.node_id);

    // 3. Strand: convert genome-oriented node coords to consensus orientation, flipping the
    // feature strand on reverse-strand nodes.
    let reverse = node.strand().is_reverse();
    let (node_start, node_end, strand_on_consensus) = if reverse {
      (
        len_node - seg.node_b,
        len_node - seg.node_a,
        feature.strand.map(|s| s.reverse()),
      )
    } else {
      (seg.node_a, seg.node_b, feature.strand)
    };

    // 4. Node → consensus, keeping the per-endpoint in-insertion flags.
    let ((cons_start, start_in_insertion), (cons_end, end_in_insertion)) =
      consensus_coords_from_node_flagged((node_start, node_end), edits, block_l);

    // Terminus flags, in the consensus-endpoint frame. The genome-low feature end (f_s) lives in
    // the segment whose genome overlap starts at f_s; the genome-high end (f_e) in the segment
    // whose overlap ends at f_e. Reverse-strand nodes swap which consensus endpoint each maps to.
    let holds_feature_start = seg.genome_start == f_s;
    let holds_feature_end = seg.genome_end == f_e;
    let (start_is_terminus, end_is_terminus) = if reverse {
      (holds_feature_end, holds_feature_start)
    } else {
      (holds_feature_start, holds_feature_end)
    };

    let frac_covered = (seg.genome_end - seg.genome_start) as f64 / feature_len;
    let feature_id = make_feature_id(parent, feature, segment_idx, n_segments);

    out.push(LiftedAnnotation {
      feature_id,
      parent_feature_id: feature.id.clone(),
      segment_idx,
      n_segments,
      genome: genome.clone(),
      block_id,
      node_id: seg.node_id,
      strand_on_consensus,
      node_start,
      node_end,
      cons_start,
      cons_end,
      start_is_terminus,
      end_is_terminus,
      start_in_insertion,
      end_in_insertion,
      frac_covered,
      feature_type: feature.feature_type.clone(),
      name: feature.name.clone(),
      attributes: feature.attributes.clone(),
    });
  }

  Ok(out)
}

/// Lift every matched feature onto the graph, producing the flat node-level annotation table.
///
/// `grouped` is the output of [`match_features_to_paths`](crate::annotation::matching::match_features_to_paths):
/// features already keyed by the path they belong to.
pub fn lift_features(
  grouped: &BTreeMap<PathId, Vec<Feature>>,
  graph: &Pangraph,
) -> Result<Vec<LiftedAnnotation>, Report> {
  let mut out = Vec::new();
  for (path_id, features) in grouped {
    let path = graph
      .paths
      .get(path_id)
      .ok_or_else(|| make_internal_report!("When lifting features: path {path_id} not found in graph"))?;
    for feature in features {
      out.extend(lift_feature(feature, path, graph)?);
    }
  }
  Ok(out)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::{Del, Edit, Ins, Sub};
  use crate::pangraph::pangraph_block::PangraphBlock;
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use crate::utils::interval::Interval;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  /// Spec for one node (its own block) when building a single-path test graph.
  struct NodeSpec {
    consensus: &'static str,
    edits: Edit,
    strand: Strand,
    position: (usize, usize),
  }

  /// Build a one-path graph from node specs. Ids are deterministic: `NodeId(i)`, `BlockId(i)`,
  /// `PathId(0)`. The caller is responsible for making each node's `position` span equal its
  /// reconstructed length (consensus_len + insertions − deletions).
  fn build_graph(name: &str, tot_len: usize, circular: bool, specs: Vec<NodeSpec>) -> (Pangraph, PangraphPath) {
    let mut nodes = BTreeMap::new();
    let mut blocks = BTreeMap::new();
    let mut node_ids = Vec::new();
    for (i, spec) in specs.into_iter().enumerate() {
      let nid = NodeId(i);
      let bid = BlockId(i);
      let node = PangraphNode::new(Some(nid), bid, PathId(0), spec.strand, spec.position);
      let block = PangraphBlock::new(bid, spec.consensus, btreemap! { nid => spec.edits });
      nodes.insert(nid, node);
      blocks.insert(bid, block);
      node_ids.push(nid);
    }
    let path = PangraphPath::new(
      Some(PathId(0)),
      node_ids,
      tot_len,
      circular,
      Some(name.to_owned()),
      None,
    );
    let graph = Pangraph {
      paths: btreemap! { PathId(0) => path.clone() },
      blocks,
      nodes,
    };
    (graph, path)
  }

  fn feat(start: usize, end: usize, strand: Option<Strand>) -> Feature {
    Feature {
      seqid: "p".to_owned(),
      source: None,
      feature_type: "CDS".to_owned(),
      interval: Interval::new(start, end),
      strand,
      id: Some("g1".to_owned()),
      name: Some("geneA".to_owned()),
      attributes: vec![("ID".to_owned(), "g1".to_owned())],
    }
  }

  /// Assert the coordinate-bearing fields of a lifted segment in one place.
  fn assert_seg(
    a: &LiftedAnnotation,
    node: (usize, usize),
    cons: (usize, usize),
    strand: Option<Strand>,
    termini: (bool, bool),
    in_ins: (bool, bool),
  ) {
    assert_eq!((a.node_start, a.node_end), node, "node coords");
    assert_eq!((a.cons_start, a.cons_end), cons, "consensus coords");
    assert_eq!(a.strand_on_consensus, strand, "strand_on_consensus");
    assert_eq!((a.start_is_terminus, a.end_is_terminus), termini, "termini");
    assert_eq!((a.start_in_insertion, a.end_in_insertion), in_ins, "in_insertion");
  }

  #[test]
  fn test_lift_forward_single_node_substitution_only() {
    // Substitutions never move coordinates, so the lift is the identity.
    let (graph, path) = build_graph(
      "p",
      10,
      false,
      vec![NodeSpec {
        consensus: "ACGTACGTAC",
        edits: Edit::new(vec![], vec![], vec![Sub::new(3, 'T')]),
        strand: Forward,
        position: (0, 10),
      }],
    );
    let lifted = lift_feature(&feat(2, 7, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 1);
    let a = &lifted[0];
    assert_seg(a, (2, 7), (2, 7), Some(Forward), (true, true), (false, false));
    assert_eq!(a.feature_id, "g1");
    assert_eq!(a.parent_feature_id.as_deref(), Some("g1"));
    assert_eq!((a.segment_idx, a.n_segments), (0, 1));
    assert!((a.frac_covered - 1.0).abs() < 1e-9);
    assert_eq!(a.block_id, BlockId(0));
    assert_eq!(a.node_id, NodeId(0));
  }

  #[test]
  fn test_lift_forward_deletion_inside_span_shifts_cons_end() {
    // A deletion of 3 inside the feature span stretches the consensus interval by 3.
    let (graph, path) = build_graph(
      "p",
      9,
      false,
      vec![NodeSpec {
        consensus: "ACGTACGTACGT",
        edits: Edit::new(vec![], vec![Del::new(4, 3)], vec![]),
        strand: Forward,
        position: (0, 9),
      }],
    );
    let lifted = lift_feature(&feat(1, 8, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 1);
    assert_seg(&lifted[0], (1, 8), (1, 11), Some(Forward), (true, true), (false, false));
  }

  #[test]
  fn test_lift_forward_insertion_endpoint_inside_insertion() {
    // Feature end lands inside an inserted region -> snaps to the insertion anchor and flags it.
    let (graph, path) = build_graph(
      "p",
      11,
      false,
      vec![NodeSpec {
        consensus: "ACGTACGT",
        edits: Edit::new(vec![Ins::new(4, "TTT")], vec![], vec![]),
        strand: Forward,
        position: (0, 11),
      }],
    );
    let lifted = lift_feature(&feat(3, 6, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 1);
    assert_seg(&lifted[0], (3, 6), (3, 4), Some(Forward), (true, true), (false, true));
  }

  #[test]
  fn test_lift_forward_feature_wholly_inside_insertion() {
    // The whole feature sits in inserted bases -> empty consensus interval, both flags set.
    let (graph, path) = build_graph(
      "p",
      11,
      false,
      vec![NodeSpec {
        consensus: "ACGTACGT",
        edits: Edit::new(vec![Ins::new(4, "TTT")], vec![], vec![]),
        strand: Forward,
        position: (0, 11),
      }],
    );
    let lifted = lift_feature(&feat(5, 6, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 1);
    assert_seg(&lifted[0], (5, 6), (4, 4), Some(Forward), (true, true), (true, true));
  }

  #[test]
  fn test_lift_reverse_single_node_flips_offsets_and_strand() {
    // On a reverse node the genome-oriented offsets flip and the feature strand is reversed.
    let (graph, path) = build_graph(
      "p",
      10,
      false,
      vec![NodeSpec {
        consensus: "AAAACCCCGG",
        edits: Edit::empty(),
        strand: Reverse,
        position: (0, 10),
      }],
    );
    let lifted = lift_feature(&feat(2, 5, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 1);
    assert_seg(&lifted[0], (5, 8), (5, 8), Some(Reverse), (true, true), (false, false));
  }

  #[test]
  fn test_lift_multi_node_segments_termini_and_frac() {
    // A feature spanning two forward nodes splits into two segments; only the outer endpoints are
    // termini and frac_covered sums to 1.
    let (graph, path) = build_graph(
      "p",
      20,
      false,
      vec![
        NodeSpec {
          consensus: "ACGTACGTAC",
          edits: Edit::empty(),
          strand: Forward,
          position: (0, 10),
        },
        NodeSpec {
          consensus: "TGCATGCATG",
          edits: Edit::empty(),
          strand: Forward,
          position: (10, 20),
        },
      ],
    );
    let lifted = lift_feature(&feat(5, 15, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 2);

    let s0 = &lifted[0];
    assert_eq!((s0.segment_idx, s0.n_segments), (0, 2));
    assert_eq!(s0.node_id, NodeId(0));
    assert_eq!(s0.feature_id, "g1.seg0");
    assert_eq!(s0.parent_feature_id.as_deref(), Some("g1"));
    assert_seg(s0, (5, 10), (5, 10), Some(Forward), (true, false), (false, false));
    assert!((s0.frac_covered - 0.5).abs() < 1e-9);

    let s1 = &lifted[1];
    assert_eq!((s1.segment_idx, s1.n_segments), (1, 2));
    assert_eq!(s1.node_id, NodeId(1));
    assert_eq!(s1.feature_id, "g1.seg1");
    assert_seg(s1, (0, 5), (0, 5), Some(Forward), (false, true), (false, false));
    assert!((s1.frac_covered - 0.5).abs() < 1e-9);
  }

  #[test]
  fn test_lift_reverse_multi_node_terminus_mapping() {
    // Two reverse nodes: the genome-low feature end (f_s) maps to the cons_end of its segment, so
    // terminus flags land on the opposite consensus endpoints from the forward case.
    let (graph, path) = build_graph(
      "p",
      20,
      false,
      vec![
        NodeSpec {
          consensus: "ACGTACGTAC",
          edits: Edit::empty(),
          strand: Reverse,
          position: (0, 10),
        },
        NodeSpec {
          consensus: "TGCATGCATG",
          edits: Edit::empty(),
          strand: Reverse,
          position: (10, 20),
        },
      ],
    );
    let lifted = lift_feature(&feat(5, 15, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 2);
    // Segment over node 0 holds the feature's genome start (f_s=5) -> its cons_end is the terminus.
    assert_seg(&lifted[0], (0, 5), (0, 5), Some(Reverse), (false, true), (false, false));
    // Segment over node 1 holds the feature's genome end (f_e=15) -> its cons_start is the terminus.
    assert_seg(
      &lifted[1],
      (5, 10),
      (5, 10),
      Some(Reverse),
      (true, false),
      (false, false),
    );
  }

  #[test]
  fn test_lift_circular_wrapping_node_high_side() {
    // Origin-wrapping node [15,5) on a circular path of length 20; a feature in the high part
    // [15,20) gets node-local offsets measured from p0.
    let (graph, path) = wrapping_graph();
    let lifted = lift_feature(&feat(16, 19, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 1);
    assert_eq!(lifted[0].node_id, NodeId(0));
    assert_seg(&lifted[0], (1, 4), (1, 4), Some(Forward), (true, true), (false, false));
  }

  #[test]
  fn test_lift_circular_wrapping_node_low_side() {
    // The same wrapping node also covers [0,5) after the origin: genome position 1 is node-local 6.
    let (graph, path) = wrapping_graph();
    let lifted = lift_feature(&feat(1, 4, Some(Forward)), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 1);
    assert_eq!(lifted[0].node_id, NodeId(0));
    assert_seg(&lifted[0], (6, 9), (6, 9), Some(Forward), (true, true), (false, false));
  }

  /// Circular path of length 20 with an origin-wrapping node 0 (position `(15,5)`, covering
  /// `[15,20)∪[0,5)`) and a normal node 1 (`(5,15)`).
  fn wrapping_graph() -> (Pangraph, PangraphPath) {
    build_graph(
      "p",
      20,
      true,
      vec![
        NodeSpec {
          consensus: "ACGTACGTAC",
          edits: Edit::empty(),
          strand: Forward,
          position: (15, 5),
        },
        NodeSpec {
          consensus: "TGCATGCATG",
          edits: Edit::empty(),
          strand: Forward,
          position: (5, 15),
        },
      ],
    )
  }

  #[test]
  fn test_lift_unstranded_feature_stays_unstranded() {
    // An unstranded feature stays unstranded on the consensus, even on a reverse node.
    let (graph, path) = build_graph(
      "p",
      10,
      false,
      vec![NodeSpec {
        consensus: "AAAACCCCGG",
        edits: Edit::empty(),
        strand: Reverse,
        position: (0, 10),
      }],
    );
    let lifted = lift_feature(&feat(2, 5, None), &path, &graph).unwrap();
    assert_eq!(lifted.len(), 1);
    assert_eq!(lifted[0].strand_on_consensus, None);
  }

  #[test]
  fn test_lift_feature_without_id_uses_coordinate_fallback() {
    let (graph, path) = build_graph(
      "p",
      10,
      false,
      vec![NodeSpec {
        consensus: "ACGTACGTAC",
        edits: Edit::empty(),
        strand: Forward,
        position: (0, 10),
      }],
    );
    let mut f = feat(2, 7, Some(Forward));
    f.id = None;
    let lifted = lift_feature(&f, &path, &graph).unwrap();
    assert_eq!(lifted[0].parent_feature_id, None);
    assert_eq!(lifted[0].feature_id, "p:2-7");
  }

  #[test]
  fn test_lift_features_groups_over_paths() {
    let (graph, path) = build_graph(
      "p",
      10,
      false,
      vec![NodeSpec {
        consensus: "ACGTACGTAC",
        edits: Edit::empty(),
        strand: Forward,
        position: (0, 10),
      }],
    );
    let grouped = btreemap! { path.id() => vec![feat(1, 4, Some(Forward)), feat(5, 9, Some(Forward))] };
    let lifted = lift_features(&grouped, &graph).unwrap();
    assert_eq!(lifted.len(), 2);
  }
}
