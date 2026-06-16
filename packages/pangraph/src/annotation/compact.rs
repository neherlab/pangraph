use crate::annotation::lift::LiftedAnnotation;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_path::PangraphPath;
use crate::pangraph::strand::Strand;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// One segment of a feature's block-consensus crossing: a block, its consensus coordinates
/// `[cons_start, cons_end)`, and the strand relative to that block's consensus.
type Segment = (BlockId, usize, usize, Option<Strand>);

/// Identity of a block-level cluster: the feature type plus the feature's whole crossing as a
/// canonical, 5'→3'-ordered multiset of [`Segment`]s. Two per-genome instances are "the same"
/// consensus feature exactly when these match — every block, coordinate and per-block strand along
/// the crossing agrees (and, for a duplicated block crossed twice, the multiplicity agrees too).
type ClusterKey = (String, Vec<Segment>);

/// Identity of one placement on a single block, used to tally per-segment support across the whole
/// node-level table (independent of clustering): `(feature_type, block, cons_start, cons_end, strand)`.
type SegmentKey = (String, BlockId, usize, usize, Option<Strand>);

/// Per-genome metadata kept for one supporter of a cluster: its feature `name` and `attributes`.
type SupporterMeta = (Option<String>, Vec<(String, String)>);

/// A block-level consensus annotation: **one segment** of a feature crossing that recurs, at the
/// *same* block-consensus coordinates, across enough of the genomes that carry it. This is the
/// compacted, opinionated view built on top of the lossless node-level table (it never replaces it).
///
/// A multi-block feature emits one row per segment, all sharing a `cluster_id` and ordered 5'→3' by
/// `segment_idx` (a duplicated block crossed twice yields two rows, distinguished by `segment_idx`).
/// `n_support`/`n_total` are **crossing-level** (identical across a cluster's rows);
/// `n_support_segment`/`n_total_segment` are **this segment's** reproducibility on its own block.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct BlockAnnotation {
  /// Feature type shared by the cluster, e.g. `"CDS"` or `"gene"`.
  pub feature_type: String,

  /// Identifier of the consensus feature (crossing) this segment belongs to; shared by the rows of
  /// one cluster, assigned `0..` in deterministic key order over the emitted clusters.
  pub cluster_id: usize,

  /// 0-based position of this segment within the feature, ordered 5'→3' (or, for unstranded
  /// features, by the canonical multiset order).
  pub segment_idx: usize,

  /// Total number of segments in the crossing.
  pub n_segments: usize,

  /// Block whose consensus this segment lies on.
  pub block_id: BlockId,

  /// This segment's block-consensus coordinates `[cons_start, cons_end)`.
  pub cons_start: usize,
  pub cons_end: usize,

  /// Strand of this segment relative to its block consensus (`None` if unstranded).
  pub strand_on_consensus: Option<Strand>,

  /// Majority feature name across the supporters, when it clears the property threshold; else `None`.
  pub consensus_name: Option<String>,

  /// Per-attribute majority values that clear the property threshold (key-sorted).
  pub consensus_attributes: Vec<(String, String)>,

  /// `M` — genomes sharing this exact crossing (crossing-level; same on every row of the cluster).
  pub n_support: usize,

  /// `N` — genomes structurally capable of the crossing, i.e. traversing every block in it with the
  /// required multiplicity (crossing-level; same on every row of the cluster).
  pub n_total: usize,

  /// `M_seg` — genomes placing a feature-segment at exactly this `(block, coords, strand)`, counted
  /// across the whole node-level table (so a body segment shared by several crossings reports the
  /// same value in each).
  pub n_support_segment: usize,

  /// `N_seg` — distinct genomes traversing this segment's block (its depth).
  pub n_total_segment: usize,
}

/// A pluggable strategy for compacting node-level annotations to the block level.
///
/// Compaction is an opinionated, configurable summary built on top of the node-level lift; making
/// it a trait keeps the policy swappable (coordinate-consensus today; name- or ortholog-based
/// strategies could be added later) without touching the lift or the writers.
pub trait BlockCompactionStrategy {
  /// Compact node-level [`LiftedAnnotation`]s into block-level [`BlockAnnotation`]s.
  fn compact(&self, node_annotations: &[LiftedAnnotation], graph: &Pangraph) -> Result<Vec<BlockAnnotation>, Report>;
}

/// Compact by **block-consensus coordinate agreement**.
///
/// Each genome's feature instance is reduced to its whole crossing — the 5'→3'-ordered multiset of
/// per-segment `(block, cons_start, cons_end, strand)` placements; genomes whose `(feature_type,
/// crossing)` match exactly form a cluster. A cluster is emitted when its support `M` reaches
/// `ceil(min_frequency * N)` (at least 1), where `N` is the number of genomes structurally capable
/// of the crossing (traversing every block in it with the required multiplicity). Each emitted
/// cluster produces one row per segment. Consensus `name`/attributes are the per-value majorities
/// among the supporters that clear `property_threshold`.
pub struct CoordinateConsensusStrategy {
  /// Minimum supporter fraction `M / N` (of genomes capable of the crossing) for a cluster to be emitted.
  pub min_frequency: f64,

  /// Minimum supporter fraction for a metadata value (name / attribute) to be promoted to consensus.
  pub property_threshold: f64,
}

impl Default for CoordinateConsensusStrategy {
  /// Provisional defaults; the user-facing defaults are owned by the CLI (P5.2).
  fn default() -> Self {
    Self {
      min_frequency: 0.9,
      property_threshold: 0.5,
    }
  }
}

impl BlockCompactionStrategy for CoordinateConsensusStrategy {
  fn compact(&self, node_annotations: &[LiftedAnnotation], graph: &Pangraph) -> Result<Vec<BlockAnnotation>, Report> {
    // Aggregate per cluster key: dedup supporters by genome (a genome counts once) and reduce the
    // per-block node requirement to its element-wise minimum over them.
    let mut clusters: BTreeMap<ClusterKey, ClusterAgg> = BTreeMap::new();
    for inst in feature_instances(node_annotations) {
      let agg = clusters.entry(inst.key).or_default();
      agg
        .supporters
        .entry(inst.genome)
        .or_insert((inst.name, inst.attributes));
      for (block_id, req) in inst.node_req {
        agg
          .node_req
          .entry(block_id)
          .and_modify(|cur| *cur = (*cur).min(req))
          .or_insert(req);
      }
    }

    // Per-block genome multiplicities (for crossing-level `N` and per-segment `N_seg`) and the
    // per-segment support tally (for `M_seg`), both independent of the clustering above.
    let block_depths = block_genome_counts(graph);
    let segment_support = segment_support_sets(node_annotations);

    // Iterating the BTreeMap yields clusters in key order, so `cluster_id` (assigned only to emitted
    // clusters) and the row order are deterministic.
    let mut out = Vec::new();
    let mut cluster_id = 0;
    for (key, agg) in clusters {
      let (feature_type, segments) = key;
      let ClusterAgg { supporters, node_req } = agg;
      let m = supporters.len();
      let n = crossing_n_total(&block_depths, &node_req);
      let min_support = min_count(self.min_frequency, n).max(1);
      if m < min_support {
        continue; // below threshold: stays only in the node-level table
      }
      let consensus_name = majority_name(&supporters, m, self.property_threshold);
      let consensus_attributes = majority_attributes(&supporters, m, self.property_threshold);
      let n_segments = segments.len();
      for (segment_idx, &(block_id, cons_start, cons_end, strand)) in segments.iter().enumerate() {
        let seg_key = (feature_type.clone(), block_id, cons_start, cons_end, strand);
        let n_support_segment = segment_support.get(&seg_key).map_or(0, BTreeSet::len);
        let n_total_segment = block_depths.get(&block_id).map_or(0, BTreeMap::len);
        out.push(BlockAnnotation {
          feature_type: feature_type.clone(),
          cluster_id,
          segment_idx,
          n_segments,
          block_id,
          cons_start,
          cons_end,
          strand_on_consensus: strand,
          consensus_name: consensus_name.clone(),
          consensus_attributes: consensus_attributes.clone(),
          n_support: m,
          n_total: n,
          n_support_segment,
          n_total_segment,
        });
      }
      cluster_id += 1;
    }
    Ok(out)
  }
}

/// One genome's instance of a feature, reduced to its cluster key (the whole crossing) plus the
/// metadata compaction needs (representative `name`/`attributes`, identical across a feature's
/// segments) and the **distinct nodes used per block** (for the structural-capability denominator).
struct FeatureInstance {
  genome: String,
  key: ClusterKey,
  /// Distinct nodes the crossing uses on each block — *not* segment count, so an origin-spanning
  /// node split into two coverage pieces still requires only one node of its block.
  node_req: BTreeMap<BlockId, usize>,
  name: Option<String>,
  attributes: Vec<(String, String)>,
}

/// Aggregated supporters of one cluster, plus the per-block node requirement reduced across them.
#[derive(Default)]
struct ClusterAgg {
  /// Supporters keyed by genome (a genome counts once), each with its representative metadata.
  supporters: BTreeMap<String, SupporterMeta>,
  /// Per-block minimum distinct-node requirement over the supporters — the loosest hosting, so
  /// every supporter satisfies it and `N >= M` always holds.
  node_req: BTreeMap<BlockId, usize>,
}

/// Group node-level rows into per-genome feature instances and reduce each to its cluster key.
///
/// Rows are grouped by `(genome, feature_type, base feature id)`; instances that are not bounded by
/// exactly two real termini (e.g. partial / truncated features) are dropped from compaction and
/// remain in the node-level table.
fn feature_instances(node_annotations: &[LiftedAnnotation]) -> Vec<FeatureInstance> {
  let mut groups: BTreeMap<(String, String, String), Vec<&LiftedAnnotation>> = BTreeMap::new();
  for a in node_annotations {
    let group_key = (a.genome.clone(), a.feature_type.clone(), base_feature_id(a));
    groups.entry(group_key).or_default().push(a);
  }

  groups
    .into_iter()
    .filter_map(|((genome, feature_type, _base), rows)| reduce_instance(genome, feature_type, &rows))
    .collect()
}

/// Reduce one feature instance's node-level rows to its [`FeatureInstance`], or `None` when the
/// crossing is not bounded by exactly two real termini (a partial / truncated feature).
///
/// The cluster key is the feature's whole crossing as a 5'→3'-ordered multiset of
/// `(block, cons_start, cons_end, strand_on_consensus)`. Node-level `segment_idx` is genome
/// (low→high coordinate) order, which is already 5'→3' for a forward feature; a reverse feature is
/// reversed, and an unstranded feature (no reading direction) is canonicalized by the
/// lexicographically smaller of the two orders so homologous instances still collapse.
fn reduce_instance(genome: String, feature_type: String, rows: &[&LiftedAnnotation]) -> Option<FeatureInstance> {
  // A clean feature is bounded by exactly two real termini (its 5' and 3' ends); the internal
  // boundaries where it crosses node/block splits are not termini. Anything else is partial.
  let terminus_count: usize = rows
    .iter()
    .map(|a| usize::from(a.start_is_terminus) + usize::from(a.end_is_terminus))
    .sum();
  if terminus_count != 2 {
    return None;
  }

  // Segments in genome (arc) order, then oriented 5'→3'.
  let mut rows_sorted: Vec<&&LiftedAnnotation> = rows.iter().collect();
  rows_sorted.sort_by_key(|a| a.segment_idx);
  let mut segments: Vec<Segment> = rows_sorted
    .iter()
    .map(|a| (a.block_id, a.cons_start, a.cons_end, a.strand_on_consensus))
    .collect();
  match rows_sorted[0].feature_strand {
    Some(Strand::Reverse) => segments.reverse(),
    Some(Strand::Forward) => {}, // genome order is already 5'→3'
    None => {
      let mut reversed = segments.clone();
      reversed.reverse();
      if reversed < segments {
        segments = reversed;
      }
    },
  }

  // Distinct nodes used per block (not segments): an origin-spanning node split into two coverage
  // pieces still needs only one node of its block, while a genuinely duplicated block needs as many.
  let mut node_sets: BTreeMap<BlockId, BTreeSet<_>> = BTreeMap::new();
  for a in rows {
    node_sets.entry(a.block_id).or_default().insert(a.node_id);
  }
  let node_req = node_sets.into_iter().map(|(b, nodes)| (b, nodes.len())).collect();

  // Name/attributes are identical across a feature's segments; take them from the 5'-most.
  let rep = rows_sorted[0];
  Some(FeatureInstance {
    genome,
    key: (feature_type, segments),
    node_req,
    name: rep.name.clone(),
    attributes: rep.attributes.clone(),
  })
}

/// Recover the base feature id shared across a feature's segments: its `parent_feature_id` when
/// present, else the per-row `feature_id` with any `.seg{idx}` suffix stripped.
fn base_feature_id(a: &LiftedAnnotation) -> String {
  if let Some(parent) = &a.parent_feature_id {
    return parent.clone();
  }
  if a.n_segments > 1 {
    let suffix = format!(".seg{}", a.segment_idx);
    if let Some(base) = a.feature_id.strip_suffix(&suffix) {
      return base.to_owned();
    }
  }
  a.feature_id.clone()
}

/// The genome label for a path, matching the `genome` field the node-level lift produces (the path
/// name, or its id when unnamed).
fn genome_label(path: &PangraphPath) -> String {
  path.name().clone().unwrap_or_else(|| path.id().to_string())
}

/// Per block, the number of nodes (block instances) each genome carries of it. A genome appears with
/// a count > 1 exactly when the block is duplicated on its path; `map.len()` is the block's depth
/// (distinct genomes traversing it). Drives both crossing-level `N` and per-segment `N_seg`.
fn block_genome_counts(graph: &Pangraph) -> BTreeMap<BlockId, BTreeMap<String, usize>> {
  graph
    .blocks
    .iter()
    .map(|(&bid, block)| {
      let mut counts: BTreeMap<String, usize> = BTreeMap::new();
      for pid in block.isolates(graph) {
        if let Some(path) = graph.paths.get(&pid) {
          *counts.entry(genome_label(path)).or_default() += 1;
        }
      }
      (bid, counts)
    })
    .collect()
}

/// `N` for a crossing: genomes **structurally capable** of it — those traversing every block in the
/// crossing at least `required` times (the crossing's per-block distinct-node count, so an
/// origin-split node is not double-counted). With all requirements 1 this is the intersection of the
/// per-block genome sets. `required` is reduced to its element-wise min over the supporters, so every
/// supporter is capable and `N >= M` holds.
fn crossing_n_total(
  block_depths: &BTreeMap<BlockId, BTreeMap<String, usize>>,
  required: &BTreeMap<BlockId, usize>,
) -> usize {
  let empty = BTreeMap::new();
  let mut blocks = required.iter();
  let Some((first_block, &first_mult)) = blocks.next() else {
    return 0;
  };
  // Genomes carrying the first block enough times, then narrowed by every remaining block.
  let mut candidates: BTreeSet<&String> = block_depths
    .get(first_block)
    .unwrap_or(&empty)
    .iter()
    .filter(|&(_, &count)| count >= first_mult)
    .map(|(genome, _)| genome)
    .collect();
  for (block_id, &mult) in blocks {
    let counts = block_depths.get(block_id).unwrap_or(&empty);
    candidates.retain(|genome| counts.get(*genome).is_some_and(|&count| count >= mult));
  }
  candidates.len()
}

/// Per-segment support across the whole node-level table: for each distinct placement
/// `(feature_type, block, cons_start, cons_end, strand)`, the set of genomes that place a
/// feature-segment there. Independent of clustering, so a body segment shared by several crossings
/// reports the same support in each.
fn segment_support_sets(node_annotations: &[LiftedAnnotation]) -> BTreeMap<SegmentKey, BTreeSet<String>> {
  let mut map: BTreeMap<SegmentKey, BTreeSet<String>> = BTreeMap::new();
  for a in node_annotations {
    let key = (
      a.feature_type.clone(),
      a.block_id,
      a.cons_start,
      a.cons_end,
      a.strand_on_consensus,
    );
    map.entry(key).or_default().insert(a.genome.clone());
  }
  map
}

/// Minimum supporter count to clear a fraction `f` of `n`: `ceil(f * n)`.
fn min_count(fraction: f64, n: usize) -> usize {
  (fraction * n as f64).ceil() as usize
}

/// The majority feature name across the supporters, when its support clears `threshold * M`.
/// Ties break to the lexicographically smallest name for determinism.
fn majority_name(per_genome: &BTreeMap<String, SupporterMeta>, m: usize, threshold: f64) -> Option<String> {
  let mut counts: BTreeMap<String, usize> = BTreeMap::new();
  for (name, _) in per_genome.values() {
    if let Some(name) = name {
      *counts.entry(name.clone()).or_default() += 1;
    }
  }
  let needed = min_count(threshold, m);
  counts
    .into_iter()
    .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.cmp(&a.0)))
    .filter(|(_, c)| *c >= needed)
    .map(|(name, _)| name)
}

/// Per-attribute majority values across the supporters that clear `threshold * M`, key-sorted.
/// A supporter contributes each of its `(key, value)` pairs once (duplicates within a supporter are
/// collapsed). Per key, the value with the most support wins; ties break to the smaller value.
fn majority_attributes(
  per_genome: &BTreeMap<String, SupporterMeta>,
  m: usize,
  threshold: f64,
) -> Vec<(String, String)> {
  let mut pair_counts: BTreeMap<(String, String), usize> = BTreeMap::new();
  for (_, attrs) in per_genome.values() {
    let distinct: BTreeSet<&(String, String)> = attrs.iter().collect();
    for kv in distinct {
      *pair_counts.entry(kv.clone()).or_default() += 1;
    }
  }

  let mut by_key: BTreeMap<String, (String, usize)> = BTreeMap::new();
  for ((k, v), c) in pair_counts {
    let best = by_key.entry(k).or_insert_with(|| (v.clone(), c));
    if c > best.1 || (c == best.1 && v < best.0) {
      *best = (v, c);
    }
  }

  let needed = min_count(threshold, m);
  by_key
    .into_iter()
    .filter(|(_, (_, c))| *c >= needed)
    .map(|(k, (v, _))| (k, v))
    .collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pangraph::edits::Edit;
  use crate::pangraph::pangraph_block::PangraphBlock;
  use crate::pangraph::pangraph_node::{NodeId, PangraphNode};
  use crate::pangraph::pangraph_path::PathId;
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use pretty_assertions::assert_eq;

  /// Build a minimal graph: `path_names[i]` becomes `PathId(i)`; each `(block_idx, paths)` entry
  /// makes `BlockId(block_idx)` traversed by those path indices (one node per listed index, so a
  /// path index listed twice gives that genome two nodes of the block — a duplication). Only the
  /// block→path traversal and path names matter to compaction, so positions/consensus are placeholders.
  fn graph_with(path_names: &[&str], blocks: &[(usize, &[usize])]) -> Pangraph {
    let mut nodes = BTreeMap::new();
    let mut blocks_map = BTreeMap::new();
    for &(bidx, path_idxs) in blocks {
      let mut alignments = BTreeMap::new();
      for (occ, &pidx) in path_idxs.iter().enumerate() {
        let nid = NodeId(bidx * 1000 + pidx * 10 + occ);
        nodes.insert(
          nid,
          PangraphNode::new(Some(nid), BlockId(bidx), PathId(pidx), Forward, (0, 0)),
        );
        alignments.insert(nid, Edit::empty());
      }
      blocks_map.insert(BlockId(bidx), PangraphBlock::new(BlockId(bidx), "A", alignments));
    }
    let paths = path_names
      .iter()
      .enumerate()
      .map(|(i, name)| {
        (
          PathId(i),
          PangraphPath::new(
            Some(PathId(i)),
            Vec::<NodeId>::new(),
            0,
            false,
            Some((*name).to_owned()),
            None,
          ),
        )
      })
      .collect::<BTreeMap<_, _>>();
    Pangraph {
      paths,
      blocks: blocks_map,
      nodes,
    }
  }

  /// A single-segment lifted annotation (both endpoints are termini).
  fn lifted(
    genome: &str,
    block: usize,
    cons: (usize, usize),
    strand: Option<Strand>,
    ftype: &str,
    name: Option<&str>,
    id: &str,
    attrs: &[(&str, &str)],
  ) -> LiftedAnnotation {
    LiftedAnnotation {
      feature_id: id.to_owned(),
      parent_feature_id: Some(id.to_owned()),
      segment_idx: 0,
      n_segments: 1,
      genome: genome.to_owned(),
      block_id: BlockId(block),
      node_id: NodeId(0),
      strand_on_consensus: strand,
      // Single segment: ordering is irrelevant, so the feature strand can mirror the consensus strand.
      feature_strand: strand,
      node_start: cons.0,
      node_end: cons.1,
      cons_start: cons.0,
      cons_end: cons.1,
      start_is_terminus: true,
      end_is_terminus: true,
      start_in_insertion: false,
      end_in_insertion: false,
      frac_covered: 1.0,
      feature_type: ftype.to_owned(),
      name: name.map(str::to_owned),
      attributes: attrs.iter().map(|(k, v)| ((*k).to_owned(), (*v).to_owned())).collect(),
    }
  }

  /// One segment of a multi-segment feature: `idx` is its node-level (genome/arc-order) index,
  /// `soc` its strand on this block's consensus, and `fstrand` the feature's genome strand (the same
  /// across all its segments). A reverse `fstrand` means the gene reads 5'→3' against genome order.
  #[allow(clippy::too_many_arguments)]
  fn seg(
    genome: &str,
    base: &str,
    idx: usize,
    n: usize,
    block: usize,
    cons: (usize, usize),
    termini: (bool, bool),
    soc: Option<Strand>,
    fstrand: Option<Strand>,
  ) -> LiftedAnnotation {
    LiftedAnnotation {
      feature_id: if n > 1 {
        format!("{base}.seg{idx}")
      } else {
        base.to_owned()
      },
      parent_feature_id: Some(base.to_owned()),
      segment_idx: idx,
      n_segments: n,
      genome: genome.to_owned(),
      block_id: BlockId(block),
      // Distinct per segment so a block crossed at two segments counts as two nodes (a duplication),
      // matching how the lift assigns a node per crossed block instance.
      node_id: NodeId(idx),
      strand_on_consensus: soc,
      feature_strand: fstrand,
      node_start: cons.0,
      node_end: cons.1,
      cons_start: cons.0,
      cons_end: cons.1,
      start_is_terminus: termini.0,
      end_is_terminus: termini.1,
      start_in_insertion: false,
      end_in_insertion: false,
      frac_covered: 0.5,
      feature_type: "CDS".to_owned(),
      name: Some("geneA".to_owned()),
      attributes: vec![],
    }
  }

  fn strat(min_frequency: f64, property_threshold: f64) -> CoordinateConsensusStrategy {
    CoordinateConsensusStrategy {
      min_frequency,
      property_threshold,
    }
  }

  #[test]
  fn test_all_agree_single_block() {
    let graph = graph_with(&["g0", "g1", "g2"], &[(1, &[0, 1, 2])]);
    let anns = vec![
      lifted("g0", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f0", &[]),
      lifted("g1", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f1", &[]),
      lifted("g2", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f2", &[]),
    ];
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    assert_eq!(out.len(), 1);
    let b = &out[0];
    assert_eq!((b.cluster_id, b.segment_idx, b.n_segments), (0, 0, 1));
    assert_eq!((b.block_id, b.cons_start, b.cons_end), (BlockId(1), 10, 200));
    assert_eq!((b.n_support, b.n_total), (3, 3));
    assert_eq!((b.n_support_segment, b.n_total_segment), (3, 3));
    assert_eq!(b.consensus_name.as_deref(), Some("geneA"));
    assert_eq!(b.strand_on_consensus, Some(Forward));
  }

  #[test]
  fn test_below_threshold_dropped_but_kept_when_lenient() {
    let graph = graph_with(&["g0", "g1", "g2", "g3"], &[(1, &[0, 1, 2, 3])]);
    let anns = vec![
      lifted("g0", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f0", &[]),
      lifted("g1", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f1", &[]),
    ];
    // 2 of 4 -> ceil(0.9*4)=4 required -> dropped.
    assert!(strat(0.9, 0.5).compact(&anns, &graph).unwrap().is_empty());
    // ceil(0.5*4)=2 required -> kept.
    let out = strat(0.5, 0.5).compact(&anns, &graph).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!((out[0].n_support, out[0].n_total), (2, 4));
  }

  #[test]
  fn test_feature_type_splits_clusters() {
    let graph = graph_with(&["g0", "g1"], &[(1, &[0, 1])]);
    let anns = vec![
      lifted("g0", 1, (10, 200), Some(Forward), "gene", Some("geneA"), "gene0", &[]),
      lifted("g0", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "cds0", &[]),
      lifted("g1", 1, (10, 200), Some(Forward), "gene", Some("geneA"), "gene1", &[]),
      lifted("g1", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "cds1", &[]),
    ];
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    assert_eq!(out.len(), 2);
    let types: Vec<&str> = out.iter().map(|b| b.feature_type.as_str()).collect();
    assert_eq!(types, vec!["CDS", "gene"]); // same coords -> sorted by type
    assert!(out.iter().all(|b| (b.n_support, b.n_total) == (2, 2)));
  }

  #[test]
  fn test_strand_splits_clusters() {
    let graph = graph_with(&["g0", "g1"], &[(1, &[0, 1])]);
    let anns = vec![
      lifted("g0", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f0", &[]),
      lifted("g1", 1, (10, 200), Some(Reverse), "CDS", Some("geneA"), "f1", &[]),
    ];
    let out = strat(0.0, 0.5).compact(&anns, &graph).unwrap();
    assert_eq!(out.len(), 2);
    assert_eq!(out[0].strand_on_consensus, Some(Forward)); // Forward < Reverse
    assert_eq!(out[1].strand_on_consensus, Some(Reverse));
    assert_eq!((out[0].n_support, out[0].n_total), (1, 2));
  }

  #[test]
  fn test_name_and_attribute_majority_threshold() {
    let graph = graph_with(&["g0", "g1", "g2"], &[(1, &[0, 1, 2])]);
    let anns = vec![
      lifted(
        "g0",
        1,
        (10, 200),
        Some(Forward),
        "CDS",
        Some("geneA"),
        "f0",
        &[("product", "widget")],
      ),
      lifted(
        "g1",
        1,
        (10, 200),
        Some(Forward),
        "CDS",
        Some("geneA"),
        "f1",
        &[("product", "widget")],
      ),
      lifted(
        "g2",
        1,
        (10, 200),
        Some(Forward),
        "CDS",
        Some("geneB"),
        "f2",
        &[("product", "gadget")],
      ),
    ];
    // 2/3 majority clears ceil(0.5*3)=2.
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].consensus_name.as_deref(), Some("geneA"));
    assert_eq!(
      out[0].consensus_attributes,
      vec![("product".to_owned(), "widget".to_owned())]
    );
    // 2/3 does not clear ceil(0.7*3)=3.
    let out = strat(1.0, 0.7).compact(&anns, &graph).unwrap();
    assert_eq!(out[0].consensus_name, None);
    assert!(out[0].consensus_attributes.is_empty());
  }

  #[test]
  fn test_multi_block_feature_emits_one_row_per_segment() {
    let graph = graph_with(&["g0", "g1"], &[(1, &[0, 1]), (2, &[0, 1])]);
    let mut anns = Vec::new();
    for (g, base) in [("g0", "f0"), ("g1", "f1")] {
      anns.push(seg(
        g,
        base,
        0,
        2,
        1,
        (150, 200),
        (true, false),
        Some(Forward),
        Some(Forward),
      ));
      anns.push(seg(
        g,
        base,
        1,
        2,
        2,
        (0, 50),
        (false, true),
        Some(Forward),
        Some(Forward),
      ));
    }
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    // One cluster, two segments ordered 5'→3'.
    assert_eq!(out.len(), 2);
    assert!(out.iter().all(|b| b.cluster_id == 0 && b.n_segments == 2));
    assert_eq!(
      (out[0].segment_idx, out[0].block_id, out[0].cons_start, out[0].cons_end),
      (0, BlockId(1), 150, 200)
    );
    assert_eq!(
      (out[1].segment_idx, out[1].block_id, out[1].cons_start, out[1].cons_end),
      (1, BlockId(2), 0, 50)
    );
    assert!(out.iter().all(|b| (b.n_support, b.n_total) == (2, 2)));
  }

  #[test]
  fn test_n_total_is_paths_traversing_block_not_all_paths() {
    // Block 1 is accessory: present only in g0, g1 (not g2, g3).
    let graph = graph_with(&["g0", "g1", "g2", "g3"], &[(1, &[0, 1])]);
    let anns = vec![
      lifted("g0", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f0", &[]),
      lifted("g1", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f1", &[]),
    ];
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!((out[0].n_support, out[0].n_total), (2, 2)); // N is 2, not 4
  }

  #[test]
  fn test_partial_feature_without_two_termini_is_excluded() {
    let graph = graph_with(&["g0"], &[(1, &[0])]);
    let mut a = lifted("g0", 1, (10, 200), Some(Forward), "CDS", Some("geneA"), "f0", &[]);
    a.end_is_terminus = false; // only one real terminus remains
    assert!(strat(0.0, 0.5).compact(&[a], &graph).unwrap().is_empty());
  }

  /// Regression for the yjeM bug: the same gene in two genomes of opposite global orientation —
  /// body on block 1 (`+` on its consensus), inverted tail on block 3 (`-` on its consensus). The
  /// forward genome traverses them in arc order [B1, B3]; the reverse genome in [B3, B1] with the
  /// per-block strands unchanged. Ordering 5'→3' must collapse both into one cluster (not split by a
  /// strand sampled from the genome-lowest segment), giving M = N = 2.
  #[test]
  fn test_inversion_crossing_collapses_to_one_cluster() {
    let graph = graph_with(&["gplus", "gminus"], &[(1, &[0, 1]), (3, &[0, 1])]);
    let anns = vec![
      // Forward genome: gene reads 5'→3' with genome order, B1 then B3.
      seg(
        "gplus",
        "fp",
        0,
        2,
        1,
        (9228, 10694),
        (true, false),
        Some(Forward),
        Some(Forward),
      ),
      seg(
        "gplus",
        "fp",
        1,
        2,
        3,
        (705, 742),
        (false, true),
        Some(Reverse),
        Some(Forward),
      ),
      // Reverse genome: arc order is B3 then B1; feature_strand Reverse flips it back to 5'→3'.
      seg(
        "gminus",
        "fm",
        0,
        2,
        3,
        (705, 742),
        (true, false),
        Some(Reverse),
        Some(Reverse),
      ),
      seg(
        "gminus",
        "fm",
        1,
        2,
        1,
        (9228, 10694),
        (false, true),
        Some(Forward),
        Some(Reverse),
      ),
    ];
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    // A single cluster with two segments, body (B1,+) before tail (B3,-).
    assert_eq!(out.len(), 2);
    assert!(out.iter().all(|b| b.cluster_id == 0 && b.n_segments == 2));
    assert_eq!(
      (
        out[0].segment_idx,
        out[0].block_id,
        out[0].cons_start,
        out[0].cons_end,
        out[0].strand_on_consensus
      ),
      (0, BlockId(1), 9228, 10694, Some(Forward))
    );
    assert_eq!(
      (
        out[1].segment_idx,
        out[1].block_id,
        out[1].cons_start,
        out[1].cons_end,
        out[1].strand_on_consensus
      ),
      (1, BlockId(3), 705, 742, Some(Reverse))
    );
    assert!(out.iter().all(|b| (b.n_support, b.n_total) == (2, 2)));
    // Both genomes place each segment at the same coords, so per-segment support is also 2 of 2.
    assert!(out.iter().all(|b| (b.n_support_segment, b.n_total_segment) == (2, 2)));
  }

  /// A feature crossing a duplicated block twice (B5, B6, B5) keeps both occurrences as separate
  /// segment rows, and `N` requires genomes to carry block 5 at least twice: g2, which has it only
  /// once, is excluded from the denominator.
  #[test]
  fn test_duplicated_block_crossed_twice() {
    let graph = graph_with(&["g0", "g1", "g2"], &[(5, &[0, 0, 1, 1, 2]), (6, &[0, 1, 2])]);
    let mut anns = Vec::new();
    for (g, base) in [("g0", "f0"), ("g1", "f1")] {
      anns.push(seg(
        g,
        base,
        0,
        3,
        5,
        (0, 50),
        (true, false),
        Some(Forward),
        Some(Forward),
      ));
      anns.push(seg(
        g,
        base,
        1,
        3,
        6,
        (0, 30),
        (false, false),
        Some(Forward),
        Some(Forward),
      ));
      anns.push(seg(
        g,
        base,
        2,
        3,
        5,
        (70, 100),
        (false, true),
        Some(Forward),
        Some(Forward),
      ));
    }
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    // Three segment rows; block 5 appears twice, at its two distinct consensus coordinates.
    assert_eq!(out.len(), 3);
    assert!(out.iter().all(|b| b.cluster_id == 0 && b.n_segments == 3));
    let blocks: Vec<_> = out
      .iter()
      .map(|b| (b.segment_idx, b.block_id, b.cons_start, b.cons_end))
      .collect();
    assert_eq!(
      blocks,
      vec![(0, BlockId(5), 0, 50), (1, BlockId(6), 0, 30), (2, BlockId(5), 70, 100),]
    );
    // N counts only genomes carrying block 5 twice (g0, g1) — not g2, which has a single copy.
    assert!(out.iter().all(|b| (b.n_support, b.n_total) == (2, 2)));
  }

  /// Two segments of one feature landing on the **same node** (an origin-spanning node split into
  /// two coverage pieces) need only one node of the block — so `N` must still count the genome and
  /// `N >= M` must hold. Regression for the `M > N` bug on whole-genome `region` features, where
  /// counting required multiplicity by segment (2) instead of distinct node (1) gave `N = 0`.
  #[test]
  fn test_origin_split_same_node_not_double_counted() {
    let graph = graph_with(&["g0"], &[(7, &[0])]); // block 7 present once on g0
    let mut s0 = seg(
      "g0",
      "f",
      0,
      2,
      7,
      (100, 150),
      (true, false),
      Some(Forward),
      Some(Forward),
    );
    let mut s1 = seg("g0", "f", 1, 2, 7, (0, 40), (false, true), Some(Forward), Some(Forward));
    s0.node_id = NodeId(0);
    s1.node_id = NodeId(0); // same node as s0: origin-split, not a duplication
    let out = strat(1.0, 0.5).compact(&[s0, s1], &graph).unwrap();
    assert_eq!(out.len(), 2);
    assert!(out.iter().all(|b| b.block_id == BlockId(7) && b.n_segments == 2));
    assert!(out.iter().all(|b| (b.n_support, b.n_total) == (1, 1)));
  }

  /// Mini-yjeM: a body shared on block 1 with three different tails (blocks 2/3/4) across genomes.
  /// Each tail variant is its own cluster, every one at M/N = 1.0, and the shared body reports the
  /// same per-segment support (6 of 6) in all three clusters.
  #[test]
  fn test_shared_body_three_tails_each_confident() {
    let graph = graph_with(
      &["g0", "g1", "g2", "g3", "g4", "g5"],
      &[(1, &[0, 1, 2, 3, 4, 5]), (2, &[0, 1, 2]), (3, &[3, 4]), (4, &[5])],
    );
    let mut anns = Vec::new();
    for (g, base, tail) in [
      ("g0", "a0", 2),
      ("g1", "a1", 2),
      ("g2", "a2", 2),
      ("g3", "a3", 3),
      ("g4", "a4", 3),
      ("g5", "a5", 4),
    ] {
      anns.push(seg(
        g,
        base,
        0,
        2,
        1,
        (9228, 10694),
        (true, false),
        Some(Forward),
        Some(Forward),
      ));
      anns.push(seg(
        g,
        base,
        1,
        2,
        tail,
        (0, 37),
        (false, true),
        Some(Forward),
        Some(Forward),
      ));
    }
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    // Three clusters (B1+B2, B1+B3, B1+B4), two rows each.
    assert_eq!(out.len(), 6);
    let cluster_ids: BTreeSet<_> = out.iter().map(|b| b.cluster_id).collect();
    assert_eq!(cluster_ids, BTreeSet::from([0, 1, 2]));
    // Every cluster is fully supported among the genomes capable of it.
    assert!(out.iter().all(|b| b.n_support == b.n_total));
    // The shared body segment reports the same reproducibility (6 of 6) in every cluster.
    assert!(
      out
        .iter()
        .filter(|b| b.block_id == BlockId(1))
        .all(|b| (b.n_support_segment, b.n_total_segment) == (6, 6))
    );
  }
}
