use crate::annotation::lift::LiftedAnnotation;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::BlockId;
use crate::pangraph::pangraph_path::{PangraphPath, PathId};
use crate::pangraph::strand::Strand;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// A consensus endpoint: a `(block, block-consensus coordinate)` pair.
type Endpoint = (BlockId, usize);

/// Identity of a block-level cluster. Two per-genome feature instances are "the same" consensus
/// feature exactly when these match: the two block-consensus terminus endpoints (in canonical,
/// sorted order) together with the feature type and the consensus strand.
type ClusterKey = (String, Option<Strand>, Endpoint, Endpoint);

/// Per-genome metadata kept for one supporter of a cluster: its feature `name` and `attributes`.
type SupporterMeta = (Option<String>, Vec<(String, String)>);

/// A block-level consensus annotation: a placement that recurs, at the *same* block-consensus
/// coordinates, across enough of the genomes that carry it. This is the compacted, opinionated
/// view built on top of the lossless node-level table (it never replaces it).
///
/// The two endpoints are stored in **canonical (sorted) order**, so `(start_block_id, cons_start)`
/// is the lower endpoint and `(end_block_id, cons_end)` the higher — *not* necessarily the genome
/// 5'/3' ends. `strand_on_consensus` disambiguates orientation. For the common single-block feature
/// `start_block_id == end_block_id`.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct BlockAnnotation {
  /// Feature type shared by the cluster, e.g. `"CDS"` or `"gene"`.
  pub feature_type: String,

  /// Strand relative to the block consensus shared by the cluster (`None` if unstranded).
  pub strand_on_consensus: Option<Strand>,

  /// Lower (canonical) terminus endpoint: its block and block-consensus coordinate.
  pub start_block_id: BlockId,
  pub cons_start: usize,

  /// Higher (canonical) terminus endpoint: its block and block-consensus coordinate.
  pub end_block_id: BlockId,
  pub cons_end: usize,

  /// Majority feature name across the supporters, when it clears the property threshold; else `None`.
  pub consensus_name: Option<String>,

  /// Per-attribute majority values that clear the property threshold (key-sorted).
  pub consensus_attributes: Vec<(String, String)>,

  /// `M` — number of genomes sharing this exact placement.
  pub n_support: usize,

  /// `N` — number of genomes traversing the cluster's block(s) (the supporters' denominator).
  pub n_total: usize,
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
/// Each genome's feature instance is reduced to its two block-consensus terminus endpoints; genomes
/// whose `(feature_type, strand, endpoints)` match exactly form a cluster. A cluster is emitted when
/// its support `M` reaches `ceil(min_frequency * N)` (at least 1), where `N` is the number of
/// genomes traversing the cluster's block(s). Consensus `name`/attributes are the per-value
/// majorities among the supporters that clear `property_threshold`.
pub struct CoordinateConsensusStrategy {
  /// Minimum supporter fraction `M / N` (of genomes traversing the block) for a cluster to be emitted.
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
    // Aggregate per cluster key, deduplicating supporters by genome (a genome counts once).
    let mut clusters: BTreeMap<ClusterKey, BTreeMap<String, SupporterMeta>> = BTreeMap::new();
    for inst in feature_instances(node_annotations) {
      clusters
        .entry(inst.key)
        .or_default()
        .entry(inst.genome)
        .or_insert((inst.name, inst.attributes));
    }

    let block_paths = block_path_sets(graph);

    let mut out = Vec::new();
    for (key, per_genome) in clusters {
      let (feature_type, strand, lo, hi) = key;
      let m = per_genome.len();
      let n = cluster_n_total(graph, &block_paths, lo.0, hi.0);
      let min_support = min_count(self.min_frequency, n).max(1);
      if m < min_support {
        continue; // below threshold: stays only in the node-level table
      }
      out.push(BlockAnnotation {
        feature_type,
        strand_on_consensus: strand,
        start_block_id: lo.0,
        cons_start: lo.1,
        end_block_id: hi.0,
        cons_end: hi.1,
        consensus_name: majority_name(&per_genome, m, self.property_threshold),
        consensus_attributes: majority_attributes(&per_genome, m, self.property_threshold),
        n_support: m,
        n_total: n,
      });
    }

    // Deterministic order: by placement, then type, then strand.
    out.sort_by(|a, b| {
      (a.start_block_id, a.cons_start, a.end_block_id, a.cons_end)
        .cmp(&(b.start_block_id, b.cons_start, b.end_block_id, b.cons_end))
        .then_with(|| a.feature_type.cmp(&b.feature_type))
        .then_with(|| a.strand_on_consensus.cmp(&b.strand_on_consensus))
    });
    Ok(out)
  }
}

/// One genome's instance of a feature, reduced to its cluster key plus the metadata compaction
/// needs (representative `name`/`attributes`, taken from the 5'-most segment — they are identical
/// across a feature's segments).
struct FeatureInstance {
  genome: String,
  key: ClusterKey,
  name: Option<String>,
  attributes: Vec<(String, String)>,
}

/// Group node-level rows into per-genome feature instances and reduce each to its cluster key.
///
/// Rows are grouped by `(genome, feature_type, base feature id)`; instances that do not resolve to
/// exactly two terminus endpoints (e.g. partial features) are dropped from compaction and remain in
/// the node-level table.
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

/// Reduce one feature instance's node-level rows to its [`FeatureInstance`], or `None` when it does
/// not have exactly two terminus endpoints.
fn reduce_instance(genome: String, feature_type: String, rows: &[&LiftedAnnotation]) -> Option<FeatureInstance> {
  let mut endpoints: Vec<Endpoint> = Vec::with_capacity(2);
  for a in rows {
    if a.start_is_terminus {
      endpoints.push((a.block_id, a.cons_start));
    }
    if a.end_is_terminus {
      endpoints.push((a.block_id, a.cons_end));
    }
  }
  if endpoints.len() != 2 {
    return None;
  }
  endpoints.sort_unstable();
  let (lo, hi) = (endpoints[0], endpoints[1]);

  // The 5'-most segment (segment_idx 0) provides the representative strand/name/attributes.
  let rep = rows.iter().min_by_key(|a| a.segment_idx).copied()?;
  Some(FeatureInstance {
    genome,
    key: (feature_type, rep.strand_on_consensus, lo, hi),
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

/// Map every block to the set of paths that traverse it.
fn block_path_sets(graph: &Pangraph) -> BTreeMap<BlockId, BTreeSet<PathId>> {
  graph
    .blocks
    .iter()
    .map(|(&bid, block)| (bid, block.isolates(graph).collect::<BTreeSet<_>>()))
    .collect()
}

/// `N` for a cluster: the number of distinct genomes traversing **all** of the cluster's block(s)
/// (the intersection of the per-block path sets; just the one set for a single-block cluster).
fn cluster_n_total(
  graph: &Pangraph,
  block_paths: &BTreeMap<BlockId, BTreeSet<PathId>>,
  start_block: BlockId,
  end_block: BlockId,
) -> usize {
  let empty = BTreeSet::new();
  let start = block_paths.get(&start_block).unwrap_or(&empty);
  let paths: BTreeSet<PathId> = if start_block == end_block {
    start.clone()
  } else {
    let end = block_paths.get(&end_block).unwrap_or(&empty);
    start.intersection(end).copied().collect()
  };
  paths
    .iter()
    .filter_map(|pid| graph.paths.get(pid))
    .map(genome_label)
    .collect::<BTreeSet<_>>()
    .len()
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
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use pretty_assertions::assert_eq;

  /// Build a minimal graph: `path_names[i]` becomes `PathId(i)`; each `(block_idx, paths)` entry
  /// makes `BlockId(block_idx)` traversed by those paths (one node each). Only the block→path
  /// traversal and path names matter to compaction, so positions/consensus are placeholders.
  fn graph_with(path_names: &[&str], blocks: &[(usize, &[usize])]) -> Pangraph {
    let mut nodes = BTreeMap::new();
    let mut blocks_map = BTreeMap::new();
    for &(bidx, path_idxs) in blocks {
      let mut alignments = BTreeMap::new();
      for &pidx in path_idxs {
        let nid = NodeId(bidx * 100 + pidx);
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

  /// One segment of a multi-segment forward feature, with explicit terminus flags.
  fn seg(
    genome: &str,
    base: &str,
    idx: usize,
    n: usize,
    block: usize,
    cons: (usize, usize),
    termini: (bool, bool),
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
      node_id: NodeId(0),
      strand_on_consensus: Some(Forward),
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
    assert_eq!(
      (b.start_block_id, b.cons_start, b.end_block_id, b.cons_end),
      (BlockId(1), 10, BlockId(1), 200)
    );
    assert_eq!((b.n_support, b.n_total), (3, 3));
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
  fn test_multi_block_feature_uses_both_endpoints() {
    let graph = graph_with(&["g0", "g1"], &[(1, &[0, 1]), (2, &[0, 1])]);
    let mut anns = Vec::new();
    for (g, base) in [("g0", "f0"), ("g1", "f1")] {
      anns.push(seg(g, base, 0, 2, 1, (150, 200), (true, false)));
      anns.push(seg(g, base, 1, 2, 2, (0, 50), (false, true)));
    }
    let out = strat(1.0, 0.5).compact(&anns, &graph).unwrap();
    assert_eq!(out.len(), 1);
    let b = &out[0];
    assert_eq!((b.start_block_id, b.cons_start), (BlockId(1), 150));
    assert_eq!((b.end_block_id, b.cons_end), (BlockId(2), 50));
    assert_eq!((b.n_support, b.n_total), (2, 2));
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
}
