mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use pangraph::annotation::compact::{BlockCompactionStrategy, CoordinateConsensusStrategy};
  use pangraph::annotation::feature::Feature;
  use pangraph::annotation::lift::lift_features;
  use pangraph::annotation::matching::match_features_to_paths;
  use pangraph::annotation::writer::{AnnotationWriter, CsvAnnotationWriter};
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::pangraph::strand::Strand;
  use pangraph::utils::interval::Interval;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::fs::read_to_string;
  use tempfile::tempdir;

  const GRAPH: &str = "../../data/test_graph.json";

  /// A `CDS` feature covering the full genome extent `[p0, p1)` of a node.
  fn whole_node_feature(seqid: &str, p0: usize, p1: usize) -> Feature {
    let id = format!("{seqid}_core");
    Feature {
      seqid: seqid.to_owned(),
      source: None,
      feature_type: "CDS".to_owned(),
      interval: Interval::new(p0, p1),
      strand: Some(Strand::Forward),
      id: Some(id.clone()),
      name: Some("core_gene".to_owned()),
      attributes: vec![("ID".to_owned(), id), ("product".to_owned(), "core_product".to_owned())],
    }
  }

  /// End-to-end: placing a feature over the full extent of a core block's node on every genome lifts
  /// each to the same consensus interval `[0, L)`, so compaction must collapse them into consensus
  /// feature(s) for that block whose supporters sum to the number of genomes carrying it.
  #[test]
  fn itest_compact_core_block_collapses_across_genomes() -> Result<(), Report> {
    let graph = Pangraph::from_path(&Some(GRAPH))?;
    let n_paths = graph.paths.len();

    // A core block whose nodes are all non-wrapping, so a whole-node feature is a clean interval.
    let bid = graph
      .core_block_ids()
      .find(|bid| {
        graph.blocks[bid].alignment_keys().iter().all(|nid| {
          let (p0, p1) = graph.nodes[nid].position();
          p0 < p1
        })
      })
      .expect("a core block with only non-wrapping nodes");
    let block_len = graph.blocks[&bid].consensus_len();

    // One whole-node feature per genome over that block.
    let mut features = Vec::new();
    for path in graph.paths.values() {
      let name = path.name().as_deref().unwrap();
      let nid = path
        .nodes()
        .iter()
        .copied()
        .find(|nid| graph.nodes[nid].block_id() == bid)
        .expect("core block present on every path");
      let (p0, p1) = graph.nodes[&nid].position();
      features.push(whole_node_feature(name, p0, p1));
    }

    let grouped = match_features_to_paths(features, &graph, &BTreeMap::new())?;
    let lifted = lift_features(&grouped, &graph)?;

    let strategy = CoordinateConsensusStrategy {
      min_frequency: 0.0,
      property_threshold: 0.5,
    };
    let blocks = strategy.compact(&lifted, &graph)?;

    // Every emitted annotation is coherent.
    for b in &blocks {
      assert!(b.n_support >= 1 && b.n_support <= b.n_total, "support within total");
    }

    // The core-block placement(s): one per traversal strand, together supported by all genomes.
    let for_block: Vec<_> = blocks
      .iter()
      .filter(|b| b.start_block_id == bid && b.end_block_id == bid)
      .collect();
    assert!(!for_block.is_empty(), "core block produced a consensus annotation");
    let total_support: usize = for_block.iter().map(|b| b.n_support).sum();
    assert_eq!(total_support, n_paths, "every genome supports the core-block placement");
    for b in &for_block {
      assert_eq!(
        (b.cons_start, b.cons_end),
        (0, block_len),
        "whole node maps to whole consensus"
      );
      assert_eq!(b.n_total, n_paths, "N is all genomes traversing the core block");
      assert_eq!(b.feature_type, "CDS");
      assert_eq!(b.consensus_name.as_deref(), Some("core_gene"));
      assert!(
        b.consensus_attributes
          .contains(&("product".to_owned(), "core_product".to_owned())),
        "shared product is promoted to consensus"
      );
    }

    // The block-level writer round-trips the real output.
    let dir = tempdir()?;
    let out = dir.path().join("block_annotations.csv");
    {
      let mut writer = CsvAnnotationWriter::new(&out, b',')?;
      writer.write_block_annotations(&blocks)?;
    }
    let contents = read_to_string(&out)?;
    assert!(
      contents.starts_with("type,strand_on_consensus,start_block_id"),
      "block CSV header present"
    );
    Ok(())
  }
}
