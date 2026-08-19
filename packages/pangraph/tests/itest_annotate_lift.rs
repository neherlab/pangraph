mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use pangraph::annotation::feature::Feature;
  use pangraph::annotation::lift::{LiftedAnnotation, lift_feature, lift_features};
  use pangraph::annotation::matching::match_features_to_paths;
  use pangraph::annotation::writer::{AnnotationWriter, CsvAnnotationWriter};
  use pangraph::commands::reconstruct::reconstruct_run::reconstruct;
  use pangraph::io::csv::parse_csv;
  use pangraph::io::seq::reverse_complement;
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::pangraph::pangraph_block::BlockId;
  use pangraph::pangraph::pangraph_path::PangraphPath;
  use pangraph::pangraph::strand::Strand;
  use pangraph::representation::seq::Seq;
  use pangraph::utils::interval::Interval;
  use pretty_assertions::assert_eq;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use std::fs::read_to_string;
  use tempfile::tempdir;

  const GRAPH: &str = "../../data/test_graph.json";

  /// Subset of the CSV columns asserted in the writer test; other columns are ignored on
  /// deserialization (some fields only drive deserialization).
  #[allow(dead_code)]
  #[derive(Deserialize)]
  struct Row {
    genome: String,
    block_id: usize,
    cons_start: usize,
    cons_end: usize,
  }

  /// Load the test graph together with each path's genome reconstructed from the graph (the
  /// independent oracle for round-trip checks). Genomes are keyed by path name and laid out in
  /// absolute genome coordinates `[0, tot_len)`.
  fn load() -> Result<(Pangraph, BTreeMap<String, Seq>), Report> {
    let graph = Pangraph::from_path(&Some(GRAPH))?;
    let mut genomes = BTreeMap::new();
    for rec in reconstruct(&graph) {
      let rec = rec?;
      genomes.insert(rec.seq_name.clone(), rec.seq);
    }
    Ok((graph, genomes))
  }

  fn feature_at(seqid: &str, start: usize, end: usize, strand: Strand) -> Feature {
    let id = format!("{seqid}_{start}_{end}");
    Feature {
      seqid: seqid.to_owned(),
      source: None,
      feature_type: "CDS".to_owned(),
      interval: Interval::new(start, end),
      strand: Some(strand),
      id: Some(id.clone()),
      name: Some("test_gene".to_owned()),
      attributes: vec![("ID".to_owned(), id)],
    }
  }

  /// Reassemble the genome bases of a lifted feature from its segments: each node's
  /// consensus-oriented sequence sliced at the node-local coordinates, reverse-complemented on
  /// reverse-strand nodes, concatenated in segment order.
  fn assemble_from_segments(graph: &Pangraph, lifted: &[LiftedAnnotation]) -> Seq {
    let mut out = Seq::new();
    for seg in lifted {
      let block = &graph.blocks[&seg.block_id];
      let node = &graph.nodes[&seg.node_id];
      let node_seq = block.alignment(seg.node_id).apply(block.consensus()).unwrap();
      let piece = Seq::from(&node_seq[seg.node_start..seg.node_end]);
      let piece = if node.strand().is_reverse() {
        reverse_complement(&piece).unwrap()
      } else {
        piece
      };
      out.extend_seq(&piece);
    }
    out
  }

  /// Lift a feature and check structural invariants plus the strong round-trip oracle: the
  /// reassembled segment bases must equal the genome substring `[f_s, f_e)`.
  fn check_round_trip(
    graph: &Pangraph,
    genome: &Seq,
    path: &PangraphPath,
    f_s: usize,
    f_e: usize,
    strand: Strand,
  ) -> Result<(), Report> {
    let name = path.name().as_deref().unwrap();
    let lifted = lift_feature(&feature_at(name, f_s, f_e, strand), path, graph)?;
    assert!(
      !lifted.is_empty(),
      "feature [{f_s},{f_e}) on {name} produced no segments"
    );

    for (i, seg) in lifted.iter().enumerate() {
      assert_eq!(seg.segment_idx, i, "segment_idx order");
      assert_eq!(seg.n_segments, lifted.len(), "n_segments consistency");
      assert!(seg.node_start <= seg.node_end, "node coords ordered");
      assert!(seg.cons_start <= seg.cons_end, "cons coords ordered");
    }

    // A fully-covered feature has exactly two real termini (its genome start and end).
    let termini: usize = lifted
      .iter()
      .map(|s| usize::from(s.start_is_terminus) + usize::from(s.end_is_terminus))
      .sum();
    assert_eq!(termini, 2, "feature [{f_s},{f_e}) termini count");

    let frac: f64 = lifted.iter().map(|s| s.frac_covered).sum();
    assert!(
      (frac - 1.0).abs() < 1e-9,
      "feature [{f_s},{f_e}) frac_covered sum {frac}"
    );

    let assembled = assemble_from_segments(graph, &lifted);
    let expected = Seq::from(&genome[f_s..f_e]);
    assert_eq!(
      assembled, expected,
      "feature [{f_s},{f_e}) strand {strand:?} round-trip"
    );
    Ok(())
  }

  #[test]
  fn itest_lift_round_trip_across_genomes() -> Result<(), Report> {
    let (graph, genomes) = load()?;
    for path in graph.paths.values() {
      let name = path.name().as_deref().unwrap();
      let genome = &genomes[name];
      let n = path.tot_len();
      // A spread of intervals: short ones, ones spanning many blocks, and (for the longer paths)
      // ones overlapping the origin-wrapping node, in both strands.
      let cases: [(usize, usize); 6] = [
        (1, 50),
        (1000, 1200),
        (5000, 5100),
        (100, n - 100),
        (n * 3 / 4, n),
        (n - 1500, n - 10),
      ];
      for (f_s, f_e) in cases {
        assert!(f_s < f_e && f_e <= n);
        for strand in [Strand::Forward, Strand::Reverse] {
          check_round_trip(&graph, genome, path, f_s, f_e, strand)?;
        }
      }
    }
    Ok(())
  }

  #[test]
  fn itest_lift_single_node_interior_strand_mapping() -> Result<(), Report> {
    let (graph, genomes) = load()?;
    for path in graph.paths.values() {
      let name = path.name().as_deref().unwrap();
      let genome = &genomes[name];
      // First non-wrapping node with room for an interior feature.
      let node = path
        .nodes()
        .iter()
        .map(|nid| &graph.nodes[nid])
        .find(|node| {
          let (s, e) = node.position();
          s < e && e - s >= 20
        })
        .expect("a non-wrapping node");
      let (s, e) = node.position();
      let (f_s, f_e) = (s + 5, e - 5);

      let lifted = lift_feature(&feature_at(name, f_s, f_e, Strand::Forward), path, &graph)?;
      assert_eq!(lifted.len(), 1, "interior feature is a single segment");
      let seg = &lifted[0];
      assert_eq!(seg.node_id, node.id());
      let expected_strand = if node.strand().is_reverse() {
        Strand::Reverse
      } else {
        Strand::Forward
      };
      assert_eq!(
        seg.strand_on_consensus,
        Some(expected_strand),
        "strand maps through node strand"
      );
      assert!(
        seg.start_is_terminus && seg.end_is_terminus,
        "single segment: both termini"
      );
      check_round_trip(&graph, genome, path, f_s, f_e, Strand::Forward)?;
    }
    Ok(())
  }

  #[test]
  fn itest_lift_cross_block_boundary_segments() -> Result<(), Report> {
    let (graph, genomes) = load()?;
    let path = graph.paths.values().next().expect("at least one path");
    let name = path.name().as_deref().unwrap();
    let genome = &genomes[name];

    // Find an internal boundary shared by two non-wrapping neighbours, with room on each side.
    let nodes: Vec<_> = path.nodes().iter().map(|nid| &graph.nodes[nid]).collect();
    let mut found = false;
    for pair in nodes.windows(2) {
      let [n0, n1] = pair else { continue };
      let (s0, e0) = n0.position();
      let (s1, e1) = n1.position();
      if s0 < e0 && s1 < e1 && e0 == s1 && e0 - s0 >= 80 && e1 - s1 >= 80 {
        let boundary = e0;
        let (f_s, f_e) = (boundary - 60, boundary + 60);

        let lifted = lift_feature(&feature_at(name, f_s, f_e, Strand::Forward), path, &graph)?;
        assert_eq!(
          lifted.len(),
          2,
          "feature crossing one boundary splits into two segments"
        );
        assert_eq!(lifted[0].node_id, n0.id());
        assert_eq!(lifted[1].node_id, n1.id());
        assert_eq!(lifted[0].parent_feature_id, lifted[1].parent_feature_id);
        assert_ne!(lifted[0].feature_id, lifted[1].feature_id);
        check_round_trip(&graph, genome, path, f_s, f_e, Strand::Forward)?;
        found = true;
        break;
      }
    }
    assert!(
      found,
      "expected an internal block boundary between two non-wrapping nodes"
    );
    Ok(())
  }

  #[test]
  fn itest_lift_to_csv_writer() -> Result<(), Report> {
    let (graph, _genomes) = load()?;
    let path = graph.paths.values().next().expect("at least one path");
    let name = path.name().as_deref().unwrap().to_owned();

    let features = vec![
      feature_at(&name, 1000, 1200, Strand::Forward),
      feature_at(&name, 100, 5000, Strand::Reverse),
    ];
    let grouped = match_features_to_paths(features, &graph, &BTreeMap::new())?;
    let lifted = lift_features(&grouped, &graph)?;
    assert!(!lifted.is_empty());

    let dir = tempdir()?;
    let out = dir.path().join("annotations.csv");
    {
      let mut writer = CsvAnnotationWriter::new(&out, b',')?;
      writer.write_node_annotations(&lifted)?;
    }

    let contents = read_to_string(&out)?;
    assert!(
      contents.starts_with("feature_id,parent_feature_id,segment_idx"),
      "header row present"
    );

    let rows: Vec<Row> = parse_csv(&contents)?;
    assert_eq!(rows.len(), lifted.len(), "one CSV row per lifted segment");
    for r in &rows {
      assert_eq!(r.genome, name, "genome column is the path name");
      let block = graph.blocks.get(&BlockId(r.block_id)).expect("block exists");
      assert!(r.cons_start <= r.cons_end, "consensus coords ordered");
      assert!(r.cons_end <= block.consensus_len(), "consensus end within block");
    }
    Ok(())
  }
}
