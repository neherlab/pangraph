mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use pangraph::annotation::lift::lift_features;
  use pangraph::annotation::matching::match_features_to_paths;
  use pangraph::io::gff::GffReader;
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::pangraph::strand::Strand;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};

  // Real Klebsiella annotations (downloaded from NCBI) for two of the genomes in
  // `data/klebs_graph.json.gz`. The GFF seqid column has been normalized to the bare accession (the
  // FASTA record id / pangraph path name) by stripping the `.N` version, so annotations match the
  // graph paths by exact string equality — no seqid map needed.
  //                     bare accession      gff path
  const KLEBS: &[(&str, &str)] = &[
    ("NZ_CP013711", "../../data/klebs_annotations/NZ_CP013711.gff.gz"),
    ("NC_017540", "../../data/klebs_annotations/NC_017540.gff.gz"),
  ];

  #[rstest]
  #[case(0)]
  #[case(1)]
  fn smoke_parse_real_klebs_gff(#[case] i: usize) -> Result<(), Report> {
    let (name, gff) = KLEBS[i];
    let features = GffReader::from_path(gff)?.read_many()?;
    assert!(
      features.len() > 1000,
      "expected thousands of features, got {}",
      features.len()
    );
    assert!(features.iter().any(|f| f.feature_type == "CDS"));
    assert!(features.iter().any(|f| f.strand == Some(Strand::Reverse)));
    assert!(features.iter().any(|f| f.name.is_some()));
    // Seqids are the bare accession (version stripped), matching the graph path name exactly.
    assert!(
      features.iter().all(|f| f.seqid == name),
      "all seqids should be the bare accession {name}"
    );
    Ok(())
  }

  #[test]
  fn smoke_lift_klebs_annotations_to_graph() -> Result<(), Report> {
    let graph = Pangraph::from_path(&Some("../../data/klebs_graph.json.gz"))?;

    let mut all_features = Vec::new();
    for (_, gff) in KLEBS {
      all_features.extend(GffReader::from_path(gff)?.read_many()?);
    }

    // Exact seqid -> path matching (empty seqid map): both annotated genomes group cleanly.
    let grouped = match_features_to_paths(all_features, &graph, &BTreeMap::new())?;
    assert_eq!(grouped.len(), 2);
    for (name, _) in KLEBS {
      let pid = graph.path_id_by_name(name)?;
      assert!(grouped.get(&pid).is_some_and(|features| !features.is_empty()));
    }

    // Lifting thousands of real features onto the real graph succeeds and stays in-bounds.
    let lifted = lift_features(&grouped, &graph)?;
    assert!(
      lifted.len() > 1000,
      "expected many lifted segments, got {}",
      lifted.len()
    );
    let names: BTreeSet<&str> = KLEBS.iter().map(|(n, _)| *n).collect();
    for ann in &lifted {
      assert!(names.contains(ann.genome.as_str()), "unexpected genome {}", ann.genome);
      assert!(ann.cons_start <= ann.cons_end, "consensus coords ordered");
      let block = graph.blocks.get(&ann.block_id).expect("block exists");
      assert!(ann.cons_end <= block.consensus_len(), "consensus end within block");
    }
    Ok(())
  }
}
