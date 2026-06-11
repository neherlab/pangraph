mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use pangraph::annotation::matching::match_features_to_paths;
  use pangraph::io::gff::GffReader;
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::pangraph::strand::Strand;
  use rstest::rstest;
  use std::collections::BTreeMap;

  // Real Klebsiella annotations (downloaded from NCBI) for two of the genomes in
  // `data/klebs_graph.json.gz`. The path name is the bare accession (the FASTA record id),
  // while the annotation files carry the versioned accession as their seqid.
  //                     path name (bare)   gff path
  const KLEBS: &[(&str, &str)] = &[
    ("NZ_CP013711", "../../data/klebs_annotations/NZ_CP013711.gff.gz"),
    ("NC_017540", "../../data/klebs_annotations/NC_017540.gff.gz"),
  ];

  #[rstest]
  #[case(0)]
  #[case(1)]
  fn smoke_parse_real_klebs_gff(#[case] i: usize) -> Result<(), Report> {
    let (_, gff) = KLEBS[i];
    let features = GffReader::from_path(gff)?.read_many()?;
    assert!(
      features.len() > 1000,
      "expected thousands of features, got {}",
      features.len()
    );
    assert!(features.iter().any(|f| f.feature_type == "CDS"));
    assert!(features.iter().any(|f| f.strand == Some(Strand::Reverse)));
    assert!(features.iter().any(|f| f.name.is_some()));
    Ok(())
  }

  #[test]
  fn smoke_match_klebs_annotations_to_graph() -> Result<(), Report> {
    let graph = Pangraph::from_path(&Some("../../data/klebs_graph.json.gz"))?;

    let mut all_features = Vec::new();
    let mut seqid_map = BTreeMap::new();
    for (name, gff) in KLEBS {
      let features = GffReader::from_path(gff)?.read_many()?;
      // The annotation seqid is the versioned accession (e.g. `NZ_CP013711.1`); map it onto
      // the bare-accession path name used when the graph was built.
      let seqid = features.first().expect("annotation has features").seqid.clone();
      seqid_map.insert(seqid, (*name).to_owned());
      all_features.extend(features);
    }

    let grouped = match_features_to_paths(all_features, &graph, &seqid_map)?;
    assert_eq!(grouped.len(), 2);
    for (name, _) in KLEBS {
      let pid = graph.path_id_by_name(name)?;
      assert!(grouped.get(&pid).is_some_and(|features| !features.is_empty()));
    }
    Ok(())
  }
}
