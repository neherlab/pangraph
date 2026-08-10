mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use itertools::Itertools;
  use pangraph::align::alignment_args::GraphMergeParams;
  use pangraph::commands::build::build_args::PangraphBuildArgs;
  use pangraph::commands::build::build_run::build;
  use pangraph::commands::merge::merge_args::PangraphMergeArgs;
  use pangraph::commands::merge::merge_run::merge_run;
  use pangraph::io::fasta::{FastaReader, FastaRecord};
  use pangraph::io::json::{JsonPretty, json_write_file};
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::pangraph::reconstruct::reconstruct;
  use pangraph::representation::seq::Seq;
  use pangraph::utils::error::report_to_string;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use std::path::{Path, PathBuf};
  use tempfile::{TempDir, tempdir};

  /// Reads the first `n_total` records of a FASTA file and splits them in two groups, so that the
  /// two resulting graphs are built from disjoint but homologous sets of genomes.
  fn read_and_split(path: &str, n_left: usize, n_total: usize) -> Result<(Vec<FastaRecord>, Vec<FastaRecord>), Report> {
    let mut fastas = FastaReader::from_paths(&[PathBuf::from(path)])?.read_many()?;
    fastas.truncate(n_total);
    let right = fastas.split_off(n_left);
    Ok((fastas, right))
  }

  /// Builds a graph out of the given records and writes it to a JSON file in `dir`.
  fn build_graph_file(dir: &TempDir, name: &str, fastas: Vec<FastaRecord>) -> Result<PathBuf, Report> {
    let args = PangraphBuildArgs {
      circular: false,
      ..PangraphBuildArgs::default()
    };
    let graph = build(fastas, &args, true)?;
    let path = dir.path().join(name);
    json_write_file(&path, &graph, JsonPretty(false))?;
    Ok(path)
  }

  /// Reconstructs the genomes of a graph, keyed by genome name.
  fn sequences(graph: &Pangraph) -> Result<BTreeMap<String, Seq>, Report> {
    reconstruct(graph).map(|r| r.map(|r| (r.seq_name, r.seq))).collect()
  }

  fn read_graph(path: &Path) -> Result<Pangraph, Report> {
    Pangraph::from_path(&Some(path))
  }

  fn merge_args(left: PathBuf, right: PathBuf, output: PathBuf) -> PangraphMergeArgs {
    PangraphMergeArgs {
      left_graph: left,
      right_graph: right,
      output_json: output,
      merge_params: GraphMergeParams::default(),
      verify: true,
    }
  }

  /// Two graphs built independently from homologous genomes merge into a graph that reconstructs
  /// every input genome exactly. This is the main end-to-end guarantee of the `merge` command.
  #[rstest]
  fn itest_merge_homologous_graphs() -> Result<(), Report> {
    let dir = tempdir()?;
    let (left_fastas, right_fastas) = read_and_split("../../data/ges-1.fa", 3, 6)?;
    let left_names: BTreeSet<String> = left_fastas.iter().map(|f| f.seq_name.clone()).collect();
    let right_names: BTreeSet<String> = right_fastas.iter().map(|f| f.seq_name.clone()).collect();
    let expected_names = left_names.union(&right_names).cloned().collect_vec();

    let left = build_graph_file(&dir, "left.json", left_fastas)?;
    let right = build_graph_file(&dir, "right.json", right_fastas)?;
    let output = dir.path().join("merged.json");

    // `--verify` is on, so this already checks that every genome round-trips
    merge_run(&merge_args(left.clone(), right.clone(), output.clone()))?;

    let merged = read_graph(&output)?;
    #[cfg(debug_assertions)]
    merged.sanity_check()?;

    // all genomes of both inputs are present, exactly once
    assert_eq!(merged.paths.len(), 6);
    assert_eq!(
      sequences(&merged)?.keys().cloned().sorted().collect_vec(),
      expected_names
    );

    // and their sequences are unchanged
    let mut expected = sequences(&read_graph(&left)?)?;
    expected.extend(sequences(&read_graph(&right)?)?);
    assert_eq!(sequences(&merged)?, expected);

    // merging actually found homology across the two graphs: some blocks are now shared between
    // genomes that came from different inputs. (Block counts alone say little: reweaving splits
    // blocks as it merges them, so the total can go either way.)
    let cross_graph_blocks = merged
      .blocks
      .values()
      .filter(|block| {
        let names = block
          .alignment_keys()
          .iter()
          .filter_map(|nid| merged.paths[&merged.nodes[nid].path_id()].name.clone())
          .collect_vec();
        names.iter().any(|n| left_names.contains(n)) && names.iter().any(|n| right_names.contains(n))
      })
      .count();
    assert!(
      cross_graph_blocks > 0,
      "no block is shared between the two input graphs"
    );

    Ok(())
  }

  /// Appending a single genome to an existing graph: the case the command is meant to serve.
  #[rstest]
  fn itest_merge_single_genome() -> Result<(), Report> {
    let dir = tempdir()?;
    let (left_fastas, right_fastas) = read_and_split("../../data/ges-1.fa", 3, 4)?;
    assert_eq!(right_fastas.len(), 1);

    let left = build_graph_file(&dir, "left.json", left_fastas)?;
    let right = build_graph_file(&dir, "right.json", right_fastas)?;
    let output = dir.path().join("merged.json");

    merge_run(&merge_args(left, right, output.clone()))?;

    let merged = read_graph(&output)?;
    #[cfg(debug_assertions)]
    merged.sanity_check()?;
    assert_eq!(merged.paths.len(), 4);

    Ok(())
  }

  /// Merging a graph with itself duplicates every genome name, and must be rejected.
  #[rstest]
  fn itest_merge_rejects_duplicate_genome_names() -> Result<(), Report> {
    let dir = tempdir()?;
    let (fastas, _) = read_and_split("../../data/ges-1.fa", 3, 3)?;
    let graph = build_graph_file(&dir, "graph.json", fastas)?;
    let output = dir.path().join("merged.json");

    let result = merge_run(&merge_args(graph.clone(), graph, output));

    let error = report_to_string(&result.unwrap_err());
    assert!(
      error.contains("Duplicate genome names"),
      "unexpected error message: {error}"
    );

    Ok(())
  }

  /// `build --verify` used to pair reconstructed genomes with input records by `FastaRecord::index`,
  /// which is only valid when the records' indices happen to be exactly `0..n-1`. Here they are
  /// `[3, 4, 5]` for a 3-record slice, which panicked with an out-of-bounds index before genomes
  /// were matched by name.
  ///
  /// Note that the intermediate-clade half of this check only runs in debug builds.
  #[rstest]
  fn itest_build_verify_with_nonzero_record_indices() -> Result<(), Report> {
    let (_, right_fastas) = read_and_split("../../data/ges-1.fa", 3, 6)?;
    assert_eq!(right_fastas.iter().map(|f| f.index).collect_vec(), vec![3, 4, 5]);

    let graph = build(right_fastas, &PangraphBuildArgs::default(), true)?;
    assert_eq!(graph.paths.len(), 3);

    Ok(())
  }

  /// `build` enforces the same uniqueness invariant on its input FASTA records, so that a graph
  /// can never carry duplicate genome names into a later merge.
  #[rstest]
  fn itest_build_rejects_duplicate_sequence_names() -> Result<(), Report> {
    let (mut fastas, _) = read_and_split("../../data/ges-1.fa", 2, 2)?;
    fastas[1].seq_name = fastas[0].seq_name.clone();

    let result = build(fastas, &PangraphBuildArgs::default(), false);

    let error = report_to_string(&result.unwrap_err());
    assert!(
      error.contains("Duplicate sequence names"),
      "unexpected error message: {error}"
    );

    Ok(())
  }
}
