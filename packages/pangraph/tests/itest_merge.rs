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
  use pangraph::pangraph::pangraph_block::{BlockId, PangraphBlock};
  use pangraph::pangraph::pangraph_node::{NodeId, PangraphNode};
  use pangraph::pangraph::pangraph_path::{PangraphPath, PathId};
  use pangraph::pangraph::reconstruct::reconstruct;
  use pangraph::pangraph::strand::Strand::Forward;
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

  /// Reads the first `n` records of a FASTA file.
  fn read_records(path: &str, n: usize) -> Result<Vec<FastaRecord>, Report> {
    let mut fastas = FastaReader::from_paths(&[PathBuf::from(path)])?.read_many()?;
    fastas.truncate(n);
    Ok(fastas)
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
    assert_eq!(merged.paths.len(), 4);

    Ok(())
  }

  /// Appending to a graph that is itself the result of an earlier merge. Identifiers minted by the
  /// first merge survive into its output whenever a block or node finds no homologue, so this is
  /// the case that used to need a varying relabeling salt to keep the third graph's identifiers off
  /// them. Deriving identifiers from genome names removes the problem at the source. The three
  /// genomes here are mutually unrelated, so nothing aligns and every identifier survives;
  /// homologous appends never exercised this.
  #[rstest]
  fn itest_merge_appends_to_an_already_merged_graph() -> Result<(), Report> {
    let dir = tempdir()?;

    let first = build_graph_file(&dir, "first.json", read_records("../../data/flu-h1.fa", 2)?)?;
    let second = build_graph_file(&dir, "second.json", read_records("../../data/sc2.fa", 1)?)?;
    let third = build_graph_file(&dir, "third.json", read_records("../../data/mpox.fa", 1)?)?;

    let merged_once = dir.path().join("merged-once.json");
    merge_run(&merge_args(first, second, merged_once.clone()))?;
    assert_eq!(read_graph(&merged_once)?.paths.len(), 3);

    let merged_twice = dir.path().join("merged-twice.json");
    merge_run(&merge_args(merged_once, third, merged_twice.clone()))?;

    let merged = read_graph(&merged_twice)?;
    assert_eq!(merged.paths.len(), 4);

    Ok(())
  }

  /// Maps each genome name to the ids of the blocks its path walks through, in order. Comparing the
  /// id *sets* of two graphs would not do: `build` used to label singleton blocks `0, 1, 2, ...`
  /// whatever the genome, so the sets matched even when every genome held a different block.
  fn blocks_by_genome(graph: &Pangraph) -> BTreeMap<String, Vec<BlockId>> {
    graph
      .paths
      .values()
      .map(|path| {
        let blocks = path.nodes().iter().map(|nid| graph.nodes[nid].block_id()).collect_vec();
        (path.name().clone().unwrap_or_default(), blocks)
      })
      .collect()
  }

  /// Block and node identifiers are derived from genome names rather than from the order the input
  /// sequences were read in, so the same genomes under the same guide tree produce the same graph
  /// whichever order they arrive in. Only path ids, which deliberately record the input order, are
  /// expected to differ.
  #[rstest]
  fn itest_build_ids_do_not_depend_on_input_order() -> Result<(), Report> {
    let dir = tempdir()?;
    let fastas = read_records("../../data/ges-1.fa", 4)?;
    let names = fastas.iter().map(|f| f.seq_name.clone()).collect_vec();

    // Pinned, so that the input order is the only thing that varies between the two builds.
    let newick = dir.path().join("guide.nwk");
    std::fs::write(
      &newick,
      format!("(({},{}),({},{}));", names[0], names[1], names[2], names[3]),
    )?;

    let build_in_order = |fastas: Vec<FastaRecord>| -> Result<Pangraph, Report> {
      let args = PangraphBuildArgs {
        circular: false,
        guide_tree: Some(newick.clone()),
        ..PangraphBuildArgs::default()
      };
      build(fastas, &args, true)
    };

    // Reversed *and* re-indexed, which is what reading the same genomes from a reordered file does.
    let mut backwards = fastas.clone();
    backwards.reverse();
    for (index, record) in backwards.iter_mut().enumerate() {
      record.index = index;
    }

    let forward = build_in_order(fastas)?;
    let reversed = build_in_order(backwards)?;

    assert_eq!(blocks_by_genome(&forward), blocks_by_genome(&reversed));
    assert_eq!(
      forward.node_ids().collect::<BTreeSet<_>>(),
      reversed.node_ids().collect::<BTreeSet<_>>()
    );
    assert_eq!(sequences(&forward)?, sequences(&reversed)?);

    // path ids still record the input order, so the genomes come back in the order they were read
    assert_eq!(forward.path_names().flatten().collect_vec(), names);
    assert_eq!(
      reversed.path_names().flatten().collect_vec(),
      names.iter().rev().collect_vec()
    );

    Ok(())
  }

  /// A graph written by pangraph 1.3 or earlier labels its first genome with block, node and path
  /// id `0`, whatever that genome is called, so two of them collide even though their names differ.
  /// `merge` has to say so: without the check, `graph_join` panics on the conflicting key.
  #[rstest]
  fn itest_merge_rejects_graphs_with_colliding_identifiers() -> Result<(), Report> {
    let dir = tempdir()?;

    // A singleton graph in the pre-1.4 identifier scheme, where ids came from the record index.
    let legacy_graph = |name: &str, seq: &str| {
      let (bid, nid, pid) = (BlockId(0), NodeId(0), PathId(0));
      Pangraph {
        blocks: BTreeMap::from([(bid, PangraphBlock::from_consensus(seq, bid, nid))]),
        nodes: BTreeMap::from([(nid, PangraphNode::new(Some(nid), bid, pid, Forward, (0, seq.len())))]),
        paths: BTreeMap::from([(
          pid,
          PangraphPath::new(Some(pid), [nid], seq.len(), false, Some(name.to_owned()), None),
        )]),
      }
    };

    let left = dir.path().join("left.json");
    let right = dir.path().join("right.json");
    json_write_file(&left, &legacy_graph("genome_a", "ACGTACGTAC"), JsonPretty(false))?;
    json_write_file(&right, &legacy_graph("genome_b", "TTTTGGGGCC"), JsonPretty(false))?;

    let result = merge_run(&merge_args(left, right, dir.path().join("out.json")));

    let error = report_to_string(&result.unwrap_err());
    assert!(
      error.contains("share block or node identifiers"),
      "unexpected error message: {error}"
    );

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

  /// Unnamed paths used to slip past the duplicate-name guard entirely: it scanned
  /// `path_names().flatten()`, which drops `None`, and the check that would have caught them ran
  /// only under `--verify`. Merging such a graph with itself therefore succeeded and silently
  /// emitted every genome twice. `verify: false` is the whole point of this test.
  #[rstest]
  fn itest_merge_rejects_unnamed_genomes_without_verify() -> Result<(), Report> {
    let dir = tempdir()?;
    let (fastas, _) = read_and_split("../../data/ges-1.fa", 3, 3)?;

    let mut graph = read_graph(&build_graph_file(&dir, "named.json", fastas)?)?;
    for path in graph.paths.values_mut() {
      path.name = None;
    }
    let anonymous = dir.path().join("anonymous.json");
    json_write_file(&anonymous, &graph, JsonPretty(false))?;

    let output = dir.path().join("merged.json");
    let result = merge_run(&PangraphMergeArgs {
      verify: false,
      ..merge_args(anonymous.clone(), anonymous, output.clone())
    });

    let error = report_to_string(&result.unwrap_err());
    assert!(
      error.contains("no name or an empty name"),
      "unexpected error message: {error}"
    );
    assert!(!output.exists(), "a rejected merge must not write an output graph");

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

  /// A FASTA header of the form `> id` leaves the record unnamed, with the identifier in the
  /// description. Such a genome could not be addressed by name afterwards, so `build` rejects it.
  #[rstest]
  fn itest_build_rejects_empty_sequence_name() -> Result<(), Report> {
    let (mut fastas, _) = read_and_split("../../data/ges-1.fa", 2, 2)?;
    fastas[1].desc = Some(fastas[1].seq_name.clone());
    fastas[1].seq_name = String::new();

    let result = build(fastas, &PangraphBuildArgs::default(), false);

    let error = report_to_string(&result.unwrap_err());
    assert!(error.contains("empty name"), "unexpected error message: {error}");

    Ok(())
  }
}
