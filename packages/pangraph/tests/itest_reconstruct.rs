mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pangraph::align::alignment_args::GraphMergeParams;
  use pangraph::commands::build::build_args::PangraphBuildArgs;
  use pangraph::commands::build::build_run::build;
  use pangraph::commands::merge::merge_args::PangraphMergeArgs;
  use pangraph::commands::merge::merge_run::merge_run;
  use pangraph::commands::reconstruct::reconstruct_args::PangraphReconstructArgs;
  use pangraph::commands::reconstruct::reconstruct_run::reconstruct_run;
  use pangraph::io::fasta::{FastaReader, FastaRecord, FastaWriter};
  use pangraph::io::json::{JsonPretty, json_write_file};
  use pangraph::pangraph::edits::Edit;
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::pangraph::pangraph_block::{BlockId, PangraphBlock};
  use pangraph::pangraph::pangraph_node::{NodeId, PangraphNode};
  use pangraph::pangraph::pangraph_path::{PangraphPath, PathId};
  use pangraph::pangraph::strand::Strand::Forward;
  use pangraph::representation::seq::Seq;
  use pangraph::utils::error::report_to_string;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::{Path, PathBuf};
  use tempfile::{TempDir, tempdir};

  const GES: &str = "../../data/ges-1.fa";

  /// Builds a merged graph from two disjoint halves of a FASTA file, and returns its path together
  /// with all the input records. The merger renumbers path ids, which is what makes this graph the
  /// interesting case: neither record order nor `FastaRecord::index` can pair genomes up any more.
  fn merged_graph(dir: &TempDir) -> Result<(PathBuf, Vec<FastaRecord>), Report> {
    let mut fastas = FastaReader::from_paths(&[PathBuf::from(GES)])?.read_many()?;
    fastas.truncate(6);
    let right_fastas = fastas.split_off(3);
    let all = fastas.iter().chain(right_fastas.iter()).cloned().collect_vec();

    let args = PangraphBuildArgs::default();
    let left = write_graph(dir, "left.json", &build(fastas, &args, true)?)?;
    let right = write_graph(dir, "right.json", &build(right_fastas, &args, true)?)?;

    let output = dir.path().join("merged.json");
    merge_run(&PangraphMergeArgs {
      left_graph: left,
      right_graph: right,
      output_json: output.clone(),
      merge_params: GraphMergeParams::default(),
      verify: true,
    })?;

    Ok((output, all))
  }

  fn write_graph(dir: &TempDir, name: &str, graph: &Pangraph) -> Result<PathBuf, Report> {
    let path = dir.path().join(name);
    json_write_file(&path, graph, JsonPretty(false))?;
    Ok(path)
  }

  fn write_fasta(dir: &TempDir, name: &str, records: &[FastaRecord]) -> Result<PathBuf, Report> {
    let path = dir.path().join(name);
    let mut writer = FastaWriter::from_path(&path)?;
    for record in records {
      writer.write(record.seq_name.clone(), &record.desc, &record.seq)?;
    }
    drop(writer);
    Ok(path)
  }

  fn verify_args(graph: &Path, verify: &Path) -> PangraphReconstructArgs {
    PangraphReconstructArgs {
      input_graph: Some(graph.to_owned()),
      output_fasta: PathBuf::from("-"),
      verify: Some(verify.to_owned()),
    }
  }

  /// The headline case: a merged graph verified against its genomes in a different order. Before
  /// verification was keyed by name this failed spuriously, reporting a mismatch between two
  /// sequences of identical length, because whole `FastaRecord`s were compared and `index` differs
  /// once a merge has renumbered the path ids.
  #[rstest]
  fn itest_reconstruct_verify_merged_graph_ignores_order() -> Result<(), Report> {
    let dir = tempdir()?;
    let (graph, mut records) = merged_graph(&dir)?;
    records.reverse();

    let verify = write_fasta(&dir, "verify.fa", &records)?;
    reconstruct_run(&verify_args(&graph, &verify))?;

    Ok(())
  }

  /// A record the graph does not contain used to be read past and silently ignored, so the command
  /// exited 0 while the verification file and the graph disagreed.
  #[rstest]
  fn itest_reconstruct_verify_rejects_surplus_record() -> Result<(), Report> {
    let dir = tempdir()?;
    let (graph, mut records) = merged_graph(&dir)?;
    records.push(FastaRecord {
      seq_name: "not_in_the_graph".to_owned(),
      desc: None,
      seq: Seq::from_str("ACGTACGT"),
      index: 0,
    });

    let verify = write_fasta(&dir, "verify.fa", &records)?;
    let err = report_to_string(&reconstruct_run(&verify_args(&graph, &verify)).unwrap_err());

    assert!(err.contains("not_in_the_graph"), "unexpected error: {err}");
    assert!(
      err.contains("which the graph does not contain"),
      "unexpected error: {err}"
    );

    Ok(())
  }

  #[rstest]
  fn itest_reconstruct_verify_rejects_missing_record() -> Result<(), Report> {
    let dir = tempdir()?;
    let (graph, mut records) = merged_graph(&dir)?;
    let dropped = records.remove(0).seq_name;

    let verify = write_fasta(&dir, "verify.fa", &records)?;
    let err = report_to_string(&reconstruct_run(&verify_args(&graph, &verify)).unwrap_err());

    assert!(err.contains(&dropped), "unexpected error: {err}");
    assert!(err.contains("missing"), "unexpected error: {err}");

    Ok(())
  }

  #[rstest]
  fn itest_reconstruct_verify_rejects_duplicate_record() -> Result<(), Report> {
    let dir = tempdir()?;
    let (graph, mut records) = merged_graph(&dir)?;
    records.push(records[0].clone());

    let verify = write_fasta(&dir, "verify.fa", &records)?;
    let err = report_to_string(&reconstruct_run(&verify_args(&graph, &verify)).unwrap_err());

    assert!(err.contains("more than once"), "unexpected error: {err}");

    Ok(())
  }

  /// A real sequence difference must be reported with the genome name and the position, rather than
  /// the old "expected length N but got N".
  #[rstest]
  fn itest_reconstruct_verify_reports_mutated_base() -> Result<(), Report> {
    let dir = tempdir()?;
    let (graph, mut records) = merged_graph(&dir)?;

    let mutated = records[1].seq_name.clone();
    let seq = records[1].seq.as_str().to_owned();
    let original = seq.as_bytes()[100] as char;
    let replacement = if original == 'A' { 'T' } else { 'A' };
    let mut bases = seq.into_bytes();
    bases[100] = replacement as u8;
    records[1].seq = Seq::from_str(std::str::from_utf8(&bases)?);

    let verify = write_fasta(&dir, "verify.fa", &records)?;
    let err = report_to_string(&reconstruct_run(&verify_args(&graph, &verify)).unwrap_err());

    assert!(err.contains(&mutated), "unexpected error: {err}");
    assert!(err.contains("at position 100"), "unexpected error: {err}");

    Ok(())
  }

  /// Verification matches by name, so a graph with unnamed paths cannot be verified at all.
  #[rstest]
  fn itest_reconstruct_verify_rejects_unnamed_path() -> Result<(), Report> {
    let dir = tempdir()?;
    let graph = Pangraph {
      blocks: btreemap! {
        BlockId(0) => PangraphBlock::new(BlockId(0), "ACGTACGT", btreemap!{ NodeId(0) => Edit::empty() }),
      },
      nodes: btreemap! {
        NodeId(0) => PangraphNode::new(NodeId(0), BlockId(0), PathId(0), Forward, (0, 8)),
      },
      paths: btreemap! {
        PathId(0) => PangraphPath::new(PathId(0), [NodeId(0)], 8, false, None, None),
      },
    };

    let graph = write_graph(&dir, "unnamed.json", &graph)?;
    let verify = write_fasta(
      &dir,
      "verify.fa",
      &[FastaRecord {
        seq_name: "a".to_owned(),
        desc: None,
        seq: Seq::from_str("ACGTACGT"),
        index: 0,
      }],
    )?;

    let err = report_to_string(&reconstruct_run(&verify_args(&graph, &verify)).unwrap_err());
    assert!(err.contains("no name or an empty name"), "unexpected error: {err}");

    Ok(())
  }

  /// The non-verify path is untouched: every genome is written out, and matching them by name
  /// recovers the inputs exactly. Their order is deliberately not asserted.
  #[rstest]
  fn itest_reconstruct_writes_every_genome() -> Result<(), Report> {
    let dir = tempdir()?;
    let (graph, records) = merged_graph(&dir)?;
    let output = dir.path().join("out.fa");

    reconstruct_run(&PangraphReconstructArgs {
      input_graph: Some(graph),
      output_fasta: output.clone(),
      verify: None,
    })?;

    let written: BTreeMap<String, Seq> = FastaReader::from_path(&output)?
      .read_many()?
      .into_iter()
      .map(|record| (record.seq_name, record.seq))
      .collect();
    let expected: BTreeMap<String, Seq> = records
      .into_iter()
      .map(|record| (record.seq_name, record.seq))
      .collect();

    assert_eq!(written, expected);

    Ok(())
  }

  /// A graph read from a file was not necessarily written by pangraph, so nothing guarantees its
  /// cross-references resolve. `reconstruct` used to index the node map directly, so a path naming
  /// a node the graph does not contain aborted the process with a bare "no entry found for key"
  /// panic, and in a release build, where `sanity_check` is compiled out, with no indication of
  /// which file was at fault. Validation now happens where the file name is still known.
  #[rstest]
  fn itest_reconstruct_rejects_malformed_graph_naming_the_file() -> Result<(), Report> {
    let dir = tempdir()?;

    let mut graph = Pangraph {
      blocks: btreemap! {
        BlockId(0) => PangraphBlock::new(BlockId(0), "ACGTACGT", btreemap! { NodeId(0) => Edit::empty() }),
      },
      nodes: btreemap! {
        NodeId(0) => PangraphNode::new(NodeId(0), BlockId(0), PathId(0), Forward, (0, 8)),
      },
      paths: btreemap! {
        PathId(0) => PangraphPath::new(PathId(0), [NodeId(0)], 8, false, Some("a".to_owned()), None),
      },
    };
    graph.paths.get_mut(&PathId(0)).unwrap().nodes = vec![NodeId(99)];
    let graph = write_graph(&dir, "malformed.json", &graph)?;

    let args = PangraphReconstructArgs {
      input_graph: Some(graph),
      output_fasta: PathBuf::from("-"),
      verify: None,
    };
    let err = report_to_string(&reconstruct_run(&args).unwrap_err());

    assert!(
      err.contains("Node 99 from path 0 not found in graph"),
      "unexpected error: {err}"
    );
    assert!(err.contains("malformed.json"), "unexpected error: {err}");

    Ok(())
  }
}
