#[cfg(test)]
mod tests {
  use eyre::Report;
  use pangraph::commands::annotate::annotate_args::{
    AnnotateCommonArgs, PangraphAnnotateArgs, PangraphAnnotateBlocksArgs, PangraphAnnotateNodesArgs,
  };
  use pangraph::commands::annotate::annotate_run::annotate_run;
  use pangraph::io::csv::parse_csv;
  use pangraph::io::file::open_file_or_stdin;
  use pangraph::pangraph::pangraph::Pangraph;
  use serde::Deserialize;
  use std::fs::{read_to_string, write};
  use std::io::Read;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;

  const GRAPH: &str = "../../data/test_graph.json";
  const NODE_HEADER_PREFIX: &str = "feature_id,parent_feature_id,segment_idx";
  const BLOCK_HEADER_PREFIX: &str = "type,strand_on_consensus,start_block_id";

  /// Subset of the node-level CSV columns; remaining columns are ignored on deserialization.
  #[derive(Deserialize)]
  struct Row {
    genome: String,
    parent_feature_id: String,
  }

  /// Name of the first path in the test graph (a valid `seqid` to match against).
  fn first_path_name() -> Result<String, Report> {
    let graph = Pangraph::from_path(&Some(GRAPH))?;
    let path = graph.paths.values().next().expect("at least one path");
    Ok(path.name().as_deref().unwrap().to_owned())
  }

  /// Write a small GFF3 file (1-based inclusive coords) with two features on `seqid`: a short
  /// forward gene (single node) and a long reverse CDS (spans block boundaries).
  fn write_gff(dir: &Path, seqid: &str) -> Result<PathBuf, Report> {
    let gff = format!(
      "##gff-version 3\n\
       {seqid}\ttest\tgene\t1001\t1200\t.\t+\t.\tID=gene1;Name=g1\n\
       {seqid}\ttest\tCDS\t101\t5000\t.\t-\t0\tID=cds1;Name=c1\n"
    );
    let path = dir.join("ann.gff");
    write(&path, gff)?;
    Ok(path)
  }

  /// Build `annotate nodes` args against the test graph for the given GFFs and output path.
  fn nodes_args(gff: Vec<PathBuf>, output: PathBuf) -> PangraphAnnotateArgs {
    PangraphAnnotateArgs::Nodes(PangraphAnnotateNodesArgs {
      common: AnnotateCommonArgs {
        input: Some(PathBuf::from(GRAPH)),
        gff,
        output,
      },
    })
  }

  /// Build `annotate blocks` args against the test graph; `property_threshold` is fixed at 0.5.
  fn blocks_args(gff: Vec<PathBuf>, output: PathBuf, min_frequency: f64) -> PangraphAnnotateArgs {
    PangraphAnnotateArgs::Blocks(PangraphAnnotateBlocksArgs {
      common: AnnotateCommonArgs {
        input: Some(PathBuf::from(GRAPH)),
        gff,
        output,
      },
      min_frequency,
      property_threshold: 0.5,
    })
  }

  #[test]
  fn itest_annotate_cli_node_level_csv() -> Result<(), Report> {
    let name = first_path_name()?;
    let dir = tempdir()?;
    let gff = write_gff(dir.path(), &name)?;
    let out = dir.path().join("lifted.csv");

    annotate_run(nodes_args(vec![gff], out.clone()))?;

    let contents = read_to_string(&out)?;
    assert!(contents.starts_with(NODE_HEADER_PREFIX), "header row present");

    let rows: Vec<Row> = parse_csv(&contents)?;
    assert!(!rows.is_empty(), "at least one lifted segment");
    // Both source features were lifted (parent_feature_id carries the GFF ID).
    assert!(rows.iter().any(|r| r.parent_feature_id == "gene1"), "gene1 lifted");
    assert!(rows.iter().any(|r| r.parent_feature_id == "cds1"), "cds1 lifted");
    // Every row belongs to the genome we annotated.
    for r in &rows {
      assert_eq!(r.genome, name, "genome column is the path name");
    }
    Ok(())
  }

  #[test]
  fn itest_annotate_cli_blocks_csv() -> Result<(), Report> {
    let name = first_path_name()?;
    let dir = tempdir()?;
    let gff = write_gff(dir.path(), &name)?;
    let out = dir.path().join("blocks.csv");

    // `min_frequency: 0.0` keeps every cluster, so the single annotated genome still produces rows.
    annotate_run(blocks_args(vec![gff], out.clone(), 0.0))?;

    let contents = read_to_string(&out)?;
    assert!(
      contents.starts_with(BLOCK_HEADER_PREFIX),
      "block-level header row present"
    );
    // At least one consensus feature emitted (header line plus one or more data lines).
    assert!(
      contents.lines().filter(|l| !l.is_empty()).count() >= 2,
      "at least one block-level annotation row"
    );
    Ok(())
  }

  #[test]
  fn itest_annotate_cli_unmatched_seqid_errors() -> Result<(), Report> {
    // `example.gff` has seqids `chr1`/`chr2`, which match no path in the test graph: exact-match
    // policy must fail loudly and name every offending seqid.
    let dir = tempdir()?;
    let out = dir.path().join("lifted.csv");

    let err = annotate_run(nodes_args(vec![PathBuf::from("../../data/example.gff")], out)).unwrap_err();

    let msg = err.to_string();
    assert!(msg.contains("chr1"), "error names chr1: {msg}");
    assert!(msg.contains("chr2"), "error names chr2: {msg}");
    Ok(())
  }

  #[test]
  fn itest_annotate_cli_gz_output_roundtrips() -> Result<(), Report> {
    let name = first_path_name()?;
    let dir = tempdir()?;
    let gff = write_gff(dir.path(), &name)?;
    let out = dir.path().join("lifted.csv.gz");

    annotate_run(nodes_args(vec![gff], out.clone()))?;

    // The `.gz` output is real gzip: read it back through transparent decompression.
    let mut decompressed = String::new();
    open_file_or_stdin(&Some(&out))?.read_to_string(&mut decompressed)?;
    assert!(
      decompressed.starts_with(NODE_HEADER_PREFIX),
      "decompressed header present"
    );
    Ok(())
  }
}
