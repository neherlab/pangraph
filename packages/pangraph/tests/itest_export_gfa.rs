mod common;

#[cfg(test)]
mod tests {

  use eyre::Report;
  use itertools::Itertools;
  use pangraph::commands::export::export_args::PangraphExportGfaArgs;
  use pangraph::commands::export::export_gfa::export_gfa;
  use pangraph::io::fs::read_file_to_string;
  use pangraph::io::gfa::GfaWriteParams;
  use pangraph::pangraph::pangraph::Pangraph;
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;
  use tempfile::tempdir;

  #[test]
  fn itest_export_gfa() -> Result<(), Report> {
    let input_json = Some(PathBuf::from("../../data/test_graph.json"));
    let output = tempdir()?.path().join("output.gfa");
    let graph = Pangraph::from_path(&input_json)?;

    export_gfa(PangraphExportGfaArgs {
      input_json,
      output: output.clone(),
      params: GfaWriteParams {
        minimum_length: Some(100),
        maximum_length: None,
        minimum_depth: None,
        maximum_depth: None,
        include_sequences: false,
        no_duplicated: false,
      },
    })?;

    let output = read_file_to_string(output)?;
    let lines = output.lines().collect_vec();

    let segment_lines = lines.iter().filter(|l| l.starts_with("S\t")).collect_vec();
    assert_eq!(segment_lines.len(), graph.blocks.len());

    Ok(())
  }
}
