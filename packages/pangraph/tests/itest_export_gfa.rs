mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use itertools::Itertools;
  use pangraph::commands::export::export_args::PangraphExportGfaArgs;
  use pangraph::commands::export::export_gfa::export_gfa;
  use pangraph::io::fs::read_file_to_string;
  use pangraph::io::gfa::GfaWriteParams;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::PathBuf;
  use tempfile::tempdir;

  #[rustfmt::skip]
  #[rstest]
  //      min_length  min_depth   export        n expected
  //                              duplicated    gfa segments
  #[case( Some(1000), Some(2),    true,         8          )]
  #[case( Some(1000), Some(2),    false,        7          )]
  #[case( None,       None,       true,         14         )]
  #[case( None,       None,       false,        13         )]
  //
  #[trace]
  fn itest_export_gfa(
    #[case] minimum_length: Option<usize>,
    #[case] minimum_depth: Option<usize>,
    #[case] export_duplicated: bool,
    #[case] n_expected_gfa_segments: usize,
  ) -> Result<(), Report> {
    let input_json = Some(PathBuf::from("../../data/test_graph.json"));
    let output = tempdir()?.path().join("output.gfa");

    export_gfa(PangraphExportGfaArgs {
      input_json,
      output: output.clone(),
      params: GfaWriteParams {
        minimum_length,
        maximum_length: None,
        minimum_depth,
        maximum_depth: None,
        include_sequences: false,
        no_duplicated: !export_duplicated,
      },
    })?;

    let output = read_file_to_string(output)?;
    let lines = output.lines().collect_vec();

    let segment_lines = lines.iter().filter(|l| l.starts_with("S\t")).collect_vec();
    assert_eq!(segment_lines.len(), n_expected_gfa_segments);

    Ok(())
  }
}
