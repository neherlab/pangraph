mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use itertools::Itertools;
  use pangraph::commands::export::export_args::PangraphExportCoreAlignmentArgs;
  use pangraph::commands::export::export_core_genome::{export_core_genome, ExportCoreAlignmentParams};
  use pangraph::io::fasta::read_many_fasta;
  use pangraph::o;
  use pangraph::pangraph::pangraph::Pangraph;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::PathBuf;
  use tempfile::tempdir;

  #[rstest]
  #[case::aligned(true)]
  #[case::unaligned(false)]
  #[trace]
  fn itest_export_core_genome(#[case] aligned: bool) -> Result<(), Report> {
    let input_json = Some(PathBuf::from("../../data/test_graph.json"));
    let output = tempdir()?.path().join("core_alignment.fa");
    let guide_strain = o!("pCAV1344-40");

    let graph = Pangraph::from_path(&input_json)?;

    export_core_genome(PangraphExportCoreAlignmentArgs {
      input_json,
      output: output.clone(),
      params: ExportCoreAlignmentParams {
        guide_strain,
        unaligned: !aligned,
      },
    })?;

    let records = read_many_fasta(&[output])?;

    let fasta_names = records.iter().map(|r| &r.seq_name).sorted().collect_vec();
    let path_names = graph.path_names().map(|n| n.unwrap()).sorted().collect_vec();
    assert_eq!(fasta_names, path_names);

    Ok(())
  }
}
