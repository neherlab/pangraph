mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use itertools::Itertools;
  use pangraph::commands::export::export_args::PangraphExportCoreAlignmentArgs;
  use pangraph::commands::export::export_core_genome::{ExportCoreAlignmentParams, export_core_genome};
  use pangraph::io::fasta::{Alphabet, FastaReader};
  use pangraph::io::json::{JsonPretty, json_write_str};
  use pangraph::o;
  use pangraph::pangraph::pangraph::Pangraph;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
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

    let alphabet = if aligned {
      Alphabet::DnaWithGap
    } else {
      Alphabet::DnaWithoutGap
    };
    let records = FastaReader::from_paths(&[output])?
      .with_alphabet(alphabet)
      .read_many()?;

    let fasta_names = records.iter().map(|r| &r.seq_name).sorted().collect_vec();
    let path_names = graph.path_names().map(|n| n.unwrap()).sorted().collect_vec();
    assert_eq!(fasta_names, path_names);

    if aligned {
      let lengths: BTreeMap<_, _> = records.iter().map(|r| (r.seq_name.clone(), r.seq.len())).collect();
      let length_0 = *lengths.values().next().unwrap();
      assert!(
        lengths.iter().all(|(_, &len)| len == length_0),
        "Aligned sequences must have the same length, but found:\n{}",
        json_write_str(&lengths, JsonPretty(true))?
      );
    }

    Ok(())
  }
}
