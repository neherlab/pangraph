mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use pangraph::commands::export::export_args::PangraphExportBlockSequencesArgs;
  use pangraph::commands::export::export_block_sequences::{ExportBlockSequencesParams, export_block_sequences};
  use pangraph::io::fasta::read_many_fasta;
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::pangraph::pangraph_node::NodeId;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::PathBuf;
  use tempfile::tempdir;

  #[rstest]
  #[case::aligned(true)]
  #[case::unaligned(false)]
  #[trace]
  fn itest_export_block_sequences(#[case] aligned: bool) -> Result<(), Report> {
    let input_json = Some(PathBuf::from("../../data/test_graph.json"));
    let output = tempdir()?.path().to_owned();
    let graph = Pangraph::from_path(&input_json)?;

    export_block_sequences(PangraphExportBlockSequencesArgs {
      input_json,
      output: output.clone(),
      params: ExportBlockSequencesParams { unaligned: !aligned },
    })?;

    for (block_id, block) in &graph.blocks {
      let block_fa = output.join(format!("block_{block_id}.fa"));
      let records = read_many_fasta(&[block_fa])?;

      assert_eq!(records.len(), block.alignments().len());

      for r in records {
        let node_id = NodeId::from_str(r.seq_name.split_whitespace().next().unwrap())?;
        assert!(block.alignment_keys().contains(&node_id));

        if aligned {
          assert_eq!(r.seq.len(), block.consensus_len());
        } else {
          let seq_len = block.unaligned_len_for_node(&node_id);
          assert_eq!(r.seq.len(), seq_len, "node_id: {}", node_id);
        }
      }
    }

    Ok(())
  }
}
