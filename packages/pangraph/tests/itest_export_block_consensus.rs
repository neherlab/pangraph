mod common;

#[cfg(test)]
mod tests {
  use eyre::Report;
  use itertools::Itertools;
  use pangraph::commands::export::export_args::PangraphExportBlockConsensusArgs;
  use pangraph::commands::export::export_block_consensus::export_block_consensus;
  use pangraph::io::fasta::read_many_fasta;
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::pangraph::pangraph_block::BlockId;
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;
  use tempfile::tempdir;

  #[test]
  fn itest_export_block_consensus() -> Result<(), Report> {
    let input_json = Some(PathBuf::from("../../data/test_graph.json"));
    let output = tempdir()?.path().join("block_consensus.fa");
    let graph = Pangraph::from_path(&input_json)?;

    export_block_consensus(PangraphExportBlockConsensusArgs {
      input_json,
      output: output.clone(),
    })?;

    let records = read_many_fasta(&[output])?;

    let fasta_bids = records
      .iter()
      .map(|r| BlockId::from_str(&r.seq_name).unwrap())
      .sorted()
      .collect_vec();

    let graph_bids = graph.blocks.keys().copied().sorted().collect_vec();

    assert_eq!(fasta_bids, graph_bids);

    for r in records {
      let bid = BlockId::from_str(&r.seq_name)?;
      assert!(graph.blocks.contains_key(&bid));
      assert_eq!(r.seq, graph.blocks[&bid].consensus());
    }

    Ok(())
  }
}
