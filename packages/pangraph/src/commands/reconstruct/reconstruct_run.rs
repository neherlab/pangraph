use crate::commands::reconstruct::reconstruct_args::PangraphReconstructArgs;
use crate::io::fasta::{FastaRecord, FastaWriter};
use crate::io::json::json_read_file;
use crate::io::seq::reverse_complement;
use crate::make_internal_report;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_path::PangraphPath;
use eyre::Report;
use itertools::Itertools;

pub fn reconstruct_run(args: &PangraphReconstructArgs) -> Result<(), Report> {
  let PangraphReconstructArgs {
    input_graph,
    output_fasta,
  } = &args;

  let graph: Pangraph = json_read_file(input_graph)?;

  let mut writer = FastaWriter::from_path(output_fasta)?;

  reconstruct(&graph).try_for_each(|fasta| {
    let fasta = fasta?;
    writer.write(fasta.seq_name, fasta.seq)
  })?;

  Ok(())
}

// Note: this function could belong to the path class.
fn reconstruct_path_sequence(path: &PangraphPath, graph: &Pangraph) -> Result<String, Report> {
  // collect here the sequence block by block
  let mut seq = Vec::new();

  for node_id in &path.nodes {
    let node = graph
      .nodes
      .get(node_id)
      .ok_or_else(|| make_internal_report!("Node {node_id} not found in graph"))?;

    let block_id = node.block_id();
    let block = graph
      .blocks
      .get(&block_id)
      .ok_or_else(|| make_internal_report!("Block {block_id} not found in graph"))?;

    // Get edits and apply them to the consensus sequence
    let edits = block.alignment(*node_id);

    let mut s = edits.apply(block.consensus())?;

    // Reverse-complement if on opposite strand
    if node.strand().is_reverse() {
      s = reverse_complement(s)?;
    }

    seq.push(s);
  }

  Ok(seq.join(""))
}

pub fn reconstruct(graph: &Pangraph) -> impl Iterator<Item = Result<FastaRecord, Report>> + '_ {
  graph
    .paths
    .iter()
    .sorted_by_key(|(path_id, _)| **path_id)
    .map(|(path_id, path)| {
      let index = path_id.0;
      let seq = reconstruct_path_sequence(path, graph)?;
      let seq_name = path
        .name()
        .clone()
        .unwrap_or_else(|| format!("Unknown sequence #{path_id}"));
      Ok(FastaRecord { seq_name, seq, index })
    })
}
