use crate::commands::reconstruct::reconstruct_args::PangraphReconstructArgs;
use crate::io::fasta::{FastaReader, FastaRecord, FastaWriter};
use crate::io::json::json_read_file;
use crate::io::seq::reverse_complement;
use crate::make_internal_report;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_node::NodeId;
use crate::pangraph::pangraph_path::PangraphPath;
use crate::representation::seq::Seq;
use eyre::Report;
use itertools::Itertools;
use log::info;

pub fn reconstruct_run(args: &PangraphReconstructArgs) -> Result<(), Report> {
  let PangraphReconstructArgs {
    input_graph,
    output_fasta,
    verify,
  } = &args;

  let graph: Pangraph = json_read_file(input_graph)?;
  let mut results = reconstruct(&graph);

  if let Some(verify) = verify {
    info!("Verifying sequences reconstructed from pangenome graph");
    let mut reader = FastaReader::from_path(verify)?;
    results.try_for_each(|actual| -> Result<(), Report> {
      let actual = actual?;
      let mut expected = FastaRecord::new();
      reader.read(&mut expected)?;
      compare_sequences(&expected, &actual);
      Ok(())
    })?;
  } else {
    let mut writer = FastaWriter::from_path(output_fasta)?;
    results.try_for_each(|fasta| {
      let fasta = fasta?;
      writer.write(fasta.seq_name, &fasta.desc, &fasta.seq)
    })?;
  }

  Ok(())
}

pub fn compare_sequences(left: &FastaRecord, right: &FastaRecord) -> bool {
  pretty_assertions::assert_eq!(left, right);
  true
}

pub fn reconstruct(graph: &Pangraph) -> impl Iterator<Item = Result<FastaRecord, Report>> + '_ {
  graph
    .paths
    .iter()
    .sorted_by_key(|(path_id, _)| **path_id)
    .map(|(path_id, path)| {
      let index = path_id.0;
      let seq = reconstruct_path_sequence(graph, path)?;
      let seq_name = path
        .name()
        .clone()
        .unwrap_or_else(|| format!("Unknown sequence #{path_id}"));
      let desc = path.desc().clone();
      Ok(FastaRecord {
        seq_name,
        desc,
        seq,
        index,
      })
    })
}

fn reconstruct_path_sequence(graph: &Pangraph, path: &PangraphPath) -> Result<Seq, Report> {
  if let Some(first_node_id) = path.nodes.first() {
    let first_node_pos = graph.nodes[first_node_id].position().0;

    let mut genome: Seq = path
      .nodes
      .iter()
      .map(|node_id| reconstruct_block_sequence(graph, *node_id))
      .collect::<Result<Seq, Report>>()?;

    let genome_len = path.tot_len();
    assert_eq!(genome.len(), genome_len);

    genome.rotate_right(first_node_pos);

    Ok(genome)
  } else {
    Ok(Seq::new())
  }
}

fn reconstruct_block_sequence(graph: &Pangraph, node_id: NodeId) -> Result<Seq, Report> {
  let node = graph
    .nodes
    .get(&node_id)
    .ok_or_else(|| make_internal_report!("Node {node_id} not found in graph"))?;

  let block_id = node.block_id();
  let block = graph
    .blocks
    .get(&block_id)
    .ok_or_else(|| make_internal_report!("Block {block_id} not found in graph"))?;

  // Get edits and apply them to the consensus sequence
  let edits = block.alignment(node_id);

  let mut s = edits.apply(block.consensus())?;

  // Reverse-complement if on opposite strand
  if node.strand().is_reverse() {
    s = reverse_complement(&s)?;
  }
  Ok(s)
}
