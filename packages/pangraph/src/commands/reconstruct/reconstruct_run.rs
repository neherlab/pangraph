use crate::commands::reconstruct::reconstruct_args::PangraphReconstructArgs;
use crate::io::fasta::FastaRecord;
use crate::io::seq::reverse_complement;
use crate::pangraph::pangraph::Pangraph;
use crate::{make_internal_error, make_internal_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;

pub fn reconstruct_run(args: &PangraphReconstructArgs) -> Result<(), Report> {
  let PangraphReconstructArgs { input_graph, .. } = &args;

  // Load input graph from json file
  // let graph = ...

  // Reconstruct sequences from the graph
  // let fasta_output = reconstruct(graph, &args)?;

  // write fasta to stdout
  // println!("{}", fasta_output);

  Ok(())
}

// this function could belong to the path class.
// fn reconstruct_path_sequence(path: &PangraphPath, graph: &Pangraph) -> Result<String, Report> {
//   let mut seq = Vec::new();
//   for node_id in &path.nodes {
//     let node = graph
//       .nodes
//       .get(node_id)
//       .ok_or_else(|| make_internal_report!("Node {} not found in graph", node_id))?;
//     let block_id = n.block_id();
//     let block = graph
//       .blocks
//       .get(&block_id)
//       .ok_or_else(|| make_internal_report!("Block {} not found in graph", block_id))?;
//     let edits = block.alignment(node_id);
//     let mut s = edits.apply(block.consensus())?;
//     if node.strand() == Strand::Reverse {
//       // TODO: reverse complement
//       s = reverse_complement(s)
//     }
//     seq.push(s);
//   }
//   Ok(seq.join(""))
// }

// pub fn reconstruct(graph: &Pangraph, args: &PangraphReconstructArgs) -> Result<Vec<FastaRecord>, Report> {
//   let seqs = graph
//     .paths
//     .iter()
//     .enumerate()
//     .map(|(i, p)| (i, reconstruct_path_sequence(p, graph), p.name.clone()))
//     .map(|(i, s, n)| FastaRecord {
//       index: i,
//       seq: s,
//       seq_name: n,
//     });
//   Ok(seqs)
// }
