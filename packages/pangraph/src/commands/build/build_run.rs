use crate::commands::build::build_args::{AlignmentBackend, PangraphBuildArgs};
use crate::commands::reconstruct::reconstruct_run::{compare_sequences, reconstruct};
use crate::io::fasta::{FastaRecord, read_many_fasta};
use crate::io::json::{JsonPretty, json_write_file};
use crate::pangraph::graph_merging::merge_graphs;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::strand::Strand::Forward;
use crate::tree::clade::postorder;
use crate::tree::neighbor_joining::build_tree_using_neighbor_joining;
use crate::utils::progress_bar::ProgressBar;
use crate::{make_error, make_internal_error, make_internal_report};
use color_eyre::owo_colors::{AnsiColors, OwoColorize};
use color_eyre::{Help, SectionExt};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::info;

pub fn build_cmd_preliminary_checks(args: &PangraphBuildArgs) -> Result<(), Report> {
  // alignment kernel checks
  if args.alignment_kernel == AlignmentBackend::Mmseqs {
    // check that mmseqs is available in PATH
    std::process::Command::new("mmseqs")
      .arg("--help")
      .output()
      .wrap_err("When executing `mmseqs --help`")
      .section(
        "Please make sure that `mmseqs` is installed, available in PATH and is functional outside of pangraph. For more details, refer to mmseqs documentation at https://github.com/soedinglab/MMseqs2"
          .color(AnsiColors::Cyan)
          .header("Suggestion:"),
      )?;
  }

  Ok(())
}

pub fn reconstruct_and_compare_graph_seqs(graph: &Pangraph, fastas: &[FastaRecord]) -> Result<(), Report> {
  // Reconstruct sequences from the given graph.
  let mut results = reconstruct(graph);

  // Check that the reconstructed sequences match the original FASTA records.
  results.try_for_each(|actual| -> Result<(), Report> {
    let actual = actual?;
    let expected = &fastas[actual.index];
    compare_sequences(expected, &actual)?;
    Ok(())
  })?;

  Ok(())
}

pub fn graph_sanity_checks(graph: &Pangraph, fastas: &[FastaRecord]) -> Result<(), Report> {
  // check that graph internal structure (blocks, paths, nodes, edits...) is valid
  #[cfg(debug_assertions)]
  graph
    .sanity_check()
    .wrap_err("When performing sanity check on the pangraph")?;

  // Reconstruct sequences from the graph and compare them with the original FASTA records.
  reconstruct_and_compare_graph_seqs(graph, fastas)
    .wrap_err("When comparing reconstructed sequences with original FASTA records")?;

  Ok(())
}

pub fn build_run(args: &PangraphBuildArgs) -> Result<(), Report> {
  let input_fastas = &args.input_fastas;

  let fastas = read_many_fasta(input_fastas)?;

  // Ensure that no input sequences contain gap characters
  ensure_no_gap_characters(&fastas)?;

  // TODO: adjust fasta letter case if `upper_case` is set
  // TODO: check for duplicate fasta names

  build_cmd_preliminary_checks(args).wrap_err("When performing preliminary checks before building the pangraph.")?;

  let pangraph = build(fastas, args, args.verify)?;

  json_write_file(&args.output_json, &pangraph, JsonPretty(true))?;

  Ok(())
}

/// Ensure that no input sequences contain gap characters (`-`).
/// If any sequences with gaps are found, an error is returned listing the offending sequences.
fn ensure_no_gap_characters(fastas: &[FastaRecord]) -> Result<(), Report> {
  let sequences_with_gaps: Vec<&str> = fastas
    .iter()
    .filter_map(|record| record.seq.contains_str("-").then_some(record.seq_name.as_str()))
    .collect();

  if sequences_with_gaps.is_empty() {
    return Ok(());
  }

  let names_listing = sequences_with_gaps.iter().join("\n");

  make_error!(
    "Detected gap character '-' in {} input sequences: {}. Pangraph does not support gap characters in the input.",
    sequences_with_gaps.len(),
    names_listing
  )
  .with_section(|| names_listing.header("Sequences containing gaps:"))
  .with_section(|| {
    "Remove gap characters from the input FASTA before running `pangraph build`."
      .color(AnsiColors::Cyan)
      .header("Suggestion:")
  })
}

pub fn build(fastas: Vec<FastaRecord>, args: &PangraphBuildArgs, verify: bool) -> Result<Pangraph, Report> {
  // If verification is requested, we need to keep a copy of the original FASTA records
  // to compare them with the sequences reconstructed from the graph.
  let fasta_copy = verify.then(|| fastas.clone());

  // Build singleton graphs from input sequences
  // TODO: initial graphs can potentially be constructed when initializing tree clades. This could avoid a lot of boilerplate code.
  let n_paths = fastas.len();
  let graphs = fastas
    .into_iter()
    .map(|fasta| Pangraph::singleton(fasta, Forward, args.circular)) // FIXME: strand hardcoded
    .collect_vec();

  // Build guide tree
  let tree = build_tree_using_neighbor_joining(graphs)?;

  // Instantiate the progress bar
  let pb = ProgressBar::new(n_paths - 1, args.no_progress_bar)?;

  // Main loop: traverse the tree starting from leaf nodes and build the graphs bottom-up all the way to the root node.
  // The graph of the root node is the graph we are looking for.
  postorder(&tree, |clade| {
    match (&clade.left, &clade.right) {
      (None, None) => {
        // Case: leaf node. Action: nothing to do.
        Ok(())
      },
      (Some(left), Some(right)) => {
        // Case: internal node with two children. Action: produce graph for this node based on the graphs of its children.
        // Assumption: Child nodes are assumed to be already visited at this point.
        match (&left.read().data, &right.read().data) {
          (Some(left), Some(right)) => {
            info!(
              "=== Graph merging start:     clades sizes {} + {}",
              left.paths.len(),
              right.paths.len()
            );

            clade.data = Some(merge_graphs(left, right, args).wrap_err("When merging graphs")?);

            // increase progress bar
            pb.inc(1);

            info!(
              "=== Graph merging completed: clades sizes {} + {} -> {}",
              left.paths.len(),
              right.paths.len(),
              clade.data.as_ref().unwrap().paths.len()
            );

            // perform checks only in debug mode and if requested
            #[cfg(debug_assertions)]
            {
              if verify {
                // verify the graph if requested
                graph_sanity_checks(clade.data.as_ref().unwrap(), fasta_copy.as_ref().unwrap())
                  .wrap_err("When performing sanity checks on the merged graph")?;
              }
            }

            Ok(())
          },
          _ => {
            make_internal_error!("Found internal clade with two children, of which one or both have no graph attached.")
          },
        }
      },
      (None, Some(_)) | (Some(_), None) => {
        // Case: internal node with one child. Action: ???
        unimplemented!("What to do if there's only one child?");
      },
    }
  })
  .wrap_err("When traversing guide tree")?;

  // Finish progress bar
  pb.finish_with_message("Graph merging completed");

  let graph = tree
    .write()
    .data
    .take()
    .ok_or_else(|| make_internal_report!("Root clade of the guide tree contains no graph after graph alignment"))?;

  // verify the final graph if requested
  if verify {
    graph_sanity_checks(&graph, fasta_copy.as_ref().unwrap())
      .wrap_err("When performing sanity checks on the final pangraph")?;
  }

  Ok(graph)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::representation::seq::Seq;
  use crate::utils::error::report_to_string;

  #[test]
  fn ensure_no_gap_characters_accepts_clean_sequences() {
    let mut record = FastaRecord::default();
    record.seq_name = "clean_seq".to_string();
    record.seq = Seq::from_str("AACCGGTT");

    assert!(ensure_no_gap_characters(std::slice::from_ref(&record)).is_ok());
  }

  #[test]
  fn ensure_no_gap_characters_rejects_sequences_with_gaps() {
    let mut clean = FastaRecord::default();
    clean.seq_name = "clean".to_string();
    clean.seq = Seq::from_str("ACGT");

    let mut gappy = FastaRecord::default();
    gappy.seq_name = "gappy".to_string();
    gappy.seq = Seq::from_str("AC-GT");

    let fastas = vec![clean, gappy];

    let err = ensure_no_gap_characters(&fastas).unwrap_err();
    let err_string = report_to_string(&err);

    assert!(err_string.contains("Detected gap character '-'"));
    assert!(err_string.contains("gappy"));
  }
}
