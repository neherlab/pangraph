use crate::commands::build::build_args::{AlignmentBackend, PangraphBuildArgs};
use crate::commands::reconstruct::reconstruct_run::{compare_sequences, reconstruct};
use crate::io::fasta::{read_many_fasta, FastaRecord};
use crate::io::json::{json_write_file, JsonPretty};
use crate::pangraph::graph_merging::merge_graphs;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::strand::Strand::Forward;
use crate::tree::clade::postorder;
use crate::tree::neighbor_joining::build_tree_using_neighbor_joining;
use crate::{make_internal_error, make_internal_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::info;

pub fn build_cmd_preliminary_checks(args: &PangraphBuildArgs) -> Result<(), Report> {
  // alignment kernel checks
  match args.alignment_kernel {
    AlignmentBackend::Mmseqs => {
      // check that mmseqs is available in PATH
      std::process::Command::new("mmseqs")
        .arg("--help")
        .output()
        .map_err(|_| {
          eyre::eyre!("mmseqs command not found in PATH. Please install mmseqs and make sure it is available in PATH.")
        })?;
    }
    _ => {}
  }

  Ok(())
}

pub fn build_run(args: &PangraphBuildArgs) -> Result<(), Report> {
  let input_fastas = &args.input_fastas;

  let fastas = read_many_fasta(input_fastas)?;

  // TODO: adjust fasta letter case if `upper_case` is set
  // TODO: check for duplicate fasta names

  build_cmd_preliminary_checks(args).wrap_err("When performing preliminary checks before building the pangraph.")?;

  let pangraph = if args.verify {
    let pangraph = build(fastas.clone(), args)?;
    {
      let mut results = reconstruct(&pangraph);
      results.try_for_each(|actual| -> Result<(), Report> {
        let actual = actual?;
        let expected = &fastas[actual.index];
        compare_sequences(expected, &actual);
        Ok(())
      })?;
    }
    Ok(pangraph)
  } else {
    build(fastas, args)
  }?;

  json_write_file(&args.output_json, &pangraph, JsonPretty(true))?;

  Ok(())
}

pub fn build(fastas: Vec<FastaRecord>, args: &PangraphBuildArgs) -> Result<Pangraph, Report> {
  // Build singleton graphs from input sequences
  // TODO: initial graphs can potentially be constructed when initializing tree clades. This could avoid a lot of boilerplate code.
  let graphs = fastas
    .into_iter()
    .map(|fasta| Pangraph::singleton(fasta, Forward, args.circular)) // FIXME: strand hardcoded
    .collect_vec();

  // Build guide tree
  let tree = build_tree_using_neighbor_joining(graphs)?;

  // Main loop: traverse the tree starting from leaf nodes and build the graphs bottom-up all the way to the root node.
  // The graph of the root node is the graph we are looking for.
  postorder(&tree, |clade| {
    match (&clade.left, &clade.right) {
      (None, None) => {
        // Case: leaf node. Action: nothing to do.
        Ok(())
      }
      (Some(left), Some(right)) => {
        // Case: internal node with two children. Action: produce graph for this node based on the graphs of its children.
        // Assumption: Child nodes are assumed to be already visited at this point.
        if let (Some(left), Some(right)) = (&left.read().data, &right.read().data) {
          info!(
            "=== Graph merging start:     clades sizes {} + {}",
            left.paths.len(),
            right.paths.len()
          );

          clade.data = Some(merge_graphs(left, right, args).wrap_err("When merging graphs")?);

          info!(
            "=== Graph merging completed: clades sizes {} + {} -> {}",
            left.paths.len(),
            right.paths.len(),
            clade.data.as_ref().unwrap().paths.len()
          );

          #[cfg(debug_assertions)]
          clade
            .data
            .as_ref()
            .unwrap()
            .sanity_check()
            .wrap_err("failed sanity check after merging graphs.")?;

          Ok(())
        } else {
          make_internal_error!("Found internal clade with two children, of which one or both have no graph attached.")
        }
      }
      (None, Some(child)) | (Some(child), None) => {
        // Case: internal node with one child. Action: ???
        unimplemented!("What to do if there's only one child?");
      }
    }
  })
  .into_iter()
  .collect::<Result<Vec<_>, Report>>()
  .wrap_err("When traversing guide tree")?;

  let graph = tree
    .write()
    .data
    .take()
    .ok_or_else(|| make_internal_report!("Root clade of the guide tree contains no graph after graph alignment"))?;

  Ok(graph)
}
