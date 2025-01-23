use crate::align::alignment::Alignment;
use crate::align::alignment_args::AlignmentArgs;
use crate::align::energy::alignment_energy2;
use crate::align::minimap2_lib::align_with_minimap2_lib::align_with_minimap2_lib;
use crate::align::mmseqs::align_with_mmseqs::align_with_mmseqs;
use crate::circularize::circularize::remove_transitive_edges;
use crate::commands::build::build_args::{AlignmentBackend, PangraphBuildArgs};
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::reweave::reweave;
use crate::pangraph::split_matches::split_matches;
use crate::reconsensus::reconsensus::reconsensus_graph;
use crate::utils::interval::{have_no_overlap, Interval};
use crate::utils::map_merge::{map_merge, ConflictResolution};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::{debug, trace, warn};
use maplit::btreemap;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use std::collections::BTreeMap;

/// This is the function that is called when a node of the guide tree is visited.
/// Take two graphs and merge them into a single one. The merged graph is passed to the parent node.
pub fn merge_graphs(
  left_graph: &Pangraph,
  right_graph: &Pangraph,
  args: &PangraphBuildArgs,
) -> Result<Pangraph, Report> {
  // put the two graphs in a single one, by simply joining
  // the two sets of blocks and paths. No merging is performed
  let mut graph = graph_join(left_graph, right_graph);

  // iteratively try to merge homologous regions in blocks.
  // We map the consensus sequences of blocks to each other and merge
  // blocks in which we find matches. We iterate this until no more merging
  // is possible.
  let mut i = 0;
  loop {
    // print keys of left/right graph paths
    let left_keys = left_graph.paths.keys().map(|k| k.to_string()).join(", ");
    let right_keys = right_graph.paths.keys().map(|k| k.to_string()).join(", ");
    debug!("Self-merge iteration {i} of left/right graph {left_keys} <---> {right_keys}");

    let (graph_new, has_changed) =
      self_merge(graph, args).wrap_err_with(|| format!("During self-merge iteration {i}"))?;
    graph = graph_new;

    // stop when no more mergers are possible
    // or when the maximum number of iterations is reached
    if !has_changed {
      debug!("Graph merge {left_keys} <---> {right_keys} complete.");
      break;
    } else if i >= args.max_self_map {
      warn!(
        "Reached maximum number of self-merge iterations at graph merging {left_keys} <---> {right_keys}, consider increasing the current limit -x {}",
        args.max_self_map
      );
      break;
    }
    i += 1;
  }

  debug!("Removing transitive edges");
  remove_transitive_edges(&mut graph).wrap_err("When removing transitive edges")?;

  #[cfg(any(test, debug_assertions))]
  graph.sanity_check().wrap_err("After merging graphs")?;

  Ok(graph)
}

pub fn graph_join(left_graph: &Pangraph, right_graph: &Pangraph) -> Pangraph {
  // simply join the two sets of blocks and paths
  Pangraph {
    blocks: map_merge(
      &left_graph.blocks,
      &right_graph.blocks,
      ConflictResolution::Custom(|(kl, vl), (kr, vr)| panic!("Conflicting key: '{kl}'")),
    ),
    paths: map_merge(
      &left_graph.paths,
      &right_graph.paths,
      ConflictResolution::Custom(|(kl, vl), (kr, vr)| panic!("Conflicting key: '{kl}'")),
    ),
    nodes: map_merge(
      &left_graph.nodes,
      &right_graph.nodes,
      ConflictResolution::Custom(|(kl, vl), (kr, vr)| panic!("Conflicting key: '{kl}'")),
    ),
  }
}

pub fn self_merge(graph: Pangraph, args: &PangraphBuildArgs) -> Result<(Pangraph, bool), Report> {
  // use minimap2 or other aligners to find matches between the consensus
  // sequences of the blocks
  let matches = find_matches(&graph.blocks, args)?;
  debug!("Found matches: {}", matches.len());
  trace!("{matches:#?}");

  // exclude self-alignments (a block with itself)
  // and split matches:
  // - whenever an alignment contains an in/del longer than the threshold length
  //   (parameter - default 100 bp) we want to split the alignment in two)
  let matches = matches
    .iter()
    .filter(|m| m.qry.name != m.reff.name)
    .map(|m| split_matches(m, &args.aln_args))
    .collect::<Result<Vec<Vec<Alignment>>, Report>>()?
    .into_iter()
    .flatten()
    .collect_vec();
  debug!("Matches after splitting: {}", matches.len());
  trace!("{matches:#?}");

  // filter matches:
  // - calculate energy and keep only matches with E < 0
  // - sort them by energy
  // - discard incompatible matches (the ones that have overlapping regions)
  let mut matches = filter_matches(&matches, &args.aln_args);
  debug!("Matches after filtering: {}", matches.len());
  trace!("{matches:#?}");

  // If there's no changes, then break out of the loop
  if matches.is_empty() {
    return Ok((graph, false));
  }

  // complex function: takes the list of desired matches and the two
  // graphs. It splits blocks and re-weaves paths through them. Paths
  // will not be updated after this.
  // The function does not perform any merging, but returns a list of
  // blocks that should be merged.
  // Nb: this function should already take care of:
  // - creating new nodes
  // - substituting the old nodes with the new ones in paths
  // - splitting blocks and updating the node ids.
  // - adding the blocks that do not need merging to the preliminary graph
  // - return the set of blocks that should be merged
  let (mut graph, mergers) =
    reweave(&mut matches, graph, args.aln_args.indel_len_threshold).wrap_err("During reweave")?;

  let merged_blocks = mergers
    .into_par_iter()
    .map(|mut merge_promise| {
      merge_promise
        .solve_promise()
        .wrap_err_with(|| format!("When solving merge promise: {merge_promise:#?}"))
    })
    .collect::<Result<Vec<_>, Report>>()?
    .into_iter()
    .map(|block| (block.id(), block))
    .collect::<BTreeMap<_, _>>();

  // add the new blocks to the graph
  graph.blocks = map_merge(
    &graph.blocks,
    &merged_blocks,
    ConflictResolution::Custom(|(k1, _), (k2, _)| panic!("Conflicting key '{k1}'")),
  );

  // update consensus and alignment of merged blocks.
  let merge_block_ids = merged_blocks.keys().copied().collect_vec();
  reconsensus_graph(&mut graph, merge_block_ids).wrap_err("During reconsensus")?;

  Ok((graph, true))
}

// Call an aligner (default: minimap2) to find matches between the consensus sequences of the blocks.
// Returns a list of alignment objects.
pub fn find_matches(
  blocks: &BTreeMap<BlockId, PangraphBlock>,
  args: &PangraphBuildArgs,
) -> Result<Vec<Alignment>, Report> {
  match args.alignment_kernel {
    AlignmentBackend::Minimap2 => align_with_minimap2_lib(blocks, &args.aln_args),
    AlignmentBackend::Mmseqs => align_with_mmseqs(blocks, &args.aln_args),
  }
  .wrap_err_with(|| format!("When trying to align sequences using {}", &args.alignment_kernel))
}

pub fn filter_matches(alns: &[Alignment], args: &AlignmentArgs) -> Vec<Alignment> {
  // - evaluates the energy of the alignments
  // - keeps only matches with E < 0
  // - sorts them by energy
  // - discards incompatible matches (the ones that have overlapping regions)

  // TODO: energy is calculated for each alignment.
  // Consider calculating it earlier and making it a property to simplify filtering and sorting.
  let alns = alns
    .iter()
    .map(|aln| (aln, alignment_energy2(aln, args)))
    .filter(|(_, energy)| energy < &0.0)
    .sorted_by_key(|(_, energy)| OrderedFloat(*energy))
    .map(|(aln, _)| aln)
    .collect_vec();

  // Iteratively accept alignments if they do not overlap with previously accepted ones
  let mut accepted_alns = vec![];
  let mut accepted_intervals = btreemap![];

  for aln in alns {
    debug_assert!(aln.qry.name != aln.reff.name);
    if is_match_compatible(aln, &accepted_intervals) {
      accepted_alns.push(aln.clone());
      update_intervals(aln, &mut accepted_intervals);
    }
  }

  accepted_alns
}

pub fn is_match_compatible(aln: &Alignment, accepted_intervals: &BTreeMap<BlockId, Vec<Interval>>) -> bool {
  let ref_compatible = have_no_overlap(
    accepted_intervals.get(&aln.reff.name).unwrap_or(&vec![]),
    &aln.reff.interval,
  );

  let qry_compatible = have_no_overlap(
    accepted_intervals.get(&aln.qry.name).unwrap_or(&vec![]),
    &aln.qry.interval,
  );

  ref_compatible && qry_compatible
}

pub fn update_intervals(aln: &Alignment, accepted_intervals: &mut BTreeMap<BlockId, Vec<Interval>>) {
  accepted_intervals
    .entry(aln.reff.name)
    .or_default()
    .push(aln.reff.interval.clone());

  accepted_intervals
    .entry(aln.qry.name)
    .or_default()
    .push(aln.qry.interval.clone());
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::align::alignment::Hit;
  use crate::align::bam::cigar::parse_cigar_str;
  use crate::pangraph::strand::Strand;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_is_match_compatible() -> Result<(), Report> {
    let accepted_intervals = btreemap! {
      BlockId(0) => vec![Interval::new(100,200), Interval::new(300,400)],
      BlockId(1) => vec![Interval::new(200,300), Interval::new(400,500)],
    };

    let aln = Alignment {
      qry: Hit::new(BlockId(0), 1000, (210, 290)),
      reff: Hit::new(BlockId(1), 1000, (310, 390)),
      matches: 80,
      length: 80,
      quality: 10,
      orientation: Strand::Reverse,
      new_block_id: None, // FIXME
      anchor_block: None, // FIXME
      cigar: parse_cigar_str("90M").unwrap(),
      divergence: Some(0.05),
      align: None,
    };

    assert!(is_match_compatible(&aln, &accepted_intervals));

    Ok(())
  }

  #[rstest]
  fn test_is_match_compatible_not() -> Result<(), Report> {
    let accepted_intervals = btreemap! {
      BlockId(0) => vec![Interval::new(100,200), Interval::new(300,400)],
      BlockId(1) => vec![Interval::new(200,300), Interval::new(400,500)],
    };

    let aln = Alignment {
      qry: Hit::new(BlockId(0), 1000, (310, 390)),
      reff: Hit::new(BlockId(1), 1000, (310, 390)),
      matches: 80,
      length: 80,
      quality: 10,
      orientation: Strand::Reverse,
      new_block_id: None, // FIXME
      anchor_block: None, // FIXME
      cigar: parse_cigar_str("90M").unwrap(),
      divergence: Some(0.05),
      align: None,
    };

    assert!(!is_match_compatible(&aln, &accepted_intervals));

    Ok(())
  }

  #[rstest]
  fn test_filter_matches() -> Result<(), Report> {
    let aln_0 = Alignment {
      qry: Hit::new(BlockId(0), 500, (100, 200)),
      reff: Hit::new(BlockId(1), 500, (200, 300)),
      matches: 100,
      length: 0,
      quality: 0,
      orientation: Strand::default(),
      new_block_id: None, // FIXME
      anchor_block: None, // FIXME
      cigar: parse_cigar_str("100M").unwrap(),
      divergence: Some(0.05),
      align: None,
    };

    let aln_1 = Alignment {
      qry: Hit::new(BlockId(2), 500, (100, 200)),
      reff: Hit::new(BlockId(3), 500, (200, 300)),
      matches: 100,
      length: 0,
      quality: 0,
      orientation: Strand::default(),
      new_block_id: None, // FIXME
      anchor_block: None, // FIXME
      cigar: parse_cigar_str("100M").unwrap(),
      divergence: Some(0.02),
      align: None,
    };

    let aln_2 = Alignment {
      qry: Hit::new(BlockId(2), 500, (150, 250)),
      reff: Hit::new(BlockId(4), 500, (200, 300)),
      matches: 100,
      length: 0,
      quality: 0,
      orientation: Strand::default(),
      new_block_id: None, // FIXME
      anchor_block: None, // FIXME
      cigar: parse_cigar_str("100M").unwrap(),
      divergence: Some(0.05),
      align: None,
    };

    let aln_3 = Alignment {
      qry: Hit::new(BlockId(5), 500, (100, 200)),
      reff: Hit::new(BlockId(6), 500, (200, 300)),
      matches: 100,
      length: 0,
      quality: 0,
      orientation: Strand::default(),
      new_block_id: None, // FIXME
      anchor_block: None, // FIXME
      cigar: parse_cigar_str("100M").unwrap(),
      divergence: Some(0.1),
      align: None,
    };

    let args = AlignmentArgs {
      alpha: 10.0,
      beta: 10.0,
      ..Default::default()
    };

    let alns = [aln_0.clone(), aln_1.clone(), aln_2, aln_3];

    assert_eq!(filter_matches(&alns, &args), vec![aln_1, aln_0]);

    Ok(())
  }
}
