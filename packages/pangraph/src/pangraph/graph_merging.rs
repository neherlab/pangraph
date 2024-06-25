use crate::align::alignment::Alignment;
use crate::align::alignment_args::AlignmentArgs;
use crate::align::energy::alignment_energy2;
use crate::align::minimap2::align_with_minimap2::align_with_minimap2;
use crate::align::mmseqs::align_with_mmseqs::align_with_mmseqs;
use crate::commands::build::build_args::{AlignmentBackend, PangraphBuildArgs};
use crate::o;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::PangraphBlock;
use crate::pangraph::split_matches::split_matches;
use crate::utils::id::Id;
use crate::utils::interval::{have_no_overlap, Interval};
use crate::utils::map_merge::{map_merge, ConflictResolution};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use ordered_float::OrderedFloat;
use std::borrow::Cow;
use std::collections::BTreeMap;

/// This is the function that is called when a node of the guide tree is visited.
/// Take two graphs and merge them into a single one. The merged graph is passed to the parent node.
fn merge_graphs(left_graph: &Pangraph, right_graph: &Pangraph, args: &PangraphBuildArgs) -> Result<Pangraph, Report> {
  // put the two graphs in a single one, by simply joining
  // the two sets of blocks and paths. No merging is performed
  let graph = Cow::Owned(graph_join(left_graph, right_graph));

  // iteratively try to merge homologous regions in blocks.
  // We map the consensus sequences of blocks to each other and merge
  // blocks in which we find matches. We iterate this until no more merging
  // is possible.
  loop {
    let (graph, has_changed) = self_merge(&graph, args)?;
    // stop when no more mergers are possible
    if !has_changed {
      break Ok(graph.into_owned());
    }
  }
}

fn graph_join(left_graph: &Pangraph, right_graph: &Pangraph) -> Pangraph {
  // simply join the two sets of blocks and paths
  Pangraph {
    blocks: map_merge(&left_graph.blocks, &right_graph.blocks, ConflictResolution::Left),
    paths: map_merge(&left_graph.paths, &right_graph.paths, ConflictResolution::Left),
    nodes: map_merge(&left_graph.nodes, &right_graph.nodes, ConflictResolution::Left),
  }
}

fn self_merge<'a>(graph: &'a Cow<'a, Pangraph>, args: &PangraphBuildArgs) -> Result<(Cow<'a, Pangraph>, bool), Report> {
  // use minimap2 or other aligners to find matches between the consensus
  // sequences of the blocks
  let matches = find_matches(graph.blocks.values(), args)?;

  // split matches:
  // - whenever an alignment contains an in/del longer than the threshold length
  //   (parameter - default 100 bp) we want to split the alignment in two)
  let matches = matches
    .iter()
    .map(|m| split_matches(m, &args.aln_args))
    .collect::<Result<Vec<Vec<Alignment>>, Report>>()?
    .into_iter()
    .flatten()
    .collect_vec();

  // filter matches:
  // - calculate energy and keep only matches with E < 0
  // - sort them by energy
  // - discard incompatible matches (the ones that have overlapping regions)
  let matches = filter_matches(&matches, &args.aln_args);

  // If there's no changes, then break out of the loop
  if matches.is_empty() {
    return Ok((Cow::clone(graph), false));
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
  let (mut graph, mergers) = reweave_graph(graph, &matches);

  // this can be parallelized
  let merged_blocks = mergers
    .into_iter()
    .map(|merger| perform_merger(&merger))
    .map(|block| (block.id(), block))
    .collect();

  // add the new blocks to the graph
  graph.blocks = map_merge(&graph.blocks, &merged_blocks, ConflictResolution::Right);

  // we might need this step for some final updates and consistency checks.
  // TODO: here we could also take care of transitive edges, which is useful
  // in the case of circular paths.
  let graph = consolidate(graph);

  Ok((Cow::Owned(graph), true))
}

fn find_matches<'a>(
  blocks: impl IntoIterator<Item = &'a PangraphBlock>,
  args: &PangraphBuildArgs,
) -> Result<Vec<Alignment>, Report> {
  // This function calls an aligner (default: minimap2) to find matches
  // between the consensus sequences of the blocks. It returns a list of
  // alignment objects.
  let seqs = blocks.into_iter().map(|block| block.consensus.as_str()).collect_vec();

  match args.alignment_kernel {
    AlignmentBackend::Minimap2 => align_with_minimap2(&seqs, &seqs, &args.aln_args),
    AlignmentBackend::Mmseqs => align_with_mmseqs(&seqs, &seqs, &args.aln_args),
  }
}

fn filter_matches(alns: &[Alignment], args: &AlignmentArgs) -> Vec<Alignment> {
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
    if is_match_compatible(aln, &accepted_intervals) {
      accepted_alns.push(aln.clone());
      update_intervals(aln, &mut accepted_intervals);
    }
  }

  accepted_alns
}

fn is_match_compatible(aln: &Alignment, accepted_intervals: &BTreeMap<String, Vec<Interval>>) -> bool {
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

fn update_intervals(aln: &Alignment, accepted_intervals: &mut BTreeMap<String, Vec<Interval>>) {
  accepted_intervals
    .entry(aln.reff.name.clone())
    .or_default()
    .push(aln.reff.interval.clone());

  accepted_intervals
    .entry(aln.qry.name.clone())
    .or_default()
    .push(aln.qry.interval.clone());
}

fn reweave_graph(graph: &Pangraph, alns: &[Alignment]) -> (Pangraph, Vec<Merger>) {
  // TODO: complex function. I will expand more on this, but it should:
  // - create a new graph with a copy of the paths, and only the blocks
  //   that do not undergo any merging.
  // - for each of the suggested matches:
  //   - split the block and create new nodes (including ones for mergers)
  //   - add the blocks that should not be processed further to the new graph
  //   - update the path with the new nodes
  //   - for blocks that should be merged, create a merger object and add it
  //     to the list of mergers.
  // TODO: this function should also update the position of all of the nodes!

  let graph = Pangraph {
    paths: btreemap! {},
    blocks: btreemap! {},
    nodes: btreemap! {},
  };

  let mergers = vec![];

  (graph, mergers)
}

#[derive(Clone)]
struct Merger {
  merged_block_id: String,
  deep_block: PangraphBlock,
  shallow_block: PangraphBlock,
  alignment: Alignment,
}

fn perform_merger(merger: &Merger) -> PangraphBlock {
  // TODO: this part is basically contained in the block_operations notebook.
  // in the merger we add all sequences from the shallow block to the deep
  // block, and return the updated block.

  // TODO: potentially here later we can also take care of consolidating the
  // consensus sequence, updating it with indels or substitutions that might
  // have become present in the majority of sequences.
  PangraphBlock {
    consensus: o!(""),
    alignments: btreemap! {},
  }
}

fn consolidate(graph: Pangraph) -> Pangraph {
  // TODO: final updates and consistency checks.
  // - we can take care of removing _transitive edges_, i.e.
  //   pairs of blocks that are always adjacent and connected in the same
  //   way. We have an algorithm to quickly check for this. This is useful when
  //   dealing with circular paths.
  // - we can also add optional consistency checks, to make sure that sequence
  //   lengths and path lengths are conserved, and optionally check that we can
  //   reconstruct the full genomes exactly.
  graph
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
      o!("block_0") => vec![Interval::new(100,200), Interval::new(300,400)],
      o!("block_1") => vec![Interval::new(200,300), Interval::new(400,500)],
    };

    let aln = Alignment {
      qry: Hit::new("block_0", 1000, (210, 290)),
      reff: Hit::new("block_1", 1000, (310, 390)),
      matches: 80,
      length: 80,
      quality: 10,
      orientation: Strand::Reverse,
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
      o!("block_0") => vec![Interval::new(100,200), Interval::new(300,400)],
      o!("block_1") => vec![Interval::new(200,300), Interval::new(400,500)],
    };

    let aln = Alignment {
      qry: Hit::new("block_0", 1000, (310, 390)),
      reff: Hit::new("block_1", 1000, (310, 390)),
      matches: 80,
      length: 80,
      quality: 10,
      orientation: Strand::Reverse,
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
      qry: Hit::new("bl0", 500, (100, 200)),
      reff: Hit::new("bl1", 500, (200, 300)),
      matches: 100,
      length: 0,
      quality: 0,
      orientation: Strand::default(),
      cigar: parse_cigar_str("100M").unwrap(),
      divergence: Some(0.05),
      align: None,
    };

    let aln_1 = Alignment {
      qry: Hit::new("bl2", 500, (100, 200)),
      reff: Hit::new("bl3", 500, (200, 300)),
      matches: 100,
      length: 0,
      quality: 0,
      orientation: Strand::default(),
      cigar: parse_cigar_str("100M").unwrap(),
      divergence: Some(0.02),
      align: None,
    };

    let aln_2 = Alignment {
      qry: Hit::new("bl2", 500, (150, 250)),
      reff: Hit::new("bl4", 500, (200, 300)),
      matches: 100,
      length: 0,
      quality: 0,
      orientation: Strand::default(),
      cigar: parse_cigar_str("100M").unwrap(),
      divergence: Some(0.05),
      align: None,
    };

    let aln_3 = Alignment {
      qry: Hit::new("bl5", 500, (100, 200)),
      reff: Hit::new("bl6", 500, (200, 300)),
      matches: 100,
      length: 0,
      quality: 0,
      orientation: Strand::default(),
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
