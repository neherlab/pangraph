use crate::align::alignment::Alignment;
use crate::align::minimap2::align_with_minimap2::{align_with_minimap2, Minimap2Params};
use crate::align::mmseqs::align_with_mmseqs::{align_with_mmseqs, MmseqsParams};
use crate::commands::build::build_args::{AlignmentBackend, PangraphBuildArgs};
use crate::o;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::PangraphBlock;
use crate::pangraph::split_matches::split_matches;
use eyre::Report;
use itertools::{chain, Itertools};
use maplit::btreemap;
use ordered_float::OrderedFloat;

/// This is the function that is called when a node of the guide tree is visited.
/// Take two graphs and merge them into a single one. The merged graph is passed to the parent node.
fn merge_graphs(left_graph: &Pangraph, right_graph: &Pangraph, args: &PangraphBuildArgs) -> Result<Pangraph, Report> {
  // put the two graphs in a single one, by simply joining
  // the two sets of blocks and paths. No merging is performed
  let graph = graph_join(left_graph, right_graph);

  // iteratively try to merge homologous regions in blocks.
  // We map the consensus sequences of blocks to each other and merge
  // blocks in which we find matches. We iterate this until no more merging
  // is possible.
  loop {
    let (graph, has_changed) = self_merge(&graph, args)?;
    // stop when no more mergers are possible
    if !has_changed {
      break Ok(graph);
    }
  }
}

fn graph_join(left_graph: &Pangraph, right_graph: &Pangraph) -> Pangraph {
  // simply join the two sets of blocks and paths
  Pangraph {
    // TODO: restructure code to avoid copying
    blocks: chain!(left_graph.blocks.iter().cloned(), right_graph.blocks.iter().cloned()).collect_vec(),
    paths: chain!(left_graph.paths.iter().cloned(), right_graph.paths.iter().cloned()).collect_vec(),
  }
}

fn self_merge(graph: &Pangraph, args: &PangraphBuildArgs) -> Result<(Pangraph, bool), Report> {
  // use minimap2 or other aligners to find matches between the consensus
  // sequences of the blocks
  let matches = find_matches(&graph.blocks, args);

  // split matches:
  // - whenever an alignment contains an in/del longer than the threshold length
  //   (parameter - default 100 bp) we want to split the alignment in two)
  let matches = matches
    .iter()
    .flatten()
    .map(|m| split_matches(m, &args.split_matches_args))
    .collect::<Result<Vec<Vec<Alignment>>, Report>>()?
    .into_iter()
    .flatten()
    .collect_vec();

  // filter matches:
  // - calculate energy and keep only matches with E < 0
  // - sort them by energy
  // - discard incompatible matches (the ones that have overlapping regions)
  let matches = filter_matches(&matches);

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
  let merged_blocks = mergers.into_iter().map(|merger| perform_merger(&merger)).collect_vec();

  // add the new blocks to the graph
  graph.blocks.extend(merged_blocks);

  // we might need this step for some final updates and consistency checks.
  // TODO: here we could also take care of transitive edges, which is useful
  // in the case of circular paths.
  Ok((consolidate(graph), true))
}

fn find_matches(blocks: &[PangraphBlock], args: &PangraphBuildArgs) -> Result<Vec<Alignment>, Report> {
  // This function calls an aligner (default: minimap2) to find matches
  // between the consensus sequences of the blocks. It returns a list of
  // alignment objects.
  let seqs = blocks.iter().map(|block| block.consensus.as_str()).collect_vec();

  match args.alignment_kernel {
    AlignmentBackend::Minimap2 => align_with_minimap2(
      &seqs,
      &seqs,
      &Minimap2Params {
        kmersize: args.kmer_length,
        preset: Some(format!("asm{}", args.sensitivity)),
      },
    ),
    AlignmentBackend::Mmseqs => align_with_mmseqs(
      &seqs,
      &seqs,
      &MmseqsParams {
        kmer_len: args.kmer_length,
      },
    ),
  }
}

fn filter_matches(alns: &[Alignment]) -> Vec<Alignment> {
  // - evaluates the energy of the alignments
  // - keeps only matches with E < 0
  // - sorts them by energy
  // - discards incompatible matches (the ones that have overlapping regions)

  // TODO: energy is calculated for each alignment.
  // Consider calculating it earlier and making it a property to simplify filtering and sorting.
  let alns = alns
    .iter()
    .map(|aln| (aln, alignment_energy(aln)))
    .filter(|(_, energy)| energy < &0.0)
    .sorted_by_key(|(_, energy)| OrderedFloat(*energy))
    .map(|(aln, _)| aln)
    .collect_vec();

  // discard incompatible matches
  let mut alns_keep = vec![];
  for aln in alns {
    if is_match_compatible(aln, &alns_keep) {
      alns_keep.push(aln.clone());
    }
  }

  alns_keep
}

fn alignment_energy(aln: &Alignment) -> f64 {
  // TODO: calculate the energy of the alignment
  // this is a function of alignment length, identity, and whether it will
  // create a block split.
  0.0
}

fn is_match_compatible(aln: &Alignment, alns_keep: &[Alignment]) -> bool {
  // TODO: check that the new match does not overlap with any of the
  // matches that we already have.
  true
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
    paths: vec![],
    blocks: vec![],
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
    id: 0,
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
