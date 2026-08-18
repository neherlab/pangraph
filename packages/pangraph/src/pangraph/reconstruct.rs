use crate::io::fasta::FastaRecord;
use crate::io::seq::reverse_complement;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_node::NodeId;
use crate::pangraph::pangraph_path::{PangraphPath, PathId};
use crate::representation::seq::Seq;
use crate::utils::collections::find_duplicates;
use crate::utils::string::str_slice_safe;
use crate::{make_error, make_internal_error, make_internal_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};

/// Number of genome names listed in full in an error message before the rest are elided.
const MAX_NAMES_IN_ERROR: usize = 10;

/// Number of nucleotides shown on each side of the first difference between two genomes.
const MISMATCH_CONTEXT: usize = 10;

/// How much of the expected genome set the graph is required to reconstruct.
///
/// In both cases a genome that the graph contains but `expected` does not is an error: the
/// variants differ only in whether the graph is allowed to be missing expected genomes.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum GenomeCoverage {
  /// The graph must reconstruct every expected genome, and nothing else.
  Complete,
  /// The graph may reconstruct only some of the expected genomes. Used for the intermediate graphs
  /// of a build, which hold the genomes of one clade of the guide tree.
  Partial,
}

/// Reconstructs every genome of the graph as a FASTA record, ordered by path id.
///
/// The ordering comes for free: `paths` is a `BTreeMap` keyed by [`PathId`], so iterating it
/// already yields ascending ids and the records stay lazy.
///
/// Record order and `index` reproduce the order of the original input FASTA only for graphs
/// produced directly by `build`. A merged graph renumbers its path ids, so consumers must match
/// records by genome name.
pub fn reconstruct(graph: &Pangraph) -> impl Iterator<Item = Result<FastaRecord, Report>> + use<'_> {
  graph.paths.iter().map(|(path_id, path)| {
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

/// Checks that every path of `graphs` carries a non-empty name, and that no name occurs more than
/// once across all `graphs` taken together.
///
/// The genome name is the only identifier that survives a merge unchanged: it is the key that
/// verification matches on, and what `simplify` resolves genomes by. A graph whose genomes cannot
/// be told apart by name is therefore rejected rather than processed. A name that is empty, or made
/// of whitespace only, identifies a genome no better than a missing one, and is reported the same
/// way. A name padded with whitespace is rejected for the same reason: every comparison pangraph
/// makes on names is exact, so " a" would be a genome that cannot be addressed by the name it
/// appears to have, and that does not collide with the "a" of another graph. Names are only tested
/// here, never trimmed: they are stored exactly as given.
///
/// This is the single implementation of that invariant for graphs; [`check_sequence_names`] is its
/// counterpart for FASTA records. Note that it inspects the paths directly rather than going
/// through [`reconstruct`], which masks unnamed paths behind a placeholder name.
///
/// Path ids are only unique within one graph, so when several graphs are passed the reported ids
/// are ambiguous on their own. Callers that pass more than one graph are expected to name the
/// graphs in an error section.
pub fn check_genome_names(graphs: &[&Pangraph]) -> Result<(), Report> {
  let unnamed = graphs
    .iter()
    .flat_map(|graph| graph.paths.iter())
    .filter(|(_, path)| path.name.as_deref().is_none_or(|name| name.trim().is_empty()))
    .map(|(path_id, _)| path_id.to_string())
    .collect_vec();

  if !unnamed.is_empty() {
    return make_error!(
      "Found {} genome(s) with no name or an empty name (path ids: {}). Genomes are identified by name, so every path must be named.",
      unnamed.len(),
      format_names(&unnamed)
    );
  }

  let padded = graphs
    .iter()
    .flat_map(|graph| graph.path_names().flatten())
    .filter(|name| name.trim() != *name)
    .map(|name| format!("{name:?}"))
    .collect_vec();

  if !padded.is_empty() {
    return make_error!(
      "Found {} genome name(s) with leading or trailing whitespace: {}. Genome names are compared exactly, so surrounding whitespace makes a genome impossible to tell apart from, and impossible to address as, the name it appears to have.",
      padded.len(),
      format_names(&padded)
    );
  }

  let duplicates = find_duplicates(graphs.iter().flat_map(|graph| graph.path_names().flatten()));
  if !duplicates.is_empty() {
    return make_error!(
      "Duplicate genome names found: {}. Genome names must be unique, because they identify genomes.",
      format_names(&duplicates)
    );
  }

  Ok(())
}

/// Maps each genome name of the graph to the id of the path holding it.
///
/// Errors if any path is unnamed or if two paths share a name, via [`check_genome_names`].
pub fn path_ids_by_name(graph: &Pangraph) -> Result<BTreeMap<&str, PathId>, Report> {
  check_genome_names(&[graph])?;

  Ok(
    graph
      .paths
      .iter()
      .filter_map(|(path_id, path)| path.name.as_deref().map(|name| (name, *path_id)))
      .collect(),
  )
}

/// Reconstructs the genome of a single path.
pub fn reconstruct_genome(graph: &Pangraph, path_id: PathId) -> Result<Seq, Report> {
  let path = graph
    .paths
    .get(&path_id)
    .ok_or_else(|| make_internal_report!("Path {path_id} not found in graph"))?;
  reconstruct_path_sequence(graph, path)
}

/// Checks that every FASTA record carries a non-empty name, and that no two records share one.
///
/// The FASTA counterpart of [`check_genome_names`]: sequence names become genome names when the
/// records are built into a graph, so `build` enforces the invariant on its input to guarantee that
/// a graph can never carry unusable genome names into a later merge.
///
/// Empty names are reported before duplicates, so that a file of headers that are all empty is
/// reported for the reason it actually has. Names padded with whitespace are rejected too, so that
/// two records whose headers look alike cannot become two genomes.
pub fn check_sequence_names(fastas: &[FastaRecord]) -> Result<(), Report> {
  let empty = fastas
    .iter()
    .filter(|fasta| fasta.seq_name.trim().is_empty())
    .map(|fasta| fasta.index.to_string())
    .collect_vec();

  if !empty.is_empty() {
    return make_error!(
      "Found {} input sequence(s) with an empty name (record indices: {}). Sequences are identified by name, so every record must have one. Note that a space between '>' and the identifier makes the identifier part of the description rather than the name.",
      empty.len(),
      format_names(&empty)
    );
  }

  let padded = fastas
    .iter()
    .filter(|fasta| fasta.seq_name.trim() != fasta.seq_name)
    .map(|fasta| format!("{:?}", fasta.seq_name))
    .collect_vec();

  if !padded.is_empty() {
    return make_error!(
      "Found {} input sequence(s) whose name has leading or trailing whitespace: {}. Sequence names are compared exactly, so surrounding whitespace makes a genome impossible to tell apart from the name it appears to have. Note that whitespace between '>' and the identifier becomes part of the name, unless it is a plain space, which makes the identifier part of the description instead.",
      padded.len(),
      format_names(&padded)
    );
  }

  let duplicates = find_duplicates(fastas.iter().map(|fasta| fasta.seq_name.as_str()));
  if !duplicates.is_empty() {
    return make_error!(
      "Duplicate sequence names found: {}. Sequence names must be unique, because they identify genomes.",
      format_names(&duplicates)
    );
  }

  Ok(())
}

/// Collects FASTA records into genome sequences keyed by name.
///
/// Names are assumed to be unique already: pass the records through [`check_sequence_names`]
/// first, or records sharing a name will silently collapse into one entry.
pub(crate) fn sequences_by_name(fastas: &[FastaRecord]) -> BTreeMap<String, Seq> {
  fastas
    .iter()
    .map(|fasta| (fasta.seq_name.clone(), fasta.seq.clone()))
    .collect()
}

/// Compares one reconstructed genome against the sequence it is expected to have, reporting the
/// genome name and the first position at which the two differ.
pub fn verify_genome(name: &str, expected: &Seq, actual: &Seq) -> Result<(), Report> {
  if expected == actual {
    return Ok(());
  }

  if expected.len() != actual.len() {
    return make_error!(
      "Sequence mismatch for genome '{name}': expected length {} but got {}",
      expected.len(),
      actual.len()
    );
  }

  let pos = expected
    .iter()
    .zip(actual.iter())
    .position(|(e, a)| e != a)
    .ok_or_else(|| make_internal_report!("Genomes of '{name}' compare as different but share every character"))?;

  let (start, end) = (pos.saturating_sub(MISMATCH_CONTEXT), pos + MISMATCH_CONTEXT + 1);
  make_error!(
    "Sequence mismatch for genome '{name}' at position {pos} (length {}):\n  expected: {}\n  actual:   {}",
    expected.len(),
    str_slice_safe(expected.as_str(), start, end),
    str_slice_safe(actual.as_str(), start, end)
  )
}

/// Checks that the graph reconstructs the expected genomes, matched by name.
///
/// Genomes are never matched by position or by path id: a merge renumbers path ids, so neither the
/// order in which genomes are reconstructed nor [`FastaRecord::index`] can pair them up. Only the
/// sequence content is compared; descriptions are ignored.
///
/// Genomes are reconstructed and dropped one at a time, so this holds no more than a single genome
/// in memory beyond `expected`.
pub fn verify_graph_sequences(
  graph: &Pangraph,
  expected: &BTreeMap<String, Seq>,
  coverage: GenomeCoverage,
) -> Result<(), Report> {
  let path_ids = path_ids_by_name(graph)?;

  for (name, path_id) in &path_ids {
    let Some(expected_seq) = expected.get(*name) else {
      return make_error!("Graph contains genome '{name}', which is not among the expected genomes");
    };
    let actual =
      reconstruct_genome(graph, *path_id).wrap_err_with(|| format!("When reconstructing genome '{name}'"))?;
    verify_genome(name, expected_seq, &actual)?;
  }

  if coverage == GenomeCoverage::Complete {
    let missing = expected
      .keys()
      .filter(|name| !path_ids.contains_key(name.as_str()))
      .cloned()
      .collect_vec();

    if !missing.is_empty() {
      return make_error!(
        "Graph is missing {} expected genome(s): {}",
        missing.len(),
        format_names(&missing)
      );
    }
  }

  Ok(())
}

/// Checks that `graph` reconstructs exactly the genomes of `sources`, matched by name.
///
/// Used to verify a merged graph against the graphs it was built from. Both sides of every
/// comparison are reconstructed on demand and dropped again, so this holds two genomes at a time
/// instead of the whole sequence content of `sources`. The saving is proportional to the total
/// genome length and modest in practice — peak usage during a merge is dominated by the graphs
/// themselves — but it keeps verification consistent with `reconstruct --verify`, which streams for
/// the same reason.
///
/// Every genome of `sources` must appear in `graph`, and `graph` must contain nothing else. The
/// `sources` must have disjoint genome names: a name appearing in two of them means the two graphs
/// describe overlapping genome sets, and is reported rather than verified twice.
pub fn verify_graph_against_graphs(graph: &Pangraph, sources: &[&Pangraph]) -> Result<(), Report> {
  check_genome_names(sources).wrap_err("When checking the genome names of the input graphs")?;

  let path_ids = path_ids_by_name(graph)?;
  let mut verified: BTreeSet<&str> = BTreeSet::new();

  for source in sources {
    for (name, source_path_id) in path_ids_by_name(source)? {
      let Some(path_id) = path_ids.get(name) else {
        return make_error!("Graph is missing genome '{name}', which is present in the input graphs");
      };

      let expected = reconstruct_genome(source, source_path_id)
        .wrap_err_with(|| format!("When reconstructing genome '{name}' from the input graphs"))?;
      let actual =
        reconstruct_genome(graph, *path_id).wrap_err_with(|| format!("When reconstructing genome '{name}'"))?;

      verify_genome(name, &expected, &actual)?;
      verified.insert(name);
    }
  }

  // Every genome of `sources` has now been found in `graph` and checked, so anything left over is a
  // genome the merged graph invented.
  let extra = path_ids
    .keys()
    .filter(|name| !verified.contains(*name))
    .copied()
    .collect_vec();

  if !extra.is_empty() {
    return make_error!(
      "Graph contains {} genome(s) that are not present in the input graphs: {}",
      extra.len(),
      format_names(&extra)
    );
  }

  Ok(())
}

/// Formats a list of genome names for an error message, eliding all but the first few.
pub(crate) fn format_names<S: AsRef<str>>(names: &[S]) -> String {
  let shown = names.iter().take(MAX_NAMES_IN_ERROR).map(AsRef::as_ref).join(", ");
  if names.len() > MAX_NAMES_IN_ERROR {
    format!("[{shown}, ... and {} more]", names.len() - MAX_NAMES_IN_ERROR)
  } else {
    format!("[{shown}]")
  }
}

/// Reconstructs the genome of a path, by concatenating the sequences of its nodes and rotating the
/// result so that it starts where the genome does.
///
/// Both lookups here are guaranteed by [`Pangraph::validate`], which every graph read from a file
/// passes and every graph pangraph builds satisfies, so a failure is a bug rather than bad input.
/// They are still checked, because reporting beats aborting: indexing the node map and rotating by
/// an out-of-range offset both panic.
fn reconstruct_path_sequence(graph: &Pangraph, path: &PangraphPath) -> Result<Seq, Report> {
  if let Some(first_node_id) = path.nodes.first() {
    let first_node = graph
      .nodes
      .get(first_node_id)
      .ok_or_else(|| make_internal_report!("Node {first_node_id} not found in graph"))?;
    let first_node_pos = first_node.position().0;

    let mut genome: Seq = path
      .nodes
      .iter()
      .map(|node_id| reconstruct_block_sequence(graph, *node_id))
      .collect::<Result<Seq, Report>>()?;

    let genome_len = path.tot_len();
    if genome.len() != genome_len {
      return make_error!(
        "When reconstructing sequences, genome length mismatch: computed length {} expected {}",
        genome.len(),
        genome_len
      );
    }

    if first_node_pos > genome.len() {
      return make_internal_error!(
        "When reconstructing sequences, the first node of the genome starts at position {first_node_pos}, past the end of a genome of length {}",
        genome.len()
      );
    }
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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::pangraph::edits::Edit;
  use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
  use crate::pangraph::pangraph_node::PangraphNode;
  use crate::pangraph::strand::Strand;
  use crate::pangraph::strand::Strand::{Forward, Reverse};
  use crate::utils::error::report_to_string;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  /// A two-genome graph: `a` is a forward single-block path, `b` a reverse one, so their
  /// reconstructed sequences are hand-checkable.
  fn two_genome_graph(names: [Option<&str>; 2]) -> Pangraph {
    let blocks = btreemap! {
      BlockId(0) => PangraphBlock::new(BlockId(0), "ACGTACGT", btreemap!{ NodeId(0) => Edit::empty() }),
      BlockId(1) => PangraphBlock::new(BlockId(1), "TTTTGGGG", btreemap!{ NodeId(1) => Edit::empty() }),
    };
    let nodes = btreemap! {
      NodeId(0) => PangraphNode::new(NodeId(0), BlockId(0), PathId(0), Forward, (0, 8)),
      NodeId(1) => PangraphNode::new(NodeId(1), BlockId(1), PathId(1), Reverse, (0, 8)),
    };
    let paths = btreemap! {
      PathId(0) => PangraphPath::new(Some(PathId(0)), [NodeId(0)], 8, false, names[0].map(String::from), None),
      PathId(1) => PangraphPath::new(Some(PathId(1)), [NodeId(1)], 8, false, names[1].map(String::from), None),
    };
    Pangraph { paths, blocks, nodes }
  }

  /// A one-genome graph, so that a pair of them stands in for the two inputs of a merge. Every id
  /// is `0`, exactly as `Pangraph::singleton` assigns them, which is also what makes two of these
  /// indistinguishable by id.
  fn one_genome_graph(name: &str, consensus: &str, strand: Strand) -> Pangraph {
    let len = consensus.len();
    let blocks = btreemap! {
      BlockId(0) => PangraphBlock::new(BlockId(0), consensus, btreemap!{ NodeId(0) => Edit::empty() }),
    };
    let nodes = btreemap! {
      NodeId(0) => PangraphNode::new(NodeId(0), BlockId(0), PathId(0), strand, (0, len)),
    };
    let paths = btreemap! {
      PathId(0) => PangraphPath::new(Some(PathId(0)), [NodeId(0)], len, false, Some(name.to_owned()), None),
    };
    Pangraph { paths, blocks, nodes }
  }

  /// The two single-genome graphs whose merger `graph()` stands for.
  fn sources() -> (Pangraph, Pangraph) {
    (
      one_genome_graph("a", "ACGTACGT", Forward),
      one_genome_graph("b", "TTTTGGGG", Reverse),
    )
  }

  fn graph() -> Pangraph {
    two_genome_graph([Some("a"), Some("b")])
  }

  fn expected_genomes() -> BTreeMap<String, Seq> {
    btreemap! { o!("a") => Seq::from_str("ACGTACGT"), o!("b") => Seq::from_str("CCCCAAAA") }
  }

  #[rstest]
  fn test_path_ids_by_name_rejects_unnamed_path() {
    let graph = two_genome_graph([Some("a"), None]);
    assert!(report_to_string(&path_ids_by_name(&graph).unwrap_err()).contains("no name or an empty name"));
  }

  /// A name that is empty, or made of whitespace only, identifies a genome no better than a missing
  /// one: `--strains` cannot address it, and it round-trips to FASTA as a nameless record.
  #[rstest]
  #[case("")]
  #[case(" ")]
  fn test_path_ids_by_name_rejects_empty_path_name(#[case] name: &str) {
    let graph = two_genome_graph([Some("a"), Some(name)]);
    let err = report_to_string(&path_ids_by_name(&graph).unwrap_err());
    assert!(err.contains("no name or an empty name"), "unexpected error: {err}");
  }

  /// Every comparison on genome names is exact, so " a" would be a genome that `--strains a`
  /// cannot address, that does not collide with the "a" of the graph it is merged with, and that
  /// round-trips to a FASTA header with an invisible space in it.
  #[rstest]
  #[case(" a")]
  #[case("a ")]
  #[case("\ta")]
  fn test_path_ids_by_name_rejects_padded_path_name(#[case] name: &str) {
    let graph = two_genome_graph([Some("b"), Some(name)]);
    let err = report_to_string(&path_ids_by_name(&graph).unwrap_err());
    assert!(
      err.contains("leading or trailing whitespace"),
      "unexpected error: {err}"
    );
  }

  /// The pair a padded name is most likely to be mistaken for: without the check they are two
  /// distinct genomes, since duplicates are detected on the name as given.
  #[rstest]
  fn test_path_ids_by_name_rejects_a_name_that_differs_only_by_padding() {
    let graph = two_genome_graph([Some("a"), Some("a ")]);
    let err = report_to_string(&path_ids_by_name(&graph).unwrap_err());
    assert!(
      err.contains("leading or trailing whitespace"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_path_ids_by_name_rejects_duplicate_names() {
    let graph = two_genome_graph([Some("a"), Some("a")]);
    assert!(report_to_string(&path_ids_by_name(&graph).unwrap_err()).contains("Duplicate genome names"));
  }

  /// A name may be unique within each graph and still collide across them. This is what makes
  /// `merge` reject merging a graph with itself, and what stops two overlapping source graphs from
  /// being "verified" against a merged graph that holds their genomes only once.
  #[rstest]
  fn test_check_genome_names_rejects_name_shared_across_graphs() {
    let (left, _) = sources();
    let other = one_genome_graph("a", "GGGGCCCC", Forward);

    check_genome_names(&[&left]).unwrap();
    let err = report_to_string(&check_genome_names(&[&left, &other]).unwrap_err());
    assert!(err.contains("Duplicate genome names"), "unexpected error: {err}");
    assert!(err.contains('a'), "unexpected error: {err}");
  }

  #[rstest]
  fn test_check_sequence_names_rejects_duplicate_names() {
    let record = |name: &str| FastaRecord {
      seq_name: name.to_owned(),
      desc: None,
      seq: Seq::from_str("ACGT"),
      index: 0,
    };
    let fastas = [record("a"), record("a")];
    assert!(report_to_string(&check_sequence_names(&fastas).unwrap_err()).contains("Duplicate sequence names"));
  }

  /// `> id` — a space between the '>' and the identifier — parses into an empty name with the
  /// identifier in the description, so this common header style used to yield a nameless genome.
  /// The record index is part of the message, since the name cannot point at the offending record.
  #[rstest]
  #[case("")]
  #[case(" ")]
  fn test_check_sequence_names_rejects_empty_name(#[case] name: &str) {
    let fastas = [FastaRecord {
      seq_name: name.to_owned(),
      desc: Some(o!("NC_000913.3 Escherichia coli")),
      seq: Seq::from_str("ACGT"),
      index: 3,
    }];

    let err = report_to_string(&check_sequence_names(&fastas).unwrap_err());
    assert!(err.contains("empty name"), "unexpected error: {err}");
    assert!(err.contains('3'), "unexpected error: {err}");
  }

  /// A header indented with anything other than a plain space — `>\tid` — keeps the whitespace in
  /// the name, which is invisible in every place the name is later shown.
  #[rstest]
  #[case(" a")]
  #[case("a ")]
  #[case("\ta")]
  fn test_check_sequence_names_rejects_padded_name(#[case] name: &str) {
    let fastas = [FastaRecord {
      seq_name: name.to_owned(),
      desc: None,
      seq: Seq::from_str("ACGT"),
      index: 3,
    }];

    let err = report_to_string(&check_sequence_names(&fastas).unwrap_err());
    assert!(
      err.contains("leading or trailing whitespace"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_verify_graph_sequences_accepts_exact_match() {
    verify_graph_sequences(&graph(), &expected_genomes(), GenomeCoverage::Complete).unwrap();
  }

  /// The whole point of keying by name: `index` and `desc` must not take part in the comparison.
  /// The previous whole-`FastaRecord` comparison failed this, which is why verifying a merged
  /// graph reported a mismatch between two sequences of identical length.
  #[rstest]
  fn test_verify_graph_sequences_ignores_index_and_desc() {
    let fastas = [
      FastaRecord {
        seq_name: o!("a"),
        desc: Some(o!("some description")),
        seq: Seq::from_str("ACGTACGT"),
        index: 41,
      },
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: Seq::from_str("CCCCAAAA"),
        index: 42,
      },
    ];
    let expected = sequences_by_name(&fastas);
    verify_graph_sequences(&graph(), &expected, GenomeCoverage::Complete).unwrap();
  }

  #[rstest]
  fn test_verify_graph_sequences_detects_missing_genome() {
    let mut expected = expected_genomes();
    expected.insert(o!("c"), Seq::from_str("GGGG"));

    let err = verify_graph_sequences(&graph(), &expected, GenomeCoverage::Complete).unwrap_err();
    let err = report_to_string(&err);
    assert!(err.contains("missing"), "unexpected error: {err}");
    assert!(err.contains('c'), "unexpected error: {err}");
  }

  #[rstest]
  fn test_verify_graph_sequences_detects_extra_genome() {
    let mut expected = expected_genomes();
    expected.remove("b");

    let err = report_to_string(&verify_graph_sequences(&graph(), &expected, GenomeCoverage::Complete).unwrap_err());
    assert!(
      err.contains("not among the expected genomes"),
      "unexpected error: {err}"
    );
  }

  /// `Partial` relaxes only one direction: the graph may hold a subset of the expected genomes,
  /// but a genome the expected set does not know about is still an error.
  #[rstest]
  fn test_verify_graph_sequences_partial_accepts_superset() {
    let mut expected = expected_genomes();
    expected.insert(o!("c"), Seq::from_str("GGGG"));

    verify_graph_sequences(&graph(), &expected, GenomeCoverage::Partial).unwrap();
  }

  #[rstest]
  fn test_verify_graph_sequences_partial_still_rejects_extra_genome() {
    let mut expected = expected_genomes();
    expected.remove("b");

    let err = report_to_string(&verify_graph_sequences(&graph(), &expected, GenomeCoverage::Partial).unwrap_err());
    assert!(
      err.contains("not among the expected genomes"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_verify_graph_against_graphs_accepts_exact_match() {
    let (left, right) = sources();
    verify_graph_against_graphs(&graph(), &[&left, &right]).unwrap();
  }

  #[rstest]
  fn test_verify_graph_against_graphs_detects_missing_genome() {
    let (left, right) = sources();
    let extra = one_genome_graph("c", "GGGGCCCC", Forward);

    let err = report_to_string(&verify_graph_against_graphs(&graph(), &[&left, &right, &extra]).unwrap_err());
    assert!(err.contains("missing genome 'c'"), "unexpected error: {err}");
  }

  #[rstest]
  fn test_verify_graph_against_graphs_detects_extra_genome() {
    let (left, _) = sources();

    let err = report_to_string(&verify_graph_against_graphs(&graph(), &[&left]).unwrap_err());
    assert!(
      err.contains("not present in the input graphs"),
      "unexpected error: {err}"
    );
    assert!(err.contains('b'), "unexpected error: {err}");
  }

  /// Two sources holding the same genome name describe overlapping genome sets. Counting one
  /// expectation per (source, name) pair used to make the totals disagree while no genome was
  /// actually extra, reporting the nonsensical "contains 0 genome(s) that are not present ...: []".
  #[rstest]
  fn test_verify_graph_against_graphs_detects_name_shared_across_sources() {
    let (left, right) = sources();
    let duplicate = one_genome_graph("a", "ACGTACGT", Forward);

    let err = report_to_string(&verify_graph_against_graphs(&graph(), &[&left, &right, &duplicate]).unwrap_err());
    assert!(err.contains("Duplicate genome names"), "unexpected error: {err}");
    assert!(!err.contains("0 genome(s)"), "unexpected error: {err}");
  }

  /// The name matches but the sequence does not: the mismatch is reported against the source graph
  /// the genome came from, with the position of the first difference.
  #[rstest]
  fn test_verify_graph_against_graphs_detects_mutated_genome() {
    let (_, right) = sources();
    let mutated = one_genome_graph("a", "ACGTTCGT", Forward);

    let err = report_to_string(&verify_graph_against_graphs(&graph(), &[&mutated, &right]).unwrap_err());
    assert!(
      err.contains("Sequence mismatch for genome 'a'"),
      "unexpected error: {err}"
    );
  }

  #[rstest]
  fn test_verify_genome_reports_first_difference() {
    let expected = Seq::from_str("ACGTACGT");
    let actual = Seq::from_str("ACGTTCGT");

    let err = report_to_string(&verify_genome("a", &expected, &actual).unwrap_err());
    assert!(err.contains("at position 4"), "unexpected error: {err}");
    assert!(err.contains("ACGTACGT"), "unexpected error: {err}");
    assert!(err.contains("ACGTTCGT"), "unexpected error: {err}");
  }

  #[rstest]
  fn test_verify_genome_reports_length_mismatch() {
    let err = report_to_string(&verify_genome("a", &Seq::from_str("ACGT"), &Seq::from_str("ACG")).unwrap_err());
    assert!(err.contains("expected length 4 but got 3"), "unexpected error: {err}");
  }

  #[rstest]
  fn test_format_names_elides_long_lists() {
    let names = (0..12).map(|i| format!("g{i}")).collect_vec();
    assert_eq!(
      format_names(&names),
      "[g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, ... and 2 more]"
    );
    assert_eq!(format_names(&names[..2]), "[g0, g1]");
  }

  /// `Pangraph::validate` rejects both of these before reconstruction is ever reached, so they
  /// stand for a graph assembled in-process rather than read from a file. They are still worth
  /// pinning: both used to be indexing panics, which abort the process instead of being reported,
  /// and `sanity_check` is compiled out of release builds.
  #[rstest]
  fn test_reconstruct_reports_path_referring_to_a_missing_node() {
    let mut graph = graph();
    graph.paths.get_mut(&PathId(0)).unwrap().nodes = vec![NodeId(99)];

    let err = report_to_string(&reconstruct_genome(&graph, PathId(0)).unwrap_err());
    assert!(err.contains("Node 99 not found in graph"), "unexpected error: {err}");
  }

  #[rstest]
  fn test_reconstruct_reports_first_node_starting_past_the_end_of_the_genome() {
    let mut graph = graph();
    let node = &graph.nodes[&NodeId(0)];
    graph.nodes.insert(
      NodeId(0),
      PangraphNode::new(NodeId(0), node.block_id(), node.path_id(), node.strand(), (100, 8)),
    );

    let err = report_to_string(&reconstruct_genome(&graph, PathId(0)).unwrap_err());
    assert!(
      err.contains("starts at position 100, past the end of a genome of length 8"),
      "unexpected error: {err}"
    );
  }
}
