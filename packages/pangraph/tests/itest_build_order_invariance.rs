#[cfg(test)]
mod tests {
  use eyre::Report;
  use itertools::Itertools;
  use pangraph::commands::build::build_args::PangraphBuildArgs;
  use pangraph::commands::build::build_run::build;
  use pangraph::io::fasta::FastaRecord;
  use pangraph::pangraph::pangraph::Pangraph;
  use pangraph::representation::seq::Seq;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::io::Write;
  use tempfile::NamedTempFile;

  /// Deterministic pseudo-random nucleotide sequence, so the test does not depend on an RNG.
  fn random_seq(len: usize, seed: u64) -> String {
    let mut state = seed;
    std::iter::repeat_with(|| {
      state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
      b"ACGT"[((state >> 33) % 4) as usize] as char
    })
    .take(len)
    .collect()
  }

  /// Applies a substitution every `period` bases, to make sequences homologous but not identical.
  fn mutate(seq: &str, period: usize) -> String {
    seq
      .chars()
      .enumerate()
      .map(|(i, c)| {
        if i % period == 0 {
          match c {
            'A' => 'G',
            'G' => 'A',
            'C' => 'T',
            _ => 'C',
          }
        } else {
          c
        }
      })
      .collect()
  }

  /// Four genomes sharing two homologous core segments, each with its own accessory segment.
  ///
  /// Four genomes (rather than three) matter: with a balanced topology the final merge joins two
  /// clades of *equal* depth, and with no ambiguous bases the anchor choice is a full tie. That is
  /// the case the old reference-first tie-break resolved using the input order, so a three-genome
  /// fixture - where depth always decides - cannot detect the defect.
  fn genomes() -> Vec<(String, String)> {
    let core_a = random_seq(2500, 42);
    let core_b = random_seq(2500, 7);
    let mosaic =
      |a: usize, b: usize, acc: u64| format!("{}{}{}", mutate(&core_a, a), random_seq(700, acc), mutate(&core_b, b));
    vec![
      ("g1".to_owned(), mosaic(50, 70, 1)),
      ("g2".to_owned(), mosaic(55, 75, 2)),
      ("g3".to_owned(), mosaic(60, 80, 3)),
      ("g4".to_owned(), mosaic(65, 85, 4)),
    ]
  }

  /// Builds a graph from the given genomes, in the given order, optionally with a guide tree.
  fn build_in_order(order: &[usize], newick: Option<&str>) -> Result<Pangraph, Report> {
    let all = genomes();
    let fastas = order
      .iter()
      .enumerate()
      .map(|(index, &g)| FastaRecord {
        seq_name: all[g].0.clone(),
        desc: None,
        seq: Seq::from_str(&all[g].1),
        index,
      })
      .collect_vec();

    // The tree file has to outlive the build call.
    let tree_file = newick
      .map(|nwk| -> Result<NamedTempFile, Report> {
        let mut f = NamedTempFile::new()?;
        write!(f, "{nwk}")?;
        f.flush()?;
        Ok(f)
      })
      .transpose()?;

    let args = PangraphBuildArgs {
      guide_tree: tree_file.as_ref().map(|f| f.path().to_owned()),
      ..PangraphBuildArgs::default()
    };

    build(fastas, &args, true)
  }

  /// Multiset of block consensus sequences: the graph structure, independent of block identifiers.
  fn consensus_multiset(graph: &Pangraph) -> Vec<String> {
    graph
      .blocks
      .values()
      .map(|b| b.consensus().as_str().to_owned())
      .sorted()
      .collect()
  }

  /// All six argument orders must produce the same graph.
  ///
  /// Block identifiers are derived from the position of a record in the input, and used to be
  /// passed to the aligner as sequence names. Since `minimap2 -X` picks the direction of each
  /// pairwise alignment by comparing those names, argument order silently decided which block was
  /// query and which was reference - and the reference's consensus was the one the merged block
  /// inherited.
  ///
  /// The `no_guide_tree` case additionally guards the neighbor-joining path: `Q.argmin()` resolves
  /// ties by matrix row order, so unless the leaves are ordered by content the inferred topology
  /// itself flips with the argument order.
  #[rstest]
  #[case::no_guide_tree(None)]
  #[case::balanced_guide_tree(Some("((g1,g2),(g3,g4));"))]
  #[case::ladder_guide_tree(Some("(((g1,g2),g3),g4);"))]
  #[trace]
  fn test_build_is_invariant_under_argument_order(#[case] newick: Option<&str>) -> Result<(), Report> {
    let permutations = [
      [0, 1, 2, 3],
      [1, 0, 2, 3],
      [3, 2, 1, 0],
      [2, 3, 0, 1],
      [1, 3, 0, 2],
      [3, 0, 2, 1],
    ];

    let expected = consensus_multiset(&build_in_order(&permutations[0], newick)?);
    assert!(
      expected.len() > 1,
      "expected a non-trivial graph, got {} block(s)",
      expected.len()
    );

    for order in &permutations[1..] {
      let actual = consensus_multiset(&build_in_order(order, newick)?);
      assert_eq!(expected, actual, "argument order {order:?} produced a different graph");
    }

    Ok(())
  }

  /// Swapping two siblings in the guide tree does not change the topology, so it must not change
  /// the graph either. `graph_join` is symmetric and `merge_graphs` uses left/right only for
  /// logging, so this holds once the aligner input is canonical.
  #[rstest]
  fn test_build_is_invariant_under_guide_tree_sibling_swap() -> Result<(), Report> {
    let straight = build_in_order(&[0, 1, 2, 3], Some("((g1,g2),(g3,g4));"))?;
    let swapped = build_in_order(&[0, 1, 2, 3], Some("((g2,g1),(g4,g3));"))?;

    assert_eq!(consensus_multiset(&straight), consensus_multiset(&swapped));

    Ok(())
  }
}
