use crate::commands::export::export_args::PangraphExportCoreAlignmentArgs;
use crate::io::fasta::FastaWriter;
use crate::io::file::create_file_or_stdout;
use crate::io::seq::reverse_complement;
use crate::make_internal_report;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::RecordNaming;
use clap::Parser;
use eyre::{Context, Report};
use itertools::{izip, Itertools};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Parser, Debug, Default, Clone)]
pub struct ExportCoreAlignmentParams {
  /// Specify the strain to use as a reference for the alignment.
  /// Core blocks are ordered and oriented (forward or reverse) according to the reference strain.
  #[clap(long)]
  pub guide_strain: String,

  /// If set, then the full core sequences are exported but not aligned.
  ///
  /// They should be linearly alignable and can be fed to an external aligner.
  #[clap(long)]
  pub unaligned: bool,
}

#[allow(clippy::needless_pass_by_value)]
pub fn export_core_genome(args: PangraphExportCoreAlignmentArgs) -> Result<(), Report> {
  let PangraphExportCoreAlignmentArgs {
    input_json,
    output,
    params,
  } = &args;

  let pangraph = Pangraph::from_path(input_json)?;
  let mut output_fasta = FastaWriter::new(create_file_or_stdout(output)?);

  core_block_aln(&pangraph, params)?
    .into_iter()
    .enumerate()
    .try_for_each(|(index, (id, seq))| {
      output_fasta
        .write(&id, &seq)
        .wrap_err_with(|| format!("When writing sequence #{index} '{id}'"))
    })?;

  Ok(())
}

/// Given a block, returns a list of FastaRecord objects containing a sequence per node.
/// If aligned is True, it returns aligned sequences, with gaps for deletions and no insertions.
/// If aligned is False, it returns the full unaligned sequences.
fn core_block_aln(graph: &Pangraph, params: &ExportCoreAlignmentParams) -> Result<Vec<(String, String)>, Report> {
  let core_block_ids = graph.core_block_ids();
  let guide_path_id = graph.path_id_by_name(&params.guide_strain)?;
  let guide_path = &graph.paths[&guide_path_id];

  let guide_block_info: BTreeMap<_, _> = guide_path
    .nodes
    .iter()
    .map(|node_id| {
      let node = &graph.nodes[node_id];
      (node.block_id(), node.strand())
    })
    .collect();

  let guide_block_ids: BTreeSet<_> = guide_block_info.keys().collect();
  let guide_block_values: Vec<_> = guide_block_info.values().collect();

  // Extract sequences for all core blocks
  let mut records = vec![];
  for (bid, &strand) in izip!(core_block_ids, guide_block_values) {
    if !guide_block_ids.contains(&bid) {
      continue; // Only consider core blocks
    }

    let block = &graph.blocks[&bid];
    let mut block_records: Vec<_> = block
      .sequences(graph, !params.unaligned, RecordNaming::Path)
      .map(|(id, seq)| Ok((id, seq?)))
      .collect::<Result<_, Report>>()?; // TODO(perf): avoid allocations

    // Reverse-complement if reverse-complemented on guide strain
    if strand.is_reverse() {
      block_records = block_records
        .into_iter()
        .map(|(id, seq)| Ok((id, reverse_complement(seq)?)))
        .collect::<Result<_, Report>>()?; // TODO(perf): avoid allocations
    }

    records.push(block_records); // TODO(perf): avoid allocations
  }

  let records = if records.is_empty() {
    // If no record: return empty records
    graph
      .path_names()
      .enumerate()
      .map(|(i, path)| {
        let id = path.map_or_else(|| i.to_string(), |p| p.to_owned());
        (id, String::new())
      })
      .collect()
  } else {
    // Else concatenate them in a single record set
    concatenate_records(&records)?
  };

  Ok(records)
}

/// Given a list of lists of FastaRecord objects, concatenates all records
/// with the same name in the same sequence, and returns a single list of FastaRecord objects.
fn concatenate_records(records_list: &[Vec<(String, String)>]) -> Result<Vec<(String, String)>, Report> {
  if records_list.is_empty() {
    return Ok(vec![]);
  }

  let mut records: BTreeMap<String, String> = records_list[0]
    .iter()
    .map(|(id, _)| (id.clone(), String::new()))
    .collect();

  for entry in records_list {
    for (id, seq) in entry {
      records
        .get_mut(id)
        .ok_or_else(|| make_internal_report!("Sequence name '{id}' not found in the initial set"))?
        .push_str(seq); // TODO(perf): avoid allocations
    }
  }

  let records = records
    .into_iter()
    .enumerate()
    .map(|(idx, (name, seq))| (name, seq))
    .collect_vec();

  Ok(records)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::utils::error::report_to_string;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::str::FromStr;

  #[test]
  fn test_core_block_aln_general_case() {
    let graph = Pangraph::from_str(indoc! {
    // language=json
    r#"
      {
        "paths": {
          "0": {
            "id": 0,
            "nodes": [1, 2],
            "tot_len": 100,
            "circular": false,
            "name": "Path A"
          },
          "1": {
            "id": 1,
            "nodes": [3, 4],
            "tot_len": 100,
            "circular": false,
            "name": "Path B"
          }
        },
        "blocks": {
          "1": {
            "id": 1,
            "consensus": "ACGTACGT",
            "alignments": {
              "1": {
                "subs": [],
                "dels": [],
                "inss": []
              },
              "2": {
                "subs": [{"pos": 3, "alt": "T"}],
                "dels": [],
                "inss": []
              }
            }
          },
          "2": {
            "id": 2,
            "consensus": "TTGCA",
            "alignments": {
              "3": {
                "subs": [],
                "dels": [{"pos": 2, "len": 1}],
                "inss": []
              },
              "4": {
                "subs": [],
                "dels": [],
                "inss": [{"pos": 0, "seq": "AAA"}]
              }
            }
          }
        },
        "nodes": {
          "1": {
            "id": 1,
            "block_id": 1,
            "path_id": 0,
            "strand": "+",
            "position": [0, 8]
          },
          "2": {
            "id": 2,
            "block_id": 1,
            "path_id": 0,
            "strand": "-",
            "position": [0, 8]
          },
          "3": {
            "id": 3,
            "block_id": 2,
            "path_id": 1,
            "strand": "+",
            "position": [0, 5]
          },
          "4": {
            "id": 4,
            "block_id": 2,
            "path_id": 1,
            "strand": "+",
            "position": [0, 5]
          }
        }
      }
      "#})
    .unwrap();

    let params = ExportCoreAlignmentParams {
      guide_strain: o!("Path A"),
      unaligned: true,
    };

    let expected = vec![(o!("Path A"), o!("ACGTTGTCA--")), (o!("Path B"), o!("TT--CGAATTTGCA"))];
    let actual = core_block_aln(&graph, &params).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_concatenate_records_general_case() {
    let input = vec![
      vec![(o!("name1"), o!("ATG")), (o!("name2"), o!("CCC"))],
      vec![(o!("name1"), o!("TGA")), (o!("name2"), o!("GGG"))],
    ];
    let expected = vec![(o!("name1"), o!("ATGTGA")), (o!("name2"), o!("CCCGGG"))];
    let actual = concatenate_records(&input).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_concatenate_records_single_entry() {
    let input = vec![vec![(o!("name1"), o!("ATG"))]];
    let expected = vec![(o!("name1"), o!("ATG"))];
    let actual = concatenate_records(&input).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_concatenate_records_multiple_entries_same_name() {
    let input = vec![
      vec![(o!("name1"), o!("ATG")), (o!("name1"), o!("TGA"))],
      vec![(o!("name1"), o!("CCC"))],
    ];
    let expected = vec![(o!("name1"), o!("ATGTGACCC"))];
    let actual = concatenate_records(&input).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_concatenate_records_missing_name_in_initial_set() {
    let input = vec![
      vec![(o!("name1"), o!("ATG")), (o!("name2"), o!("CCC"))],
      vec![(o!("name3"), o!("TGA"))],
    ];
    let actual = report_to_string(&concatenate_records(&input).unwrap_err());
    let expected =
      "Sequence name 'name3' not found in the initial set. This is an internal error. Please report it to developers.";
    assert_eq!(expected, actual);
  }
}
