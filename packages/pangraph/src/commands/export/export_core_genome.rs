use crate::commands::export::export_args::PangraphExportCoreAlignmentArgs;
use crate::io::fasta::FastaWriter;
use crate::io::file::create_file_or_stdout;
use crate::io::seq::reverse_complement;
use crate::make_internal_report;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::RecordNaming;
use clap::Parser;
use eyre::{Context, Report};
use itertools::Itertools;
use std::collections::BTreeMap;

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
  let core_block_ids: Vec<_> = graph.core_block_ids().collect();
  let guide_path_id = graph.path_id_by_name(&params.guide_strain)?;
  let guide_path = &graph.paths[&guide_path_id];

  // Extract sequences for all core blocks
  let mut records = vec![];
  for (bid, strand) in guide_path
    .nodes
    .iter()
    .map(|node_id| {
      let node = &graph.nodes[node_id];
      (node.block_id(), node.strand())
    })
    .filter(|(bid, _)| core_block_ids.contains(bid))
  // Filter out non-core blocks
  {
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
    // path A: n1+, n2-
    // path B: n4+, n5+, n3-

    // block 1:
    // cons: ACCTATCGTGATCGTTCGAT
    // A n1: ACCTATCGT---CGTTCGAT
    // B n3: ACTTATCGTGATCGTTCGAT
    // reverse-complement:
    // cons: ATCGAACGATCACGATAGGT
    // A n1: ATCGAACG---ACGATAGGT
    // B n3: ATCGAACGATCACGATAAGT

    // block 2:
    // cons: CTGCAA   GTCTGATCTAGTTA
    // A n2: CTGCAATTTGTCTGATGTAGTTA
    // B n4: CT--AA   GTCTGATCTAGTTA
    // reverse-complement
    // cons: TAACTAGATCAGAC   TTGCAG
    // A n2: TAACTACATCAGACAAATTGCAG
    // B n4: TAACTAGATCAGAC   TT--AG

    // core (guide A):
    // A: ACCTATCGT---CGTTCGATTAACTACATCAGACAAATTGCAG
    // B: ACTTATCGTGATCGTTCGATTAACTAGATCAGAC   TT--AG

    // core (guide B):
    // A: CTGCAATTTGTCTGATGTAGTTAATCGAACG---ACGATAGGT
    // B: CT--AA   GTCTGATCTAGTTAATCGAACGATCACGATAAGT

    let graph = Pangraph::from_str(indoc! {
    // language=json
    r#"
      {
        "paths": {
          "0": {
            "id": 0,
            "nodes": [1, 2],
            "tot_len": 40,
            "circular": false,
            "name": "Path A"
          },
          "1": {
            "id": 1,
            "nodes": [4, 5, 3],
            "tot_len": 48,
            "circular": false,
            "name": "Path B"
          }
        },
        "blocks": {
          "1": {
            "id": 1,
            "consensus": "ACCTATCGTGATCGTTCGAT",
            "alignments": {
              "1": {
                "subs": [],
                "dels": [{"pos": 9, "len": 3}],
                "inss": []
              },
              "3": {
                "subs": [{"pos": 2, "alt": "T"}],
                "dels": [],
                "inss": []
              }
            }
          },
          "2": {
            "id": 2,
            "consensus": "CTGCAAGTCTGATCTAGTTA",
            "alignments": {
              "2": {
                "subs": [{"pos": 13, "alt": "G"}],
                "dels": [],
                "inss": [{"pos": 6, "seq": "TTT"}]
              },
              "4": {
                "subs": [],
                "dels": [{"pos": 2, "len": 2}],
                "inss": []
              }
            }
          },
          "3": {
            "id": 3,
            "consensus": "AGGCTACGAT",
            "alignments": {
              "5": {
                "subs": [],
                "dels": [],
                "inss": []
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
            "position": [0, 17]
          },
          "2": {
            "id": 2,
            "block_id": 2,
            "path_id": 0,
            "strand": "-",
            "position": [17, 40]
          },
          "3": {
            "id": 3,
            "block_id": 1,
            "path_id": 1,
            "strand": "-",
            "position": [28, 48]
          },
          "4": {
            "id": 4,
            "block_id": 2,
            "path_id": 1,
            "strand": "+",
            "position": [0, 18]
          },
          "5": {
            "id": 5,
            "block_id": 3,
            "path_id": 1,
            "strand": "+",
            "position": [18, 28]
          }
        }
      }
      "#})
    .unwrap();

    let params = ExportCoreAlignmentParams {
      guide_strain: o!("Path A"),
      unaligned: false,
    };

    let expected = vec![
      (o!("Path A"), o!("ACCTATCGT---CGTTCGATTAACTACATCAGACTTGCAG")),
      (o!("Path B"), o!("ACTTATCGTGATCGTTCGATTAACTAGATCAGACTT--AG")),
    ];
    let actual = core_block_aln(&graph, &params).unwrap();

    assert_eq!(expected, actual);

    let params = ExportCoreAlignmentParams {
      guide_strain: o!("Path A"),
      unaligned: true,
    };

    let expected = vec![
      (o!("Path A"), o!("ACCTATCGTCGTTCGATTAACTACATCAGACAAATTGCAG")),
      (o!("Path B"), o!("ACTTATCGTGATCGTTCGATTAACTAGATCAGACTTAG")),
    ];

    let actual = core_block_aln(&graph, &params).unwrap();

    assert_eq!(expected, actual);

    let params = ExportCoreAlignmentParams {
      guide_strain: o!("Path B"),
      unaligned: false,
    };

    let expected = vec![
      (o!("Path A"), o!("CTGCAAGTCTGATGTAGTTAATCGAACG---ACGATAGGT")),
      (o!("Path B"), o!("CT--AAGTCTGATCTAGTTAATCGAACGATCACGATAAGT")),
    ];

    let actual = core_block_aln(&graph, &params).unwrap();

    assert_eq!(expected, actual);

    let params = ExportCoreAlignmentParams {
      guide_strain: o!("Path B"),
      unaligned: true,
    };

    let expected = vec![
      (o!("Path A"), o!("CTGCAATTTGTCTGATGTAGTTAATCGAACGACGATAGGT")),
      (o!("Path B"), o!("CTAAGTCTGATCTAGTTAATCGAACGATCACGATAAGT")),
    ];

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
