use crate::commands::export::export_args::PangraphExportCoreAlignmentArgs;
use crate::io::fasta::{FastaRecord, FastaWriter};
use crate::io::file::create_file_or_stdout;
use crate::io::seq::reverse_complement;
use crate::make_internal_report;
use crate::pangraph::pangraph::Pangraph;
use crate::pangraph::pangraph_block::RecordNaming;
use crate::representation::seq::Seq;
use clap::Parser;
use eyre::{Context, Report};
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
    .try_for_each(|(index, record)| {
      output_fasta
        .write(&record.seq_name, &record.desc, &record.seq)
        .wrap_err_with(|| format!("When writing sequence #{index} '{}'", record.seq_name))
    })?;

  Ok(())
}

/// Given a block, returns a list of FastaRecord objects containing a sequence per node.
/// If aligned is True, it returns aligned sequences, with gaps for deletions and no insertions.
/// If aligned is False, it returns the full unaligned sequences.
fn core_block_aln(graph: &Pangraph, params: &ExportCoreAlignmentParams) -> Result<Vec<FastaRecord>, Report> {
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
    let mut block_records: Vec<FastaRecord> = block
      .sequences(graph, !params.unaligned, RecordNaming::Path)
      .collect::<Result<Vec<_>, Report>>()?;

    // Reverse-complement if reverse-complemented on guide strain
    if strand.is_reverse() {
      block_records.iter_mut().try_for_each(|record| -> Result<(), Report> {
        record.seq = reverse_complement(&record.seq)?;
        Ok(())
      })?;
    }

    records.push(block_records); // TODO(perf): avoid allocations
  }

  let records = if records.is_empty() {
    // If no record: return empty records
    graph
      .paths()
      .enumerate()
      .map(|(i, path)| {
        let seq_name = path.name.as_ref().map_or_else(|| i.to_string(), |p| p.to_owned());
        FastaRecord {
          seq_name,
          desc: path.desc.clone(),
          seq: Seq::new(),
          index: i,
        }
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
fn concatenate_records(records_list: &[Vec<FastaRecord>]) -> Result<Vec<FastaRecord>, Report> {
  if records_list.is_empty() {
    return Ok(vec![]);
  }

  let mut records: BTreeMap<String, FastaRecord> = records_list[0]
    .iter()
    .map(|record| {
      (
        record.seq_name.clone(),
        FastaRecord {
          index: record.index,
          seq_name: record.seq_name.clone(),
          desc: record.desc.clone(),
          seq: Seq::new(),
        },
      )
    })
    .collect();

  for entries in records_list {
    for entry in entries {
      let record = records
        .get_mut(&entry.seq_name)
        .ok_or_else(|| make_internal_report!("Sequence name '{}' not found in the initial set", &entry.seq_name))?;
      record.seq.extend_seq(&entry.seq);
    }
  }

  Ok(records.into_values().collect())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::o;
  use crate::utils::error::report_to_string;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::str::FromStr;

  fn fasta(id: impl AsRef<str>, seq: impl AsRef<str>) -> FastaRecord {
    FastaRecord {
      seq_name: id.as_ref().to_owned(),
      desc: None,
      seq: Seq::from_str(seq.as_ref()),
      index: 0,
    }
  }

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
      fasta("Path A", "ACCTATCGT---CGTTCGATTAACTACATCAGACTTGCAG"),
      fasta("Path B", "ACTTATCGTGATCGTTCGATTAACTAGATCAGACTT--AG"),
    ];
    let actual = core_block_aln(&graph, &params).unwrap();

    assert_eq!(expected, actual);

    let params = ExportCoreAlignmentParams {
      guide_strain: o!("Path A"),
      unaligned: true,
    };

    let expected = vec![
      fasta("Path A", "ACCTATCGTCGTTCGATTAACTACATCAGACAAATTGCAG"),
      fasta("Path B", "ACTTATCGTGATCGTTCGATTAACTAGATCAGACTTAG"),
    ];

    let actual = core_block_aln(&graph, &params).unwrap();

    assert_eq!(expected, actual);

    let params = ExportCoreAlignmentParams {
      guide_strain: o!("Path B"),
      unaligned: false,
    };

    let expected = vec![
      fasta("Path A", "CTGCAAGTCTGATGTAGTTAATCGAACG---ACGATAGGT"),
      fasta("Path B", "CT--AAGTCTGATCTAGTTAATCGAACGATCACGATAAGT"),
    ];

    let actual = core_block_aln(&graph, &params).unwrap();

    assert_eq!(expected, actual);

    let params = ExportCoreAlignmentParams {
      guide_strain: o!("Path B"),
      unaligned: true,
    };

    let expected = vec![
      fasta("Path A", "CTGCAATTTGTCTGATGTAGTTAATCGAACGACGATAGGT"),
      fasta("Path B", "CTAAGTCTGATCTAGTTAATCGAACGATCACGATAAGT"),
    ];

    let actual = core_block_aln(&graph, &params).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_concatenate_records_general_case() {
    let input = vec![
      vec![fasta("name1", "ATG"), fasta("name2", "CCC")],
      vec![fasta("name1", "TGA"), fasta("name2", "GGG")],
    ];
    let expected = vec![fasta("name1", "ATGTGA"), fasta("name2", "CCCGGG")];
    let actual = concatenate_records(&input).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_concatenate_records_single_entry() {
    let input = vec![vec![fasta("name1", "ATG")]];
    let expected = vec![fasta("name1", "ATG")];
    let actual = concatenate_records(&input).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_concatenate_records_multiple_entries_same_name() {
    let input = vec![
      vec![fasta("name1", "ATG"), fasta("name1", "TGA")],
      vec![fasta("name1", "CCC")],
    ];
    let expected = vec![fasta("name1", "ATGTGACCC")];
    let actual = concatenate_records(&input).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_concatenate_records_missing_name_in_initial_set() {
    let input = vec![
      vec![fasta("name1", "ATG"), fasta("name2", "CCC")],
      vec![fasta("name3", "TGA")],
    ];
    let actual = report_to_string(&concatenate_records(&input).unwrap_err());
    let expected =
      "Sequence name 'name3' not found in the initial set. This is an internal error. Please report it to developers.";
    assert_eq!(expected, actual);
  }
}
