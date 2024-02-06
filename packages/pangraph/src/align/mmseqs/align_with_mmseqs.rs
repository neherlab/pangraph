use crate::align::alignment::Alignment;
use crate::align::mmseqs::paf::PafTsvRecord;
use crate::io::fasta::{write_one_fasta, FastaRecord};
use crate::io::file::open_file_or_stdin;
use crate::io::fs::path_to_str;
use crate::o;
use crate::utils::subprocess::{subprocess, subprocess_with_args};
use color_eyre::{Section, SectionExt};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::io::Read;
use std::str::FromStr;
use tempfile::Builder as TempDirBuilder;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct MmseqsParams {
  kmer_len: Option<usize>,
}

pub fn align_with_mmseqs(
  ref_seq: &FastaRecord,
  qry_seq: &FastaRecord,
  params: &MmseqsParams,
) -> Result<Alignment, Report> {
  let temp_dir = TempDirBuilder::new().tempdir()?;
  let temp_dir_path = temp_dir.path();
  let qry_path = path_to_str(&temp_dir_path.join("qry.fa"))?;
  let ref_path = path_to_str(&temp_dir_path.join("reff.fa"))?;
  let res_path = path_to_str(&temp_dir_path.join("res.paf"))?;
  let tmp_path = path_to_str(&temp_dir_path.join("tmp"))?;

  write_one_fasta(&qry_path, &qry_seq.seq_name, &qry_seq.seq)?;
  write_one_fasta(&ref_path, &ref_seq.seq_name, &ref_seq.seq)?;

  #[rustfmt::skip]
  let mut args: Vec<String> = vec![
    "easy-search",
    &qry_path, &ref_path, &res_path, &tmp_path,
    "--threads", "1",
    "--max-seq-len", "10000",
    "-a",
    "--search-type", "3",
    "--format-output", &PafTsvRecord::fields_names().join(","),
  ].into_iter().map(ToOwned::to_owned).collect_vec();

  if let Some(kmer_len) = params.kmer_len {
    if kmer_len > 0 {
      args.extend(vec![o!("-k"), kmer_len.to_string()]);
    }
  }

  let args_ref: Vec<&str> = args.iter().map(AsRef::as_ref).collect_vec();

  let res = subprocess_with_args("mmseqs", args_ref)
    .wrap_err_with(|| {
      format!(
        "When trying to align sequences using mmseqs:\n  ref: '{}'\n  qry: '{}'",
        ref_seq.seq_name, qry_seq.seq_name
      )
    })
    .with_section(|| ref_seq.seq_name.clone().header("ref"))
    .with_section(|| qry_seq.seq_name.clone().header("qry"))?;

  let mut paf_str = String::new();
  open_file_or_stdin(&Some(res_path))?.read_to_string(&mut paf_str)?;

  Alignment::from_paf_str(&paf_str)
}
