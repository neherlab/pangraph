use crate::align::alignment::{Alignment, Hit};
use crate::align::alignment_args::AlignmentArgs;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::strand::Strand;
use crate::{make_error, make_internal_error};
use eyre::{Report, WrapErr};
use itertools::{Itertools, izip};
use minimap2::{Minimap2Args, Minimap2Index, Minimap2Mapper, Minimap2Preset, Minimap2Result};
use noodles::sam::record::Cigar;
use num_traits::clamp_min;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::str::FromStr;

pub fn align_with_minimap2_lib(
  blocks: &BTreeMap<BlockId, PangraphBlock>,
  params: &AlignmentArgs,
) -> Result<Vec<Alignment>, Report> {
  let (names, seqs): (Vec<String>, Vec<&str>) = blocks
    .iter()
    .map(|(id, block)| (id.to_string(), block.consensus().as_str()))
    .unzip();

  let alns: Vec<Alignment> = align_with_minimap2_lib_impl(&seqs, &names, params)?;

  Ok(alns)
}

fn align_with_minimap2_lib_impl(
  seqs: &[impl AsRef<str>],
  names: &[impl AsRef<str>],
  params: &AlignmentArgs,
) -> Result<Vec<Alignment>, Report> {
  if names.len() != seqs.len() {
    return make_internal_error!(
      "Number of sequences and number of sequence names is expected to be the same, but found: {} sequences and {} names",
      seqs.len(),
      names.len()
    );
  }

  let preset = match params.sensitivity {
    5 => Ok(Minimap2Preset::Asm5),
    10 => Ok(Minimap2Preset::Asm10),
    20 => Ok(Minimap2Preset::Asm20),
    _ => make_error!("Unknown sensitivity preset: {}", params.sensitivity),
  }?;

  let args = Minimap2Args {
    x: Some(preset),
    k: params.kmer_length.map(|v| v as i32),
    c: true,
    X: true,
    s: Some(clamp_min(params.indel_len_threshold - 10, 5) as i32),
    bucket_bits: Some(14),
    ..Default::default()
  };

  let seqs = seqs.iter().map(AsRef::as_ref).collect_vec();
  let names = names.iter().map(AsRef::as_ref).collect_vec();

  let idx = Minimap2Index::new(&seqs, &names, &args)?;

  let results: Vec<Minimap2Result> = izip!(&seqs, &names)
    .par_bridge()
    .map_init(
      || Minimap2Mapper::new(&idx).unwrap(),
      move |mapper, (seq, name)| {
        mapper
          .run_map(seq, name)
          .wrap_err_with(|| format!("When aligning sequence '{name}'"))
      },
    )
    .collect::<Result<Vec<_>, Report>>()?;

  let alns = results
    .into_iter()
    .map(Alignment::from_minimap_paf_obj)
    .collect::<Result<Vec<Vec<_>>, Report>>()?
    .into_iter()
    .flatten()
    .collect_vec();

  Ok(alns)
}

#[allow(clippy::multiple_inherent_impl)]
impl Alignment {
  pub fn from_minimap_paf_obj(res: Minimap2Result) -> Result<Vec<Self>, Report> {
    let Minimap2Result { pafs, .. } = res;
    pafs
      .into_iter()
      .map(|paf| {
        if let Some(cg) = &paf.cg {
          Ok(Alignment {
            qry: Hit::new(
              BlockId::from_str(&paf.q.name)?,
              paf.q.len,
              (paf.q.start as usize, paf.q.end as usize),
            ),
            reff: Hit::new(
              BlockId::from_str(&paf.t.name)?,
              paf.t.len,
              (paf.t.start as usize, paf.t.end as usize),
            ),
            matches: paf.mlen as usize,
            length: paf.blen as usize,
            quality: paf.mapq as usize,
            orientation: Strand::from_char(paf.strand)?,
            new_block_id: None, // FIXME: initialize?
            anchor_block: None, // FIXME: initialize?
            cigar: Cigar::from_str(cg)?,
            divergence: paf.de.map(|x| x.0),
            align: paf.AS.map(|x| x as f64),
          })
        } else {
          make_internal_error!("Unable to find CIGAR string in the result")
        }
      })
      .collect()
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  use crate::align::alignment::Hit;
  use crate::align::bam::cigar::parse_cigar_str;
  use crate::io::fasta::read_many_fasta_str;
  use crate::pangraph::strand::Strand;
  use pretty_assertions::assert_eq;

  // TODO: add test cases for alignment of multiple sequences at once
  #[test]
  fn test_align_with_minimap2_lib_one_general_case() -> Result<(), Report> {
    let fasta = r#"
      >1
      GTAGTTTTTGTACCCCCCGACTGTACGTCCTCCGTCGATAGAAGCAATAAAGGTGACGTC
      TGACTACTTTTGGGTGTATATGAAATTCAACACACAAGGTAGCCAAGGCAACGATTTGTT
      CAGCACTCCTATGGTGCGGTCCGGTGGCAAACGATTTCTGACTAACATCCGCATCCGGCT
      AACAGCAAGGAAGCCAGGTTTCGGTTCTGTGCACTTCAGGGGTCGCGCATCATGTTGACT
      TGCTTCGCCAGAAATAAATCTTTACAGTCGTAGCGAAGCATGGTATGGCTCGTCCTACTT
      CATCATTTGGATACTTATACTAGCTTGGGCGCTGACTAGTAGTGGCTGCTTGAGCCCCCT
      CACACGCAATCAACCAAGTCCATCTGTCATCATGACTATCGCTGAGCTGGATGAGCGCCC
      TACACACCGTCTGATTTCCACATTCCTGCAGGAGCTTACCCGGACCACCGTCAATACACG
      GGATAATAAGGCATTTGATCTGTCTTAAACCTGTTTGCGAATAACTATTCACTATACCAC
      CATCGTTTCATACATGCAAAAGGTTTGGGTGTACCTTATGCTAGGCAGGACCTTTTTAGG
      TGTATAGATAGGCCCCAGTAAATAAATGCAATATGGAGATACAACCAATAAACCAAATAT
      CGTCTACTATAGGTAAATAGTCCGTTATACATCTCAATTGGAGCGGTTGATGAGCTGACA
      TGTTGGTTACTGTCGACGTCTAAGATGGCCGCGCAGAATATCCCGCACTTCTAATCATGA
      AAGATAAAGCTGCTGGTCTGGTGGGATTCTGGCGATCTCACCACTGATGCGAGGTCCGTG
      CAATCCAGATGACGAAAGCGCTCCCGCCAACGATGACGCAAATAATCGTGCGGTAGGGAG
      TCCTCGTCCCGCACTTTCGGGGACACGTACCCGAAGGGTGTAACGGATGCCCTATGTGAG
      GTGGCGCACATTTGATGGTGACTAAGCTGCCAAACTGATT
      >0
      GTAGTTCTTGTACCCCCCGACTGTACGTACTCCGTCGATTGAATCAATAAAGGTGACGTC
      TGACTACTATTGGGTGTATATGAAATTCAACACACAAGGTAGCCAAGGCAACGATTTGTT
      CAGCACTCCTATGGTGCGGTCCGGTGGCAAACGATTTCTGACTTACATCCGCATCCGGCT
      AACAGCAAGGAAGCCAGGTTTCCGTGCTGTGCACTTCAGGGGTCGCGCATCATGTTGACT
      TGCTTCGCCAGAAATAAATCTTTACAGTCGTAGCGAAGCATGGTATGGCTCGTCCCACTT
      CATCATTTGGATACTTATACTAGATTGGGCGCTGACTAGTGGTGGCTGCTTGAGCCCCCT
      CACACGCAATCAACCAAGTCCATCTGTCATCATGACTATCGCTGAGCTGGATGCGCGCCC
      TACACACCGTCTGACTTCCATATTGCAGCAGGAGCTTACCCGGACCACCGTCAATACACG
      GGATAAGAAGGCATTTGATCTGTCTTAAACCTGTTTGCGAATAACTATTCACTATACCAC
      CATCGCTCATACCTGCAAAAGGTTTGAGTGTACCTTATGCTAGGCAGGGCCTTTTTAGGT
      GTATAGATAGGCCCCAGTAAATAAATGCAATATGGAGATACAACCAATAAACCAAATATC
      GTCTACCATAGGTAAATAGTCCGTTATACATCTTAATTGGAGCGGTTGATGAGCTGACAT
      GTTGGTTACTGTTGACGTCTAAGATGGCCGCGCAGAATATCCCGCACTTCAATCATGAAA
      GACAAAGCTGCTGGTCTGGTGGGATTCTGGCGATCTCACCACTGATGCGAGGTCCGTGCA
      ATCCAGATGACGAAAGCGCTCCCGCCCACGATGACGCAAATAATCGTGCGGTAGGGAGTC
      CTCGTCCCGCACTTTCGGGGACACGTACCCGAAGGGTGTAACGGATGCCCTATGTGAGGT
      GGCGCAGATTTGATGGTGACTAAGCTGCCAAACTGAGT
      "#;

    let fasta = read_many_fasta_str(fasta)?;

    let (seqs, names): (Vec<String>, Vec<String>) = fasta.into_iter().map(|f| (f.seq.to_string(), f.seq_name)).unzip();

    let params = AlignmentArgs {
      kmer_length: Some(10),
      sensitivity: 20,
      ..AlignmentArgs::default()
    };

    let actual = align_with_minimap2_lib_impl(&seqs, &names, &params)?;

    let expected = vec![Alignment {
      qry: Hit::new(BlockId(0), 998, (0, 996)),
      reff: Hit::new(BlockId(1), 1000, (0, 998)),
      matches: 969,
      length: 998,
      quality: 0,
      orientation: Strand::Forward,
      new_block_id: None, // FIXME
      anchor_block: None, // FIXME
      cigar: parse_cigar_str("545M1D225M1D226M")?,
      divergence: Some(0.029058116232464903),
      align: Some(845.0),
    }];

    assert_eq!(expected, actual);
    Ok(())
  }
}
