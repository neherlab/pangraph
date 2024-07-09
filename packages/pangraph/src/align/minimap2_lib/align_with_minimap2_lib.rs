use crate::align::alignment::{Alignment, Hit};
use crate::align::alignment_args::AlignmentArgs;
use crate::make_internal_error;
use crate::pangraph::pangraph_block::{BlockId, PangraphBlock};
use crate::pangraph::strand::Strand;
use eyre::{Report, WrapErr};
use minimap2::{Minimap2Index, Minimap2Mapper, Minimap2Options, Minimap2Preset, Minimap2Result};
use noodles::sam::record::Cigar;
use std::collections::BTreeMap;
use std::str::FromStr;

pub fn align_with_minimap2_lib(
  blocks: &BTreeMap<BlockId, PangraphBlock>,
  params: &AlignmentArgs,
) -> Result<Vec<Alignment>, Report> {
  // let kmer_size = create_arg_optional("-k", &params.kmer_length);
  // let preset = create_arg_optional("-x", &Some(format!("asm{}", &params.sensitivity)));
  let options = Minimap2Options::with_preset(Minimap2Preset::Asm20)?;

  let (names, seqs): (Vec<String>, Vec<&str>) = blocks
    .iter()
    .map(|(id, block)| (id.to_string(), block.consensus()))
    .unzip();

  let idx =
    Minimap2Index::new(&seqs, &names, options).wrap_err("When initializing alignment index using minimap2 library")?;

  let mut mapper = Minimap2Mapper::new(&idx).wrap_err("When initializing alignment mapper using minimap2 library")?;

  let results = blocks
    .iter()
    .map(|(id, block)| (id.to_string(), block.consensus()))
    .map(|(name, seq)| {
      let res = mapper
        .run_map(&name, seq)
        .wrap_err_with(|| format!("When trying to align sequence using minimap2 library: '{name}'"))?;

      Ok((name, seq, res))
    })
    .collect::<Result<Vec<_>, Report>>()?
    .into_iter()
    .flat_map(|(name, seq, res)| {
      let Minimap2Result { regs, pafs } = res;
      pafs.into_iter().map(move |paf| {
        if let Some(cg) = paf.cg {
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
            cigar: Cigar::from_str(&cg)?,
            divergence: paf.de.map(|x| x.0),
            align: paf.AS.map(|x| x as f64),
          })
        } else {
          make_internal_error!("Unable to find CIGAR string in the result")
            .wrap_err_with(|| format!("When trying to align sequence using minimap2 library: '{name}'"))
        }
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;

  dbg!(&results);

  Ok(results)
}

// #[cfg(test)]
// mod tests {
//   use super::*;
//   use super::*;
//   use crate::align::alignment::Hit;
//   use crate::pangraph::strand::Strand;
//   use pretty_assertions::assert_eq;
//   use rstest::rstest;
//
//   // TODO: add test cases for alignment of multiple sequences at once
//   #[rstest]
//   fn test_align_with_minimap2_lib_one_general_case() -> Result<(), Report> {
//     let ref_seq = "GTAGTTTTTGTACCCCCCGACTGTACGTCCTCCGTCGATAGAAGCAATAAAGGTGACGTC
//     TGACTACTTTTGGGTGTATATGAAATTCAACACACAAGGTAGCCAAGGCAACGATTTGTT
//     CAGCACTCCTATGGTGCGGTCCGGTGGCAAACGATTTCTGACTAACATCCGCATCCGGCT
//     AACAGCAAGGAAGCCAGGTTTCGGTTCTGTGCACTTCAGGGGTCGCGCATCATGTTGACT
//     TGCTTCGCCAGAAATAAATCTTTACAGTCGTAGCGAAGCATGGTATGGCTCGTCCTACTT
//     CATCATTTGGATACTTATACTAGCTTGGGCGCTGACTAGTAGTGGCTGCTTGAGCCCCCT
//     CACACGCAATCAACCAAGTCCATCTGTCATCATGACTATCGCTGAGCTGGATGAGCGCCC
//     TACACACCGTCTGATTTCCACATTCCTGCAGGAGCTTACCCGGACCACCGTCAATACACG
//     GGATAATAAGGCATTTGATCTGTCTTAAACCTGTTTGCGAATAACTATTCACTATACCAC
//     CATCGTTTCATACATGCAAAAGGTTTGGGTGTACCTTATGCTAGGCAGGACCTTTTTAGG
//     TGTATAGATAGGCCCCAGTAAATAAATGCAATATGGAGATACAACCAATAAACCAAATAT
//     CGTCTACTATAGGTAAATAGTCCGTTATACATCTCAATTGGAGCGGTTGATGAGCTGACA
//     TGTTGGTTACTGTCGACGTCTAAGATGGCCGCGCAGAATATCCCGCACTTCTAATCATGA
//     AAGATAAAGCTGCTGGTCTGGTGGGATTCTGGCGATCTCACCACTGATGCGAGGTCCGTG
//     CAATCCAGATGACGAAAGCGCTCCCGCCAACGATGACGCAAATAATCGTGCGGTAGGGAG
//     TCCTCGTCCCGCACTTTCGGGGACACGTACCCGAAGGGTGTAACGGATGCCCTATGTGAG
//     GTGGCGCACATTTGATGGTGACTAAGCTGCCAAACTGATT";
//     let qry_seq = "GTAGTTCTTGTACCCCCCGACTGTACGTACTCCGTCGATTGAATCAATAAAGGTGACGTC
//     TGACTACTATTGGGTGTATATGAAATTCAACACACAAGGTAGCCAAGGCAACGATTTGTT
//     CAGCACTCCTATGGTGCGGTCCGGTGGCAAACGATTTCTGACTTACATCCGCATCCGGCT
//     AACAGCAAGGAAGCCAGGTTTCCGTGCTGTGCACTTCAGGGGTCGCGCATCATGTTGACT
//     TGCTTCGCCAGAAATAAATCTTTACAGTCGTAGCGAAGCATGGTATGGCTCGTCCCACTT
//     CATCATTTGGATACTTATACTAGATTGGGCGCTGACTAGTGGTGGCTGCTTGAGCCCCCT
//     CACACGCAATCAACCAAGTCCATCTGTCATCATGACTATCGCTGAGCTGGATGCGCGCCC
//     TACACACCGTCTGACTTCCATATTGCAGCAGGAGCTTACCCGGACCACCGTCAATACACG
//     GGATAAGAAGGCATTTGATCTGTCTTAAACCTGTTTGCGAATAACTATTCACTATACCAC
//     CATCGCTCATACCTGCAAAAGGTTTGAGTGTACCTTATGCTAGGCAGGGCCTTTTTAGGT
//     GTATAGATAGGCCCCAGTAAATAAATGCAATATGGAGATACAACCAATAAACCAAATATC
//     GTCTACCATAGGTAAATAGTCCGTTATACATCTTAATTGGAGCGGTTGATGAGCTGACAT
//     GTTGGTTACTGTTGACGTCTAAGATGGCCGCGCAGAATATCCCGCACTTCAATCATGAAA
//     GACAAAGCTGCTGGTCTGGTGGGATTCTGGCGATCTCACCACTGATGCGAGGTCCGTGCA
//     ATCCAGATGACGAAAGCGCTCCCGCCCACGATGACGCAAATAATCGTGCGGTAGGGAGTC
//     CTCGTCCCGCACTTTCGGGGACACGTACCCGAAGGGTGTAACGGATGCCCTATGTGAGGT
//     GGCGCAGATTTGATGGTGACTAAGCTGCCAAACTGAGT";
//
//     let params = AlignmentArgs {
//       kmer_length: Some(10),
//       sensitivity: 20,
//       ..AlignmentArgs::default()
//     };
//
//     let expected = vec![Alignment {
//       qry: Hit::new(BlockId(1), 998, (0, 996)),
//       reff: Hit::new(BlockId(0), 1000, (0, 998)),
//       matches: 969,
//       length: 998,
//       quality: 60,
//       orientation: Strand::Forward,
//       new_block_id: None, // FIXME
//       anchor_block: None, // FIXME
//       cigar: parse_cigar_str("545M1D225M1D226M").unwrap(),
//       divergence: Some(0.0291),
//       align: Some(845.0),
//     }];
//
//     assert_eq!(expected, actual);
//     Ok(())
//   }
// }
