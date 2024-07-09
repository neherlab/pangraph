use clap::{Parser, ValueHint};
use ctor::ctor;
use eyre::Report;
use minimap2::{Minimap2Index, Minimap2Mapper, Minimap2Options, Minimap2Preset};
use pangraph::align::alignment::{Alignment, Hit};
use pangraph::align::alignment_args::AlignmentArgs;
use pangraph::align::bam::cigar::parse_cigar_str;
use pangraph::io::fasta::read_many_fasta;
use pangraph::pangraph::pangraph_block::BlockId;
use pangraph::pangraph::strand::Strand;
use pangraph::utils::global_init::global_init;
use pretty_assertions::assert_eq;
use std::path::PathBuf;
use std::str::FromStr;

#[ctor]
fn init() {
  global_init();
}

#[derive(Parser, Debug)]
struct Args {
  #[clap(value_hint = ValueHint::FilePath)]
  pub input_query_fastas: Vec<PathBuf>,

  #[clap(flatten)]
  pub params: AlignmentArgs,
}

fn run_with_fasta() -> Result<(), Report> {
  let (names, seqs): (Vec<String>, Vec<String>) = read_many_fasta(&[&PathBuf::from_str("data/flu.fa")?])?
    .into_iter()
    .map(|f| (f.seq_name, f.seq))
    .unzip();

  let options = Minimap2Options::with_preset(Minimap2Preset::Asm20)?;
  let idx = Minimap2Index::new(&seqs[..2], &names[..2], options)?;
  let mut mapper = Minimap2Mapper::new(&idx)?;

  let result = mapper.run_map(&seqs[2], &names[2])?;

  dbg!(&result);

  Ok(())
}

fn run_with_hardcoded() -> Result<(), Report> {
  let ref_seq = "GTAGTTTTTGTACCCCCCGACTGTACGTCCTCCGTCGATAGAAGCAATAAAGGTGACGTC
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
    GTGGCGCACATTTGATGGTGACTAAGCTGCCAAACTGATT";
  let qry_seq = "GTAGTTCTTGTACCCCCCGACTGTACGTACTCCGTCGATTGAATCAATAAAGGTGACGTC
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
    GGCGCAGATTTGATGGTGACTAAGCTGCCAAACTGAGT";

  let expected = vec![Alignment {
    qry: Hit::new(BlockId(1), 998, (0, 996)),
    reff: Hit::new(BlockId(0), 1000, (0, 998)),
    matches: 969,
    length: 998,
    quality: 60,
    orientation: Strand::Forward,
    new_block_id: None, // FIXME
    anchor_block: None, // FIXME
    cigar: parse_cigar_str("545M1D225M1D226M").unwrap(),
    divergence: Some(0.0291),
    align: Some(845.0),
  }];

  let options = Minimap2Options::with_preset(Minimap2Preset::Asm20)?;
  let idx = Minimap2Index::new(&[ref_seq], &["0"], options)?;
  let mut mapper = Minimap2Mapper::new(&idx)?;
  let result = mapper.run_map(qry_seq, "1")?;

  dbg!(&result);
  let actual = Alignment::from_minimap_paf_obj(result)?;

  assert_eq!(expected, actual);

  Ok(())
}

fn main() -> Result<(), Report> {
  // run_with_fasta()?;
  run_with_hardcoded()?;
  Ok(())
}
