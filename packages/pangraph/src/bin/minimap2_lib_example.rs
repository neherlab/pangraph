#![allow(non_snake_case)]
#![allow(clippy::struct_excessive_bools)]

use clap::{AppSettings, ArgEnum, Parser};
use ctor::ctor;
use eyre::{Report, WrapErr};
use itertools::{izip, Itertools};
use minimap2::{Minimap2Args, Minimap2Index, Minimap2Mapper, Minimap2Preset, Minimap2Result};
use pangraph::io::fasta::read_many_fasta;
use pangraph::io::json::json_write;
use pangraph::utils::global_init::global_init;
use rayon::prelude::*;
use std::path::PathBuf;
use std::str::FromStr;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let cli = Minimap2CliArgs::parse();

  let (names, seqs): (Vec<String>, Vec<String>) = read_many_fasta(&cli.input_fastas)?
    .into_iter()
    .map(|f| (f.seq_name, f.seq))
    .unzip();

  let args = minimap2_from_clap(&cli);
  let idx = Minimap2Index::new(&seqs, &names, &args)?;

  // FIXME: output actual PAF
  // let output_paf = create_file_or_stdout(&cli.output_paf)?;
  // let mut csv = CsvStructWriter::new(output_paf, b'\t')?;
  // izip!(seqs, names).try_for_each(move |(seq, name)| -> Result<(), Report> {
  //   let result = mapper.run_map(&seq, &name)?;
  //   result.pafs.iter().try_for_each(|paf| csv.write(paf))?;
  //   Ok(())
  // })?;

  let results: Vec<Minimap2Result> = izip!(seqs, names)
    .par_bridge()
    .map_init(
      || Minimap2Mapper::new(&idx).unwrap(),
      move |mapper, (seq, name)| {
        mapper
          .run_map(&seq, &name)
          .wrap_err_with(|| format!("When aligning sequence '{name}'"))
      },
    )
    .collect::<Result<Vec<_>, Report>>()?;

  // Outputs PAF-like JSON
  let pafs = results.iter().flat_map(|result| &result.pafs).collect_vec();
  json_write(&cli.output_paf, &pafs)?;

  // // Outputs alignment JSON
  // let alns = results
  //   .into_iter()
  //   .map(Alignment::from_minimap_paf_obj)
  //   .collect::<Result<Vec<Vec<_>>, Report>>()?
  //   .into_iter()
  //   .flatten()
  //   .collect_vec();
  // json_write(&cli.output_aln, &alns)?;

  Ok(())
}

#[derive(Debug, Parser)]
#[clap(name = "minimap2_lib_example", trailing_var_arg = true)]
#[clap(author, version)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(verbatim_doc_comment)]
struct Minimap2CliArgs {
  #[clap(display_order = 1)]
  pub input_fastas: Vec<PathBuf>,

  /// Path to output PAF file, or "-" for stdout
  #[clap(short = 'o', long, default_value = "-", value_parser)]
  #[clap(display_order = 1)]
  pub output_paf: PathBuf,

  // /// Path to output Alignment JSON file, or "-" for stdout
  // #[clap(short = 'l', long, default_value = "-", value_parser)]
  // #[clap(display_order = 1)]
  // pub output_aln: PathBuf,
  //
  /// Preset (always applied before other options)
  #[clap(arg_enum, short = 'x')]
  pub preset: Option<Minimap2CliPreset>,

  /// k-mer size (no larger than 28)
  #[clap(short = 'k')]
  pub k: Option<i32>,

  /// Minimizer window size
  #[clap(short = 'w')]
  pub w: Option<i32>,

  /// Skip self and dual mappings (for the all-vs-all mode)
  #[clap(short = 'X')]
  pub X: bool,

  /// Use homopolymer-compressed k-mer
  #[clap(short = 'H')]
  pub H: bool,

  /// Stop chain elongation if there are no minimizers in INT-bp
  #[clap(short = 'g')]
  pub g: Option<i32>,

  /// Max intron length (effective with -xsplice; changing -r)
  #[clap(short = 'G')]
  pub max_intron_len: Option<i32>,

  /// Max fragment length (effective with -xsr or in the fragment mode)
  #[clap(short = 'F')]
  pub F: Option<i32>,

  /// Chaining/alignment bandwidth and long-join bandwidth
  #[clap(short = 'r', value_parser = parse_i32_tuple)]
  pub r: Option<(i32, i32)>,

  /// Min secondary-to-primary score ratio
  #[clap(short = 'p')]
  pub p: Option<f32>,

  /// Filter out top FLOAT fraction of repetitive minimizers
  #[clap(short = 'M')]
  pub frac_rep_min: Option<f32>,

  /// Retain at most INT secondary alignments
  #[clap(short = 'N')]
  pub sec_align: Option<i32>,

  /// Output CIGAR in PAF
  #[clap(short = 'c')]
  pub c: bool,

  /// Matching score
  #[clap(short = 'a')]
  pub a: Option<i32>,

  /// Mismatch penalty (larger value for lower divergence)
  #[clap(short = 'b')]
  pub b: Option<i32>,

  /// Minimal peak DP alignment score
  #[clap(short = 's')]
  pub s: Option<i32>,

  /// Minimal number of minimizers on a chain
  #[clap(short = 'n')]
  pub n: Option<i32>,

  /// Minimal chaining score (matching bases minus log gap penalty)
  #[clap(short = 'm')]
  pub m: Option<i32>,

  /// Gap open penalty
  #[clap(short = 'O', value_parser = parse_i32_tuple)]
  pub gap_extend: Option<(i32, i32)>,

  /// Gap extension penalty
  #[clap(short = 'E', value_parser = parse_i32_tuple)]
  pub gap_open: Option<(i32, i32)>,

  /// Z-drop score and inversion Z-drop score
  #[clap(short = 'z', value_parser = parse_i32_tuple)]
  pub z: Option<(i32, i32)>,

  /// Do not output base quality in SAM
  #[clap(short = 'Q')]
  pub Q: bool,

  /// Use soft clipping for supplementary alignments
  #[clap(short = 'Y')]
  pub Y: bool,

  /// Write CIGAR with >65535 ops at the CG tag
  #[clap(short = 'L')]
  pub L: bool,

  /// Output copy comment in SAM
  #[clap(short = 'y')]
  pub sam_out_copy: bool,

  /// SDUST threshold; 0 to disable SDUST
  #[clap(short = 'T')]
  pub T: Option<i32>,

  /// Split index for every ~NUM input bases
  #[clap(short = 'I')]
  pub I: Option<u64>,

  /// Minibatch size for mapping
  #[clap(short = 'K')]
  pub minibatch_size: Option<i64>,

  /// Occurrence distance
  #[clap(short = 'e')]
  pub e: Option<i32>,

  /// Splice mode. 0: original minimap2 model; 1: miniprot model
  #[clap(short = 'J')]
  pub J: Option<i32>,

  /// max-chain-skip
  #[clap(long)]
  pub max_chain_skip: Option<i32>,

  /// max-chain-iter
  #[clap(long)]
  pub max_chain_iter: Option<i32>,

  /// cap-sw-mat
  #[clap(long)]
  pub max_sw_mat: Option<i64>,

  /// max-qlen
  #[clap(long)]
  pub max_qlen: Option<i32>,

  /// junc-bonus
  #[clap(long)]
  pub junc_bonus: Option<i32>,

  /// chain-gap-scale
  #[clap(long)]
  pub chain_gap_scale: Option<f32>,

  /// chain-skip-scale
  #[clap(long)]
  pub chain_skip_scale: Option<f32>,

  /// alt-drop
  #[clap(long)]
  pub alt_drop: Option<f32>,

  /// mask-len
  #[clap(long)]
  pub mask_len: Option<i32>,

  /// q-occ-frac
  #[clap(long)]
  pub q_occ_frac: Option<f32>,

  /// min-occ-floor
  #[clap(long)]
  pub min_mid_occ: Option<i32>,

  /// max-occ-floor
  #[clap(long)]
  pub max_mid_occ: Option<i32>,

  /// end-seed-pen
  #[clap(long)]
  pub anchor_ext_shift: Option<i32>,

  /// split-prefix
  #[clap(long)]
  pub split_prefix: Option<String>,

  /// cap-kalloc
  #[clap(long)]
  pub cap_kalloc: Option<i64>,

  /// bucket-bits
  #[clap(long)]
  pub bucket_bits: Option<i16>,

  /// seed
  #[clap(long)]
  pub seed: Option<i32>,

  /// rmq-size-cap
  #[clap(long)]
  pub rmq_size_cap: Option<i32>,

  /// rmq-inner-dist
  #[clap(long)]
  pub rmq_inner_dist: Option<i32>,

  /// rmq-rescue-size
  #[clap(long)]
  pub rmq_rescue_size: Option<i32>,

  /// rmq-rescue-ratio
  #[clap(long)]
  pub rmq_rescue_ratio: Option<f32>,

  /// min-dp-len
  #[clap(long)]
  pub max_ksw_len: Option<i32>,

  /// Enable splice mode
  #[clap(long)]
  pub splice: bool,

  /// Disable long join
  #[clap(long)]
  pub no_ljoin: bool,

  /// Enable short read mode
  #[clap(long)]
  pub sr: bool,

  /// End bonus
  #[clap(long)]
  pub end_bonus: Option<i32>,

  /// Disable pairing
  #[clap(long)]
  pub independ_seg: bool,

  /// Index with no sequence
  #[clap(long)]
  pub idx_no_seq: bool,

  /// Forward-only mode
  #[clap(long)]
  pub for_only: bool,

  /// Reverse-only mode
  #[clap(long)]
  pub rev_only: bool,

  /// Maximum clip ratio
  #[clap(long)]
  pub max_clip_ratio: Option<f32>,

  /// Output MD in SAM
  #[clap(long)]
  pub out_md: bool,

  /// Score ambiguous bases
  #[clap(long)]
  pub sc_ambi: Option<i32>,

  /// Extended CIGAR format
  #[clap(long)]
  pub eqx: bool,

  /// PAF output without hits
  #[clap(long)]
  pub paf_no_hit: bool,

  /// No end filter
  #[clap(long)]
  pub no_end_flt: bool,

  /// Hard mask level
  #[clap(long)]
  pub hard_mlevel: bool,

  /// SAM output hits only
  #[clap(long)]
  pub sam_hit_only: bool,

  /// Query strand-specific alignments
  #[clap(long)]
  pub qstrand: bool,

  /// Do not hash names
  #[clap(long)]
  pub no_hash_name: bool,

  /// Include secondary sequences
  #[clap(long)]
  pub secondary_seq: bool,

  /// Disable diagonal chaining
  #[clap(short = 'D')]
  pub D: bool,

  /// Keep all chains
  #[clap(short = 'P')]
  pub all_chain: bool,

  /// Noncanonical penalty
  #[clap(long)]
  pub noncan: Option<i32>,

  /// Min-occ-floor
  #[clap(long)]
  pub min_occ_floor: Option<i32>,

  /// Output in color space, short format
  #[clap(long)]
  pub out_cs: Option<bool>,

  /// Output in color space, long format
  #[clap(long)]
  pub out_cs_long: Option<bool>,

  /// Enable read mapping quality control
  #[clap(long)]
  pub rmq: Option<bool>,

  /// Enable splice for forward strand
  #[clap(long)]
  pub splice_for: Option<bool>,

  /// Enable splice for reverse strand
  #[clap(long)]
  pub splice_rev: Option<bool>,

  /// Enable heap sort
  #[clap(long)]
  pub heap_sort: Option<bool>,

  /// Disable dual mode
  #[clap(long)]
  pub no_dual: Option<bool>,

  /// Suppress printing of secondary alignments
  #[clap(long)]
  pub no_print_2nd: Option<bool>,
}

#[allow(clippy::map_err_ignore)]
fn parse_i32_tuple(s: &str) -> Result<(i32, i32), &'static str> {
  let parts: Vec<&str> = s.split(',').collect();
  if parts.len() == 2 {
    let first = parts[0]
      .trim()
      .parse::<i32>()
      .map_err(|_| "Parse error for first integer")?;
    let second = parts[1]
      .trim()
      .parse::<i32>()
      .map_err(|_| "Parse error for second integer")?;
    Ok((first, second))
  } else {
    Err("Input should be two integers separated by a comma")
  }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
enum Minimap2CliPreset {
  LrHqae,
  LrHq,
  Splice,
  SpliceHq,
  Asm,
  Asm5,
  Asm10,
  Asm20,
  Sr,
  MapPb,
  MapHifi,
  MapOnt,
  AvaPb,
  AvaOnt,
  Short,
  Map10k,
  Cdna,
}

fn minimap2_from_clap_preset(preset: Minimap2CliPreset) -> Minimap2Preset {
  match preset {
    Minimap2CliPreset::LrHqae => Minimap2Preset::LrHqae,
    Minimap2CliPreset::LrHq => Minimap2Preset::LrHq,
    Minimap2CliPreset::Splice => Minimap2Preset::Splice,
    Minimap2CliPreset::SpliceHq => Minimap2Preset::SpliceHq,
    Minimap2CliPreset::Asm => Minimap2Preset::Asm,
    Minimap2CliPreset::Asm5 => Minimap2Preset::Asm5,
    Minimap2CliPreset::Asm10 => Minimap2Preset::Asm10,
    Minimap2CliPreset::Asm20 => Minimap2Preset::Asm20,
    Minimap2CliPreset::Sr => Minimap2Preset::Sr,
    Minimap2CliPreset::MapPb => Minimap2Preset::MapPb,
    Minimap2CliPreset::MapHifi => Minimap2Preset::MapHifi,
    Minimap2CliPreset::MapOnt => Minimap2Preset::MapOnt,
    Minimap2CliPreset::AvaPb => Minimap2Preset::AvaPb,
    Minimap2CliPreset::AvaOnt => Minimap2Preset::AvaOnt,
    Minimap2CliPreset::Short => Minimap2Preset::Short,
    Minimap2CliPreset::Map10k => Minimap2Preset::Map10k,
    Minimap2CliPreset::Cdna => Minimap2Preset::Cdna,
  }
}

fn minimap2_from_clap(args: &Minimap2CliArgs) -> Minimap2Args {
  Minimap2Args {
    x: args.preset.as_ref().map(|x| minimap2_from_clap_preset(*x)),
    k: args.k,
    w: args.w,
    X: args.X,
    max_chain_skip: args.max_chain_skip,
    max_chain_iter: args.max_chain_iter,
    H: args.H,
    g: args.g,
    G: args.max_intron_len,
    F: args.F,
    r: args.r,
    p: args.p,
    M: args.frac_rep_min,
    N: args.sec_align,
    c: args.c,
    a: args.a,
    b: args.b,
    s: args.s,
    n: args.n,
    m: args.m,
    A: args.a,
    B: args.b,
    O: args.gap_extend,
    E: args.gap_open,
    z: args.z,
    Q: args.Q,
    Y: args.Y,
    L: args.L,
    y: args.sam_out_copy,
    T: args.T,
    I: args.I,
    K: args.minibatch_size,
    e: args.e,
    J: args.J,
    max_sw_mat: args.max_sw_mat,
    max_qlen: args.max_qlen,
    junc_bonus: args.junc_bonus,
    chain_gap_scale: args.chain_gap_scale,
    chain_skip_scale: args.chain_skip_scale,
    alt_drop: args.alt_drop,
    mask_len: args.mask_len,
    q_occ_frac: args.q_occ_frac,
    min_mid_occ: args.min_mid_occ,
    max_mid_occ: args.max_mid_occ,
    anchor_ext_shift: args.anchor_ext_shift,
    split_prefix: args.split_prefix.as_ref().map(|s| s.to_owned()),
    cap_kalloc: args.cap_kalloc,
    bucket_bits: args.bucket_bits,
    seed: args.seed,
    rmq_size_cap: args.rmq_size_cap,
    rmq_inner_dist: args.rmq_inner_dist,
    rmq_rescue_size: args.rmq_rescue_size,
    rmq_rescue_ratio: args.rmq_rescue_ratio,
    min_dp_len: args.max_ksw_len,
    splice: args.splice,
    no_ljoin: args.no_ljoin,
    sr: args.sr,
    end_bonus: args.end_bonus,
    independ_seg: args.independ_seg,
    idx_no_seq: args.idx_no_seq,
    for_only: args.for_only,
    rev_only: args.rev_only,
    max_clip_ratio: args.max_clip_ratio,
    out_md: args.out_md,
    sc_ambi: args.sc_ambi,
    eqx: args.eqx,
    paf_no_hit: args.paf_no_hit,
    no_end_flt: args.no_end_flt,
    hard_mlevel: args.hard_mlevel,
    sam_hit_only: args.sam_hit_only,
    qstrand: args.qstrand,
    no_hash_name: args.no_hash_name,
    secondary_seq: args.secondary_seq,
    D: args.D,
    P: args.all_chain,
    noncan: args.noncan,
    min_occ_floor: args.min_occ_floor,
    out_cs: args.out_cs,
    out_cs_long: args.out_cs_long,
    rmq: args.rmq,
    splice_for: args.splice_for,
    splice_rev: args.splice_rev,
    heap_sort: args.heap_sort,
    no_dual: args.no_dual,
    no_print_2nd: args.no_print_2nd,
  }
}
