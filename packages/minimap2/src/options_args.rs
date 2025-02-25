#![allow(non_snake_case)]
#![allow(clippy::struct_excessive_bools)]

use crate::Minimap2Preset;
use minimap2_sys::{
  MM_F_ALL_CHAINS, MM_F_CIGAR, MM_F_COPY_COMMENT, MM_F_EQX, MM_F_FOR_ONLY, MM_F_HARD_MLEVEL, MM_F_HEAP_SORT,
  MM_F_INDEPEND_SEG, MM_F_LONG_CIGAR, MM_F_NO_DIAG, MM_F_NO_DUAL, MM_F_NO_END_FLT, MM_F_NO_HASH_NAME, MM_F_NO_INV,
  MM_F_NO_LJOIN, MM_F_NO_PRINT_2ND, MM_F_NO_QUAL, MM_F_OUT_CG, MM_F_OUT_CS, MM_F_OUT_CS_LONG, MM_F_OUT_MD,
  MM_F_PAF_NO_HIT, MM_F_QSTRAND, MM_F_REV_ONLY, MM_F_RMQ, MM_F_SAM_HIT_ONLY, MM_F_SECONDARY_SEQ, MM_F_SOFTCLIP,
  MM_F_SPLICE, MM_F_SPLICE_FOR, MM_F_SPLICE_REV, MM_F_SR, MM_I_HPC, MM_I_NO_SEQ, mm_idxopt_t, mm_mapopt_t,
};
use std::os::raw::c_short;

#[derive(Debug, Default)]
pub struct Minimap2Args {
  /// preset (always applied before other options)
  pub x: Option<Minimap2Preset>,

  /// "k-mer size (no larger than 28)"
  pub k: Option<i32>,

  /// "minimizer window size"
  pub w: Option<i32>,

  /// "skip self and dual mappings (for the all-vs-all mode)"
  pub X: bool,

  /// "--max-chain-skip"
  pub max_chain_skip: Option<i32>,

  /// "--max-chain-iter"
  pub max_chain_iter: Option<i32>,

  /// "use homopolymer-compressed k-mer"
  pub H: bool,

  /// "stop chain enlongation if there are no minimizers in INT-bp"
  pub g: Option<i32>,

  /// "max intron length (effective with -xsplice; changing -r)"
  pub G: Option<i32>,

  /// "max fragment length (effective with -xsr or in the fragment mode)"
  pub F: Option<i32>,

  /// "chaining/alignment bandwidth and long-join bandwidth"
  pub r: Option<(i32, i32)>,

  /// "min secondary-to-primary score ratio"
  pub p: Option<f32>,

  /// "filter out top FLOAT fraction of repetitive minimizers"
  pub M: Option<f32>,

  /// "retain at most INT secondary alignments"
  pub N: Option<i32>,

  /// "output CIGAR in PAF"
  pub c: bool,

  /// "matching score"
  pub a: Option<i32>,

  /// "mismatch penalty (larger value for lower divergence)"
  pub b: Option<i32>,

  /// "minimal peak DP alignment score"
  pub s: Option<i32>,

  /// "minimal number of minimizers on a chain"
  pub n: Option<i32>,

  /// "minimal chaining score (matching bases minus log gap penalty)"
  pub m: Option<i32>,

  /// "matching score"
  pub A: Option<i32>,

  /// "mismatch penalty"
  pub B: Option<i32>,

  /// "gap open penalty"
  pub O: Option<(i32, i32)>,

  /// "gap extension penalty"
  pub E: Option<(i32, i32)>,

  /// "Z-drop score and inversion Z-drop score"
  pub z: Option<(i32, i32)>,

  /// "do not output base quality in SAM"
  pub Q: bool,

  /// "use soft clipping for supplementary alignments"
  pub Y: bool,

  /// "write CIGAR with >65535 ops at the CG tag"
  pub L: bool,

  /// "output copy comment in SAM"
  pub y: bool,

  /// "SDUST threshold; 0 to disable SDUST"
  pub T: Option<i32>,

  /// "split index for every ~NUM input bases"
  pub I: Option<u64>,

  /// "minibatch size for mapping"
  pub K: Option<i64>,

  /// "occurance distance"
  pub e: Option<i32>,

  /// "splice mode. 0: original minimap2 model; 1: miniprot model"
  pub J: Option<i32>,

  /// "disable diagonal chaining"
  pub D: bool,

  /// "keep all chains"
  pub P: bool,

  /// "noncanonical penalty"
  pub noncan: Option<i32>,

  /// "--cap-sw-mat"
  pub max_sw_mat: Option<i64>,

  /// "--max-qlen"
  pub max_qlen: Option<i32>,

  /// "--junc-bonus"
  pub junc_bonus: Option<i32>,

  /// "--chain-gap-scale"
  pub chain_gap_scale: Option<f32>,

  /// "--chain-skip-scale"
  pub chain_skip_scale: Option<f32>,

  /// "--alt-drop"
  pub alt_drop: Option<f32>,

  /// "--mask-len"
  pub mask_len: Option<i32>,

  /// "--q-occ-frac"
  pub q_occ_frac: Option<f32>,

  /// "--min-occ-floor"
  pub min_mid_occ: Option<i32>,

  /// "--max-occ-floor"
  pub max_mid_occ: Option<i32>,

  /// "--end-seed-pen"
  pub anchor_ext_shift: Option<i32>,

  /// "--split-prefix"
  pub split_prefix: Option<String>,

  /// "--cap-kalloc"
  pub cap_kalloc: Option<i64>,

  /// "--bucket-bits"
  pub bucket_bits: Option<i16>,

  /// "--seed"
  pub seed: Option<i32>,

  /// "--rmq-size-cap"
  pub rmq_size_cap: Option<i32>,

  /// "--rmq-inner-dist"
  pub rmq_inner_dist: Option<i32>,

  /// "--rmq-rescue-size"
  pub rmq_rescue_size: Option<i32>,

  /// "--rmq-rescue-ratio"
  pub rmq_rescue_ratio: Option<f32>,

  /// "--min-dp-len"
  pub min_dp_len: Option<i32>,

  /// "--splice"
  pub splice: bool,

  /// "--no-long-join"
  pub no_ljoin: bool,

  /// "--sr"
  pub sr: bool,

  /// "--end-bonus"
  pub end_bonus: Option<i32>,

  /// "--no-pairing"
  pub independ_seg: bool,

  /// "--idx-no-seq"
  pub idx_no_seq: bool,

  /// "--min-occ-floor"
  pub min_occ_floor: Option<i32>,

  /// "--for-only"
  pub for_only: bool,

  /// "--rev-only"
  pub rev_only: bool,

  /// "--max-clip-ratio"
  pub max_clip_ratio: Option<f32>,

  /// "--MD"
  pub out_md: bool,

  /// "--score-N"
  pub sc_ambi: Option<i32>,

  /// "--eqx"
  pub eqx: bool,

  /// "--paf-no-hit"
  pub paf_no_hit: bool,

  /// "--no-end-flt"
  pub no_end_flt: bool,

  /// "--hard-mask-level"
  pub hard_mlevel: bool,

  /// "--sam-hit-only"
  pub sam_hit_only: bool,

  /// "--qstrand"
  pub qstrand: bool,

  /// "--no-hash-name"
  pub no_hash_name: bool,

  /// "--secondary-seq"
  pub secondary_seq: bool,

  /// "Output in color space, short format"
  pub out_cs: Option<bool>,

  /// "Output in color space, long format"
  pub out_cs_long: Option<bool>,

  /// "Enable read mapping quality control"
  pub rmq: Option<bool>,

  /// "Enable splice for forward strand"
  pub splice_for: Option<bool>,

  /// "Enable splice for reverse strand"
  pub splice_rev: Option<bool>,

  /// "Enable heap sort"
  pub heap_sort: Option<bool>,

  /// "Disable dual mode"
  pub no_dual: Option<bool>,

  /// "Suppress printing of secondary alignments"
  pub no_print_2nd: Option<bool>,
}

#[allow(clippy::redundant_pattern_matching)]
pub fn init_opts(args: &Minimap2Args, idx_opt: &mut mm_idxopt_t, map_opt: &mut mm_mapopt_t) {
  if let Some(k) = args.k {
    idx_opt.k = k as i16;
  }
  if let Some(w) = args.w {
    idx_opt.w = w as i16;
  }
  if args.H {
    idx_opt.flag |= MM_I_HPC as c_short;
  }
  if let Some(g) = args.g {
    map_opt.max_gap = g;
  }
  if let Some(_) = args.G {
    unimplemented!("-G")
  }
  if let Some(F) = args.F {
    map_opt.max_frag_len = F;
  }
  if let Some((r1, r2)) = args.r {
    map_opt.bw = r1;
    map_opt.bw_long = r2;
  }
  if let Some(p) = args.p {
    map_opt.pri_ratio = p;
  }
  if let Some(M) = args.M {
    map_opt.mask_level = M;
  }
  if let Some(N) = args.N {
    map_opt.best_n = N;
  }
  if args.c {
    map_opt.flag |= (MM_F_OUT_CG | MM_F_CIGAR) as i64;
  }
  if let Some(a) = args.a {
    map_opt.a = a;
  }
  if let Some(b) = args.b {
    map_opt.b = b;
  }
  if let Some(s) = args.s {
    map_opt.min_dp_max = s;
  }
  if let Some(n) = args.n {
    map_opt.min_cnt = n;
  }
  if let Some(m) = args.m {
    map_opt.min_chain_score = m;
  }
  if args.X {
    map_opt.flag |= (MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN) as i64;
  }
  if let Some(O) = args.O {
    map_opt.q = O.0;
    map_opt.q2 = O.1;
  }
  if let Some(E) = args.E {
    map_opt.e = E.0;
    map_opt.e2 = E.1;
  }
  if let Some(z) = args.z {
    map_opt.zdrop = z.0;
    map_opt.zdrop_inv = z.1;
  }
  if args.Q {
    map_opt.flag |= MM_F_NO_QUAL as i64;
  }
  if args.Y {
    map_opt.flag |= MM_F_SOFTCLIP as i64;
  }
  if args.L {
    map_opt.flag |= MM_F_LONG_CIGAR as i64;
  }
  if args.y {
    map_opt.flag |= MM_F_COPY_COMMENT as i64;
  }
  if let Some(T) = args.T {
    map_opt.sdust_thres = T;
  }
  if let Some(I) = args.I {
    idx_opt.batch_size = I;
  }
  if let Some(K) = args.K {
    map_opt.mini_batch_size = K;
  }
  if let Some(e) = args.e {
    map_opt.occ_dist = e;
  }
  if let Some(_) = args.J {
    unimplemented!("-J")
  }
  if args.D {
    map_opt.flag |= MM_F_NO_DIAG as i64;
  }
  if args.P {
    map_opt.flag |= MM_F_ALL_CHAINS as i64;
  }
  if let Some(noncan) = args.noncan {
    map_opt.noncan = noncan;
  }
  if let Some(min_occ_floor) = args.min_occ_floor {
    map_opt.min_mid_occ = min_occ_floor;
  }
  if let Some(max_sw_mat) = args.max_sw_mat {
    map_opt.max_sw_mat = max_sw_mat;
  }
  if let Some(max_qlen) = args.max_qlen {
    map_opt.max_qlen = max_qlen;
  }
  if let Some(junc_bonus) = args.junc_bonus {
    map_opt.junc_bonus = junc_bonus;
  }
  if let Some(chain_gap_scale) = args.chain_gap_scale {
    map_opt.chain_gap_scale = chain_gap_scale;
  }
  if let Some(chain_skip_scale) = args.chain_skip_scale {
    map_opt.chain_skip_scale = chain_skip_scale;
  }
  if let Some(alt_drop) = args.alt_drop {
    map_opt.alt_drop = alt_drop;
  }
  if let Some(mask_len) = args.mask_len {
    map_opt.mask_len = mask_len;
  }
  if let Some(q_occ_frac) = args.q_occ_frac {
    map_opt.q_occ_frac = q_occ_frac;
  }
  if let Some(min_mid_occ) = args.min_mid_occ {
    map_opt.min_mid_occ = min_mid_occ;
  }
  if let Some(max_mid_occ) = args.max_mid_occ {
    map_opt.max_mid_occ = max_mid_occ;
  }
  if let Some(anchor_ext_shift) = args.anchor_ext_shift {
    map_opt.anchor_ext_shift = anchor_ext_shift;
  }
  if let Some(_) = &args.split_prefix {
    unimplemented!("--split-prefix")
  }
  if let Some(cap_kalloc) = args.cap_kalloc {
    map_opt.cap_kalloc = cap_kalloc;
  }
  if let Some(bucket_bits) = args.bucket_bits {
    idx_opt.bucket_bits = bucket_bits;
  }
  if let Some(seed) = args.seed {
    map_opt.seed = seed;
  }
  if let Some(rmq_size_cap) = args.rmq_size_cap {
    map_opt.rmq_size_cap = rmq_size_cap;
  }
  if let Some(rmq_inner_dist) = args.rmq_inner_dist {
    map_opt.rmq_inner_dist = rmq_inner_dist;
  }
  if let Some(rmq_rescue_size) = args.rmq_rescue_size {
    map_opt.rmq_rescue_size = rmq_rescue_size;
  }
  if let Some(rmq_rescue_ratio) = args.rmq_rescue_ratio {
    map_opt.rmq_rescue_ratio = rmq_rescue_ratio;
  }
  if let Some(min_ksw_len) = args.min_dp_len {
    map_opt.min_ksw_len = min_ksw_len;
  }
  if args.splice {
    map_opt.flag |= MM_F_SPLICE as i64;
  }
  if args.no_ljoin {
    map_opt.flag |= MM_F_NO_LJOIN as i64;
  }
  if args.sr {
    map_opt.flag |= MM_F_SR as i64;
  }
  if let Some(end_bonus) = args.end_bonus {
    map_opt.end_bonus = end_bonus;
  }
  if args.independ_seg {
    map_opt.flag |= MM_F_INDEPEND_SEG as i64;
  }
  if args.idx_no_seq {
    idx_opt.flag |= MM_I_NO_SEQ as c_short;
  }
  if args.for_only {
    map_opt.flag |= MM_F_FOR_ONLY as i64;
  }
  if args.rev_only {
    map_opt.flag |= MM_F_REV_ONLY as i64;
  }
  if let Some(max_clip_ratio) = args.max_clip_ratio {
    map_opt.max_clip_ratio = max_clip_ratio;
  }
  if args.out_md {
    map_opt.flag |= MM_F_OUT_MD as i64;
  }
  if let Some(sc_ambi) = args.sc_ambi {
    map_opt.sc_ambi = sc_ambi;
  }
  if args.eqx {
    map_opt.flag |= MM_F_EQX as i64;
  }
  if args.paf_no_hit {
    map_opt.flag |= MM_F_PAF_NO_HIT as i64;
  }
  if args.no_end_flt {
    map_opt.flag |= MM_F_NO_END_FLT as i64;
  }
  if args.hard_mlevel {
    map_opt.flag |= MM_F_HARD_MLEVEL as i64;
  }
  if args.sam_hit_only {
    map_opt.flag |= MM_F_SAM_HIT_ONLY as i64;
  }
  if args.qstrand {
    map_opt.flag |= (MM_F_QSTRAND | MM_F_NO_INV) as i64;
  }
  if args.no_hash_name {
    map_opt.flag |= MM_F_NO_HASH_NAME as i64;
  }
  if args.secondary_seq {
    map_opt.flag |= MM_F_SECONDARY_SEQ as i64;
  }
  if let Some(out_cs) = args.out_cs {
    if out_cs {
      map_opt.flag |= MM_F_OUT_CS as i64;
    } else {
      map_opt.flag &= !(MM_F_OUT_CS as i64);
    }
  }
  if let Some(out_cs_long) = args.out_cs_long {
    if out_cs_long {
      map_opt.flag |= MM_F_OUT_CS_LONG as i64;
    } else {
      map_opt.flag &= !(MM_F_OUT_CS_LONG as i64);
    }
  }
  if let Some(rmq) = args.rmq {
    if rmq {
      map_opt.flag |= MM_F_RMQ as i64;
    } else {
      map_opt.flag &= !(MM_F_RMQ as i64);
    }
  }
  if let Some(splice_for) = args.splice_for {
    if splice_for {
      map_opt.flag |= MM_F_SPLICE_FOR as i64;
    } else {
      map_opt.flag &= !(MM_F_SPLICE_FOR as i64);
    }
  }
  if let Some(splice_rev) = args.splice_rev {
    if splice_rev {
      map_opt.flag |= MM_F_SPLICE_REV as i64;
    } else {
      map_opt.flag &= !(MM_F_SPLICE_REV as i64);
    }
  }

  if let Some(heap_sort) = args.heap_sort {
    if heap_sort {
      map_opt.flag |= MM_F_HEAP_SORT as i64;
    } else {
      map_opt.flag &= !(MM_F_HEAP_SORT as i64);
    }
  }
  if let Some(no_dual) = args.no_dual {
    if no_dual {
      map_opt.flag |= MM_F_NO_DUAL as i64;
    } else {
      map_opt.flag &= !(MM_F_NO_DUAL as i64);
    }
  }
  if let Some(no_print_2nd) = args.no_print_2nd {
    if no_print_2nd {
      map_opt.flag |= MM_F_NO_PRINT_2ND as i64;
    } else {
      map_opt.flag &= !(MM_F_NO_PRINT_2ND as i64);
    }
  }
}
