use clap::{Parser, ValueHint};
use ctor::ctor;
use eyre::{Report, WrapErr};
use log::info;
use pangraph::commands::build::build_args::PangraphBuildArgs;
use pangraph::commands::reconstruct::reconstruct_run::{compare_sequences, reconstruct};
use pangraph::io::json::{JsonPretty, json_write_file};
use pangraph::pangraph::graph_merging::merge_graphs;
use pangraph::pangraph::pangraph::Pangraph;
use pangraph::utils::global_init::global_init;
use pangraph::utils::global_init::setup_logger;
use std::collections::HashMap;
use std::path::PathBuf;

#[ctor]
fn init() {
  global_init();
  setup_logger(log::LevelFilter::Debug);
}

#[derive(Parser, Debug)]
struct Args {
  #[clap(flatten)]
  pub build_args: PangraphBuildArgs,

  #[clap(long, short = 'L')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub left_graph: Option<PathBuf>,

  #[clap(long, short = 'R')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub right_graph: Option<PathBuf>,

  // binary flag: if set, verify the sequences in the merged graph did not change
  #[clap(long, short = 'V')]
  pub merge_verify: bool,
}

fn main() -> Result<(), Report> {
  let mut args = Args::parse();

  rayon::ThreadPoolBuilder::new()
    .num_threads(num_cpus::get())
    .build_global()?;

  // Manually empty the input_fastas field in PangraphBuildArgs
  args.build_args.input_fastas = vec![];

  info!("Reading left graph from {:?}", args.left_graph);
  info!("Reading right graph from {:?}", args.right_graph);

  let left = Pangraph::from_path(&args.left_graph).wrap_err("When reading left graph")?;
  let right = Pangraph::from_path(&args.right_graph).wrap_err("When reading right graph")?;

  info!("left graph sanity check");
  #[cfg(debug_assertions)]
  left.sanity_check()?;

  info!("right graph sanity check");
  #[cfg(debug_assertions)]
  right.sanity_check()?;

  info!("Merging graphs");

  let merged = merge_graphs(&left, &right, &args.build_args)?;

  #[cfg(debug_assertions)]
  merged.sanity_check()?;

  info!("Writing merged graph to {}", args.build_args.output_json.display());

  json_write_file(&args.build_args.output_json, &merged, JsonPretty(true))?;

  if args.merge_verify {
    let initial_left_seqs = reconstruct(&left);
    let initial_right_seqs = reconstruct(&right);

    // create a map of {seq.idx -> seq}
    let mut initial_seqs = HashMap::new();
    initial_left_seqs
      .chain(initial_right_seqs)
      .try_for_each(|seq| -> Result<(), Report> {
        let seq = seq?;
        initial_seqs.insert(seq.index, seq);
        Ok(())
      })?;

    // verify that the merged graph contains the same sequences
    let mut merged_seqs = reconstruct(&merged);
    merged_seqs.try_for_each(|seq| -> Result<(), Report> {
      let seq = seq?;
      let idx = seq.index;
      if let Some(expected) = initial_seqs.get(&idx) {
        compare_sequences(expected, &seq)?;
      } else {
        return Err(eyre::eyre!("Sequence with index {idx} not found in initial sequences"));
      }
      Ok(())
    })?;
    info!("Merged graph contains the same sequences as the original graphs");
  }

  Ok(())
}
