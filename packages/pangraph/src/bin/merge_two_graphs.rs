use clap::{Parser, ValueHint};
use ctor::ctor;
use eyre::{Report, WrapErr};
use log::info;
use pangraph::commands::build::build_args::PangraphBuildArgs;
use pangraph::io::json::{JsonPretty, json_write_file};
use pangraph::pangraph::graph_merging::merge_graphs;
use pangraph::pangraph::pangraph::Pangraph;
use pangraph::utils::global_init::global_init;
use pangraph::utils::global_init::setup_logger;
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

  info!("Writing merged graph to {:?}", args.build_args.output_json);

  json_write_file(&args.build_args.output_json, &merged, JsonPretty(true))?;

  Ok(())
}
