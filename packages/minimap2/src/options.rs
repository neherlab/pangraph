#![allow(unsafe_code)]
use crate::options_args::{init_opts, Minimap2Args};
use crate::ptr::is_null;
use eyre::{eyre, Report};
use minimap2_sys::{mm_check_opt, mm_idxopt_t, mm_mapopt_t, mm_set_opt};
use std::ffi::CString;
use std::mem::MaybeUninit;
use std::ptr::null;

#[derive(Copy, Clone, Debug)]
pub enum Minimap2Preset {
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

impl Minimap2Preset {
  pub fn as_str(&self) -> &str {
    match self {
      Minimap2Preset::LrHqae => "lr:hqae",
      Minimap2Preset::LrHq => "lr:hq",
      Minimap2Preset::Splice => "splice",
      Minimap2Preset::SpliceHq => "splice:hq",
      Minimap2Preset::Asm => "asm",
      Minimap2Preset::Asm5 => "asm5",
      Minimap2Preset::Asm10 => "asm10",
      Minimap2Preset::Asm20 => "asm20",
      Minimap2Preset::Sr => "sr",
      Minimap2Preset::MapPb => "map-pb",
      Minimap2Preset::MapHifi => "map-hifi",
      Minimap2Preset::MapOnt => "map-ont",
      Minimap2Preset::AvaPb => "ava-pb",
      Minimap2Preset::AvaOnt => "ava-ont",
      Minimap2Preset::Short => "short",
      Minimap2Preset::Map10k => "map10k",
      Minimap2Preset::Cdna => "cdna",
    }
  }
}

#[must_use]
#[derive(Clone, Debug)]
pub struct Minimap2Options {
  idx_opt: mm_idxopt_t,
  map_opt: mm_mapopt_t,
}

impl Minimap2Options {
  pub fn new(args: &Minimap2Args) -> Result<Self, Report> {
    init_options(args)
  }

  pub(crate) fn idx_opt(&self) -> &mm_idxopt_t {
    &self.idx_opt
  }

  pub(crate) fn idx_opt_mut(&mut self) -> *mut mm_idxopt_t {
    &mut self.idx_opt
  }

  pub(crate) fn map_opt(&self) -> &mm_mapopt_t {
    &self.map_opt
  }

  pub(crate) fn map_opt_mut(&mut self) -> *mut mm_mapopt_t {
    &mut self.map_opt
  }
}

fn init_options(args: &Minimap2Args) -> Result<Minimap2Options, Report> {
  let mut idx_opt = MaybeUninit::<mm_idxopt_t>::uninit();
  let mut map_opt = MaybeUninit::<mm_mapopt_t>::uninit();

  // SAFETY: calling unsafe function
  if unsafe { mm_set_opt(null(), idx_opt.as_mut_ptr(), map_opt.as_mut_ptr()) } != 0 {
    return Err(eyre!(
      "minimap2: mm_set_opt(null, ...): failed to set options: incorrect preset"
    ));
  }

  if let Some(preset) = args.x {
    let preset = preset.as_str();
    let preset = CString::new(preset)?;
    let preset = preset.as_ptr();

    // SAFETY: calling unsafe function
    if unsafe { mm_set_opt(preset, idx_opt.as_mut_ptr(), map_opt.as_mut_ptr()) } != 0 {
      return Err(eyre!(
        "minimap2: mm_set_opt(preset, ...): failed to set options: incorrect preset"
      ));
    }
  }

  {
    // SAFETY: calling unsafe function
    let idx_opt = unsafe { idx_opt.assume_init_mut() };

    // SAFETY: calling unsafe function
    let map_opt = unsafe { map_opt.assume_init_mut() };

    init_opts(args, idx_opt, map_opt);
  }

  if is_null(&idx_opt) || is_null(&map_opt) {
    return Err(eyre!("minimap2: mm_set_opt(): null pointer returned"));
  }

  // SAFETY: calling unsafe function
  if unsafe { mm_check_opt(idx_opt.as_ptr(), map_opt.as_ptr()) } != 0 {
    return Err(eyre!("minimap2: mm_check_opt(): options are invalid"));
  }

  if is_null(&idx_opt) || is_null(&map_opt) {
    return Err(eyre!("minimap2: mm_check_opt(): null pointer returned"));
  }

  // SAFETY: calling unsafe function
  let idx_opt = unsafe { idx_opt.assume_init() };

  // SAFETY: calling unsafe function
  let map_opt = unsafe { map_opt.assume_init() };

  Ok(Minimap2Options { idx_opt, map_opt })
}
