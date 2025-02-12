#![allow(warnings)]
#![allow(clippy::all)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

unsafe impl Send for mm_mapopt_t {}
unsafe impl Sync for mm_mapopt_t {}

use std::mem::MaybeUninit;

pub const MM_CIGAR_CHARS: &[u8] = MM_CIGAR_STR.to_bytes();

impl Default for mm_mapopt_t {
  fn default() -> Self {
    unsafe {
      let mut opt = MaybeUninit::uninit();
      mm_mapopt_init(opt.as_mut_ptr());
      opt.assume_init()
    }
  }
}

impl Default for mm_idxopt_t {
  fn default() -> Self {
    unsafe {
      let mut opt = MaybeUninit::uninit();
      mm_idxopt_init(opt.as_mut_ptr());
      opt.assume_init()
    }
  }
}

// TODO: Add more tests!
#[cfg(test)]
mod tests {
  use super::*;
  use std::mem::MaybeUninit;

  #[test]
  fn set_index_and_mapping_opts() {
    let mut mm_idxopt = MaybeUninit::uninit();
    let mut mm_mapopt = MaybeUninit::uninit();

    unsafe { mm_set_opt(std::ptr::null(), mm_idxopt.as_mut_ptr(), mm_mapopt.as_mut_ptr()) };
    println!("{:#?}", unsafe { mm_idxopt.assume_init() });
    println!("{:#?}", unsafe { mm_mapopt.assume_init() }); // Run tests with --nocapture to see the output
  }

  #[test]
  fn mapopt() {
    let x: mm_mapopt_t = Default::default();
  }

  #[test]
  fn idxopt() {
    let x: mm_idxopt_t = Default::default();
  }
}
