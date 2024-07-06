#![allow(unsafe_code)]
use crate::Minimap2Options;
use eyre::{eyre, Report};
use minimap2_sys::{mm_idx_destroy, mm_idx_str, mm_idx_t, mm_mapopt_update};
use std::ffi::CString;
use std::os::raw::{c_char, c_int};

#[must_use]
#[derive(Debug)]
pub struct Minimap2Index {
  idx: *mut mm_idx_t,
  options: Minimap2Options,
}

impl Minimap2Index {
  pub fn new(seqs: &[&str], names: &[&str], mut options: Minimap2Options) -> Result<Self, Report> {
    assert_eq!(seqs.len(), names.len());
    assert!(seqs.len() > 0);

    let idxopt = options.idx_opt();
    let w: c_int = idxopt.w.into();
    let k: c_int = idxopt.k.into();
    let is_hpc: c_int = (idxopt.flag & 1).into();
    let bucket_bits: c_int = idxopt.bucket_bits.into();

    let n: c_int = seqs.len().try_into()?;

    let seqs: Vec<CString> = seqs.iter().map(|&s| CString::new(s).unwrap()).collect();
    let mut seqs: Vec<*const c_char> = seqs.iter().map(|s| s.as_ptr()).collect();

    let names: Vec<CString> = names.iter().map(|&s| CString::new(s).unwrap()).collect();
    let mut names: Vec<*const c_char> = names.iter().map(|s| s.as_ptr()).collect();

    // SAFETY: calling unsafe function
    let idx: *mut mm_idx_t = unsafe { mm_idx_str(w, k, is_hpc, bucket_bits, n, seqs.as_mut_ptr(), names.as_mut_ptr()) };
    if idx.is_null() {
      return Err(eyre!("minimap2: failed to create index"));
    }

    // SAFETY: calling unsafe function
    unsafe { mm_mapopt_update(options.map_opt_mut(), idx) };

    Ok(Self { idx, options })
  }

  pub fn options(&self) -> &Minimap2Options {
    &self.options
  }

  pub(crate) fn get(&self) -> *const mm_idx_t {
    self.idx
  }

  pub(crate) fn get_mut(&self) -> *mut mm_idx_t {
    self.idx
  }
}

impl Drop for Minimap2Index {
  fn drop(&mut self) {
    // SAFETY: calling unsafe function
    unsafe {
      mm_idx_destroy(self.idx);
    }
  }
}
