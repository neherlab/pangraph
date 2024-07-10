#![allow(unsafe_code)]
use crate::options::Minimap2Options;
use crate::options_args::Minimap2Args;
use eyre::{eyre, Report};
use minimap2_sys::{mm_idx_destroy, mm_idx_str, mm_idx_t, mm_mapopt_update, MM_I_HPC};
use std::ffi::{c_short, CString};
use std::os::raw::{c_char, c_int};

#[must_use]
#[derive(Debug)]
pub struct Minimap2Index {
  idx: *mut mm_idx_t,
  options: Minimap2Options,
}

impl Minimap2Index {
  pub fn new<S, N>(seqs: &[S], names: &[N], args: &Minimap2Args) -> Result<Self, Report>
  where
    S: AsRef<str>,
    N: AsRef<str>,
  {
    assert_eq!(seqs.len(), names.len());
    assert!(seqs.len() > 0);

    let mut options = Minimap2Options::new(args)?;

    let idxopt = options.idx_opt();
    let w: c_int = idxopt.w.into();
    let k: c_int = idxopt.k.into();
    let is_hpc: c_int = (idxopt.flag & MM_I_HPC as c_short).into();
    let bucket_bits: c_int = idxopt.bucket_bits.into();

    let n: c_int = seqs.len().try_into()?;

    let seqs: Vec<CString> = seqs.iter().map(|s| CString::new(s.as_ref()).unwrap()).collect();
    let mut seqs: Vec<*const c_char> = seqs.iter().map(|s| s.as_ptr()).collect();

    let names: Vec<CString> = names.iter().map(|s| CString::new(s.as_ref()).unwrap()).collect();
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

  pub(crate) fn get_mut(&mut self) -> *mut mm_idx_t {
    self.idx
  }

  pub(crate) fn get_ref(&self) -> Result<&mm_idx_t, Report> {
    // SAFETY: dereferencing raw pointer
    unsafe { self.idx.as_ref() }.ok_or_else(|| eyre!("minimap2: index is null"))
  }

  pub(crate) fn get_ref_mut(&mut self) -> Result<&mut mm_idx_t, Report> {
    // SAFETY: dereferencing raw pointer
    unsafe { self.idx.as_mut() }.ok_or_else(|| eyre!("minimap2: index is null"))
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
