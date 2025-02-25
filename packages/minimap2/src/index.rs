use crate::options::Minimap2Options;
use crate::options_args::Minimap2Args;
use eyre::{Report, eyre};
use minimap2_sys::{MM_I_HPC, mm_idx_destroy, mm_idx_str, mm_idx_t, mm_mapopt_update};
use std::ffi::{CString, c_char};
use std::os::raw::{c_int, c_short};
use std::sync::Arc;
use std::sync::atomic::AtomicPtr;

#[derive(Debug)]
pub struct Minimap2Index {
  idx: Arc<AtomicPtr<mm_idx_t>>,
  options: Minimap2Options,
}

impl Minimap2Index {
  pub fn new<S, N>(seqs: &[S], names: &[N], args: &Minimap2Args) -> Result<Self, Report>
  where
    S: AsRef<str>,
    N: AsRef<str>,
  {
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

    // SAFETY: Explained by the correct construction of inputs and non-overlapping data.
    let idx = unsafe { mm_idx_str(w, k, is_hpc, bucket_bits, n, seqs.as_mut_ptr(), names.as_mut_ptr()) };
    if idx.is_null() {
      return Err(eyre!("minimap2: failed to create index"));
    }

    // SAFETY: Ensured that `idx` is non-null and correctly initialized.
    unsafe {
      mm_mapopt_update(options.map_opt_mut(), idx);
    }

    Ok(Self {
      idx: Arc::new(AtomicPtr::new(idx)),
      options,
    })
  }

  pub fn options(&self) -> &Minimap2Options {
    &self.options
  }

  pub fn get(&self) -> *const mm_idx_t {
    self.idx.load(std::sync::atomic::Ordering::SeqCst)
  }

  pub fn get_ref(&self) -> Result<&mm_idx_t, Report> {
    let ptr = self.idx.load(std::sync::atomic::Ordering::SeqCst);
    // SAFETY: Ensured that `ptr` points to a valid `mm_idx_t`.
    unsafe { ptr.as_ref() }.ok_or_else(|| eyre!("minimap2: index is null"))
  }
}

impl Drop for Minimap2Index {
  fn drop(&mut self) {
    let ptr = self.idx.load(std::sync::atomic::Ordering::SeqCst);
    // SAFETY: Ensured that `ptr` is a valid pointer on allocation and not used after being freed.
    unsafe { mm_idx_destroy(ptr) };
  }
}
