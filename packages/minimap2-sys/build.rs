use bindgen::Formatter;
use std::env;
use std::path::PathBuf;

// Configure for minimap2
fn configure(cc: &mut cc::Build) {
  println!("cargo:rerun-if-changed=minimap2/*.c");

  cc.include("minimap2");
  cc.opt_level(2);

  #[cfg(feature = "sse2only")]
  sse2only(cc);

  #[cfg(feature = "simde")]
  simde(cc);

  // Include ksw2.h kalloc.h
  cc.include("minimap2/");
  // cc.include("minimap2/");

  let files: Vec<_> = std::fs::read_dir("minimap2")
    .unwrap()
    .map(|f| f.unwrap().path())
    .collect();

  assert_ne!(
    files.len(),
    0,
    "No files found in minimap2 directory -- Did you forget to clone the submodule? git submodule init --recursive"
  );

  for file in files {
    // Skip "main.c" and "example.c"
    if file.file_name().unwrap() == "main.c" || file.file_name().unwrap() == "example.c" {
      continue;
    }

    // Ignore all "neon"
    if file.file_name().unwrap().to_str().unwrap().contains("neon") {
      continue;
    }

    // Ignore all "ksw"
    if file.file_name().unwrap().to_str().unwrap().contains("ksw") {
      continue;
    }

    if let Some(x) = file.extension() {
      if x == "c" {
        cc.file(file);
      }
    }
  }

  cc.file("minimap2/ksw2_ll_sse.c");

  #[cfg(not(feature = "noopt"))]
  target_specific(cc);
}

#[cfg(not(feature = "noopt"))]
fn target_specific(cc: &mut cc::Build) {
  let target = env::var("TARGET").unwrap_or_default();

  if target.contains("aarch64") | target.contains("arm") {
    cc.include("minimap2/sse2neon/");
    // For aarch64 targets with neon
    // Add the following files:
    // ksw2_extz2_neon.o ksw2_extd2_neon.o ksw2_exts2_neon.o
    cc.file("minimap2/ksw2_extz2_sse.c");
    cc.file("minimap2/ksw2_extd2_sse.c");
    cc.file("minimap2/ksw2_exts2_sse.c");
    cc.flag("-DKSW_SSE2_ONLY");

    // CFLAGS+=-D_FILE_OFFSET_BITS=64 -fsigned-char -Isse2neon -D__SSE2__
    cc.flag("-D_FILE_OFFSET_BITS=64");
    cc.flag("-fsigned-char");
    cc.flag("-Isse2neon");
    cc.flag("-D__SSE2__");
  } else if target.contains("x86_64") {
    #[cfg(all(target_feature = "sse4.1", not(feature = "simde"), not(feature = "sse2only")))]
    cc.flag("-msse4.1");

    #[cfg(all(not(target_feature = "sse4.1"), target_feature = "sse2",))]
    {
      cc.flag("-msse2");
    }

    #[cfg(all(not(target_feature = "sse4.1"), target_feature = "sse2"))]
    cc.flag("-DKSW_SSE2_ONLY");

    // #[cfg(all(not(target_feature = "sse4.1"), target_feature = "sse2"))]
    // cc.flag("-DKSW_CPU_DISPATCH");

    #[cfg(all(not(target_feature = "sse4.1"), target_feature = "sse2",))]
    cc.flag("-mno-sse4.1");

    // OBJS+=ksw2_extz2_sse41.o ksw2_extd2_sse41.o ksw2_exts2_sse41.o ksw2_extz2_sse2.o ksw2_extd2_sse2.o ksw2_exts2_sse2.o ksw2_dispatch.o
    cc.file("minimap2/ksw2_extz2_sse.c");
    cc.file("minimap2/ksw2_extd2_sse.c");
    cc.file("minimap2/ksw2_exts2_sse.c");
    // cc.file("minimap2/ksw2_dispatch.c");
  }
}

#[cfg(feature = "simde")]
fn simde(cc: &mut cc::Build) {
  cc.include("minimap2/lib/simde");
  cc.flag("-DSIMDE_ENABLE_NATIVE_ALIASES");
  cc.flag("-DUSE_SIMDE");
  cc.flag("-std=c99");
}

fn compile() {
  let out_path = PathBuf::from(env::var_os("OUT_DIR").unwrap());

  let _host = env::var("HOST").unwrap();
  let _target = env::var("TARGET").unwrap();

  println!("{}", _target);

  println!("cargo:rerun-if-env-changed=PKG_CONFIG_SYSROOT_DIR");

  println!("cargo:rustc-link-lib=m");

  let mut cc = cc::Build::new();

  cc.warnings(false);
  cc.flag("-Wc++-compat");
  cc.out_dir(&out_path);

  configure(&mut cc);

  cc.flag("-DHAVE_KALLOC");

  #[cfg(feature = "static")]
  if _target.ends_with("windows-gnu") {
    let mingw_lib_dir = env::var("MINGW_LIB_DIR").unwrap_or("/usr/x86_64-w64-mingw32/lib".to_owned());
    println!("cargo:rustc-link-search=native={mingw_lib_dir}");
    println!("cargo:rustc-link-lib=static=pthread");
  } else {
    println!("cargo:rustc-link-lib=pthread");
  }

  if cfg!(feature = "static") {
    cc.static_flag(true);

    let libz_include = env::var_os("DEP_Z_INCLUDE").unwrap();
    cc.include(libz_include);

    let libz_root = env::var_os("DEP_Z_ROOT").unwrap();
    let libz_lib = PathBuf::from(&libz_root).join("lib");
    let libz_lib = libz_lib.to_str().unwrap();
    println!("cargo:rustc-link-search={libz_lib}");
    println!("cargo:rustc-link-lib=static=z");
  } else {
    println!("cargo:rustc-link-lib=z");
  }

  cc.compile("libminimap");
}

#[cfg(feature = "sse2only")]
fn sse2only(cc: &mut cc::Build) {
  #[cfg(all(target_feature = "sse2", not(target_feature = "sse4.1")))]
  cc.flag("-DKSW_SSE2_ONLY");

  #[cfg(all(target_feature = "sse2", not(target_feature = "sse4.1")))]
  cc.flag("-mno-sse4.1");
  let target = env::var("TARGET").unwrap_or_default();
  if target.contains("x86_64") {
    #[cfg(all(
      not(target_feature = "sse4.1"),
      target_feature = "sse2",
      not(target_arch = "aarch64")
    ))]
    {
      cc.flag("-msse2");
    }
  }
}

#[cfg(feature = "bindgen")]
fn gen_bindings() {
  let out_path = PathBuf::from(env::var_os("OUT_DIR").unwrap());

  let bigings = bindgen::Builder::default()
    .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
    .formatter(Formatter::Rustfmt)
    .header("minimap2.h")
    .generate_cstr(true)
    .derive_debug(true)
    .generate_comments(true)
    .generate()
    .expect("Couldn't write bindings!");

  bigings
    .write_to_file(out_path.join("bindings.rs"))
    .expect("Unable to create bindings");

  bigings
    .write_to_file("src/bindings.rs")
    .expect("Unable to create bindings");
}

#[cfg(not(feature = "bindgen"))]
fn gen_bindings() {}

fn main() {
  compile();
  gen_bindings();
}
