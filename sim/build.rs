extern crate bindgen;
extern crate openmp_sys;
use std::process::Command;

use std::path::PathBuf;
fn main() {
    // compile the C static library for running the bit simulation in folder bit_one_step

    cc::Build::new()
        .file("bit_one_step/libbitonestep.c")
        // .opt_level(147483640)
        .compile("libbitonestep.a"); //run with optimization, no difference in output but higher speed.

    // env::var("DEP_OPENMP_FLAG").unwrap().split(" ").for_each(|f| { cc_build.flag(f); });
    // Tell cargo to tell rustc to link the library bitonestep
    println!("cargo:rustc-link-lib=static=bitonestep");
    println!("cargo:rustc-link-search=all=bit-rustc-sim/bit_one_step/");
    println!("cargo:rustc-link-lib=static=gomp");
    println!("cargo:rustc-link-search=all=/usr/lib/gcc/x86_64-linux-gnu/9/");

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    println!("cargo:rerun-if-changed=bitonestep.h");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("bit_one_step/libbitonestep.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from("src/");
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
