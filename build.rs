use std::process::Command;

fn main() {
    Command::new("make").args(&["-C", "lib/hashtree/"]).status().expect("failed to make!");
    println!("cargo:rustc-link-lib=hashtree");
    println!("cargo:rustc-link-search=native=./lib/hashtree/src");

    // Command::new("make").args(&["-C", "lib/XKCP/", "AVX512/libXKCP.a"]).status().expect("failed to make!");
    // println!("cargo:rustc-link-lib=XKCP");
    // println!("cargo:rustc-link-search=native=./lib/XKCP/bin/AVX512/");
}
