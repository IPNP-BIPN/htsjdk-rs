//! One thread, BGZF-sized blocks: this crate's deflater against the C zlib, on files given as
//! arguments, with no checks beyond the sizes. The bytes are `against_zlib`'s business; this is for
//! a cycle counter: `/usr/bin/time -l jdk_speed zlib 5 file` against `... port 5 file`.

use flate2::{Compress, Compression, FlushCompress};
use gkl_deflate::{Flavour, ReusableDeflater};

const BLOCK: usize = 65498;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args[0] == "check" {
        return check(&args[1..]);
    }
    let level: usize = args[1].parse().expect("a level");
    let data = std::fs::read(&args[2]).expect("reads");
    let data = &data[..data.len().min(32 << 20)];
    let mut out = Vec::with_capacity(1 << 17);
    let mut total = 0usize;
    if args[0] == "zlib" {
        let mut z = Compress::new(Compression::new(level as u32), false);
        for block in data.chunks(BLOCK) {
            z.reset();
            out.clear();
            z.compress_vec(block, &mut out, FlushCompress::Finish)
                .unwrap();
            total += out.len();
        }
    } else {
        let flavour = if args[0] == "gkl" {
            Flavour::Gkl { sse42: true }
        } else {
            Flavour::Jdk
        };
        let mut state = ReusableDeflater::new(level, flavour);
        for block in data.chunks(BLOCK) {
            state.deflate_into(block, &mut out);
            total += out.len();
        }
    }
    println!("{total}");
}

/// Every level against the C zlib on the given files, block by block through one reused state.
fn check(paths: &[String]) {
    for path in paths {
        let data = std::fs::read(path).expect("reads");
        let data = &data[..data.len().min(32 << 20)];
        for level in 1..=9 {
            let mut z = Compress::new(Compression::new(level as u32), false);
            let mut state = ReusableDeflater::new(level, Flavour::Jdk);
            let (mut expected, mut out) = (Vec::with_capacity(1 << 17), Vec::new());
            for block in data.chunks(BLOCK) {
                z.reset();
                expected.clear();
                z.compress_vec(block, &mut expected, FlushCompress::Finish)
                    .unwrap();
                state.deflate_into(block, &mut out);
                assert_eq!(expected, out, "{path} level {level}");
            }
            println!("{path} level {level}: identical");
        }
    }
}
