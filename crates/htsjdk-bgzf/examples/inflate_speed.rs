//! One BGZF file through `decompress_all` or through zlib block by block, for a cycle counter:
//! `RAYON_NUM_THREADS=1 /usr/bin/time -l inflate_speed port file` against `... zlib file`.

use flate2::{Decompress, FlushDecompress};
use htsjdk_bgzf::{decompress_all, BLOCK_FOOTER_LENGTH, BLOCK_HEADER_LENGTH};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let bytes = std::fs::read(&args[1]).expect("reads");
    let total = if args[0] == "zlib" {
        let mut zlib = Decompress::new(false);
        let mut out = vec![0u8; 1 << 16];
        let (mut at, mut total) = (0, 0usize);
        while at < bytes.len() {
            let length = u16::from_le_bytes([bytes[at + 16], bytes[at + 17]]) as usize + 1;
            let body = &bytes[at + BLOCK_HEADER_LENGTH..at + length - BLOCK_FOOTER_LENGTH];
            zlib.reset(false);
            zlib.decompress(body, &mut out, FlushDecompress::Finish)
                .expect("inflates");
            total += zlib.total_out() as usize;
            at += length;
        }
        total
    } else {
        decompress_all(&bytes).expect("inflates").len()
    };
    println!("{total}");
}
