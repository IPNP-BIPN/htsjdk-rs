//! Times one CRAM block codec's decoder on real blocks, for comparing two builds by hand.
//!
//! The input is blocks framed as a little-endian `u32` length and the bytes, which is how
//! `cram-bench` writes its `port_*_compress.bin` files. The checksum of what was decoded is
//! printed so that two builds can be seen to agree.
//!
//! Usage: decode_blocks nx16|tok|fqz|range|rans file reps

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let data = std::fs::read(&args[2]).expect("reads the blocks");
    let reps: usize = args[3].parse().expect("a number of repetitions");
    let mut blocks = Vec::new();
    let mut at = 0;
    while at < data.len() {
        let length = u32::from_le_bytes(data[at..at + 4].try_into().unwrap()) as usize;
        blocks.push(&data[at + 4..at + 4 + length]);
        at += 4 + length;
    }
    let (mut sum, mut bytes) = (0u64, 0usize);
    let start = std::time::Instant::now();
    for _ in 0..reps {
        for block in &blocks {
            let out = match args[1].as_str() {
                "nx16" => htsjdk_cram::rans_nx16::uncompress(block).unwrap(),
                "tok" => htsjdk_cram::name_tokeniser::uncompress(
                    block,
                    htsjdk_cram::name_tokeniser::NAME_SEPARATOR,
                )
                .unwrap(),
                "fqz" => htsjdk_cram::fqzcomp::uncompress(block).unwrap(),
                "range" => htsjdk_cram::range::uncompress(block).unwrap(),
                "rans" => htsjdk_cram::rans::uncompress(block).unwrap(),
                other => panic!("no codec {other}"),
            };
            bytes += out.len();
            sum = sum
                .wrapping_mul(31)
                .wrapping_add(out.iter().map(|&b| u64::from(b)).sum::<u64>());
        }
    }
    let seconds = start.elapsed().as_secs_f64();
    println!(
        "{} blocks, {:.1} MB out, {seconds:.3} s, {:.0} MB/s, checksum {sum}",
        blocks.len(),
        bytes as f64 / 1e6,
        bytes as f64 / 1e6 / seconds
    );
}
