use flate2::{Compress, Compression, FlushCompress};
use md5::{Digest, Md5};

fn hex(b: &[u8]) -> String { b.iter().map(|x| format!("{:02x}", x)).collect() }

fn main() {
    let mut input = vec![0u8; 65536];
    let mut s: u64 = 12345;
    for i in 0..input.len() {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        input[i] = (s >> 58) as u8;
    }
    let mut h = Md5::new(); h.update(&input);
    println!("input  len={} md5={}", input.len(), hex(&h.finalize()));

    for lvl in [1u32, 5, 6, 9] {
        // false = raw deflate, no zlib header: matches Java's new Deflater(lvl, true)
        let mut c = Compress::new(Compression::new(lvl), false);
        let mut out = Vec::with_capacity(input.len() * 2);
        c.compress_vec(&input, &mut out, FlushCompress::Finish).unwrap();
        let mut h = Md5::new(); h.update(&out);
        println!("level={} outlen={} md5={}", lvl, out.len(), hex(&h.finalize()));
    }
}
