//! `jmath::to_string` against the pinned JVM's `Double.toString` and `Float.toString`.
//!
//! The sweep comes from `tools/jmath-conformance/ToStringSweep.java`, stratified over the families
//! decision 0046 measured. It is generated in the oracle container rather than committed, as the
//! `%f` / `%e` sweep is; CI runs this with `TOSTRING_SWEEP` pointing at the fresh file.

use std::io::{BufRead, BufReader};

#[test]
#[ignore = "needs TOSTRING_SWEEP from the pinned oracle"]
fn every_line_matches_the_jvm() {
    let path = std::env::var("TOSTRING_SWEEP").expect("TOSTRING_SWEEP");
    let reader = BufReader::new(std::fs::File::open(&path).expect("sweep file"));
    let (mut doubles, mut floats) = (0u64, 0u64);
    let mut mismatches = Vec::new();
    for line in reader.lines() {
        let line = line.expect("utf-8 line");
        let mut fields = line.split('\t');
        let (Some(kind), Some(hex), Some(jvm)) = (fields.next(), fields.next(), fields.next())
        else {
            continue;
        };
        let bits = u64::from_str_radix(hex, 16).expect("hex bits");
        let ours = match kind {
            "d" => {
                doubles += 1;
                jmath::to_string::double_to_string(f64::from_bits(bits))
            }
            "f" => {
                floats += 1;
                jmath::to_string::float_to_string(f32::from_bits(bits as u32))
            }
            other => panic!("unknown kind {other:?}"),
        };
        if ours != jvm && mismatches.len() < 20 {
            mismatches.push(format!("{kind} {hex}: ours={ours} jvm={jvm}"));
        }
    }
    println!("doubles={doubles} floats={floats}");
    assert!(
        doubles > 1_500_000 && floats > 800_000,
        "the sweep is incomplete"
    );
    assert!(
        mismatches.is_empty(),
        "diverging lines:\n{}",
        mismatches.join("\n")
    );
}
