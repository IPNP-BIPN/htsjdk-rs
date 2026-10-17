//! Conformance for BCF2, written and read, against the oracle.
//!
//! Golden from `tools/vcf-conformance/Bcf2Dump.java`, which runs `VcfFormatConverter`'s three
//! conversions on each case and prints the case's VCF text beside the outputs, so this replays
//! exactly what was measured:
//!
//! ```text
//! in    <case>  <VCF text>
//! bcf   <case>  <the BCF written from it>
//! vcf   <case>  <that BCF read and written as VCF>
//! bcf2  <case>  <that BCF read and written as BCF>
//! read  <case>  <a damaged BCF>
//! err   <case>  <stage>  <class>  <message>
//! ```
//!
//! Every byte is compared: the INFO block's hash order, the GQ cap, the copied genotype blocks
//! and the refusals' messages are all in it.

use std::collections::HashMap;
use std::io::Read;

use htsjdk_vcf::bcf2::codec::read_bcf;
use htsjdk_vcf::bcf2::writer::{Bcf2Writer, Origin};
use htsjdk_vcf::bcf2::Bcf2Error;
use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::reader::read_vcf;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bcf2.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).expect("hex"))
        .collect()
}

/// What one stage produced: its bytes, or the exception.
#[derive(Debug, PartialEq, Eq)]
enum Outcome {
    Bytes(Vec<u8>),
    Refused(String, String),
}

fn refused(e: Bcf2Error) -> Outcome {
    Outcome::Refused(e.class.to_string(), e.message)
}

fn to_bcf(text: &str) -> Outcome {
    let file = match read_vcf(text) {
        Ok(file) => file,
        Err(f) => return Outcome::Refused(f.error.class().into(), f.error.message()),
    };
    let mut writer = match Bcf2Writer::new(&file.header, false) {
        Ok(w) => w,
        Err(e) => return refused(e),
    };
    let mut out = writer.header_bytes();
    for vc in &file.records {
        match writer.encode(vc, Origin::VcfCodec) {
            Ok(bytes) => out.extend(bytes),
            Err(e) => return refused(e),
        }
    }
    Outcome::Bytes(out)
}

/// Read a BCF and write it as VCF (`bcf2 == false`) or BCF again.
fn from_bcf(bytes: &[u8], bcf2: bool, source: &str) -> Outcome {
    let file = match read_bcf(bytes) {
        Ok(file) => file,
        Err(failure) => return refused(failure.reader_exception(source)),
    };
    if bcf2 {
        let mut writer = match Bcf2Writer::new(&file.header, false) {
            Ok(w) => w,
            Err(e) => return refused(e),
        };
        let mut out = writer.header_bytes();
        for vc in &file.records {
            match writer.encode(vc, Origin::Bcf2Codec) {
                Ok(b) => out.extend(b),
                Err(e) => return refused(e),
            }
        }
        Outcome::Bytes(out)
    } else {
        let encoder = VcfEncoder::new(&file.header);
        let mut out = file.header.write();
        for vc in &file.records {
            encoder.encode_into(vc, &mut out).expect("encodes");
            out.push('\n');
        }
        Outcome::Bytes(out.into_bytes())
    }
}

#[test]
fn bcf2_matches_the_oracle() {
    let text = golden();
    let mut expected: HashMap<(String, String), Outcome> = HashMap::new();
    let mut inputs: Vec<(String, String, Vec<u8>)> = Vec::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let parts: Vec<&str> = line.split('\t').collect();
        match parts[0] {
            "in" | "read" => inputs.push((parts[0].into(), parts[1].into(), unhex(parts[2]))),
            "err" => {
                expected.insert(
                    (parts[1].into(), parts[2].into()),
                    Outcome::Refused(parts[3].into(), parts[4].into()),
                );
            }
            stage => {
                expected.insert(
                    (parts[1].into(), stage.into()),
                    Outcome::Bytes(unhex(parts[2])),
                );
            }
        }
    }
    assert!(inputs.len() >= 20, "{} inputs", inputs.len());

    // The reader names the file it failed on, and the path is the run's, not the file's.
    let source_of = |case: &str, stage: &str| -> String {
        match expected.get(&(case.to_string(), stage.to_string())) {
            Some(Outcome::Refused(_, message)) => message
                .rsplit_once(", for input source: ")
                .map(|(_, s)| s.to_string())
                .unwrap_or_default(),
            _ => String::new(),
        }
    };

    let mut failures = Vec::new();
    let mut check = |case: &str, stage: &str, got: Outcome| {
        let want = expected.get(&(case.to_string(), stage.to_string()));
        if want != Some(&got) {
            let show = |o: Option<&Outcome>| match o {
                Some(Outcome::Bytes(b)) => format!("{} bytes", b.len()),
                Some(Outcome::Refused(c, m)) => format!("{c}: {m}"),
                None => "nothing".to_string(),
            };
            failures.push(format!(
                "{case}/{stage}: oracle {} port {}",
                show(want),
                show(Some(&got))
            ));
        }
    };

    for (kind, case, bytes) in &inputs {
        if kind == "read" {
            check(case, "vcf", from_bcf(bytes, false, &source_of(case, "vcf")));
            continue;
        }
        let text = String::from_utf8(bytes.clone()).expect("utf-8");
        let bcf = to_bcf(&text);
        let Outcome::Bytes(bcf_bytes) = &bcf else {
            check(case, "bcf", bcf);
            continue;
        };
        let bcf_bytes = bcf_bytes.clone();
        check(case, "bcf", bcf);
        check(case, "vcf", from_bcf(&bcf_bytes, false, ""));
        check(case, "bcf2", from_bcf(&bcf_bytes, true, ""));
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
