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

/// Every BCF the golden reads, eager and deferred (decision 0048): the same records and the same
/// bytes for a file that reads, and for one that does not, the same refusal, at the read or from
/// the first record whose genotypes are looked at.
#[test]
fn deferred_bcf_reads_refuse_what_eager_reads_refuse() {
    use htsjdk_vcf::bcf2::codec::read_bcf_with;
    use htsjdk_vcf::reader::GenotypeParsing;

    let text = golden();
    let mut files: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let parts: Vec<&str> = line.split('\t').collect();
        match parts[0] {
            "read" => files.push((parts[1].into(), unhex(parts[2]))),
            "in" => {
                let vcf = String::from_utf8(unhex(parts[2])).expect("utf-8");
                if let Outcome::Bytes(bcf) = to_bcf(&vcf) {
                    files.push((parts[1].into(), bcf));
                }
            }
            _ => {}
        }
    }
    // None of the golden's damaged files breaks a genotype block, so one is made here: a BCF with
    // sorted samples whose first record's first FORMAT key is typed CHAR rather than an integer
    // offset, which `getDictionaryString` cannot cast.
    let vcf = "##fileformat=VCFv4.2\n\
               ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
               ##contig=<ID=chr1,length=1000>\n\
               #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
               chr1\t10\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1\n\
               chr1\t20\t.\tC\tG\t50\tPASS\t.\tGT\t0/0\t0/1\n";
    let Outcome::Bytes(mut damaged) = to_bcf(vcf) else {
        panic!("the made-up VCF converts")
    };
    let le = |at: usize| i32::from_le_bytes(damaged[at..at + 4].try_into().unwrap()) as usize;
    let first_record = 5 + 4 + le(5);
    let genotypes_at = first_record + 8 + le(first_record);
    damaged[genotypes_at] = 0x17;
    files.push(("damaged-genotype-key".into(), damaged));

    let mut refused_late = 0;
    for (case, bytes) in &files {
        let eager = read_bcf(bytes);
        let deferred = read_bcf_with(bytes, GenotypeParsing::Deferred);
        match (eager, deferred) {
            (Ok(eager), Ok(deferred)) => {
                assert_eq!(eager.records, deferred.records, "{case}");
                for (e, d) in eager.records.iter().zip(&deferred.records) {
                    let mut we = Bcf2Writer::new(&eager.header, false).expect("a writer");
                    let mut wd = Bcf2Writer::new(&deferred.header, false).expect("a writer");
                    assert_eq!(
                        we.encode(e, Origin::Bcf2Codec),
                        wd.encode(d, Origin::Bcf2Codec),
                        "{case}"
                    );
                }
            }
            (Err(eager), Err(deferred)) => {
                assert_eq!(eager.error, deferred.error, "{case}");
                assert_eq!(eager.records, deferred.records, "{case}");
            }
            (Err(eager), Ok(deferred)) => {
                let at = eager.records.len();
                assert_eq!(&deferred.records[..at], &eager.records[..], "{case}");
                let error = deferred.records[at]
                    .genotypes
                    .try_get()
                    .expect_err("the record the eager read refused");
                assert_eq!(
                    (error.class, error.message.as_str()),
                    (eager.error.class, eager.error.message.as_str()),
                    "{case}"
                );
                refused_late += 1;
            }
            (Ok(_), Err(deferred)) => {
                panic!("{case}: only the deferred read refused: {deferred:?}")
            }
        }
    }
    assert!(files.len() >= 15, "{} files", files.len());
    eprintln!(
        "{} files, {refused_late} refused on first access",
        files.len()
    );
}
