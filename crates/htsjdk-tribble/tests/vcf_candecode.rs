//! Conformance for `VCFCodec.canDecode` and `VCF3Codec.canDecode`, against the oracle.
//!
//! Golden from `tools/tribble-conformance/VcfCanDecodeDump.java`.
//!
//! The rows that carry the claim, and the reason `-L` has two kinds of codec:
//!
//! ```text
//! candecode  vcf4-bed-extension        .bed            true   false
//! candecode  bed-body                  .bed            false  false
//! candecode  vcf4-bgzf-bed-extension   .bed            true   false
//! candecode  magic-then-junk           .vcf            true   false
//! ```
//!
//! Same extension, opposite answers, decided by the first eighteen bytes. And `magic-then-junk`
//! is `true`: nothing after the magic is read, so a file that is not a VCF at all past its first
//! line is still a Feature file to `FeatureManager`.

use std::io::Read;

use htsjdk_tribble::vcf::{can_decode, can_decode_v3};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_candecode.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

const VCF4: &str = "##fileformat=VCFv4.2\n\
                    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                    chr1\t100\t.\tA\tC\t.\t.\t.\n";
const VCF3: &str = "##fileformat=VCFv3.3\n\
                    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn gzip(body: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(body).expect("gzip");
    encoder.finish().expect("gzip")
}

/// A one-block BGZF member: a gzip member carrying the `BC` extra field.
///
/// The port only needs it to be gzip, which is the point of the row: the reference's second
/// attempt already answers for a block-compressed file.
fn bgzf(body: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut writer = htsjdk_bgzf::BgzfWriter::new(Vec::new());
    writer.write_all(body).expect("bgzf");
    writer.into_inner().expect("bgzf")
}

/// The bytes each label was measured over. `None` is a path the reference could not read at all,
/// which is the absent file and the directory.
fn body(label: &str) -> Option<Vec<u8>> {
    Some(match label {
        "vcf4-plain" | "vcf4-list-extension" | "vcf4-bed-extension" | "vcf4-no-extension" => {
            VCF4.as_bytes().to_vec()
        }
        "magic-only" => b"##fileformat=VCFv4".to_vec(),
        "magic-then-junk" => b"##fileformat=VCFv4NONSENSE".to_vec(),
        "leading-space" => format!(" {VCF4}").into_bytes(),
        "leading-newline" => format!("\n{VCF4}").into_bytes(),
        "vcf40" => b"##fileformat=VCFv4.0\n".to_vec(),
        "vcf43" => b"##fileformat=VCFv4.3\n".to_vec(),
        "vcf3" => VCF3.as_bytes().to_vec(),
        "vcf5" => b"##fileformat=VCFv5.0\n".to_vec(),
        "truncated-magic" => b"##fileformat=V".to_vec(),
        "one-byte" => b"#".to_vec(),
        "empty" => Vec::new(),
        "vcf4-gzip" => gzip(VCF4.as_bytes()),
        "vcf4-bgzf" | "vcf4-bgzf-bed-extension" => bgzf(VCF4.as_bytes()),
        "gzip-not-vcf" => gzip(b"hello there\n"),
        "broken-gzip" => vec![0x1f, 0x8b, 0x08, 0, 0, 0, 0, 0],
        "bed-body" => b"chr1\t0\t10\n".to_vec(),
        "interval-list-body" => b"@HD\tVN:1.6\nchr1\t1\t10\t+\t.\n".to_vec(),
        // A path that cannot be read is an IOException the reference catches, which is the same
        // answer as an empty file by a different route. The port takes bytes, so the caller is
        // what decides this; the rows are kept so the answer is recorded.
        "absent" | "directory" => return None,
        other => panic!("{other} has no body"),
    })
}

#[test]
fn every_file_decodes_or_does_not_as_the_reference_says() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("candecode\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let label = fields.next().expect("a label");
        let _extension = fields.next().expect("an extension");
        let vcf4 = fields.next().expect("the VCF4 answer");
        let vcf3 = fields.next().expect("the VCF3 answer");
        let Some(bytes) = body(label) else {
            // Both codecs answer false for a path they cannot read.
            assert_eq!((vcf4, vcf3), ("false", "false"), "{label}");
            count += 1;
            continue;
        };
        assert_eq!(can_decode(&bytes).to_string(), vcf4, "VCFCodec on {label}");
        assert_eq!(
            can_decode_v3(&bytes).to_string(),
            vcf3,
            "VCF3Codec on {label}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no candecode rows");
    println!("{count} canDecode answers identical");
}

/// The extension is never consulted, which is what separates this codec from the other two.
#[test]
fn the_same_extension_answers_both_ways() {
    assert!(can_decode(VCF4.as_bytes()));
    assert!(!can_decode(b"chr1\t0\t10\n"));
    // And nothing past the magic is read.
    assert!(can_decode(b"##fileformat=VCFv4NONSENSE"));
}
