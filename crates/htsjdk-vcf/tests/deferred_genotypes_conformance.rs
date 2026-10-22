//! Conformance for genotypes parsed on first access (decision 0048), against the oracle.
//!
//! Golden from `tools/vcf-conformance/VcfDeferredGenotypesDump.java`.
//!
//! What the golden pins is *whose* line a late refusal names. `createGenotypeMap` runs with the
//! codec as it stands when something first looks, and `VCFFileReader` reads one record ahead, so a
//! tool touching record N gets record N+1's line number and N+1's columns in the "too many keys"
//! message; a tool that collected every record first gets the last line's. The port reproduces both
//! through [`VcfGenotypeSource::advance`], moved here exactly as each way of reading moves the codec.
//! A record nothing looks at is never parsed, so the copy of a malformed file succeeds, as Picard's
//! `SortVcf` and `MergeVcfs` do.

use std::io::Read;
use std::sync::Arc;

use htsjdk_vcf::genotype_parse::{DeferredVcfGenotypes, VcfGenotypeSource};
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::reader::{read_vcf, read_vcf_with, GenotypeParsing};
use htsjdk_vcf::record_parse::decode_line;
use htsjdk_vcf::variant::VariantContext;
use htsjdk_vcf::vcf_file::write_vcf;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/vcf_deferred_genotypes.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

/// The dump's header, with the sample columns it is given.
fn header(samples: &str) -> String {
    format!(
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
         ##contig=<ID=chr1,length=100000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n"
    )
}

/// The dump's records: a good one after each bad one.
const BODY: &str = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ\t0/1:60\t1/1:70\n\
                    chr1\t200\t.\tC\tG\t50\tPASS\tDP=10\tGT:GQ\t0/1:60:99\t1/1:70\n\
                    chr1\t300\t.\tG\tA\t50\tPASS\tDP=10\tGT:GQ:DP\t0/0:30:5\t0/1:40:6\n\
                    chr1\t400\t.\tT\tC\t50\tPASS\tDP=10\tGT\t0/1\n\
                    chr1\t500\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0/3\t0/1\n\
                    chr1\t600\t.\tC\tT\t50\tPASS\tDP=10\tGT:GQ\t1/1:80\t0/1:90\n";

/// `size()` and `isEmpty()`, which never parse, then a look.
fn touch(out: &mut Vec<String>, scenario: &str, record: usize, vc: &VariantContext) {
    out.push(format!(
        "size\t{scenario}\t{record}\t{}\t{}",
        vc.genotypes.len(),
        vc.genotypes.is_empty()
    ));
    out.push(match vc.genotypes.try_get() {
        Ok(genotypes) => format!("access\t{scenario}\t{record}\tok\t{}", genotypes.len()),
        Err(error) => format!(
            "access\t{scenario}\t{record}\tthrow\t{}\t{}",
            error.class,
            escape(&error.message)
        ),
    });
}

/// The records of the dump's file, decoded line by line with their blocks deferred, each with the
/// counter and the line the codec held after decoding it.
fn decode_deferred(text: &str) -> (Arc<VcfGenotypeSource>, Vec<(VariantContext, usize, &str)>) {
    let header_lines = text.lines().take_while(|l| l.starts_with('#')).count();
    let header_text: String = text
        .lines()
        .take(header_lines)
        .map(|l| format!("{l}\n"))
        .collect();
    let file = read_vcf(&header_text).expect("the header reads");
    let source = VcfGenotypeSource::new(file.header.clone(), file.codec_version);
    let samples: Arc<[String]> = Arc::from(file.header.samples.clone());
    let mut line_number = header_lines;
    let mut records = Vec::new();
    for line in text.lines().skip(header_lines) {
        let decoded = decode_line(line, &file.header, line_number, file.codec_version)
            .expect("the sites decode")
            .expect("a record");
        line_number += 1;
        let mut variant = decoded.variant;
        let block = decoded.genotype_block.expect("genotype columns");
        let parser = DeferredVcfGenotypes::new(
            &source,
            variant.alleles.clone(),
            &variant.contig,
            variant.start,
        );
        variant.genotypes = GenotypesContext::deferred_vcf(block, samples.clone(), parser);
        records.push((variant, line_number, line));
    }
    (source, records)
}

fn port() -> Vec<String> {
    let mut out = Vec::new();
    let sorted = header("s0\ts1") + BODY;

    // stream: `WFIterator` has read the next record before it hands this one out.
    let (source, records) = decode_deferred(&sorted);
    for i in 0..records.len() {
        let (_, line_number, line) = records[(i + 1).min(records.len() - 1)];
        source.advance(line_number, line);
        touch(&mut out, "stream", i, &records[i].0);
    }

    // whole: every record read before any is touched.
    let file = read_vcf_with(&sorted, GenotypeParsing::Deferred).expect("the sites all read");
    for (i, vc) in file.records.iter().enumerate() {
        touch(&mut out, "whole", i, vc);
    }

    // codec: each record touched before the next line is decoded.
    let (source, records) = decode_deferred(&sorted);
    for (i, (vc, line_number, line)) in records.iter().enumerate() {
        source.advance(*line_number, line);
        touch(&mut out, "codec", i, vc);
    }

    // copy: written untouched, so nothing is parsed and the malformed columns are copied.
    let file = read_vcf_with(&sorted, GenotypeParsing::Deferred).expect("the sites all read");
    let written = write_vcf(&file.header, &file.records).expect("nothing to refuse");
    out.push(format!("write\tcopy\tok\t{}", escape(&written)));

    // unsorted: parsed at the read. The suffix is the reader's, which names the file; the library
    // does not know it.
    let unsorted = header("s1\ts0") + BODY;
    let failure = read_vcf_with(&unsorted, GenotypeParsing::Deferred)
        .expect_err("an unsorted file is parsed at the read");
    out.push(format!(
        "read\tunsorted\tthrow\t{}\t{}",
        failure.error.class(),
        escape(&format!(
            "{}, for input source: file://<dir>/unsorted.vcf",
            failure.error.message()
        ))
    ));
    out
}

#[test]
fn late_refusals_name_the_line_htsjdk_names() {
    let golden = golden();
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();
    let actual = port();
    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(got, want, "row {i}");
    }
    assert_eq!(actual.len(), expected.len(), "row count");
}

/// Deferred or not, a well-formed file reads to the same records and writes the same bytes, copied
/// or touched.
#[test]
fn deferral_changes_nothing_for_a_well_formed_file() {
    for samples in ["s0\ts1", "s1\ts0"] {
        let text = header(samples)
            + "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ:DP\t0/1:60:10\t1/1:70:20\n\
               chr1\t200\t.\tC\tG\t50\tPASS\tDP=10\tGT:GQ\t0/0:30\t./.:.\n";
        let eager = read_vcf(&text).expect("reads");
        let deferred = read_vcf_with(&text, GenotypeParsing::Deferred).expect("reads");
        assert_eq!(
            write_vcf(&eager.header, &eager.records).unwrap(),
            write_vcf(&deferred.header, &deferred.records).unwrap(),
            "copied, samples {samples}"
        );
        assert_eq!(eager.records, deferred.records, "samples {samples}");
        for (e, d) in eager.records.iter().zip(&deferred.records) {
            let _ = e.genotypes.first();
            let _ = d.genotypes.try_get().expect("well formed");
        }
        assert_eq!(
            write_vcf(&eager.header, &eager.records).unwrap(),
            write_vcf(&deferred.header, &deferred.records).unwrap(),
            "touched, samples {samples}"
        );
    }
}
