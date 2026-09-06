//! The VCF header round trip is idempotent here, which is the half of U.3 this repository owns.
//!
//! `noodles`' `Header` groups records by category on the way out, so a `FILTER` line that preceded
//! an `INFO` line in the source comes back after it (IPNP-BIPN/htsjdk-rs#112, sent upstream as
//! zaeleus/noodles#414). A port that cannot reproduce the input's line order cannot reproduce the
//! input's bytes, so this port keeps the lines in one list in the order they were read and sorts
//! only where htsjdk sorts -- which is `VCFHeader`'s own comparator, decision 0016, and not a
//! grouping by kind.
//!
//! The test is here rather than in the header suite because what it asserts is not a rendering: it
//! is that reading and writing is the identity on the bytes, for a header whose categories are
//! deliberately interleaved.

use htsjdk_vcf::reader::read_vcf;

/// A header whose FILTER, INFO and FORMAT lines are interleaved rather than grouped.
const INTERLEAVED: &str = concat!(
    "##fileformat=VCFv4.2\n",
    "##FILTER=<ID=q10,Description=\"Quality below 10\">\n",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n",
    "##contig=<ID=chr1,length=100000>\n",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
);

/// Read the header and write it back, which is what a tool that copies a file does.
fn round_trip(text: &str) -> String {
    let file = read_vcf(text).expect("the header parses");
    file.header.write()
}

/// The first pass is not the identity, and that is htsjdk's doing rather than a defect: the
/// writer emits `getMetaDataInSortedOrder`, so the lines come out sorted by their rendered text
/// (decision 0016) whatever order they were read in.
#[test]
fn the_first_pass_sorts_because_htsjdk_sorts() {
    let once = round_trip(INTERLEAVED);
    assert_ne!(
        once, INTERLEAVED,
        "the interleaved input is not already sorted"
    );
    let lines: Vec<&str> = once.lines().collect();
    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert_eq!(
        lines[lines.len() - 1],
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    );
}

/// The property U.3 is about: writing what was read reproduces it, so a tool that copies a file
/// twice gets the same bytes both times. `noodles` groups records by CATEGORY on the way out, so a
/// FILTER line that preceded an INFO line comes back after it and the trip is not idempotent --
/// which a byte-identical port cannot absorb.
#[test]
fn and_from_there_it_is_idempotent() {
    let once = round_trip(INTERLEAVED);
    let twice = round_trip(&once);
    assert_eq!(twice, once, "reading what was written reproduces it");
    assert_eq!(round_trip(&twice), twice, "and again, for a third pass");
}

/// The sort key is the whole rendered line, not the kind of line: `##FILTER` sorts before
/// `##FORMAT` and both before `##INFO` because `F` precedes `I`, and two lines of one kind end up
/// adjacent for that reason rather than because they were grouped. What distinguishes the two is
/// the ORDER, which this pins: alphabetical by rendered text, exactly as `VCFHeaderLine.compareTo`
/// does it.
#[test]
fn the_sort_is_by_the_rendered_line() {
    let once = round_trip(INTERLEAVED);
    let meta: Vec<&str> = once
        .lines()
        .filter(|line| line.starts_with("##") && !line.starts_with("##fileformat"))
        .collect();
    let mut sorted = meta.clone();
    sorted.sort();
    assert_eq!(meta, sorted, "the lines come out in the comparator's order");
}
