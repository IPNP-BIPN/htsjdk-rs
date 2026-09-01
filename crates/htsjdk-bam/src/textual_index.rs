//! `htsjdk.samtools.TextualBAMIndexWriter`: a `.bai` printed as text.
//!
//! The only place a BAM index is written as anything but bytes, and the format is a set of
//! decisions nobody would make twice the same way. Three of them decide whether a port agrees:
//!
//! * **An empty reference is spaced differently.** `writeNullContent` prints `n_bin=0` and
//!   `n_intv=0` with no space after the `=`, where every other line prints `n_bin= 4` with one. A
//!   port that formatted both the same way differs on any file with an unused contig.
//! * **The metadata bin is counted with the others and printed apart.** `n_bin` is the real bins
//!   plus one whenever the pseudo-bin is present, and the pseudo-bin is then printed after every
//!   real bin, out of numeric order, always claiming two chunks: the first pair is offsets, the
//!   second is the aligned and unaligned record counts printed in the same hexadecimal.
//! * **The numbers are `Long.toString(value, 16)`**, which is a signed hexadecimal: a negative
//!   count prints with a leading minus rather than as two's complement.
//!
//! `gatk-rs` ported this inside `PrintFileDiagnostics`, the tool that reaches it. The format is
//! htsjdk's.

use crate::bin::{bin_summary_string, MAX_BINS};
use crate::index::BamIndex;

/// `Long.toString(value, 16)`: signed hexadecimal, minus sign and all.
fn hex(value: u64) -> String {
    let signed = value as i64;
    if signed < 0 {
        format!("-{:x}", signed.unsigned_abs())
    } else {
        format!("{signed:x}")
    }
}

/// `Chunk.toString`: `blockAddress:blockOffset`, twice, joined by a dash.
fn chunk_string(pointer: u64) -> String {
    format!("{}:{}", pointer >> 16, pointer & 0xFFFF)
}

/// `BlockCompressedFilePointerUtil.asAddressOffsetString`.
fn address_offset(pointer: u64) -> String {
    format!("{}:{}", pointer >> 16, pointer & 0xFFFF)
}

/// `TextualBAMIndexWriter` over an index already read.
pub fn render(index: &BamIndex) -> String {
    let mut out = String::new();
    out.push_str(&format!("n_ref={}\n", index.references.len()));
    for (reference, content) in index.references.iter().enumerate() {
        let bins: Vec<_> = content
            .bins
            .iter()
            .filter(|bin| bin.bin_number != MAX_BINS)
            .collect();
        if bins.is_empty() {
            // `writeNullContent`, whose two lines have no space after the `=`.
            out.push_str(&format!("Reference {reference} has n_bin=0\n"));
            out.push_str(&format!("Reference {reference} has n_intv=0\n"));
            continue;
        }
        let counted = bins.len() + usize::from(content.metadata.is_some());
        out.push_str(&format!("Reference {reference} has n_bin= {counted}\n"));
        for bin in &bins {
            out.push_str(&format!(
                "  Ref {reference} bin {} ({}) has n_chunk= {}\n",
                bin.bin_number,
                bin_summary_string(bin.bin_number),
                bin.chunks.len()
            ));
            if bin.chunks.is_empty() {
                out.push('\n');
            }
            for chunk in &bin.chunks {
                out.push_str(&format!(
                    "     Chunk: {}-{} start: {} end: {}\n",
                    chunk_string(chunk.start),
                    chunk_string(chunk.end),
                    hex(chunk.start),
                    hex(chunk.end)
                ));
            }
        }

        // `writeChunkMetaData`: always bin 37450, always two chunks when the metadata is there, and
        // the "Chunk:" label is followed by two spaces because the chunk itself is not printed.
        match content.metadata {
            None => {
                out.push_str(&format!("  Ref {reference} bin 37450 has n_chunk= 0\n"));
                out.push('\n');
            }
            Some(metadata) => {
                out.push_str(&format!("  Ref {reference} bin 37450 has n_chunk= 2\n"));
                out.push_str(&format!(
                    "     Chunk:  start: {} end: {}\n",
                    hex(metadata.first_offset),
                    hex(metadata.last_offset)
                ));
                out.push_str(&format!(
                    "     Chunk:  start: {} end: {}\n",
                    hex(metadata.aligned as u64),
                    hex(metadata.unaligned as u64)
                ));
            }
        }

        if content.linear_index.is_empty() {
            out.push_str(&format!("Reference {reference} has n_intv= 0\n"));
            continue;
        }
        out.push_str(&format!(
            "Reference {reference} has n_intv= {}\n",
            content.linear_index.len()
        ));
        for (window, entry) in content.linear_index.iter().enumerate() {
            // A zero entry is skipped rather than printed, so the windows in the output are sparse.
            if *entry != 0 {
                out.push_str(&format!(
                    "  Ref {reference} ioffset for {window} is {}\n",
                    address_offset(*entry)
                ));
            }
        }
    }
    out.push_str(&format!(
        "No Coordinate Count={}\n",
        index.no_coordinate_records
    ));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{Bin, Chunk, PseudoBin, ReferenceContent};

    fn index() -> BamIndex {
        BamIndex {
            references: vec![
                ReferenceContent {
                    bins: vec![Bin {
                        bin_number: 4681,
                        chunks: vec![Chunk {
                            start: 0x1234_0010,
                            end: 0x5678_0020,
                        }],
                    }],
                    metadata: Some(PseudoBin {
                        first_offset: 0x1234_0010,
                        last_offset: 0x5678_0020,
                        aligned: 100,
                        unaligned: 3,
                    }),
                    linear_index: vec![0, 0x1234_0010],
                },
                ReferenceContent {
                    bins: Vec::new(),
                    metadata: None,
                    linear_index: Vec::new(),
                },
            ],
            no_coordinate_records: 7,
        }
    }

    #[test]
    fn an_empty_reference_is_spaced_differently() {
        let text = render(&index());
        assert!(text.contains("Reference 1 has n_bin=0\n"), "{text}");
        assert!(text.contains("Reference 0 has n_bin= 2\n"), "{text}");
    }

    #[test]
    fn the_metadata_bin_is_counted_and_printed_apart() {
        let text = render(&index());
        // One real bin plus the pseudo-bin makes two, and 37450 comes after the real ones.
        let bin_line = text.lines().position(|l| l.contains("bin 4681")).unwrap();
        let meta_line = text.lines().position(|l| l.contains("bin 37450")).unwrap();
        assert!(meta_line > bin_line);
        assert!(text.contains("     Chunk:  start: 64 end: 3\n"), "{text}");
    }

    #[test]
    fn a_zero_linear_entry_is_skipped_rather_than_printed() {
        let text = render(&index());
        assert!(!text.contains("ioffset for 0 "), "{text}");
        assert!(text.contains("ioffset for 1 is 4660:16\n"), "{text}");
    }

    #[test]
    fn the_counts_are_signed_hexadecimal() {
        assert_eq!(hex(255), "ff");
        assert_eq!(hex(u64::MAX), "-1");
    }
}
