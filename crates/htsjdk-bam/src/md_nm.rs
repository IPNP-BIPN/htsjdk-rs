//! The `MD` and `NM` alignment tags.
//!
//! Ports `htsjdk.samtools.util.SequenceUtil.calculateMdAndNmTags` at tag 4.2.0: given a read's
//! alignment and the reference bases it covers, compute the `MD` string (the mismatch/deletion
//! description) and the `NM` edit distance. `SetNmMdAndUqTags`, `ValidateSamFile`, and the alignment
//! mergers all lean on this; it is the reference-walking core they share.
//!
//! The walk is faithful to htsjdk down to its quirks: a base compares equal via the IUPAC-aware
//! [`bases_equal`](crate::sequence::bases_equal) table (a no-call read base, byte `0`, counts as a
//! match); a run past the end of the reference stops the walk mid-block; and a deletion emits `^`
//! then the deleted reference bases. `NM` counts mismatches plus every inserted and deleted base.

use crate::cigar::{Cigar, Op};
use crate::sequence::bases_equal;

/// `calculateMdAndNmTags(record, ref, true, true)`: the `MD` string and the `NM` count.
///
/// `alignment_start` is 1-based (as a record stores it); `seq` is the read bases; `ref_bases` is the
/// reference contig the read is aligned to.
pub fn calculate_md_and_nm(
    alignment_start: i32,
    cigar: &Cigar,
    seq: &[u8],
    ref_bases: &[u8],
) -> (String, i32) {
    let mut md = String::new();
    let mut match_count: i32 = 0;
    let mut nm_count: i32 = 0;

    // 0-based reference position of the current block, and 0-based read offset.
    let mut block_ref_pos = (alignment_start - 1) as isize;
    let mut block_read_start: usize = 0;

    for ce in &cigar.elements {
        let block_length = ce.length as usize;
        match ce.op {
            Op::M | Op::Eq | Op::X => {
                let mut in_block = 0;
                while in_block < block_length {
                    let read_offset = block_read_start + in_block;
                    let ref_idx = block_ref_pos + in_block as isize;
                    if ref_idx < 0 || ref_bases.len() as isize <= ref_idx {
                        break; // ran off the reference
                    }
                    let read_base = seq[read_offset];
                    let ref_base = ref_bases[ref_idx as usize];
                    if bases_equal(read_base, ref_base) || read_base == 0 {
                        match_count += 1;
                    } else {
                        md.push_str(&match_count.to_string());
                        md.push(ref_base as char);
                        match_count = 0;
                        nm_count += 1;
                    }
                    in_block += 1;
                }
                if in_block < block_length {
                    break;
                }
                block_ref_pos += block_length as isize;
                block_read_start += block_length;
            }
            Op::D => {
                md.push_str(&match_count.to_string());
                md.push('^');
                let mut in_block = 0;
                while in_block < block_length {
                    let ref_idx = block_ref_pos + in_block as isize;
                    if ref_idx < 0 || ref_idx as usize >= ref_bases.len() {
                        break;
                    }
                    let ref_base = ref_bases[ref_idx as usize];
                    if ref_base == 0 {
                        break;
                    }
                    md.push(ref_base as char);
                    in_block += 1;
                }
                match_count = 0;
                if in_block < block_length {
                    break;
                }
                block_ref_pos += block_length as isize;
                nm_count += block_length as i32;
            }
            Op::I | Op::S => {
                block_read_start += block_length;
                if ce.op == Op::I {
                    nm_count += block_length as i32;
                }
            }
            Op::N => {
                block_ref_pos += block_length as isize;
            }
            Op::H | Op::P => {} // consume neither read nor reference
        }
    }
    md.push_str(&match_count.to_string());

    (md, nm_count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::Cigar;

    const REF: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    /// Values dumped from htsjdk's `SequenceUtil.calculateMdAndNmTags` at tag 4.2.0, covering a
    /// perfect match, single and double mismatches, a deletion, an insertion, a soft clip, and a
    /// reference skip.
    #[test]
    fn md_and_nm_match_htsjdk() {
        let cases: &[(i32, &str, &[u8], &str, i32)] = &[
            (1, "8M", b"ACGTACGT", "8", 0),
            (1, "8M", b"ACCTACGT", "2G5", 1),
            (1, "4M2D4M", b"ACGTACGT", "4^AC0G0T0A0C0", 6),
            (1, "4M2I4M", b"ACGTTTACGT", "8", 2),
            (1, "2S6M", b"TTGTACGT", "0A0C0G0T0A0C0", 6),
            (1, "4M3N4M", b"ACGTACGT", "4T0A0C0G0", 4),
            (5, "8M", b"AAGTACGA", "1C5T0", 2),
        ];
        for &(start, cigar_text, seq, expect_md, expect_nm) in cases {
            let cigar = parse_cigar(cigar_text);
            let (md, nm) = calculate_md_and_nm(start, &cigar, seq, REF);
            assert_eq!(
                md,
                expect_md,
                "MD for {cigar_text} {:?}",
                std::str::from_utf8(seq)
            );
            assert_eq!(
                nm,
                expect_nm,
                "NM for {cigar_text} {:?}",
                std::str::from_utf8(seq)
            );
        }
    }

    fn parse_cigar(text: &str) -> Cigar {
        use crate::cigar::{CigarElement, Op};
        let mut elements = Vec::new();
        let mut num = String::new();
        for c in text.chars() {
            if c.is_ascii_digit() {
                num.push(c);
            } else {
                let op = match c {
                    'M' => Op::M,
                    'I' => Op::I,
                    'D' => Op::D,
                    'N' => Op::N,
                    'S' => Op::S,
                    'H' => Op::H,
                    'P' => Op::P,
                    '=' => Op::Eq,
                    'X' => Op::X,
                    _ => panic!("bad cigar op {c}"),
                };
                elements.push(CigarElement {
                    length: num.parse().unwrap(),
                    op,
                });
                num.clear();
            }
        }
        Cigar::new(elements)
    }
}
