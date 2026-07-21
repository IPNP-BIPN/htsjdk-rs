//! Reading a reference FASTA.
//!
//! Ported from `htsjdk.samtools.reference.FastaSequenceFile` and `ReferenceSequence` at
//! htsjdk 4.2.0, restricted to sequential whole-contig reads, which is what
//! `ReferenceSequenceFileWalker` gives a `SinglePassSamProgram`.
//!
//! Two behaviours that a casual FASTA reader gets wrong, both of which reach the output of any
//! tool that compares read bases to reference bases:
//!
//!  * **The bases are not uppercased.** `readSequence` copies bytes and trims trailing
//!    whitespace, and nothing in the path folds case. Soft-masked reference regions, which are
//!    conventionally lowercase, stay lowercase. It does not change a mismatch count only because
//!    `SequenceUtil.basesEqual` is itself case-insensitive; a port that uppercased here and
//!    compared bytes there would agree by luck, and stop agreeing the moment either changed.
//!  * **The name is trimmed and then truncated at the first whitespace.** `>chr1 human` is the
//!    contig `chr1`, and the description is discarded rather than being part of the name.

use std::io::{BufRead, BufReader, Read};

/// `ReferenceSequence`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceSequence {
    pub name: String,
    /// 0-based position in the file, as `FastaSequenceFile` assigns it.
    pub index: i32,
    pub bases: Vec<u8>,
}

#[derive(Debug)]
pub enum FastaError {
    Io(std::io::Error),
    /// `"Format exception reading FASTA ... Expected > but saw ..."`.
    ExpectedHeader(u8),
    /// `"Missing sequence name in FASTA"`.
    MissingName,
}

impl From<std::io::Error> for FastaError {
    fn from(e: std::io::Error) -> Self {
        FastaError::Io(e)
    }
}

/// `SAMSequenceRecord.truncateSequenceName`: everything up to the first whitespace.
fn truncate_at_whitespace(name: &str) -> &str {
    match name.find(char::is_whitespace) {
        Some(i) => &name[..i],
        None => name,
    }
}

/// Reads every sequence in a FASTA, in file order.
pub fn read_fasta<R: Read>(source: R) -> Result<Vec<ReferenceSequence>, FastaError> {
    let mut out: Vec<ReferenceSequence> = Vec::new();
    let mut index = -1;
    for line in BufReader::new(source).lines() {
        let line = line?;
        // Blank lines between records are skipped by `skipNewlines`, not treated as data.
        if line.trim().is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            let name = truncate_at_whitespace(header.trim());
            if name.is_empty() {
                return Err(FastaError::MissingName);
            }
            index += 1;
            out.push(ReferenceSequence {
                name: name.to_string(),
                index,
                bases: Vec::new(),
            });
        } else {
            match out.last_mut() {
                Some(seq) => seq.bases.extend_from_slice(line.trim_end().as_bytes()),
                None => return Err(FastaError::ExpectedHeader(line.as_bytes()[0])),
            }
        }
    }
    Ok(out)
}

/// Reads a FASTA from a path.
pub fn read_fasta_file(
    path: impl AsRef<std::path::Path>,
) -> Result<Vec<ReferenceSequence>, FastaError> {
    read_fasta(std::fs::File::open(path)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> Vec<ReferenceSequence> {
        read_fasta(s.as_bytes()).unwrap()
    }

    #[test]
    fn lines_are_concatenated_and_indices_follow_file_order() {
        let seqs = parse(">chr1\nACGT\nACGT\n>chr2\nTTTT\n");
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].name, "chr1");
        assert_eq!(seqs[0].index, 0);
        assert_eq!(seqs[0].bases, b"ACGTACGT");
        assert_eq!(seqs[1].index, 1);
    }

    /// Soft-masked regions survive as lowercase, because nothing in the read path folds case.
    #[test]
    fn bases_keep_their_case() {
        assert_eq!(parse(">c\nACgtN\n")[0].bases, b"ACgtN");
    }

    #[test]
    fn the_name_stops_at_the_first_whitespace() {
        assert_eq!(
            parse(">chr1 Homo sapiens chromosome 1\nA\n")[0].name,
            "chr1"
        );
        assert_eq!(parse(">chr1\tdescription\nA\n")[0].name, "chr1");
    }

    #[test]
    fn blank_lines_are_skipped() {
        let seqs = parse(">c\nACGT\n\n\nTTTT\n");
        assert_eq!(seqs[0].bases, b"ACGTTTTT");
    }

    #[test]
    fn bases_before_any_header_are_a_format_error() {
        assert!(matches!(
            read_fasta(&b"ACGT\n"[..]),
            Err(FastaError::ExpectedHeader(b'A'))
        ));
    }

    #[test]
    fn an_empty_name_is_rejected() {
        assert!(matches!(
            read_fasta(&b">\nACGT\n"[..]),
            Err(FastaError::MissingName)
        ));
    }

    #[test]
    fn a_sequence_may_be_empty() {
        let seqs = parse(">a\n>b\nAC\n");
        assert!(seqs[0].bases.is_empty());
        assert_eq!(seqs[1].bases, b"AC");
    }
}
