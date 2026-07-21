//! The BAM file reader.
//!
//! Ported from `htsjdk.samtools.BAMFileReader.readHeader` and its record iteration, the inverse
//! of [`crate::writer`].
//!
//! One thing here is a decision rather than a transcription, and it is worth naming because a
//! reader is where fidelity is easy to over-claim. A BAM carries its sequence dictionary
//! **twice**: once inside the SAM text header and once in the binary block after it. htsjdk
//! reads the binary one and uses the text one for everything else, and it does **not** check
//! that they agree. So a file whose two dictionaries disagree is accepted, and which one a
//! consumer sees depends on which accessor it calls. This port reproduces that rather than
//! validating, and surfaces both so a caller can compare them if it wants to.

use crate::header::{Attributes, ProgramRecord, ReadGroup, SamHeader, SequenceRecord};
use crate::record::{BamRecord, DecodeError};
use crate::writer::BAM_MAGIC;

/// Why a BAM could not be read.
#[derive(Debug, Clone, PartialEq)]
pub enum ReadError {
    /// The first four bytes are not `BAM\1`.
    NotABam([u8; 4]),
    Truncated {
        need: usize,
        have: usize,
    },
    /// A length field that cannot be honoured.
    BadLength(i32),
    Record(DecodeError),
}

/// A BAM's header, as it is actually stored: the text, and the binary dictionary separately.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BamHeader {
    /// The `@`-prefixed SAM text, parsed.
    pub text: SamHeader,
    /// The binary sequence dictionary, which htsjdk reads and never reconciles with the text.
    pub binary_sequences: Vec<(String, i32)>,
    /// The raw header text, kept so a writer can reproduce it byte for byte without
    /// re-encoding through the parser.
    pub raw_text: String,
}

/// Reads a decompressed BAM stream.
#[derive(Debug)]
pub struct BamReader<'a> {
    data: &'a [u8],
    pos: usize,
    pub header: BamHeader,
}

impl<'a> BamReader<'a> {
    /// `BAMFileReader.readHeader`.
    pub fn new(data: &'a [u8]) -> Result<Self, ReadError> {
        let need = |p: usize, n: usize, have: usize| -> Result<(), ReadError> {
            if p + n > have {
                Err(ReadError::Truncated { need: p + n, have })
            } else {
                Ok(())
            }
        };
        need(0, 4, data.len())?;
        if data[0..4] != BAM_MAGIC {
            return Err(ReadError::NotABam([data[0], data[1], data[2], data[3]]));
        }

        need(4, 4, data.len())?;
        let text_len = i32::from_le_bytes(data[4..8].try_into().unwrap());
        if text_len < 0 {
            return Err(ReadError::BadLength(text_len));
        }
        let text_len = text_len as usize;
        need(8, text_len, data.len())?;
        // One byte per UTF-16 unit on the way out, so one char per byte on the way back.
        let raw_text: String = data[8..8 + text_len].iter().map(|&b| b as char).collect();
        let mut p = 8 + text_len;

        need(p, 4, data.len())?;
        let n_ref = i32::from_le_bytes(data[p..p + 4].try_into().unwrap());
        if n_ref < 0 {
            return Err(ReadError::BadLength(n_ref));
        }
        p += 4;

        let mut binary_sequences = Vec::with_capacity(n_ref as usize);
        for _ in 0..n_ref {
            need(p, 4, data.len())?;
            let name_len = i32::from_le_bytes(data[p..p + 4].try_into().unwrap());
            if name_len < 1 {
                return Err(ReadError::BadLength(name_len));
            }
            let name_len = name_len as usize;
            p += 4;
            need(p, name_len + 4, data.len())?;
            // The stored length includes the terminator, which is not part of the name.
            let name: String = data[p..p + name_len - 1]
                .iter()
                .map(|&b| b as char)
                .collect();
            p += name_len;
            let len = i32::from_le_bytes(data[p..p + 4].try_into().unwrap());
            p += 4;
            binary_sequences.push((name, len));
        }

        Ok(BamReader {
            data,
            pos: p,
            header: BamHeader {
                text: parse_header_text(&raw_text),
                binary_sequences,
                raw_text,
            },
        })
    }

    /// Byte offset where the records begin.
    pub fn records_start(&self) -> usize {
        self.pos
    }
}

impl Iterator for BamReader<'_> {
    type Item = Result<BamRecord, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.data.len() {
            return None;
        }
        match BamRecord::decode(&self.data[self.pos..]) {
            Ok(None) => None,
            Ok(Some((record, used))) => {
                self.pos += used;
                Some(Ok(record))
            }
            Err(e) => {
                // Stop rather than spin on the same bad bytes.
                self.pos = self.data.len();
                Some(Err(ReadError::Record(e)))
            }
        }
    }
}

/// `SAMTextHeaderCodec.parseHeaderRecord`, reduced to what a well-formed header needs.
///
/// Unknown `@` line types are dropped rather than rejected, which is what htsjdk does at
/// default stringency. Attribute order is preserved, because it is the write order; see
/// decision 0009.
pub fn parse_header_text(text: &str) -> SamHeader {
    let mut header = SamHeader {
        attributes: Attributes::new(),
        sequences: Vec::new(),
        read_groups: Vec::new(),
        programs: Vec::new(),
        comments: Vec::new(),
    };

    for line in text.lines() {
        if !line.starts_with('@') {
            continue;
        }
        let mut fields = line.split('\t');
        let kind = fields.next().unwrap_or("");
        // A `@CO` line is free text after the tab, not tag:value pairs.
        if kind == "@CO" {
            header.comments.push(line.to_string());
            continue;
        }
        let pairs: Vec<(&str, &str)> = fields.filter_map(|f| f.split_once(':')).collect();
        let take = |key: &str| pairs.iter().find(|(k, _)| *k == key).map(|(_, v)| *v);

        match kind {
            "@HD" => {
                for (k, v) in &pairs {
                    header.attributes.set(k, v);
                }
            }
            "@SQ" => {
                let name = take("SN").unwrap_or("").to_string();
                let length = take("LN").and_then(|v| v.parse().ok()).unwrap_or(0);
                let mut seq = SequenceRecord::new(&name, length);
                for (k, v) in &pairs {
                    if *k != "SN" && *k != "LN" {
                        seq.attributes.set(k, v);
                    }
                }
                header.sequences.push(seq);
            }
            "@RG" => {
                let mut rg = ReadGroup::new(take("ID").unwrap_or(""));
                for (k, v) in &pairs {
                    if *k != "ID" {
                        rg.attributes.set(k, v);
                    }
                }
                header.read_groups.push(rg);
            }
            "@PG" => {
                let mut pg = ProgramRecord::new(take("ID").unwrap_or(""));
                for (k, v) in &pairs {
                    if *k != "ID" {
                        pg.attributes.set(k, v);
                    }
                }
                header.programs.push(pg);
            }
            _ => {}
        }
    }
    header
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::{Cigar, CigarElement, Op};
    use crate::writer::BamWriter;

    fn sample() -> (SamHeader, Vec<BamRecord>) {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        let mut s = SequenceRecord::new("chr1", 250_000_000);
        s.attributes.set("AS", "GRCh38");
        h.sequences.push(s);
        h.sequences.push(SequenceRecord::new("chr2", 200_000_000));
        let mut rg = ReadGroup::new("rg1");
        rg.attributes.set("SM", "sample1");
        rg.attributes.set("PL", "ILLUMINA");
        h.read_groups.push(rg);
        h.add_comment("a comment");

        let records: Vec<BamRecord> = (0..200)
            .map(|i| BamRecord {
                read_name: format!("read{i}"),
                flags: 0,
                reference_index: 0,
                alignment_start: 100 + i * 13,
                mapping_quality: 60,
                cigar: Cigar::new(vec![CigarElement {
                    length: 10,
                    op: Op::M,
                }]),
                mate_reference_index: -1,
                mate_alignment_start: 0,
                inferred_insert_size: 0,
                read_bases: b"ACGTACGTAC".to_vec(),
                base_qualities: vec![30; 10],
                tags: Default::default(),
            })
            .collect();
        (h, records)
    }

    fn write(h: &SamHeader, records: &[BamRecord]) -> Vec<u8> {
        let mut w = BamWriter::new(Vec::new(), h).unwrap();
        for r in records {
            w.write(r).unwrap();
        }
        htsjdk_bgzf::decompress_all(&w.finish().unwrap()).unwrap()
    }

    #[test]
    fn a_written_file_reads_back_identically() {
        let (h, records) = sample();
        let bytes = write(&h, &records);
        let reader = BamReader::new(&bytes).unwrap();
        assert_eq!(reader.header.text, h);
        let back: Vec<BamRecord> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(back, records);
    }

    #[test]
    fn the_raw_header_text_is_preserved_verbatim() {
        let (h, records) = sample();
        let bytes = write(&h, &records);
        let reader = BamReader::new(&bytes).unwrap();
        assert_eq!(
            reader.header.raw_text,
            h.encode(),
            "the raw text must survive so a writer can reproduce it without re-encoding"
        );
    }

    /// The dictionary is stored twice and htsjdk never reconciles the two. Both are surfaced.
    #[test]
    fn both_copies_of_the_dictionary_are_available() {
        let (h, records) = sample();
        let bytes = write(&h, &records);
        let reader = BamReader::new(&bytes).unwrap();
        assert_eq!(
            reader.header.binary_sequences,
            vec![
                ("chr1".to_string(), 250_000_000),
                ("chr2".to_string(), 200_000_000)
            ]
        );
        let from_text: Vec<(String, i32)> = reader
            .header
            .text
            .sequences
            .iter()
            .map(|s| (s.name.clone(), s.length))
            .collect();
        assert_eq!(reader.header.binary_sequences, from_text);
    }

    /// A file whose two dictionaries disagree is accepted, exactly as htsjdk accepts it. This
    /// is reproduced rather than validated, and the test exists so the choice is visible.
    #[test]
    fn disagreeing_dictionaries_are_accepted_not_rejected() {
        let (h, records) = sample();
        let mut bytes = write(&h, &records);
        // Corrupt the binary copy's first sequence length, leaving the text alone.
        let text_len = i32::from_le_bytes(bytes[4..8].try_into().unwrap()) as usize;
        let name_len_at = 8 + text_len + 4;
        let name_len =
            i32::from_le_bytes(bytes[name_len_at..name_len_at + 4].try_into().unwrap()) as usize;
        let len_at = name_len_at + 4 + name_len;
        bytes[len_at..len_at + 4].copy_from_slice(&12345i32.to_le_bytes());

        let reader = BamReader::new(&bytes).unwrap();
        assert_eq!(reader.header.binary_sequences[0].1, 12345);
        assert_eq!(
            reader.header.text.sequences[0].length, 250_000_000,
            "the two copies disagree, and htsjdk does not check"
        );
    }

    #[test]
    fn a_file_that_is_not_a_bam_is_refused() {
        assert_eq!(
            BamReader::new(b"SAM\x01rest").unwrap_err(),
            ReadError::NotABam([b'S', b'A', b'M', 1])
        );
    }

    #[test]
    fn a_truncated_header_is_refused() {
        let (h, records) = sample();
        let bytes = write(&h, &records);
        for cut in [2usize, 6, 10, 20] {
            assert!(
                matches!(
                    BamReader::new(&bytes[..cut]),
                    Err(ReadError::Truncated { .. })
                ),
                "a header cut at {cut} must be refused"
            );
        }
    }

    #[test]
    fn an_empty_file_yields_no_records() {
        let (h, _) = sample();
        let bytes = write(&h, &[]);
        let reader = BamReader::new(&bytes).unwrap();
        assert_eq!(reader.count(), 0);
    }

    #[test]
    fn header_text_parsing_keeps_attribute_order() {
        let h = parse_header_text("@HD\tVN:1.6\tSO:coordinate\n@RG\tID:x\tSM:s\tLB:l\tPL:p\n");
        let keys: Vec<&str> = h.read_groups[0].attributes.iter().map(|(k, _)| k).collect();
        assert_eq!(keys, vec!["SM", "LB", "PL"]);
        assert_eq!(h.attributes.get("SO"), Some("coordinate"));
    }

    /// A `@CO` line is free text, not tag:value pairs, and a colon in it must not be parsed.
    #[test]
    fn comments_are_kept_as_free_text() {
        let h = parse_header_text("@HD\tVN:1.6\n@CO\tnote: this has a colon\n");
        assert_eq!(h.comments, vec!["@CO\tnote: this has a colon"]);
    }

    #[test]
    fn an_unknown_header_line_is_dropped_not_fatal() {
        let h = parse_header_text("@HD\tVN:1.6\n@ZZ\tsomething\n@SQ\tSN:c\tLN:5\n");
        assert_eq!(h.sequences.len(), 1);
    }
}
