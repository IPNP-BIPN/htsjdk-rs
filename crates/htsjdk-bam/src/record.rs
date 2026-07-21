//! The BAM record codec.
//!
//! Ported from `htsjdk.samtools.BAMRecordCodec.encode` / `decode`, with the supporting
//! constants from `BAMFileConstants` and `BAMRecord`.
//!
//! The BAM record layout is specified by the format. What is not specified, and therefore what
//! this port exists for, is every choice htsjdk makes *within* the format: which bin it
//! computes, which integer width each tag takes, what order the tags come in, what a missing
//! quality string becomes, and how a CIGAR too long for the format is displaced into a tag.

use crate::bases::{self, BadBase};
use crate::bin::{self, BinError};
use crate::cigar::{Cigar, CigarElement, Op};
use crate::tag::{Tag, TagError, TagValue, Tags};

/// `BAMFileConstants.FIXED_BLOCK_SIZE`, the 8 four-byte fields before the read name.
pub const FIXED_BLOCK_SIZE: usize = 8 * 4;

/// `BAMRecord.CIGAR_SIZE_MULTIPLIER`.
pub const CIGAR_SIZE_MULTIPLIER: usize = 4;

/// `BAMRecord.MAX_CIGAR_OPERATORS`. Past this the CIGAR moves into the `CG` tag.
pub const MAX_CIGAR_OPERATORS: usize = 0xFFFF;

/// `BAMRecord.MAX_CIGAR_ELEMENT_LENGTH`, `(1 << 28) - 1`: only 28 bits hold an operator length.
pub const MAX_CIGAR_ELEMENT_LENGTH: u32 = (1 << 28) - 1;

/// `SAMFlag.READ_UNMAPPED`.
pub const READ_UNMAPPED_FLAG: u16 = 0x4;

/// The `CG` tag, which carries a CIGAR too long to encode inline.
pub fn cg_tag() -> Tag {
    Tag::new(b"CG")
}

/// A record, in the fields `BAMRecordCodec` actually reads and writes.
#[derive(Debug, Clone, PartialEq)]
pub struct BamRecord {
    pub read_name: String,
    pub flags: u16,
    /// 0-based index into the sequence dictionary, or -1 for none.
    pub reference_index: i32,
    /// 1-based, inclusive. 0 is `NO_ALIGNMENT_START`.
    pub alignment_start: i32,
    pub mapping_quality: u8,
    pub cigar: Cigar,
    pub mate_reference_index: i32,
    /// 1-based, inclusive. 0 is `NO_ALIGNMENT_START`.
    pub mate_alignment_start: i32,
    pub inferred_insert_size: i32,
    pub read_bases: Vec<u8>,
    /// Phred scores as raw bytes. Empty means "absent", which encodes as all `0xFF`.
    pub base_qualities: Vec<u8>,
    pub tags: Tags,
}

impl Default for BamRecord {
    fn default() -> Self {
        BamRecord {
            read_name: String::new(),
            flags: 0,
            reference_index: bin::NO_ALIGNMENT_REFERENCE_INDEX,
            alignment_start: bin::NO_ALIGNMENT_START,
            mapping_quality: 0,
            cigar: Cigar::default(),
            mate_reference_index: bin::NO_ALIGNMENT_REFERENCE_INDEX,
            mate_alignment_start: bin::NO_ALIGNMENT_START,
            inferred_insert_size: 0,
            read_bases: Vec::new(),
            base_qualities: Vec::new(),
            tags: Tags::new(),
        }
    }
}

/// Why a record could not be encoded.
#[derive(Debug, Clone, PartialEq)]
pub enum EncodeError {
    /// `BAMRecordCodec.encode` refuses this outright.
    QualityLengthMismatch {
        read_length: usize,
        quals: usize,
    },
    /// `BinaryCodec.writeUByte`: the name length plus its terminator must fit a byte.
    ReadNameTooLong(usize),
    /// From `makeSentinelCigar`: 28 bits cannot hold the length.
    CigarTooLongToDisplace {
        read: u32,
        reference: u32,
    },
    BadBase(BadBase),
    Tag(TagError),
    Bin(BinError),
}

impl From<TagError> for EncodeError {
    fn from(e: TagError) -> Self {
        EncodeError::Tag(e)
    }
}

/// Why a record could not be decoded.
#[derive(Debug, Clone, PartialEq)]
pub enum DecodeError {
    /// `BAMRecordCodec.decode`: "Invalid record length".
    InvalidRecordLength(i32),
    /// The variable-length part does not hold what the fixed part promised.
    Truncated {
        need: usize,
        have: usize,
    },
    UnknownCigarOperator(u32),
    /// A tag block that does not parse.
    MalformedTags(String),
}

/// `BAMRecordCodec.makeSentinelCigar`.
///
/// The sentinel is `<readLength>S<referenceLength>N`, chosen so that it has the same read and
/// reference length as the CIGAR it replaces. That is what keeps the indexing bin correct: the
/// bin is computed from the sentinel, and it lands where the real CIGAR would have put it.
pub fn make_sentinel_cigar(cigar: &Cigar) -> Result<Cigar, EncodeError> {
    let read = cigar.read_length();
    let reference = cigar.reference_length();
    if read > MAX_CIGAR_ELEMENT_LENGTH || reference > MAX_CIGAR_ELEMENT_LENGTH {
        return Err(EncodeError::CigarTooLongToDisplace { read, reference });
    }
    Ok(Cigar::new(vec![
        CigarElement {
            length: read,
            op: Op::S,
        },
        CigarElement {
            length: reference,
            op: Op::N,
        },
    ]))
}

/// `BAMRecord.isSentinelCigar`.
pub fn is_sentinel_cigar(cigar: &Cigar, read_length: u32) -> bool {
    cigar.elements.len() == 2
        && cigar.elements[1].op == Op::N
        && cigar.elements[0].op == Op::S
        && (cigar.elements[0].length == read_length || read_length == 0)
}

impl BamRecord {
    pub fn read_unmapped(&self) -> bool {
        self.flags & READ_UNMAPPED_FLAG != 0
    }

    /// `SAMRecord.getAlignmentEnd()`: 1-based inclusive, or 0 when the read is unmapped.
    pub fn alignment_end(&self) -> i32 {
        if self.read_unmapped() {
            bin::NO_ALIGNMENT_START
        } else {
            self.alignment_start + self.cigar.reference_length() as i32 - 1
        }
    }

    /// `SAMRecord.getReadLength()`.
    pub fn read_length(&self) -> usize {
        self.read_bases.len()
    }

    /// `BAMRecordCodec.encode`.
    ///
    /// Returns the record's on-disk bytes including its own leading `block_size` field, which
    /// is how records appear in a BAM stream.
    pub fn encode(&self) -> Result<Vec<u8>, EncodeError> {
        let read_length = self.read_length();

        // If cigar is too long, put into CG tag and replace with sentinel value.
        let cigar_switcharoo = self.cigar.num_elements() > MAX_CIGAR_OPERATORS;
        let (cigar_to_write, tags) = if cigar_switcharoo {
            let mut tags = self.tags.clone();
            let encoded: Vec<i32> = self.cigar.encode().into_iter().map(|v| v as i32).collect();
            // htsjdk calls setAttribute, which inserts into the sorted tag list, so CG takes
            // its sorted position rather than being appended.
            tags.insert(
                cg_tag(),
                TagValue::IntArray {
                    values: encoded,
                    unsigned: false,
                },
            );
            (make_sentinel_cigar(&self.cigar)?, tags)
        } else {
            (self.cigar.clone(), self.tags.clone())
        };

        // `getReadNameLength()` is the character count; the +1 is the null terminator.
        let read_name_len = self.read_name.encode_utf16().count();
        let mut block_size = FIXED_BLOCK_SIZE
            + read_name_len
            + 1
            + cigar_to_write.num_elements() * CIGAR_SIZE_MULTIPLIER
            + read_length.div_ceil(2) // 2 bases per byte, round up
            + read_length;
        block_size += tags.binary_size()?;

        // The sentinel has the same reference length as the real CIGAR, so displacing a long
        // CIGAR does not move the bin.
        let index_bin = if self.alignment_start != bin::NO_ALIGNMENT_START {
            // htsjdk additionally forces the bin to 0 for references longer than
            // BIN_GENOMIC_SPAN, after warning once. That needs the sequence dictionary, so it
            // belongs to the writer rather than to the record; see `encode_with_bin`.
            bin::compute_indexing_bin(self.alignment_start, self.alignment_end())
                .map_err(EncodeError::Bin)?
        } else {
            0
        };

        if read_name_len + 1 > 255 {
            return Err(EncodeError::ReadNameTooLong(read_name_len));
        }

        let mut out = Vec::with_capacity(4 + block_size);
        out.extend_from_slice(&(block_size as i32).to_le_bytes());
        out.extend_from_slice(&self.reference_index.to_le_bytes());
        // 0-based!!
        out.extend_from_slice(&(self.alignment_start - 1).to_le_bytes());
        out.push((read_name_len + 1) as u8);
        out.push(self.mapping_quality);
        out.extend_from_slice(&(index_bin as u16).to_le_bytes());
        out.extend_from_slice(&(cigar_to_write.num_elements() as u16).to_le_bytes());
        out.extend_from_slice(&self.flags.to_le_bytes());
        out.extend_from_slice(&(read_length as i32).to_le_bytes());
        out.extend_from_slice(&self.mate_reference_index.to_le_bytes());
        out.extend_from_slice(&(self.mate_alignment_start - 1).to_le_bytes());
        out.extend_from_slice(&self.inferred_insert_size.to_le_bytes());

        if read_length != self.base_qualities.len() && !self.base_qualities.is_empty() {
            return Err(EncodeError::QualityLengthMismatch {
                read_length,
                quals: self.base_qualities.len(),
            });
        }

        out.extend(self.read_name.encode_utf16().map(|u| (u & 0xFF) as u8));
        out.push(0);
        for element in cigar_to_write.encode() {
            out.extend_from_slice(&element.to_le_bytes());
        }
        out.extend_from_slice(
            &bases::bytes_to_compressed_bases(&self.read_bases).map_err(EncodeError::BadBase)?,
        );
        // An absent quality string becomes 0xFF repeated, which is how SAM's `*` is stored.
        if self.base_qualities.is_empty() {
            out.extend(std::iter::repeat_n(0xFFu8, read_length));
        } else {
            out.extend_from_slice(&self.base_qualities);
        }
        tags.write(&mut out)?;

        debug_assert_eq!(
            out.len(),
            4 + block_size,
            "block_size must describe the bytes actually written"
        );
        Ok(out)
    }

    /// `BAMRecordCodec.decode`.
    ///
    /// Returns the record and how many bytes it consumed.
    pub fn decode(input: &[u8]) -> Result<Option<(BamRecord, usize)>, DecodeError> {
        if input.len() < 4 {
            // `readInt` raising RuntimeEOFException, which `decode` turns into null.
            return Ok(None);
        }
        let record_length = i32::from_le_bytes(input[0..4].try_into().unwrap());
        if (record_length as usize) < FIXED_BLOCK_SIZE || record_length < 0 {
            return Err(DecodeError::InvalidRecordLength(record_length));
        }
        let total = 4 + record_length as usize;
        if input.len() < total {
            return Err(DecodeError::Truncated {
                need: total,
                have: input.len(),
            });
        }

        let g4 = |o: usize| i32::from_le_bytes(input[o..o + 4].try_into().unwrap());
        let reference_index = g4(4);
        let alignment_start = g4(8) + 1;
        let read_name_length = input[12] as usize;
        let mapping_quality = input[13];
        let _bin = u16::from_le_bytes(input[14..16].try_into().unwrap());
        let cigar_len = u16::from_le_bytes(input[16..18].try_into().unwrap()) as usize;
        let flags = u16::from_le_bytes(input[18..20].try_into().unwrap());
        let read_len = g4(20) as usize;
        let mate_reference_index = g4(24);
        let mate_alignment_start = g4(28) + 1;
        let inferred_insert_size = g4(32);

        let mut p = 4 + FIXED_BLOCK_SIZE;
        let need = |p: usize, n: usize| -> Result<(), DecodeError> {
            if p + n > total {
                Err(DecodeError::Truncated {
                    need: p + n,
                    have: total,
                })
            } else {
                Ok(())
            }
        };

        need(p, read_name_length)?;
        // The stored length includes the null terminator, which is not part of the name.
        let name_bytes = &input[p..p + read_name_length.saturating_sub(1)];
        let read_name: String = name_bytes.iter().map(|&b| b as char).collect();
        p += read_name_length;

        need(p, cigar_len * CIGAR_SIZE_MULTIPLIER)?;
        let mut binary_cigar = Vec::with_capacity(cigar_len);
        for i in 0..cigar_len {
            let o = p + i * 4;
            binary_cigar.push(u32::from_le_bytes(input[o..o + 4].try_into().unwrap()));
        }
        p += cigar_len * CIGAR_SIZE_MULTIPLIER;
        let cigar = Cigar::decode(&binary_cigar).ok_or_else(|| {
            DecodeError::UnknownCigarOperator(
                binary_cigar.iter().find(|v| *v & 0x0F > 8).unwrap() & 0x0F,
            )
        })?;

        let packed = read_len.div_ceil(2);
        need(p, packed)?;
        let read_bases = bases::compressed_bases_to_bytes(read_len, input, p);
        p += packed;

        need(p, read_len)?;
        let raw_quals = &input[p..p + read_len];
        // All-0xFF is the on-disk spelling of "no qualities"; htsjdk surfaces it as such.
        let base_qualities = if read_len > 0 && raw_quals.iter().all(|&b| b == 0xFF) {
            Vec::new()
        } else {
            raw_quals.to_vec()
        };
        p += read_len;

        let tags = read_tags(&input[p..total])?;

        Ok(Some((
            BamRecord {
                read_name,
                flags,
                reference_index,
                alignment_start,
                mapping_quality,
                cigar,
                mate_reference_index,
                mate_alignment_start,
                inferred_insert_size,
                read_bases,
                base_qualities,
                tags,
            },
            total,
        )))
    }
}

/// `BinaryTagCodec.readTags`.
pub fn read_tags(mut buf: &[u8]) -> Result<Tags, DecodeError> {
    let mut tags = Tags::new();
    let short = |b: &[u8]| i16::from_le_bytes([b[0], b[1]]);
    let int = |b: &[u8]| i32::from_le_bytes([b[0], b[1], b[2], b[3]]);
    let malformed = |m: &str| DecodeError::MalformedTags(m.to_string());

    while !buf.is_empty() {
        if buf.len() < 3 {
            return Err(malformed("truncated tag header"));
        }
        let tag = Tag(short(buf));
        let ty = buf[2];
        buf = &buf[3..];

        let take = |buf: &mut &[u8], n: usize| -> Result<Vec<u8>, DecodeError> {
            if buf.len() < n {
                return Err(DecodeError::MalformedTags(format!(
                    "tag value wants {n} bytes, {} left",
                    buf.len()
                )));
            }
            let (a, b) = buf.split_at(n);
            *buf = b;
            Ok(a.to_vec())
        };
        let nul_string = |buf: &mut &[u8]| -> Result<String, DecodeError> {
            let end = buf
                .iter()
                .position(|&b| b == 0)
                .ok_or_else(|| DecodeError::MalformedTags("unterminated string".into()))?;
            let s: String = buf[..end].iter().map(|&b| b as char).collect();
            *buf = &buf[end + 1..];
            Ok(s)
        };

        let value = match ty {
            b'Z' => TagValue::Str(nul_string(&mut buf)?),
            b'A' => TagValue::Char(take(&mut buf, 1)?[0]),
            // htsjdk widens 'I' into a long when it does not fit a signed int, so the
            // in-memory value is always the mathematical one.
            b'I' => TagValue::Int(int(&take(&mut buf, 4)?) as u32 as i64),
            b'i' => TagValue::Int(int(&take(&mut buf, 4)?) as i64),
            b's' => TagValue::Int(short(&take(&mut buf, 2)?) as i64),
            b'S' => TagValue::Int(short(&take(&mut buf, 2)?) as u16 as i64),
            b'c' => TagValue::Int(take(&mut buf, 1)?[0] as i8 as i64),
            b'C' => TagValue::Int(take(&mut buf, 1)?[0] as i64),
            b'f' => TagValue::Float(f32::from_le_bytes(take(&mut buf, 4)?.try_into().unwrap())),
            b'H' => {
                let hex = nul_string(&mut buf)?;
                let bytes = (0..hex.len() / 2)
                    .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|_| malformed("bad hex tag"))?;
                TagValue::Hex(bytes)
            }
            b'B' => {
                if buf.len() < 5 {
                    return Err(malformed("truncated array header"));
                }
                let element_type = buf[0];
                let n = int(&buf[1..5]) as usize;
                buf = &buf[5..];
                let unsigned = element_type.is_ascii_uppercase();
                match element_type.to_ascii_lowercase() {
                    b'c' => TagValue::ByteArray {
                        values: take(&mut buf, n)?.into_iter().map(|b| b as i8).collect(),
                        unsigned,
                    },
                    b's' => TagValue::ShortArray {
                        values: take(&mut buf, n * 2)?
                            .chunks(2)
                            .map(|c| i16::from_le_bytes([c[0], c[1]]))
                            .collect(),
                        unsigned,
                    },
                    b'i' => TagValue::IntArray {
                        values: take(&mut buf, n * 4)?.chunks(4).map(int).collect(),
                        unsigned,
                    },
                    b'f' => TagValue::FloatArray(
                        take(&mut buf, n * 4)?
                            .chunks(4)
                            .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
                            .collect(),
                    ),
                    _ => return Err(malformed("unrecognized tag array type")),
                }
            }
            _ => return Err(malformed("unrecognized tag type")),
        };
        tags.insert(tag, value);
    }
    Ok(tags)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mapped_read() -> BamRecord {
        BamRecord {
            read_name: "read1".into(),
            flags: 0,
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            cigar: Cigar::new(vec![CigarElement {
                length: 4,
                op: Op::M,
            }]),
            mate_reference_index: -1,
            mate_alignment_start: 0,
            inferred_insert_size: 0,
            read_bases: b"ACGT".to_vec(),
            base_qualities: vec![30, 30, 30, 30],
            tags: Tags::new(),
        }
    }

    #[test]
    fn block_size_describes_the_bytes_actually_written() {
        let mut rec = mapped_read();
        rec.tags.insert(Tag::new(b"NM"), TagValue::Int(1));
        rec.tags.insert(Tag::new(b"MD"), TagValue::Str("4".into()));
        let bytes = rec.encode().unwrap();
        let block_size = i32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize;
        assert_eq!(bytes.len(), 4 + block_size);
    }

    /// The fixed header, field by field, so a transposed field cannot pass.
    #[test]
    fn the_fixed_header_is_laid_out_exactly() {
        let bytes = mapped_read().encode().unwrap();
        assert_eq!(i32::from_le_bytes(bytes[4..8].try_into().unwrap()), 0);
        assert_eq!(
            i32::from_le_bytes(bytes[8..12].try_into().unwrap()),
            99,
            "alignment start is written 0-based"
        );
        assert_eq!(bytes[12], 6, "read name length includes the terminator");
        assert_eq!(bytes[13], 60);
        assert_eq!(u16::from_le_bytes(bytes[16..18].try_into().unwrap()), 1);
        assert_eq!(u16::from_le_bytes(bytes[18..20].try_into().unwrap()), 0);
        assert_eq!(i32::from_le_bytes(bytes[20..24].try_into().unwrap()), 4);
        assert_eq!(
            i32::from_le_bytes(bytes[28..32].try_into().unwrap()),
            -1,
            "an absent mate start is written as -1, i.e. 0 minus one"
        );
    }

    #[test]
    fn the_read_name_is_null_terminated() {
        let bytes = mapped_read().encode().unwrap();
        let start = 4 + FIXED_BLOCK_SIZE;
        assert_eq!(&bytes[start..start + 5], b"read1");
        assert_eq!(bytes[start + 5], 0);
    }

    /// The bin is computed from the alignment, not left at zero. A zero bin produces a file
    /// every reader accepts and every byte comparison rejects.
    #[test]
    fn the_bin_is_computed_for_a_placed_read() {
        let bytes = mapped_read().encode().unwrap();
        let b = u16::from_le_bytes(bytes[14..16].try_into().unwrap());
        assert_eq!(b, bin::compute_indexing_bin(100, 103).unwrap() as u16);
        assert_ne!(b, 0);
    }

    #[test]
    fn an_unplaced_read_gets_bin_zero() {
        let mut rec = mapped_read();
        rec.alignment_start = bin::NO_ALIGNMENT_START;
        rec.flags = READ_UNMAPPED_FLAG;
        rec.reference_index = -1;
        let bytes = rec.encode().unwrap();
        assert_eq!(u16::from_le_bytes(bytes[14..16].try_into().unwrap()), 0);
    }

    /// An unmapped read's alignment end is 0 regardless of its CIGAR, which changes the bin.
    #[test]
    fn the_unmapped_flag_overrides_the_cigar_for_the_alignment_end() {
        let mut rec = mapped_read();
        assert_eq!(rec.alignment_end(), 103);
        rec.flags |= READ_UNMAPPED_FLAG;
        assert_eq!(rec.alignment_end(), 0);
    }

    #[test]
    fn absent_qualities_become_all_ff() {
        let mut rec = mapped_read();
        rec.base_qualities = Vec::new();
        let bytes = rec.encode().unwrap();
        let quals_at = 4 + FIXED_BLOCK_SIZE + 6 + 4 + 2;
        assert_eq!(&bytes[quals_at..quals_at + 4], &[0xFF; 4]);
    }

    #[test]
    fn mismatched_quality_length_is_refused() {
        let mut rec = mapped_read();
        rec.base_qualities = vec![30, 30];
        assert_eq!(
            rec.encode(),
            Err(EncodeError::QualityLengthMismatch {
                read_length: 4,
                quals: 2
            })
        );
    }

    #[test]
    fn an_over_long_read_name_is_refused_not_truncated() {
        let mut rec = mapped_read();
        rec.read_name = "x".repeat(255);
        assert_eq!(rec.encode(), Err(EncodeError::ReadNameTooLong(255)));
        rec.read_name = "x".repeat(254);
        assert!(rec.encode().is_ok(), "254 + terminator still fits a byte");
    }

    #[test]
    fn a_record_round_trips() {
        let mut rec = mapped_read();
        rec.tags.insert(Tag::new(b"NM"), TagValue::Int(1));
        rec.tags.insert(Tag::new(b"MD"), TagValue::Str("4".into()));
        rec.tags.insert(
            Tag::new(b"BQ"),
            TagValue::ByteArray {
                values: vec![1, 2, 3],
                unsigned: true,
            },
        );
        let bytes = rec.encode().unwrap();
        let (back, used) = BamRecord::decode(&bytes).unwrap().unwrap();
        assert_eq!(used, bytes.len());
        assert_eq!(back, rec);
        assert_eq!(back.encode().unwrap(), bytes, "re-encoding must be stable");
    }

    #[test]
    fn decoding_an_empty_stream_yields_no_record() {
        assert_eq!(BamRecord::decode(&[]).unwrap(), None);
    }

    #[test]
    fn a_short_record_length_is_refused() {
        let mut bytes = mapped_read().encode().unwrap();
        bytes[0..4].copy_from_slice(&8i32.to_le_bytes());
        assert_eq!(
            BamRecord::decode(&bytes),
            Err(DecodeError::InvalidRecordLength(8))
        );
    }

    // --- the long-CIGAR displacement -------------------------------------------------

    fn long_cigar() -> Cigar {
        // 65536 elements: one past MAX_CIGAR_OPERATORS.
        let mut elements = Vec::with_capacity(MAX_CIGAR_OPERATORS + 1);
        for i in 0..=MAX_CIGAR_OPERATORS {
            elements.push(CigarElement {
                length: 1,
                op: if i % 2 == 0 { Op::M } else { Op::I },
            });
        }
        Cigar::new(elements)
    }

    #[test]
    fn a_long_cigar_moves_into_the_cg_tag_behind_a_sentinel() {
        let c = long_cigar();
        let mut rec = mapped_read();
        rec.read_bases = vec![b'A'; c.read_length() as usize];
        rec.base_qualities = vec![30; c.read_length() as usize];
        rec.cigar = c.clone();

        let bytes = rec.encode().unwrap();
        let n_cigar = u16::from_le_bytes(bytes[16..18].try_into().unwrap());
        assert_eq!(n_cigar, 2, "only the sentinel is written inline");

        let (back, _) = BamRecord::decode(&bytes).unwrap().unwrap();
        assert!(is_sentinel_cigar(&back.cigar, back.read_length() as u32));
        match back.tags.get(cg_tag()) {
            Some(TagValue::IntArray { values, unsigned }) => {
                assert!(!unsigned);
                assert_eq!(values.len(), c.num_elements());
                let decoded =
                    Cigar::decode(&values.iter().map(|&v| v as u32).collect::<Vec<_>>()).unwrap();
                assert_eq!(decoded, c, "the CG tag carries the real CIGAR intact");
            }
            other => panic!("expected an int array in CG, got {other:?}"),
        }
    }

    /// The sentinel must preserve read and reference length, because the bin is computed from
    /// it. If it did not, displacing a long CIGAR would silently move the read in the index.
    #[test]
    fn the_sentinel_preserves_both_lengths_so_the_bin_does_not_move() {
        let c = long_cigar();
        let s = make_sentinel_cigar(&c).unwrap();
        assert_eq!(s.read_length(), c.read_length());
        assert_eq!(s.reference_length(), c.reference_length());

        let mut with_long = mapped_read();
        with_long.read_bases = vec![b'A'; c.read_length() as usize];
        with_long.base_qualities = vec![30; c.read_length() as usize];
        with_long.cigar = c;
        let mut with_sentinel = with_long.clone();
        with_sentinel.cigar = s;

        let bin_long = &with_long.encode().unwrap()[14..16];
        let bin_sentinel = &with_sentinel.encode().unwrap()[14..16];
        assert_eq!(bin_long, bin_sentinel);
    }

    /// The CG tag is inserted, not appended, so it lands in packed-short order among the
    /// record's other tags. Appending would give a valid file with a different tag block.
    #[test]
    fn the_cg_tag_takes_its_sorted_position() {
        let c = long_cigar();
        let mut rec = mapped_read();
        rec.read_bases = vec![b'A'; c.read_length() as usize];
        rec.base_qualities = vec![30; c.read_length() as usize];
        rec.cigar = c;
        // "AG" packs below "CG"; "CH" packs above it.
        rec.tags.insert(Tag::new(b"AG"), TagValue::Int(1));
        rec.tags.insert(Tag::new(b"CH"), TagValue::Int(1));

        let bytes = rec.encode().unwrap();
        let (back, _) = BamRecord::decode(&bytes).unwrap().unwrap();
        let order: Vec<String> = back.tags.iter().map(|(t, _)| t.to_string()).collect();
        assert_eq!(order, vec!["AG", "CG", "CH"]);
    }

    #[test]
    fn exactly_max_operators_stays_inline() {
        let mut elements = Vec::new();
        for _ in 0..MAX_CIGAR_OPERATORS {
            elements.push(CigarElement {
                length: 1,
                op: Op::M,
            });
        }
        let c = Cigar::new(elements);
        let mut rec = mapped_read();
        rec.read_bases = vec![b'A'; c.read_length() as usize];
        rec.base_qualities = vec![30; c.read_length() as usize];
        rec.cigar = c;
        let bytes = rec.encode().unwrap();
        assert_eq!(
            u16::from_le_bytes(bytes[16..18].try_into().unwrap()) as usize,
            MAX_CIGAR_OPERATORS,
            "the switch is at > MAX_CIGAR_OPERATORS, not >="
        );
        assert!(rec.tags.get(cg_tag()).is_none());
    }

    /// Encoding must not mutate the record. htsjdk sets the CG attribute and removes it again;
    /// a port that forgot the removal would change the caller's record.
    #[test]
    fn encoding_leaves_the_record_untouched() {
        let c = long_cigar();
        let mut rec = mapped_read();
        rec.read_bases = vec![b'A'; c.read_length() as usize];
        rec.base_qualities = vec![30; c.read_length() as usize];
        rec.cigar = c;
        let before = rec.clone();
        rec.encode().unwrap();
        assert_eq!(rec, before);
    }
}
