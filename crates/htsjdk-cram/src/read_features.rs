//! How a record becomes read features: the first half of the CRAM record model.
//!
//! Ported from `htsjdk.samtools.cram.structure.CRAMRecordReadFeatures` and the twelve classes of
//! `htsjdk.samtools.cram.encoding.readfeatures` at htsjdk 4.2.0.
//!
//! The frames are pinned through [`crate::slice_header`]. A slice's records are not stored as
//! bases and a cigar: they are stored as an alignment start, a read length, and a list of read
//! features, each a one-letter operator and a payload. **Everything that matches the reference is
//! stored as nothing at all**, which is where CRAM's compression comes from.
//!
//! # The positions are one-based, and the interface says they are not
//!
//! Every construction site passes `zeroBasedPositionInRead + 1`, while `ReadFeature.getPosition`'s
//! javadoc says "zero-based position in the read". The doc is wrong about all twelve
//! implementations, and a port that believes it is off by one on every feature.
//!
//! # An insertion of n bases becomes n features, a soft clip of n becomes one
//!
//! `addInsertion` emits one `InsertBase` per base with a comment saying it should use a `Bases`
//! feature and does not, because that would need a `ByteArrayLenEncoding` and therefore a
//! frequency distribution over lengths. `addSoftClip`, five lines away, emits a single `SoftClip`
//! carrying all of them. So the `Insertion` and `Bases` features exist and the writer never emits
//! either.
//!
//! # A mismatch splits on the alphabet, not on the cigar
//!
//! An ACGTN-to-ACGTN mismatch is a [`ReadFeature::Substitution`]; anything else is a
//! [`ReadFeature::ReadBase`], which carries the quality score a second time. And `M`, `X` and `EQ`
//! all walk the same comparison, so **an `X` over bases that match emits nothing** and an `=` over
//! bases that differ emits a substitution. The cigar operator only says how far to walk.
//!
//! # Past the end of the reference every base is compared against `N`
//!
//! Which means a read base of `N` out there **matches** and is stored as nothing, while any other
//! base becomes a substitution whose reference base is `N`. Measured: `NNNN` placed so that its
//! last base is past the end produces three features, not four.
//!
//! # `SEQ="*"` manufactures `N`s, and then they mismatch
//!
//! A record with no sequence gets one `N` per read base the cigar consumes, and those `N`s are
//! compared like any other base: against an ordinary reference it produces **a substitution per
//! position**.
//!
//! # The missing-quality test is an identity test
//!
//! `baseQualities.equals(SAMRecord.NULL_QUALS)` is `Object.equals` on a `byte[]`, so it is true
//! only for that one array instance. A record whose qualities are an equal but distinct empty
//! array takes the other branch and indexes it: measured, `ArrayIndexOutOfBoundsException: Index 3
//! out of bounds for length 0`. [`Qualities`] carries that distinction rather than hiding it.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

/// `CRAMCompressionRecord.MISSING_QUALITY_SCORE`.
pub const MISSING_QUALITY_SCORE: i8 = -1;
/// `Substitution.NO_CODE`, before the substitution matrix assigns one.
pub const NO_CODE: i8 = -1;

/// `SequenceUtil.isUpperACGTN`.
///
/// Upper case only: the comparison is done on BAM bases, which are already upper case, and a
/// lower-case base would fall to the [`ReadFeature::ReadBase`] branch.
pub fn is_upper_acgtn(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T' | b'N')
}

/// One read feature: an operator and its payload.
///
/// The twelve are the reference's twelve. Two of them, [`ReadFeature::Insertion`] and
/// [`ReadFeature::Bases`], are never produced by [`create_read_features`]; they are here because
/// the reader must accept them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadFeature {
    /// `'Q'`.
    BaseQualityScore { position: i32, quality: i8 },
    /// `'b'`. Read but never written.
    Bases { position: i32, bases: Vec<u8> },
    /// `'D'`.
    Deletion { position: i32, length: i32 },
    /// `'H'`.
    HardClip { position: i32, length: i32 },
    /// `'i'`. One per inserted base.
    InsertBase { position: i32, base: u8 },
    /// `'I'`. Read but never written.
    Insertion { position: i32, sequence: Vec<u8> },
    /// `'P'`.
    Padding { position: i32, length: i32 },
    /// `'B'`, the non-ACGTN mismatch, which carries the quality score a second time.
    ReadBase {
        position: i32,
        base: u8,
        quality: i8,
    },
    /// `'N'`.
    RefSkip { position: i32, length: i32 },
    /// `'q'`.
    Scores { position: i32, scores: Vec<u8> },
    /// `'S'`, carrying every clipped base at once.
    SoftClip { position: i32, sequence: Vec<u8> },
    /// `'X'`. The code is [`NO_CODE`] until the substitution matrix assigns one.
    Substitution {
        position: i32,
        base: u8,
        reference_base: u8,
        code: i8,
    },
}

impl ReadFeature {
    /// The one-letter operator, which is how the feature is written and read.
    pub fn operator(&self) -> u8 {
        match self {
            ReadFeature::BaseQualityScore { .. } => b'Q',
            ReadFeature::Bases { .. } => b'b',
            ReadFeature::Deletion { .. } => b'D',
            ReadFeature::HardClip { .. } => b'H',
            ReadFeature::InsertBase { .. } => b'i',
            ReadFeature::Insertion { .. } => b'I',
            ReadFeature::Padding { .. } => b'P',
            ReadFeature::ReadBase { .. } => b'B',
            ReadFeature::RefSkip { .. } => b'N',
            ReadFeature::Scores { .. } => b'q',
            ReadFeature::SoftClip { .. } => b'S',
            ReadFeature::Substitution { .. } => b'X',
        }
    }

    /// The one-based position in the read, whatever the interface's javadoc says.
    pub fn position(&self) -> i32 {
        match self {
            ReadFeature::BaseQualityScore { position, .. }
            | ReadFeature::Bases { position, .. }
            | ReadFeature::Deletion { position, .. }
            | ReadFeature::HardClip { position, .. }
            | ReadFeature::InsertBase { position, .. }
            | ReadFeature::Insertion { position, .. }
            | ReadFeature::Padding { position, .. }
            | ReadFeature::ReadBase { position, .. }
            | ReadFeature::RefSkip { position, .. }
            | ReadFeature::Scores { position, .. }
            | ReadFeature::SoftClip { position, .. }
            | ReadFeature::Substitution { position, .. } => *position,
        }
    }
}

/// A record's base qualities, which carry a distinction Java makes by array identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Qualities<'a> {
    /// `SAMRecord.NULL_QUALS`, the singleton a `QUAL` of `*` leaves behind. A non-ACGTN mismatch
    /// takes [`MISSING_QUALITY_SCORE`] from it.
    Missing,
    /// Any other array, **including an empty one**. htsjdk's test is `equals` on the object, so an
    /// empty array that is not the singleton is indexed like any other and throws.
    Present(&'a [u8]),
}

/// What building read features is refused with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadFeatureError {
    /// `ArrayIndexOutOfBoundsException`, from the reference cursor, the read cursor or the
    /// qualities. The JDK's message names the index and the length and nothing else, so the three
    /// sites are indistinguishable to a caller.
    IndexOutOfBounds { index: i64, length: usize },
    /// `IllegalArgumentException: Unsupported cigar operator: <op>`.
    ///
    /// Unreachable: the switch handles all nine operators htsjdk has.
    UnsupportedCigarOperator(char),
}

impl ReadFeatureError {
    pub fn message(&self) -> String {
        match self {
            ReadFeatureError::IndexOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
            ReadFeatureError::UnsupportedCigarOperator(op) => {
                format!("Unsupported cigar operator: {op}")
            }
        }
    }
}

/// `new CRAMRecordReadFeatures(samRecord, bamReadBases, refBases)`.
///
/// `read_bases` are the BAM bases, already upper case. An empty slice is `SEQ="*"` and is replaced
/// by `N`s, one per read base the cigar consumes.
pub fn create_read_features(
    cigar: &Cigar,
    alignment_start: i32,
    read_bases: &[u8],
    base_qualities: Qualities<'_>,
    reference_bases: &[u8],
) -> Result<Vec<ReadFeature>, ReadFeatureError> {
    let manufactured;
    let bases = if read_bases.is_empty() {
        manufactured = vec![b'N'; cigar.read_length() as usize];
        &manufactured[..]
    } else {
        read_bases
    };

    let mut features = Vec::new();
    let mut zero_based_position_in_read: i32 = 0;
    let mut alignment_start_offset: i32 = 0;

    for element in &cigar.elements {
        let length = element.length as i32;
        match element.op {
            Op::D => features.push(ReadFeature::Deletion {
                position: zero_based_position_in_read + 1,
                length,
            }),
            Op::N => features.push(ReadFeature::RefSkip {
                position: zero_based_position_in_read + 1,
                length,
            }),
            Op::P => features.push(ReadFeature::Padding {
                position: zero_based_position_in_read + 1,
                length,
            }),
            Op::H => features.push(ReadFeature::HardClip {
                position: zero_based_position_in_read + 1,
                length,
            }),
            // One feature for all of the clipped bases.
            Op::S => features.push(ReadFeature::SoftClip {
                position: zero_based_position_in_read + 1,
                sequence: slice_of(bases, zero_based_position_in_read, length)?.to_vec(),
            }),
            // One feature per inserted base, which is the opposite decision five lines away.
            Op::I => {
                let inserted = slice_of(bases, zero_based_position_in_read, length)?;
                for (i, base) in inserted.iter().enumerate() {
                    features.push(ReadFeature::InsertBase {
                        position: zero_based_position_in_read + 1 + i as i32,
                        base: *base,
                    });
                }
            }
            // The cigar's own claim is not consulted: all three walk the same comparison.
            Op::M | Op::X | Op::Eq => add_mismatch_read_features(
                reference_bases,
                alignment_start,
                &mut features,
                zero_based_position_in_read,
                alignment_start_offset,
                length,
                bases,
                base_qualities,
            )?,
        }

        if element.op.consumes_read_bases() {
            zero_based_position_in_read += length;
        }
        if element.op.consumes_reference_bases() {
            alignment_start_offset += length;
        }
    }
    Ok(features)
}

fn slice_of(bases: &[u8], from: i32, length: i32) -> Result<&[u8], ReadFeatureError> {
    let start = from.max(0) as usize;
    let end = start + length.max(0) as usize;
    bases
        .get(start..end)
        .ok_or(ReadFeatureError::IndexOutOfBounds {
            index: end as i64 - 1,
            length: bases.len(),
        })
}

/// `addMismatchReadFeatures`: a stretch of match-or-mismatch, base by base.
#[allow(clippy::too_many_arguments)]
fn add_mismatch_read_features(
    reference_bases: &[u8],
    alignment_start: i32,
    features: &mut Vec<ReadFeature>,
    from_pos_in_read: i32,
    alignment_start_offset: i32,
    read_base_count: i32,
    bases: &[u8],
    base_qualities: Qualities<'_>,
) -> Result<(), ReadFeatureError> {
    // The reference walks its own cursor, which is the read's cursor shifted by the alignment: the
    // reference advances one per read base here because M, X and EQ all consume both.
    let first_position_in_read = from_pos_in_read + 1;
    let first_reference_index = i64::from(alignment_start) + i64::from(alignment_start_offset) - 1;

    for i in 0..read_base_count {
        let one_based_position_in_read = first_position_in_read + i;
        let reference_index = first_reference_index + i64::from(i);
        // Past the end is `N`; before the start is the array access, which throws. The order of
        // the two is the reference's own: only the upper bound is guarded.
        let reference_base = if reference_index >= reference_bases.len() as i64 {
            b'N'
        } else {
            *reference_bases
                .get(usize::try_from(reference_index).map_err(|_| {
                    ReadFeatureError::IndexOutOfBounds {
                        index: reference_index,
                        length: reference_bases.len(),
                    }
                })?)
                .ok_or(ReadFeatureError::IndexOutOfBounds {
                    index: reference_index,
                    length: reference_bases.len(),
                })?
        };

        let read_index = (i + from_pos_in_read) as usize;
        let read_base = *bases
            .get(read_index)
            .ok_or(ReadFeatureError::IndexOutOfBounds {
                index: read_index as i64,
                length: bases.len(),
            })?;

        if read_base != reference_base {
            if is_upper_acgtn(read_base) && is_upper_acgtn(reference_base) {
                features.push(ReadFeature::Substitution {
                    position: one_based_position_in_read,
                    base: read_base,
                    reference_base,
                    code: NO_CODE,
                });
            } else {
                let quality =
                    match base_qualities {
                        Qualities::Missing => MISSING_QUALITY_SCORE,
                        Qualities::Present(qualities) => *qualities.get(read_index).ok_or(
                            ReadFeatureError::IndexOutOfBounds {
                                index: read_index as i64,
                                length: qualities.len(),
                            },
                        )? as i8,
                    };
                features.push(ReadFeature::ReadBase {
                    position: one_based_position_in_read,
                    base: read_base,
                    quality,
                });
            }
        }
    }
    Ok(())
}

/// `getCigarForReadFeatures`: the cigar rebuilt from the features and the read length alone.
///
/// The cigar is not stored anywhere. Going forward, everything that matched became nothing;
/// coming back, **the matches are the gaps**. `gap = position - (last_op_pos + last_op_len)` is
/// the only source of `M` in the output, and nothing in the features says a base matched.
///
/// Four consequences, all measured:
///
/// - **a substitution and a `ReadBase` are both `M`**, so the rebuilt cigar never emits `X` or
///   `=`. A record written with `8X` comes back as `8M`, and that is the only shape in the corpus
///   whose round trip changes;
/// - **a feature that consumes no read bases winds the read cursor back**, `last_op_pos -=
///   length` after a `D`, `N` or `P`, because the bookkeeping is in read space;
/// - **the switch silently ignores what it does not name.** `BaseQualityScore`, `Scores` and
///   `Bases` fall through, and `Bases` carries read bases: a list holding one produces a cigar
///   that does not account for it;
/// - **the read length is what says where the read ends**, and it wins. A feature positioned past
///   it is absorbed, and a read length of 0 takes the accumulated length instead.
///
/// htsjdk guards `readFeatures == null` and `lastOperator != null`; neither can happen here, and
/// the first is unnecessary there too, since an empty list reaches the same single `M` through the
/// empty-list check at the end.
pub fn cigar_for_read_features(features: &[ReadFeature], read_length: i32) -> Cigar {
    let mut elements: Vec<CigarElement> = Vec::new();
    let mut last_operator = Op::M;
    let mut last_op_len: i32 = 0;
    let mut last_op_pos: i32 = 1;

    for feature in features {
        // Everything between the last operator and this feature matched, and is therefore an M.
        let gap = feature.position() - (last_op_pos + last_op_len);
        if gap > 0 {
            if last_operator != Op::M {
                elements.push(element(last_op_len, last_operator));
                last_op_pos += last_op_len;
                last_op_len = gap;
            } else {
                last_op_len += gap;
            }
            last_operator = Op::M;
        }

        let (operator, feature_length) = match feature {
            ReadFeature::Insertion { sequence, .. } => (Op::I, sequence.len() as i32),
            ReadFeature::SoftClip { sequence, .. } => (Op::S, sequence.len() as i32),
            ReadFeature::HardClip { length, .. } => (Op::H, *length),
            ReadFeature::InsertBase { .. } => (Op::I, 1),
            ReadFeature::Deletion { length, .. } => (Op::D, *length),
            ReadFeature::RefSkip { length, .. } => (Op::N, *length),
            ReadFeature::Padding { length, .. } => (Op::P, *length),
            // Both of the mismatch features are an ordinary M here.
            ReadFeature::Substitution { .. } | ReadFeature::ReadBase { .. } => (Op::M, 1),
            // `default: continue`, which drops the three features the switch does not name.
            ReadFeature::BaseQualityScore { .. }
            | ReadFeature::Scores { .. }
            | ReadFeature::Bases { .. } => continue,
        };

        if last_operator != operator {
            if last_op_len > 0 {
                elements.push(element(last_op_len, last_operator));
            }
            last_operator = operator;
            last_op_len = feature_length;
            last_op_pos = feature.position();
        } else {
            last_op_len += feature_length;
        }

        if !operator.consumes_read_bases() {
            last_op_pos -= feature_length;
        }
    }

    if last_operator != Op::M {
        elements.push(element(last_op_len, last_operator));
        if read_length >= last_op_pos + last_op_len {
            elements.push(element(
                read_length - (last_op_len + last_op_pos) + 1,
                Op::M,
            ));
        }
    } else if read_length == 0 || read_length > last_op_pos - 1 {
        let length = if read_length == 0 {
            last_op_len
        } else {
            read_length - last_op_pos + 1
        };
        elements.push(element(length, Op::M));
    }

    if elements.is_empty() {
        return Cigar::new(vec![element(read_length, Op::M)]);
    }
    Cigar::new(elements)
}

fn element(length: i32, op: Op) -> CigarElement {
    CigarElement {
        length: length.max(0) as u32,
        op,
    }
}

/// `getAlignmentEnd`: the span recomputed from the features rather than from the cigar.
pub fn alignment_end(features: &[ReadFeature], alignment_start: i32, read_length: i32) -> i32 {
    let mut span = read_length;
    for feature in features {
        match feature {
            ReadFeature::InsertBase { .. } => span -= 1,
            ReadFeature::Insertion { sequence, .. } => span -= sequence.len() as i32,
            ReadFeature::SoftClip { sequence, .. } => span -= sequence.len() as i32,
            ReadFeature::Deletion { length, .. } => span += length,
            ReadFeature::RefSkip { length, .. } => span += length,
            _ => {}
        }
    }
    alignment_start + span - 1
}

#[cfg(test)]
mod tests {
    use super::*;

    const REFERENCE: &[u8] = b"ACGTACGTACGTMRWSACGTACGT";

    fn cigar(text: &str) -> Cigar {
        let mut elements = Vec::new();
        let mut length = 0u32;
        for byte in text.bytes() {
            if byte.is_ascii_digit() {
                length = length * 10 + u32::from(byte - b'0');
            } else {
                let op = match byte {
                    b'M' => Op::M,
                    b'I' => Op::I,
                    b'D' => Op::D,
                    b'N' => Op::N,
                    b'S' => Op::S,
                    b'H' => Op::H,
                    b'P' => Op::P,
                    b'=' => Op::Eq,
                    b'X' => Op::X,
                    _ => panic!("cigar operator {}", byte as char),
                };
                elements.push(CigarElement { length, op });
                length = 0;
            }
        }
        Cigar::new(elements)
    }

    fn features(text: &str, bases: &[u8], start: i32) -> Vec<ReadFeature> {
        let qualities = vec![40u8; bases.len().max(8)];
        create_read_features(
            &cigar(text),
            start,
            bases,
            Qualities::Present(&qualities),
            REFERENCE,
        )
        .expect("features")
    }

    /// The twelve operators, which is what a reader dispatches on.
    #[test]
    fn the_operators_are_the_reference_operators() {
        let all = [
            (
                ReadFeature::BaseQualityScore {
                    position: 1,
                    quality: 0,
                },
                b'Q',
            ),
            (
                ReadFeature::Bases {
                    position: 1,
                    bases: vec![],
                },
                b'b',
            ),
            (
                ReadFeature::Deletion {
                    position: 1,
                    length: 1,
                },
                b'D',
            ),
            (
                ReadFeature::HardClip {
                    position: 1,
                    length: 1,
                },
                b'H',
            ),
            (
                ReadFeature::InsertBase {
                    position: 1,
                    base: b'A',
                },
                b'i',
            ),
            (
                ReadFeature::Insertion {
                    position: 1,
                    sequence: vec![],
                },
                b'I',
            ),
            (
                ReadFeature::Padding {
                    position: 1,
                    length: 1,
                },
                b'P',
            ),
            (
                ReadFeature::ReadBase {
                    position: 1,
                    base: b'A',
                    quality: 0,
                },
                b'B',
            ),
            (
                ReadFeature::RefSkip {
                    position: 1,
                    length: 1,
                },
                b'N',
            ),
            (
                ReadFeature::Scores {
                    position: 1,
                    scores: vec![],
                },
                b'q',
            ),
            (
                ReadFeature::SoftClip {
                    position: 1,
                    sequence: vec![],
                },
                b'S',
            ),
            (
                ReadFeature::Substitution {
                    position: 1,
                    base: b'A',
                    reference_base: b'C',
                    code: NO_CODE,
                },
                b'X',
            ),
        ];
        for (feature, operator) in all {
            assert_eq!(feature.operator(), operator);
            assert_eq!(feature.position(), 1);
        }
    }

    /// Everything that matches is stored as nothing.
    #[test]
    fn a_perfect_match_produces_no_features() {
        assert!(features("8M", b"ACGTACGT", 1).is_empty());
    }

    /// The positions are one-based, whatever the interface's javadoc says.
    #[test]
    fn the_first_base_of_the_read_is_position_one() {
        assert_eq!(
            features("8M", b"CCGTACGT", 1),
            vec![ReadFeature::Substitution {
                position: 1,
                base: b'C',
                reference_base: b'A',
                code: NO_CODE,
            }]
        );
    }

    /// An insertion becomes one feature per base; a soft clip becomes one feature for all.
    #[test]
    fn an_insertion_is_split_and_a_soft_clip_is_not() {
        let inserted = features("2M3I3M", b"ACTTTTAC", 1);
        let bases: Vec<_> = inserted
            .iter()
            .filter(|f| f.operator() == b'i')
            .map(|f| (f.position(), f.clone()))
            .collect();
        assert_eq!(bases.len(), 3, "one feature per inserted base");
        assert_eq!(bases.iter().map(|(p, _)| *p).collect::<Vec<_>>(), [3, 4, 5]);

        let clipped = features("3S5M", b"TTTACGTA", 4);
        assert_eq!(
            clipped[0],
            ReadFeature::SoftClip {
                position: 1,
                sequence: b"TTT".to_vec(),
            },
            "one feature for every clipped base"
        );
    }

    /// The split is on the alphabet: an ACGTN mismatch is a substitution, anything else is a
    /// ReadBase carrying the quality a second time.
    #[test]
    fn a_mismatch_splits_on_the_alphabet() {
        assert_eq!(
            features("4M", b"ACGM", 1),
            vec![ReadFeature::ReadBase {
                position: 4,
                base: b'M',
                quality: 40,
            }]
        );
        // And the same holds when it is the *reference* that is outside the alphabet.
        let against_iupac = features("4M", b"ACGT", 13);
        assert_eq!(against_iupac.len(), 4);
        assert!(against_iupac.iter().all(|f| f.operator() == b'B'));
    }

    /// The cigar's own claim about match and mismatch is not consulted.
    #[test]
    fn the_cigar_operator_only_says_how_far_to_walk() {
        assert!(
            features("8X", b"ACGTACGT", 1).is_empty(),
            "an X over bases that match emits nothing"
        );
        assert_eq!(
            features("8=", b"TTTTTTTT", 1).len(),
            6,
            "an = over bases that differ emits substitutions"
        );
    }

    /// Past the end of the reference every base is compared against N, so an N out there matches.
    #[test]
    fn past_the_end_of_the_reference_the_base_is_n() {
        let past = features("8M", b"ACGTACGT", 20);
        assert_eq!(past.len(), 8);
        let reference_bases: Vec<u8> = past
            .iter()
            .map(|f| match f {
                ReadFeature::Substitution { reference_base, .. } => *reference_base,
                other => panic!("{other:?}"),
            })
            .collect();
        assert_eq!(reference_bases, *b"TACGTNNN");

        assert_eq!(
            features("4M", b"NNNN", 22).len(),
            3,
            "the base past the end matches the N it is compared against"
        );
    }

    /// A record with no sequence gets Ns, and they mismatch like any other base.
    #[test]
    fn no_sequence_manufactures_ns_that_then_mismatch() {
        let manufactured = create_read_features(
            &cigar("4M"),
            1,
            &[],
            Qualities::Present(&[40, 40, 40, 40]),
            REFERENCE,
        )
        .expect("features");
        assert_eq!(manufactured.len(), 4);
        assert!(manufactured
            .iter()
            .all(|f| matches!(f, ReadFeature::Substitution { base: b'N', .. })));
    }

    /// The missing-quality test is an identity test, so an empty array that is not the singleton
    /// is indexed like any other.
    #[test]
    fn missing_qualities_and_an_empty_array_are_not_the_same_thing() {
        let missing = create_read_features(&cigar("4M"), 1, b"ACGM", Qualities::Missing, REFERENCE)
            .expect("features");
        assert_eq!(
            missing,
            vec![ReadFeature::ReadBase {
                position: 4,
                base: b'M',
                quality: MISSING_QUALITY_SCORE,
            }]
        );

        let empty =
            create_read_features(&cigar("4M"), 1, b"ACGM", Qualities::Present(&[]), REFERENCE);
        assert_eq!(
            empty,
            Err(ReadFeatureError::IndexOutOfBounds {
                index: 3,
                length: 0
            })
        );
        assert_eq!(
            empty.unwrap_err().message(),
            "Index 3 out of bounds for length 0"
        );
    }

    /// Only the upper bound of the reference is guarded.
    #[test]
    fn an_alignment_start_of_zero_reads_before_the_reference() {
        let error = create_read_features(
            &cigar("4M"),
            0,
            b"ACGT",
            Qualities::Present(&[40, 40, 40, 40]),
            REFERENCE,
        )
        .expect_err("refused");
        assert_eq!(
            error.message(),
            format!("Index -1 out of bounds for length {}", REFERENCE.len())
        );
    }

    fn rebuilt(features: &[ReadFeature], read_length: i32) -> String {
        cigar_for_read_features(features, read_length).to_text()
    }

    fn substitution(position: i32) -> ReadFeature {
        ReadFeature::Substitution {
            position,
            base: b'T',
            reference_base: b'A',
            code: NO_CODE,
        }
    }

    /// The matches are the gaps, and nothing else in the features says a base matched.
    #[test]
    fn a_list_of_substitutions_rebuilds_as_one_match() {
        assert_eq!(rebuilt(&[], 8), "8M");
        assert_eq!(rebuilt(&[substitution(4)], 8), "8M");
        assert_eq!(rebuilt(&[substitution(1)], 8), "8M");
        assert_eq!(rebuilt(&[substitution(8)], 8), "8M");
        assert_eq!(rebuilt(&[substitution(2), substitution(7)], 8), "8M");
        assert_eq!(
            rebuilt(
                &[ReadFeature::ReadBase {
                    position: 4,
                    base: b'M',
                    quality: 40,
                }],
                8
            ),
            "8M",
            "a ReadBase is an M like a substitution"
        );
    }

    /// A feature that consumes no read bases winds the read cursor back, which is what keeps the
    /// trailing match at the right length.
    #[test]
    fn a_reference_only_operator_winds_the_read_cursor_back() {
        assert_eq!(
            rebuilt(
                &[ReadFeature::Deletion {
                    position: 5,
                    length: 2
                }],
                8
            ),
            "4M2D4M"
        );
        assert_eq!(
            rebuilt(
                &[ReadFeature::Deletion {
                    position: 1,
                    length: 2
                }],
                8
            ),
            "2D8M",
            "and a deletion at the first position leaves the whole read after it"
        );
        assert_eq!(
            rebuilt(
                &[ReadFeature::Padding {
                    position: 5,
                    length: 2
                }],
                8
            ),
            "4M2P4M"
        );
    }

    /// The three features the switch does not name contribute nothing, including the one that
    /// carries read bases.
    #[test]
    fn the_features_the_switch_ignores_contribute_nothing() {
        assert_eq!(
            rebuilt(
                &[ReadFeature::Bases {
                    position: 1,
                    bases: b"ACGT".to_vec()
                }],
                8
            ),
            "8M",
            "a Bases feature carries read bases and is dropped anyway"
        );
        assert_eq!(
            rebuilt(
                &[
                    ReadFeature::Scores {
                        position: 1,
                        scores: b"IIII".to_vec()
                    },
                    substitution(6)
                ],
                8
            ),
            "8M"
        );
    }

    /// The read length says where the read ends, and it wins over the features.
    #[test]
    fn the_read_length_decides_where_the_read_ends() {
        assert_eq!(
            rebuilt(&[substitution(4)], 0),
            "4M",
            "a read length of 0 takes the accumulated length"
        );
        assert_eq!(
            rebuilt(&[substitution(6)], 3),
            "3M",
            "a feature past the end is absorbed"
        );
        assert_eq!(rebuilt(&[substitution(1)], 1), "1M");
    }

    #[test]
    fn a_clip_keeps_its_side() {
        assert_eq!(
            rebuilt(
                &[ReadFeature::SoftClip {
                    position: 6,
                    sequence: b"TTT".to_vec()
                }],
                8
            ),
            "5M3S"
        );
        assert_eq!(
            rebuilt(
                &[ReadFeature::HardClip {
                    position: 1,
                    length: 2
                }],
                8
            ),
            "2H8M"
        );
    }

    /// The round trip is the identity except where the cigar claimed something the features
    /// cannot carry: X and = both come back as M.
    #[test]
    fn the_round_trip_loses_only_the_x_and_the_equals() {
        for (text, bases, start) in [
            ("2M3I3M", &b"ACTTTTAC"[..], 1),
            ("3S5M", &b"TTTACGTA"[..], 4),
            ("4M2D4M", &b"ACGTCGTA"[..], 1),
            ("2M2I2M2D2M", &b"ACTTACGT"[..], 1),
        ] {
            let built = features(text, bases, start);
            let length = cigar(text).read_length() as i32;
            assert_eq!(rebuilt(&built, length), text, "{text}");
        }
        assert_eq!(rebuilt(&features("8X", b"ACGTACGT", 1), 8), "8M");
        assert_eq!(rebuilt(&features("8=", b"TTTTTTTT", 1), 8), "8M");
    }

    /// The span comes back from the features, not from the cigar.
    #[test]
    fn the_alignment_end_is_recomputed_from_the_features() {
        let inserted = features("2M3I3M", b"ACTTTTAC", 1);
        assert_eq!(alignment_end(&inserted, 1, 8), 5, "three bases inserted");
        let deleted = features("4M2D4M", b"ACGTCGTA", 1);
        assert_eq!(alignment_end(&deleted, 1, 8), 10, "two bases deleted");
        let clipped = features("2S4M2S", b"TTACGTTT", 3);
        assert_eq!(alignment_end(&clipped, 3, 8), 6, "four bases clipped");
    }
}
