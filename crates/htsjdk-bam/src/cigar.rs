//! CIGAR operators and their binary encoding.
//!
//! Ported from `htsjdk.samtools.CigarOperator`, `Cigar` and `BinaryCigarCodec`.

/// `CigarOperator`, in declaration order, which is also its binary encoding order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Match or mismatch.
    M = 0,
    /// Insertion vs. the reference.
    I = 1,
    /// Deletion vs. the reference.
    D = 2,
    /// Skipped region from the reference.
    N = 3,
    /// Soft clip.
    S = 4,
    /// Hard clip.
    H = 5,
    /// Padding.
    P = 6,
    /// Matches the reference.
    Eq = 7,
    /// Mismatches the reference.
    X = 8,
}

impl Op {
    /// `CigarOperator.enumToBinary`.
    pub fn to_binary(self) -> u32 {
        self as u32
    }

    /// `CigarOperator.binaryToEnum`.
    pub fn from_binary(b: u32) -> Option<Op> {
        Some(match b {
            0 => Op::M,
            1 => Op::I,
            2 => Op::D,
            3 => Op::N,
            4 => Op::S,
            5 => Op::H,
            6 => Op::P,
            7 => Op::Eq,
            8 => Op::X,
            _ => return None,
        })
    }

    /// `CigarOperator.enumToCharacter`. Note `Eq` prints as `=`, not `E`.
    pub fn to_char(self) -> u8 {
        match self {
            Op::M => b'M',
            Op::I => b'I',
            Op::D => b'D',
            Op::N => b'N',
            Op::S => b'S',
            Op::H => b'H',
            Op::P => b'P',
            Op::Eq => b'=',
            Op::X => b'X',
        }
    }

    /// `CigarOperator.consumesReadBases`.
    pub fn consumes_read_bases(self) -> bool {
        matches!(self, Op::M | Op::I | Op::S | Op::Eq | Op::X)
    }

    /// `CigarOperator.consumesReferenceBases`.
    pub fn consumes_reference_bases(self) -> bool {
        matches!(self, Op::M | Op::D | Op::N | Op::Eq | Op::X)
    }
}

/// `CigarElement`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarElement {
    pub length: u32,
    pub op: Op,
}

/// `Cigar`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Cigar {
    pub elements: Vec<CigarElement>,
}

impl Cigar {
    pub fn new(elements: Vec<CigarElement>) -> Self {
        Cigar { elements }
    }

    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// `Cigar.getReadLength`.
    pub fn read_length(&self) -> u32 {
        self.elements
            .iter()
            .filter(|e| e.op.consumes_read_bases())
            .map(|e| e.length)
            .sum()
    }

    /// `Cigar.getReferenceLength`.
    pub fn reference_length(&self) -> u32 {
        self.elements
            .iter()
            .filter(|e| e.op.consumes_reference_bases())
            .map(|e| e.length)
            .sum()
    }

    /// `BinaryCigarCodec.encode`: `length << 4 | op`, one `u32` per element.
    pub fn encode(&self) -> Vec<u32> {
        self.elements
            .iter()
            .map(|e| (e.length << 4) | e.op.to_binary())
            .collect()
    }

    /// `BinaryCigarCodec.decode`.
    pub fn decode(binary: &[u32]) -> Option<Cigar> {
        binary
            .iter()
            .map(|&v| Op::from_binary(v & 0x0F).map(|op| CigarElement { length: v >> 4, op }))
            .collect::<Option<Vec<_>>>()
            .map(Cigar::new)
    }

    /// The text form, as it appears in SAM.
    pub fn to_text(&self) -> String {
        if self.elements.is_empty() {
            return "*".to_string();
        }
        let mut s = String::new();
        for e in &self.elements {
            s.push_str(&e.length.to_string());
            s.push(e.op.to_char() as char);
        }
        s
    }
}

/// `CigarUtil.softClipEndOfRead(clipFrom, oldCigar)`: soft-clip the read from the 1-based position
/// `clip_from` to its end, rewriting the cigar element list.
///
/// Ported from `htsjdk.samtools.CigarUtil.softClipEndOfRead` and the `clipEndOfRead` /
/// `mergeClippingCigarElement` it delegates to, tag 4.2.0. Used by CleanSam to trim alignments that
/// hang off the end of the reference, and by the read clippers generally. This is the soft-clip
/// entry; the hard-clip path of `mergeClippingCigarElement` is reproduced faithfully but is not
/// exercised here.
pub fn soft_clip_end_of_read(clip_from: i32, old_cigar: &[CigarElement]) -> Vec<CigarElement> {
    clip_end_of_read(clip_from, old_cigar, Op::S)
}

/// `CigarUtil.clipEndOfRead`.
///
/// Walks the elements until the one that reaches or straddles `clip_from`, copying earlier elements
/// verbatim and handing the boundary element to [`merge_clipping_cigar_element`]. The clipped region
/// is always at the end, so nothing follows it.
pub fn clip_end_of_read(
    clip_from: i32,
    old_cigar: &[CigarElement],
    clipping_operator: Op,
) -> Vec<CigarElement> {
    // clippedBases = CoordMath.getLength(clipFrom, Cigar.getReadLength(oldCigar)) = readLen - clipFrom + 1.
    let read_length: i32 = old_cigar
        .iter()
        .filter(|c| c.op.consumes_read_bases())
        .map(|c| c.length as i32)
        .sum();
    let clipped_bases = read_length - clip_from + 1;

    let mut new_cigar: Vec<CigarElement> = Vec::new();
    let mut pos = 1i32;
    let last = old_cigar[old_cigar.len() - 1];
    let trailing_hard_clip_bases = if last.op == Op::H {
        last.length as i32
    } else {
        0
    };

    for c in old_cigar {
        let op = c.op;
        let length = if op.consumes_read_bases() {
            c.length as i32
        } else {
            0
        };
        let end_pos = pos + length - 1; // same as pos on the next iteration

        if end_pos < clip_from - 1 {
            // Before the clip point: copy verbatim.
            new_cigar.push(*c);
        } else {
            // Adjacent to or straddling the boundary; the rest is clipped.
            merge_clipping_cigar_element(
                &mut new_cigar,
                *c,
                (clip_from - 1) - (pos - 1),
                clipped_bases,
                clipping_operator,
                trailing_hard_clip_bases,
            );
            break;
        }
        pos = end_pos + 1;
    }
    new_cigar
}

/// `CigarUtil.mergeClippingCigarElement`.
fn merge_clipping_cigar_element(
    new_cigar: &mut Vec<CigarElement>,
    original: CigarElement,
    relative_clipped_position: i32,
    clipped_bases: i32,
    new_clipping_operator: Op,
    trailing_hard_clipped_bases: i32,
) {
    let original_operator = original.op;
    let mut clip_amount = clipped_bases;
    if new_clipping_operator == Op::H {
        clip_amount += trailing_hard_clipped_bases;
    }
    if original_operator.consumes_read_bases() {
        if (original_operator.consumes_reference_bases() || new_clipping_operator == Op::H)
            && relative_clipped_position > 0
        {
            new_cigar.push(CigarElement {
                length: relative_clipped_position as u32,
                op: original_operator,
            });
        }
        if !(original_operator.consumes_reference_bases() || new_clipping_operator == Op::H)
            || original_operator == new_clipping_operator
        {
            clip_amount = clipped_bases + relative_clipped_position;
        }
    } else if relative_clipped_position != 0 {
        panic!("Unexpected non-0 relativeClippedPosition {relative_clipped_position}");
    }
    new_cigar.push(CigarElement {
        length: clip_amount as u32,
        op: new_clipping_operator,
    });
    if new_clipping_operator == Op::S && trailing_hard_clipped_bases > 0 {
        new_cigar.push(CigarElement {
            length: trailing_hard_clipped_bases as u32,
            op: Op::H,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL: [Op; 9] = [
        Op::M,
        Op::I,
        Op::D,
        Op::N,
        Op::S,
        Op::H,
        Op::P,
        Op::Eq,
        Op::X,
    ];

    fn m(len: u32, op: Op) -> CigarElement {
        CigarElement { length: len, op }
    }

    #[test]
    fn soft_clipping_a_plain_match_splits_it() {
        // 36M clipped from base 30 -> 29M7S, verified against htsjdk in the conformance corpus.
        let out = soft_clip_end_of_read(30, &[m(36, Op::M)]);
        assert_eq!(Cigar::new(out).to_text(), "29M7S");
    }

    #[test]
    fn soft_clipping_across_an_insertion_absorbs_it() {
        // 10M5I21M clipped from base 12: the boundary lands inside the insertion, which is folded
        // into the soft clip (10M26S), not left as a separate element.
        let out = soft_clip_end_of_read(12, &[m(10, Op::M), m(5, Op::I), m(21, Op::M)]);
        assert_eq!(Cigar::new(out).to_text(), "10M26S");
    }

    #[test]
    fn the_binary_codes_are_the_declaration_order() {
        for (i, op) in ALL.iter().enumerate() {
            assert_eq!(op.to_binary(), i as u32);
            assert_eq!(Op::from_binary(i as u32), Some(*op));
        }
        assert_eq!(Op::from_binary(9), None);
    }

    /// `EQ` is spelled `=` in SAM. Emitting `E` would produce a CIGAR string that no reader
    /// accepts, which is the harmless failure; the dangerous one is emitting the wrong
    /// *binary* code, which stays readable.
    #[test]
    fn eq_prints_as_an_equals_sign() {
        assert_eq!(Op::Eq.to_char(), b'=');
        assert_eq!(Op::Eq.to_binary(), 7);
    }

    /// The two consumption predicates are what drive read length, reference length and
    /// therefore the indexing bin. Stated exhaustively so a wrong row cannot hide.
    #[test]
    fn consumption_table_is_exhaustive_and_exact() {
        let expect = [
            (Op::M, true, true),
            (Op::I, true, false),
            (Op::D, false, true),
            (Op::N, false, true),
            (Op::S, true, false),
            (Op::H, false, false),
            (Op::P, false, false),
            (Op::Eq, true, true),
            (Op::X, true, true),
        ];
        for (op, read, reference) in expect {
            assert_eq!(op.consumes_read_bases(), read, "{op:?} read bases");
            assert_eq!(
                op.consumes_reference_bases(),
                reference,
                "{op:?} reference bases"
            );
        }
    }

    #[test]
    fn lengths_ignore_the_operators_that_do_not_consume() {
        // 10S 50M 5I 20D 10H
        let c = Cigar::new(vec![
            CigarElement {
                length: 10,
                op: Op::S,
            },
            CigarElement {
                length: 50,
                op: Op::M,
            },
            CigarElement {
                length: 5,
                op: Op::I,
            },
            CigarElement {
                length: 20,
                op: Op::D,
            },
            CigarElement {
                length: 10,
                op: Op::H,
            },
        ]);
        assert_eq!(c.read_length(), 65, "S + M + I");
        assert_eq!(c.reference_length(), 70, "M + D");
        assert_eq!(c.to_text(), "10S50M5I20D10H");
    }

    #[test]
    fn binary_encoding_round_trips() {
        let c = Cigar::new(vec![
            CigarElement {
                length: 100,
                op: Op::M,
            },
            CigarElement {
                length: 3,
                op: Op::Eq,
            },
        ]);
        let bin = c.encode();
        assert_eq!(bin, vec![(100 << 4), (3 << 4) | 7]);
        assert_eq!(Cigar::decode(&bin), Some(c));
    }

    #[test]
    fn an_empty_cigar_is_a_star() {
        assert_eq!(Cigar::default().to_text(), "*");
        assert_eq!(Cigar::default().encode(), Vec::<u32>::new());
    }
}
