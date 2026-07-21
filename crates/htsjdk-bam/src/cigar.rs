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
