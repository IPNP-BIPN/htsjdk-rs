//! Canonical Huffman, the fourth codec written on the CRAM core bit stream.
//!
//! Ported from `htsjdk.samtools.cram.encoding.core.huffmanUtils.HuffmanCanoncialCodeGenerator`,
//! `CanonicalHuffmanIntegerEncoding` and `CanonicalHuffmanByteEncoding` at htsjdk 4.2.0.
//!
//! # The code words are derived, not stored
//!
//! A CRAM file carries no Huffman tree. It carries an alphabet and one code word length per
//! symbol, and both writer and reader rebuild the same code words from that pair alone: symbols
//! are grouped by length, sorted inside each group by their own natural order, and handed
//! consecutive integers, shifted left whenever the length grows. So the alphabet's order in the
//! file does not matter, only the pairing of each symbol with its length.
//!
//! The tree that produced the lengths is the writer's business and is not part of the format.
//! htsjdk has one ([`HuffmanParamsCalculator`]), but nothing in the library calls it.
//!
//! [`HuffmanParamsCalculator`]: https://github.com/samtools/htsjdk/blob/4.2.0/src/main/java/htsjdk/samtools/cram/encoding/core/huffmanUtils/HuffmanParamsCalculator.java
//!
//! # Byte symbols sort signed
//!
//! The grouping is a `TreeSet`, and Java's `Byte` compares signed, so in a byte alphabet `0x80`
//! sorts below `0x01` and takes the earlier code word. The port keeps byte symbols as [`i8`] for
//! exactly that reason.
//!
//! # A one-symbol alphabet writes nothing
//!
//! Its single length is zero, so every write emits zero bits and every read consumes zero bits and
//! returns the only symbol. Such a stream's core block is empty, and the count of symbols written
//! is not recoverable from it.
//!
//! # The overflow check counts set bits, not width
//!
//! What refuses an impossible length table is `Integer.bitCount(codeValue) > bitLength`, so it
//! fires later than a check on the code word's width would: three symbols at length 1 are accepted
//! and the third is given code word `2`, which does not fit in one bit, while four are refused.
//!
//! # An unmatched code word runs off the end of the table
//!
//! Reading walks the lengths in ascending order, consuming only the difference between one length
//! and the next, and matching against a table indexed by code word. That table is sized to the
//! largest code word, not to the largest value the bits can hold, so a truncated or foreign stream
//! comes out as an `ArrayIndexOutOfBoundsException`. The codec's own "unable to map" message is
//! reachable only with an empty alphabet.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Display;
use std::hash::Hash;

use crate::bit_stream::{BitError, BitInputStream, BitOutputStream};
use crate::varint::{read_unsigned_itf8, write_unsigned_itf8, RuntimeEof};

/// What canonical Huffman refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HuffmanError {
    /// A length table the canonical assignment cannot honour. The count is of set bits in the code
    /// word, not of its width, which is why the check fires where it does.
    BitLengthOutOfRange { bit_count: i32, symbol: i64 },
    /// More lengths than symbols. The scan is driven by the lengths, so it indexes past the
    /// alphabet; the other way round the surplus symbols are silently dropped.
    SymbolIndexOutOfBounds { index: usize, length: usize },
    /// A code word past the end of the code-word-to-symbol table.
    CodeIndexOutOfBounds { index: i32, length: usize },
    /// The codec's own message, which only an empty alphabet reaches.
    UnableToMap,
    /// A symbol that is not in the alphabet. The reference's message ends in "null" always: it
    /// prints the code word it just found to be absent.
    UnknownSymbol { symbol: i64 },
    /// The bit stream underneath refused.
    Bits(BitError),
}

impl HuffmanError {
    pub fn message(&self) -> String {
        match self {
            HuffmanError::BitLengthOutOfRange { bit_count, symbol } => {
                format!("Bit length ({bit_count}) for symbol ({symbol}) out of range")
            }
            HuffmanError::SymbolIndexOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
            HuffmanError::CodeIndexOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
            HuffmanError::UnableToMap => {
                "Unable to map huffman code from input stream to a valid symbol".to_string()
            }
            HuffmanError::UnknownSymbol { symbol } => format!(
                "Attempt to write a symbol ({symbol}) that is not in the symbol alphabet for this \
                 huffman encoder (found code word null)."
            ),
            HuffmanError::Bits(error) => error.message(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            HuffmanError::BitLengthOutOfRange { .. } => "IllegalArgumentException",
            HuffmanError::SymbolIndexOutOfBounds { .. } => "IndexOutOfBoundsException",
            HuffmanError::CodeIndexOutOfBounds { .. } => "ArrayIndexOutOfBoundsException",
            HuffmanError::UnableToMap | HuffmanError::UnknownSymbol { .. } => "RuntimeException",
            HuffmanError::Bits(error) => error.java_exception(),
        }
    }
}

impl From<BitError> for HuffmanError {
    fn from(error: BitError) -> Self {
        HuffmanError::Bits(error)
    }
}

/// A symbol of a Huffman alphabet: [`i32`] for the integer flavour, [`i8`] for the byte one.
///
/// `into_message` is how the symbol reaches an error message, where the reference formats it with
/// `%d` whichever flavour it is.
pub trait Symbol: Copy + Ord + Hash + Display {
    fn into_message(self) -> i64;
}

impl Symbol for i32 {
    fn into_message(self) -> i64 {
        i64::from(self)
    }
}

impl Symbol for i8 {
    fn into_message(self) -> i64 {
        i64::from(self)
    }
}

/// One symbol's code word and its width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HuffmanBitCode<S> {
    pub symbol: S,
    pub code_word: i32,
    pub bit_length: i32,
}

/// The code words rebuilt from an alphabet and one length per symbol.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalHuffman<S: Symbol> {
    /// Sorted by length, then by code word, which is the order reading walks.
    codes: Vec<HuffmanBitCode<S>>,
    by_symbol: HashMap<S, HuffmanBitCode<S>>,
    /// Code word to its index in `codes`, or `-1`. Sized to the largest code word, which is what
    /// makes an unmatched code word an out-of-bounds index rather than a miss.
    code_to_index: Vec<i32>,
}

impl<S: Symbol> CanonicalHuffman<S> {
    /// `new HuffmanCanoncialCodeGenerator(new HuffmanParams(symbols, lengths))`.
    pub fn new(symbols: &[S], bit_lengths: &[i32]) -> Result<Self, HuffmanError> {
        let codes = canonical_code_words(symbols, bit_lengths)?;

        let mut by_symbol = HashMap::with_capacity(codes.len());
        let mut largest = 0i32;
        for code in &codes {
            by_symbol.insert(code.symbol, *code);
            largest = largest.max(code.code_word);
        }

        // `new int[maxBitCode + 1]`, filled with -1. An empty alphabet still gets one slot.
        let mut code_to_index = vec![-1i32; largest as usize + 1];
        for (index, code) in codes.iter().enumerate() {
            code_to_index[code.code_word as usize] = index as i32;
        }

        Ok(Self {
            codes,
            by_symbol,
            code_to_index,
        })
    }

    /// The code words, in the order reading walks them.
    pub fn code_words(&self) -> &[HuffmanBitCode<S>] {
        &self.codes
    }

    /// `write`, which emits the symbol's code word and returns the number of bits it took.
    pub fn write(&self, out: &mut BitOutputStream, symbol: S) -> Result<i64, HuffmanError> {
        let code = self
            .by_symbol
            .get(&symbol)
            .ok_or(HuffmanError::UnknownSymbol {
                symbol: symbol.into_message(),
            })?;
        out.write_int_bits(code.code_word, code.bit_length)?;
        Ok(i64::from(code.bit_length))
    }

    /// `read`, which walks the lengths in ascending order and consumes only the difference between
    /// one length and the next.
    pub fn read(&self, input: &mut BitInputStream<'_>) -> Result<S, HuffmanError> {
        let mut previous_length = 0i32;
        let mut code_word = 0i32;
        let mut index = 0usize;

        while index < self.codes.len() {
            let length = self.codes[index].bit_length;
            // Java's `<<` on an int masks the count to five bits; the difference is never that
            // large here, but the arithmetic is kept rather than assumed.
            code_word = code_word.wrapping_shl((length - previous_length) as u32);
            code_word |= input.read_bits(length - previous_length)?;
            previous_length = length;

            let slot = self.code_to_index.get(code_word as usize).copied().ok_or(
                HuffmanError::CodeIndexOutOfBounds {
                    index: code_word,
                    length: self.code_to_index.len(),
                },
            )?;
            if slot > -1 && self.codes[slot as usize].bit_length == length {
                return Ok(self.codes[slot as usize].symbol);
            }

            // Advance to the end of the code words of this length. The reference reads `get(j + 1)`
            // before testing `j` against the size, so this runs off the end for a code word inside
            // the table that no length matches; measured, the table is too small for that to
            // happen before the index above does.
            let mut scan = index;
            loop {
                let next = self.codes.get(scan + 1).ok_or_else(|| {
                    HuffmanError::SymbolIndexOutOfBounds {
                        index: scan + 1,
                        length: self.codes.len(),
                    }
                })?;
                if next.bit_length != length {
                    break;
                }
                index += 1;
                scan += 1;
            }
            index += 1;
        }

        Err(HuffmanError::UnableToMap)
    }
}

/// `getCanonicalCodeWords`: group by length, sort inside the group, hand out consecutive integers.
///
/// The loop is driven by the lengths, not by the symbols. A length table shorter than the alphabet
/// therefore drops the surplus symbols in silence, and one longer indexes past the alphabet.
fn canonical_code_words<S: Symbol>(
    symbols: &[S],
    bit_lengths: &[i32],
) -> Result<Vec<HuffmanBitCode<S>>, HuffmanError> {
    let mut by_length: BTreeMap<i32, BTreeSet<S>> = BTreeMap::new();
    for (index, length) in bit_lengths.iter().enumerate() {
        let symbol = symbols
            .get(index)
            .copied()
            .ok_or(HuffmanError::SymbolIndexOutOfBounds {
                index,
                length: symbols.len(),
            })?;
        by_length.entry(*length).or_default().insert(symbol);
    }

    let mut codes = Vec::with_capacity(bit_lengths.len());
    let mut code_length = 0i32;
    let mut code_value = -1i32;

    for (bit_length, group) in by_length {
        for symbol in group {
            code_value += 1;
            let delta = bit_length - code_length;
            if delta != 0 {
                code_value = code_value.wrapping_shl(delta as u32);
                code_length += delta;
            }
            // `Integer.bitCount`, not the code word's width. On a negative code word Java counts
            // the set bits of the two's complement, which `count_ones` does too.
            let bit_count = code_value.count_ones() as i32;
            if bit_count > bit_length {
                return Err(HuffmanError::BitLengthOutOfRange {
                    bit_count,
                    symbol: symbol.into_message(),
                });
            }
            codes.push(HuffmanBitCode {
                symbol,
                code_word: code_value,
                bit_length,
            });
        }
    }

    // The reference's comparator subtracts, which the port does not need to reproduce: the values
    // it subtracts are code word lengths and code words, both non-negative here.
    codes.sort_by_key(|code| (code.bit_length, code.code_word));
    Ok(codes)
}

/// `CanonicalHuffmanIntegerEncoding.toSerializedEncodingParams`: both halves ITF8.
pub fn serialize_integer_params(symbols: &[i32], bit_lengths: &[i32]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&write_unsigned_itf8(symbols.len() as i32).0);
    for symbol in symbols {
        out.extend_from_slice(&write_unsigned_itf8(*symbol).0);
    }
    out.extend_from_slice(&write_unsigned_itf8(bit_lengths.len() as i32).0);
    for length in bit_lengths {
        out.extend_from_slice(&write_unsigned_itf8(*length).0);
    }
    out
}

/// `CanonicalHuffmanByteEncoding.toSerializedEncodingParams`: the symbols are raw bytes, the
/// lengths ITF8. The asymmetry is the format's, not an oversight of the port.
pub fn serialize_byte_params(symbols: &[i8], bit_lengths: &[i32]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&write_unsigned_itf8(symbols.len() as i32).0);
    for symbol in symbols {
        out.push(*symbol as u8);
    }
    out.extend_from_slice(&write_unsigned_itf8(bit_lengths.len() as i32).0);
    for length in bit_lengths {
        out.extend_from_slice(&write_unsigned_itf8(*length).0);
    }
    out
}

/// What a parsed encoding carries: the alphabet and its lengths, still unpaired with code words.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HuffmanParams<S> {
    pub symbols: Vec<S>,
    pub bit_lengths: Vec<i32>,
}

/// One ITF8 from `bytes` at `cursor`, advancing it.
///
/// Past the end the slice is empty and [`read_unsigned_itf8`] refuses, which is what the
/// reference's `ByteBuffer` does; a short but non-empty slice returns a silently wrong number
/// there and here alike.
fn next_itf8(bytes: &[u8], cursor: &mut usize) -> Result<i32, RuntimeEof> {
    let (value, used) = read_unsigned_itf8(bytes.get(*cursor..).unwrap_or(&[]))?;
    *cursor += used;
    Ok(value)
}

/// `CanonicalHuffmanIntegerEncoding.fromSerializedEncodingParams`.
pub fn parse_integer_params(bytes: &[u8]) -> Result<HuffmanParams<i32>, RuntimeEof> {
    let mut cursor = 0usize;

    let symbol_count = next_itf8(bytes, &mut cursor)?;
    let mut symbols = Vec::new();
    for _ in 0..symbol_count {
        symbols.push(next_itf8(bytes, &mut cursor)?);
    }
    let length_count = next_itf8(bytes, &mut cursor)?;
    let mut bit_lengths = Vec::new();
    for _ in 0..length_count {
        bit_lengths.push(next_itf8(bytes, &mut cursor)?);
    }
    Ok(HuffmanParams {
        symbols,
        bit_lengths,
    })
}

/// `CanonicalHuffmanByteEncoding.fromSerializedEncodingParams`, whose symbols are raw bytes.
pub fn parse_byte_params(bytes: &[u8]) -> Result<HuffmanParams<i8>, RuntimeEof> {
    let mut cursor = 0usize;

    let symbol_count = next_itf8(bytes, &mut cursor)?;
    let mut symbols = Vec::new();
    for _ in 0..symbol_count {
        let byte = bytes.get(cursor).copied().ok_or(RuntimeEof)?;
        symbols.push(byte as i8);
        cursor += 1;
    }

    let length_count = next_itf8(bytes, &mut cursor)?;
    let mut bit_lengths = Vec::new();
    for _ in 0..length_count {
        bit_lengths.push(next_itf8(bytes, &mut cursor)?);
    }
    Ok(HuffmanParams {
        symbols,
        bit_lengths,
    })
}
