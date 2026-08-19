//! Ported from `htsjdk.samtools.reference.FastaReferenceWriter` and the parts of
//! `FastaReferenceWriterBuilder` it validates through (htsjdk 4.2.0).
//!
//! Writes a FASTA, its `.fai` and its `.dict` together, and the three have to agree: every offset
//! in the index is a byte count of the FASTA this writer produced, so a port that wrote the
//! sequence correctly and counted differently would produce a reference no reader could seek in.
//!
//! # Where the numbers in the index come from
//!
//!  * **the offset is taken after the header line is written**, so it counts this sequence's
//!    `>name description\n` and every byte of every sequence before it;
//!  * **bytes-per-line is bases-per-line plus one**, unconditionally. It is not measured from the
//!    lines actually written, so a sequence shorter than one line still records the width it was
//!    opened with;
//!  * **the length is the bases appended**, not the bytes written, so the newlines the writer
//!    inserts are in the offsets and not in the lengths.
//!
//! # One newline per sequence, written at the end
//!
//! `appendBases` writes a separator only when a line is full and more bases are coming, and
//! `closeSequence` writes exactly one at the end. That is why a sequence whose length is a multiple
//! of the line width gets no blank line after it, and why appending in chunks that do not line up
//! with the width produces the same bytes as appending all at once: the breaks come from a running
//! count of bases on the current line rather than from the calls.
//!
//! # The md5 is of the upper-cased bases
//!
//! The digest is updated chunk by chunk with `new String(bases, next, nextLength).toUpperCase()`,
//! so a lower-case sequence keeps its case in the FASTA and hashes as though it did not. The port
//! upper-cases ASCII only, which is what a base array can hold.

use std::fmt::Write as _;

use md5::{Digest, Md5};

/// `FastaReferenceWriter.DEFAULT_BASES_PER_LINE`.
pub const DEFAULT_BASES_PER_LINE: usize = 60;

/// `FastaReferenceWriter.HEADER_START_CHAR`.
const HEADER_START: u8 = b'>';

/// `FastaReferenceWriter.HEADER_NAME_AND_DESCRIPTION_SEPARATOR`.
const NAME_AND_DESCRIPTION_SEPARATOR: u8 = b' ';

/// What the writer refuses, each of them an exception the reference throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FastaWriterError {
    /// `ValidationUtils.nonEmpty`, whose message names the argument rather than the value.
    EmptyName,
    /// A whitespace character anywhere in the name, tab included.
    BlankInName(String),
    /// An ISO control character in the name.
    ControlInName(String),
    /// An ISO control character other than tab in the description.
    ControlInDescription(String),
    /// `startSequence` or `close` with a sequence open and no base appended to it.
    NoBaseAdded,
    /// A sequence name that has already been closed.
    DuplicateName(String),
    /// A byte that is not an IUPAC code.
    InvalidBase(u8),
    /// `appendBases` before any `startSequence`.
    NoSequenceStarted,
    /// A line width of zero or less.
    BasesPerLineNotPositive,
}

impl FastaWriterError {
    /// The exception class the reference throws.
    pub fn java_class(&self) -> &'static str {
        match self {
            FastaWriterError::EmptyName
            | FastaWriterError::BlankInName(_)
            | FastaWriterError::ControlInName(_)
            | FastaWriterError::ControlInDescription(_)
            | FastaWriterError::InvalidBase(_)
            | FastaWriterError::BasesPerLineNotPositive => "java.lang.IllegalArgumentException",
            FastaWriterError::NoBaseAdded
            | FastaWriterError::DuplicateName(_)
            | FastaWriterError::NoSequenceStarted => "java.lang.IllegalStateException",
        }
    }

    /// The message the reference throws with, verbatim.
    pub fn message(&self) -> String {
        match self {
            FastaWriterError::EmptyName => "The string is empty: Sequence name".to_string(),
            FastaWriterError::BlankInName(name) => {
                format!("the input name contains blank characters: '{name}'")
            }
            FastaWriterError::ControlInName(name) => {
                format!("the input name contains control characters: '{name}'")
            }
            FastaWriterError::ControlInDescription(description) => {
                format!("the input name contains non-tab control characters: '{description}'")
            }
            FastaWriterError::NoBaseAdded => "no base was added".to_string(),
            FastaWriterError::DuplicateName(name) => {
                format!("the input sequence name '{name}' has already been added")
            }
            FastaWriterError::InvalidBase(base) => format!(
                "the input sequence contains invalid base calls like: {}",
                *base as char
            ),
            FastaWriterError::NoSequenceStarted => {
                "trying to add bases without starting a sequence".to_string()
            }
            FastaWriterError::BasesPerLineNotPositive => {
                "bases per line must be 1 or greater".to_string()
            }
        }
    }
}

/// `SequenceUtil.isIUPAC`.
///
/// The table holds `ACGT`, the eleven ambiguity codes, and `.`, in both cases. It does **not** hold
/// `X`, whatever the writer's own documentation says about it, and it holds nothing at or above
/// 127: the array is 127 long and a byte outside it is refused before it is indexed.
pub fn is_iupac(base: u8) -> bool {
    if base >= 127 {
        return false;
    }
    matches!(
        base.to_ascii_uppercase(),
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'M'
            | b'R'
            | b'W'
            | b'S'
            | b'Y'
            | b'K'
            | b'V'
            | b'H'
            | b'D'
            | b'B'
            | b'N'
    ) || base == b'.'
}

/// `Character.isISOControl`.
fn is_iso_control(c: char) -> bool {
    let code = c as u32;
    code <= 0x1f || (0x7f..=0x9f).contains(&code)
}

/// The three outputs, which the reference writes to three streams.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastaOutputs {
    /// The FASTA itself.
    pub fasta: Vec<u8>,
    /// The `.fai`.
    pub index: String,
    /// The `.dict`.
    pub dictionary: String,
}

/// `FastaReferenceWriter`.
pub struct FastaReferenceWriter {
    outputs: FastaOutputs,
    default_bases_per_line: usize,
    emit_md5: bool,
    digest: Md5,
    /// Every name already closed, which is what the duplicate check consults.
    names: Vec<String>,
    current_name: Option<String>,
    current_bases_per_line: usize,
    current_line_bases: usize,
    current_bases: u64,
    current_offset: u64,
}

impl FastaReferenceWriter {
    /// `FastaReferenceWriterBuilder.build`, reduced to the two settings that reach the output.
    ///
    /// The `@HD` line is written at construction, not at close: `encodeHeaderLine(false)` runs in
    /// the constructor, so even a writer that is closed without a sequence produces a dictionary
    /// holding one line.
    pub fn new(
        bases_per_line: usize,
        emit_md5: bool,
    ) -> Result<FastaReferenceWriter, FastaWriterError> {
        if bases_per_line == 0 {
            return Err(FastaWriterError::BasesPerLineNotPositive);
        }
        let mut outputs = FastaOutputs::default();
        // `SAMSequenceDictionaryCodec.encodeHeaderLine(false)`, which is the writer's own version
        // rather than any input's.
        outputs.dictionary.push_str("@HD\tVN:1.6\n");
        Ok(FastaReferenceWriter {
            outputs,
            default_bases_per_line: bases_per_line,
            emit_md5,
            digest: Md5::new(),
            names: Vec::new(),
            current_name: None,
            current_bases_per_line: bases_per_line,
            current_line_bases: 0,
            current_bases: 0,
            current_offset: 0,
        })
    }

    /// `startSequence(name)`, at the writer's default width and with no description.
    pub fn start_sequence(&mut self, name: &str) -> Result<(), FastaWriterError> {
        let width = self.default_bases_per_line;
        self.start_sequence_with(name, "", width)
    }

    /// `startSequence(name, description, basesPerLine)`.
    ///
    /// The order of the checks is the reference's and is observable: the name, then the
    /// description, then the width, then the previous sequence is closed, and only then is the
    /// duplicate name refused. So a run that opens a duplicate name with a bad width fails on the
    /// width, and one that opens a duplicate after an empty sequence fails on the empty sequence.
    pub fn start_sequence_with(
        &mut self,
        name: &str,
        description: &str,
        bases_per_line: usize,
    ) -> Result<(), FastaWriterError> {
        check_name(name)?;
        check_description(description)?;
        if bases_per_line == 0 {
            return Err(FastaWriterError::BasesPerLineNotPositive);
        }
        self.close_sequence()?;
        if self.names.iter().any(|existing| existing == name) {
            return Err(FastaWriterError::DuplicateName(name.to_string()));
        }

        self.current_name = Some(name.to_string());
        self.current_bases_per_line = bases_per_line;

        self.outputs.fasta.push(HEADER_START);
        self.outputs.fasta.extend_from_slice(name.as_bytes());
        if !description.is_empty() {
            self.outputs.fasta.push(NAME_AND_DESCRIPTION_SEPARATOR);
            self.outputs.fasta.extend_from_slice(description.as_bytes());
        }
        self.outputs.fasta.push(b'\n');
        // Taken here, after the header: the offset of the first base.
        self.current_offset = self.outputs.fasta.len() as u64;

        if self.emit_md5 {
            self.digest = Md5::new();
        }
        Ok(())
    }

    /// `appendBases(bases)`.
    pub fn append_bases(&mut self, bases: &[u8]) -> Result<(), FastaWriterError> {
        if self.current_name.is_none() {
            return Err(FastaWriterError::NoSequenceStarted);
        }
        for base in bases {
            if !is_iupac(*base) {
                return Err(FastaWriterError::InvalidBase(*base));
            }
        }

        let mut next = 0;
        while next < bases.len() {
            if self.current_line_bases == self.current_bases_per_line {
                self.outputs.fasta.push(b'\n');
                self.current_line_bases = 0;
            }
            let length =
                (bases.len() - next).min(self.current_bases_per_line - self.current_line_bases);
            self.outputs
                .fasta
                .extend_from_slice(&bases[next..next + length]);
            if self.emit_md5 {
                let chunk: Vec<u8> = bases[next..next + length]
                    .iter()
                    .map(|base| base.to_ascii_uppercase())
                    .collect();
                self.digest.update(&chunk);
            }
            self.current_line_bases += length;
            next += length;
        }
        self.current_bases += bases.len() as u64;
        Ok(())
    }

    /// `closeSequence`: the index row, the dictionary row, and the sequence's one newline.
    fn close_sequence(&mut self) -> Result<(), FastaWriterError> {
        let Some(name) = self.current_name.clone() else {
            return Ok(());
        };
        if self.current_bases == 0 {
            return Err(FastaWriterError::NoBaseAdded);
        }
        self.names.push(name.clone());

        let _ = writeln!(
            self.outputs.index,
            "{name}\t{}\t{}\t{}\t{}",
            self.current_bases,
            self.current_offset,
            self.current_bases_per_line,
            self.current_bases_per_line + 1
        );

        let _ = write!(
            self.outputs.dictionary,
            "@SQ\tSN:{name}\tLN:{}",
            self.current_bases
        );
        if self.emit_md5 {
            let digest = std::mem::replace(&mut self.digest, Md5::new()).finalize();
            let mut hex = String::with_capacity(32);
            for byte in digest.iter() {
                let _ = write!(hex, "{byte:02x}");
            }
            let _ = write!(self.outputs.dictionary, "\tM5:{hex}");
        }
        self.outputs.dictionary.push('\n');

        self.outputs.fasta.push(b'\n');
        self.current_bases = 0;
        self.current_line_bases = 0;
        self.current_name = None;
        Ok(())
    }

    /// `close`: closes the open sequence, if any, and hands back the three outputs.
    pub fn close(mut self) -> Result<FastaOutputs, FastaWriterError> {
        self.close_sequence()?;
        Ok(self.outputs)
    }
}

/// `checkSequenceName`, whose two refusals are in this order: blank first, then control.
fn check_name(name: &str) -> Result<(), FastaWriterError> {
    if name.is_empty() {
        return Err(FastaWriterError::EmptyName);
    }
    for c in name.chars() {
        if c.is_whitespace() {
            return Err(FastaWriterError::BlankInName(name.to_string()));
        }
        if is_iso_control(c) {
            return Err(FastaWriterError::ControlInName(name.to_string()));
        }
    }
    Ok(())
}

/// `checkDescription`: a tab is the only control character allowed.
fn check_description(description: &str) -> Result<(), FastaWriterError> {
    if description.is_empty() {
        return Ok(());
    }
    for c in description.chars() {
        if is_iso_control(c) && c != '\t' {
            return Err(FastaWriterError::ControlInDescription(
                description.to_string(),
            ));
        }
    }
    Ok(())
}
