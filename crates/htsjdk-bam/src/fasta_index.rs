//! The `.fai` and the indexed FASTA reader that uses it.
//!
//! Ported from `htsjdk.samtools.reference.FastaSequenceIndex`,
//! `FastaSequenceIndexEntry` and `AbstractIndexedFastaSequenceFile.getSubsequenceAt`
//! at htsjdk 4.2.0.
//!
//! # Why this is here rather than borrowed
//!
//! `gatk-engine` reached a FASTA through `noodles-fasta`, which is a fine reader of the format and
//! the wrong thing to depend on for the same reason `rust-htslib` is: what a GATK tool sees is
//! htsjdk's reader, and the two agree only until they do not. The arithmetic below is where they
//! could differ -- a query that spans a line boundary is answered by seeking past the newline
//! bytes, and the number of those bytes comes from the `.fai`, not from looking.
//!
//! # The five columns
//!
//! `name`, `size`, `location`, `basesPerLine`, `bytesPerLine`. The last two differ by the line
//! terminator's length, which is how a file with CRLF line endings is read correctly without ever
//! being scanned for one: `bytesPerLine - basesPerLine` is the terminator, whatever it is.
//!
//! # The bounds are htsjdk's, including the one that looks wrong
//!
//! `start > stop + 1` is the malformed-query test, so an EMPTY query -- `start == stop + 1` -- is
//! legal and answers with no bases. `stop > size` is "past end of contig". A start below one is
//! not checked at all: it is the caller's job, and `GATKTool` does it.

use std::io::{Read, Seek, SeekFrom};

/// One line of the `.fai`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaIndexEntry {
    pub name: String,
    /// The contig's length in bases.
    pub size: u64,
    /// The byte offset of the contig's first base.
    pub location: u64,
    pub bases_per_line: u64,
    pub bytes_per_line: u64,
}

/// What the reader refuses, with htsjdk's own messages where it has one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FastaIndexError {
    /// The `.fai` line did not have five usable columns.
    Malformed(String),
    /// `SAMException("Malformed query; start point %d lies after end point %d")`.
    MalformedQuery { start: i64, stop: i64 },
    /// `SAMException("Query asks for data past end of contig")`.
    PastEndOfContig,
    /// The contig is not in the index.
    UnknownContig(String),
    /// The file could not be read.
    Io(String),
}

impl FastaIndexError {
    pub fn message(&self) -> String {
        match self {
            FastaIndexError::Malformed(line) => format!("Malformed fasta index line: {line}"),
            FastaIndexError::MalformedQuery { start, stop } => {
                format!("Malformed query; start point {start} lies after end point {stop}")
            }
            FastaIndexError::PastEndOfContig => {
                "Query asks for data past end of contig".to_string()
            }
            FastaIndexError::UnknownContig(contig) => {
                format!("Unable to find entry for contig: {contig}")
            }
            FastaIndexError::Io(detail) => detail.clone(),
        }
    }
}

/// `FastaSequenceIndex`: the entries in the file's own order, which is the dictionary order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastaIndex {
    pub entries: Vec<FastaIndexEntry>,
}

impl FastaIndex {
    /// `parseIndexFile`: five whitespace-separated columns per line, blank lines skipped.
    pub fn parse(text: &str) -> Result<FastaIndex, FastaIndexError> {
        let mut entries = Vec::new();
        for line in text.lines() {
            if line.trim().is_empty() {
                continue;
            }
            let columns: Vec<&str> = line.split_whitespace().collect();
            if columns.len() < 5 {
                return Err(FastaIndexError::Malformed(line.to_string()));
            }
            let number = |value: &str| -> Result<u64, FastaIndexError> {
                value
                    .parse()
                    .map_err(|_| FastaIndexError::Malformed(line.to_string()))
            };
            entries.push(FastaIndexEntry {
                name: columns[0].to_string(),
                size: number(columns[1])?,
                location: number(columns[2])?,
                bases_per_line: number(columns[3])?,
                bytes_per_line: number(columns[4])?,
            });
        }
        Ok(FastaIndex { entries })
    }

    pub fn get(&self, contig: &str) -> Option<&FastaIndexEntry> {
        self.entries.iter().find(|entry| entry.name == contig)
    }
}

/// `IndexedFastaSequenceFile`: a FASTA and its `.fai`, queried by interval.
pub struct IndexedFasta<R: Read + Seek> {
    source: R,
    index: FastaIndex,
}

impl<R: Read + Seek> IndexedFasta<R> {
    pub fn new(source: R, index: FastaIndex) -> Self {
        Self { source, index }
    }

    pub fn index(&self) -> &FastaIndex {
        &self.index
    }

    /// `getSubsequenceAt(contig, start, stop)`: the bases of `[start, stop]`, 1-based inclusive.
    ///
    /// The loop is htsjdk's, in the form the arithmetic takes when the terminator length is known
    /// rather than searched for: read a run of at most one line's worth of bases, step over the
    /// terminator, repeat. A query inside one line reads once.
    pub fn query(
        &mut self,
        contig: &str,
        start: i64,
        stop: i64,
    ) -> Result<Vec<u8>, FastaIndexError> {
        if start > stop + 1 {
            return Err(FastaIndexError::MalformedQuery { start, stop });
        }
        let entry = self
            .index
            .get(contig)
            .ok_or_else(|| FastaIndexError::UnknownContig(contig.to_string()))?
            .clone();
        if stop as u64 > entry.size {
            return Err(FastaIndexError::PastEndOfContig);
        }
        let length = (stop - start + 1).max(0) as usize;
        let mut target = Vec::with_capacity(length);
        if length == 0 {
            return Ok(target);
        }

        let bases_per_line = entry.bases_per_line.max(1);
        let bytes_per_line = entry.bytes_per_line.max(bases_per_line);
        let terminator = bytes_per_line - bases_per_line;

        // The first base's byte offset: whole lines, then the offset within the line.
        let mut offset = ((start as u64 - 1) / bases_per_line) * bytes_per_line
            + (start as u64 - 1) % bases_per_line;

        let mut buffer = vec![0u8; bases_per_line as usize];
        while target.len() < length {
            let position_in_line = offset % bytes_per_line;
            // A `.fai` whose location lands inside a terminator: step to the next base.
            if position_in_line >= bases_per_line {
                offset += bytes_per_line - position_in_line;
                continue;
            }
            let run = std::cmp::min(
                bases_per_line - position_in_line,
                (length - target.len()) as u64,
            ) as usize;
            self.source
                .seek(SeekFrom::Start(entry.location + offset))
                .map_err(|error| FastaIndexError::Io(error.to_string()))?;
            self.source
                .read_exact(&mut buffer[..run])
                .map_err(|error| FastaIndexError::Io(error.to_string()))?;
            target.extend_from_slice(&buffer[..run]);
            offset += run as u64
                + if position_in_line + run as u64 == bases_per_line {
                    terminator
                } else {
                    0
                };
        }
        Ok(target)
    }
}

impl IndexedFasta<std::fs::File> {
    /// The FASTA at `path` with the `.fai` beside it, which is
    /// `ReferenceSequenceFileFactory.getReferenceSequenceFile` with `requireIndex`.
    pub fn open(path: &std::path::Path) -> Result<Self, FastaIndexError> {
        let index_path = {
            let mut name = path.as_os_str().to_os_string();
            name.push(".fai");
            std::path::PathBuf::from(name)
        };
        let text = std::fs::read_to_string(&index_path)
            .map_err(|error| FastaIndexError::Io(format!("{}: {error}", index_path.display())))?;
        let index = FastaIndex::parse(&text)?;
        let file = std::fs::File::open(path)
            .map_err(|error| FastaIndexError::Io(format!("{}: {error}", path.display())))?;
        Ok(IndexedFasta::new(file, index))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Two contigs, six bases per line, LF terminators.
    fn fixture() -> (Cursor<Vec<u8>>, FastaIndex) {
        let text = ">chr1\nACGTAC\nGTACGT\nAC\n>chr2\nTTTTTT\nGG\n";
        // `>chr1\n` is 6 bytes, three lines of 7, 7 and 3 follow, then `>chr2\n`: chr2's first
        // base is at 29. Getting this wrong by hand is what the `.fai` exists to prevent.
        let index =
            FastaIndex::parse("chr1\t14\t6\t6\t7\nchr2\t8\t29\t6\t7\n").expect("the index parses");
        (Cursor::new(text.as_bytes().to_vec()), index)
    }

    #[test]
    fn a_query_inside_one_line_reads_once() {
        let (source, index) = fixture();
        let mut fasta = IndexedFasta::new(source, index);
        assert_eq!(fasta.query("chr1", 1, 6).unwrap(), b"ACGTAC");
        assert_eq!(fasta.query("chr1", 2, 4).unwrap(), b"CGT");
    }

    #[test]
    fn a_query_across_lines_steps_over_the_terminator() {
        let (source, index) = fixture();
        let mut fasta = IndexedFasta::new(source, index);
        assert_eq!(fasta.query("chr1", 1, 14).unwrap(), b"ACGTACGTACGTAC");
        assert_eq!(fasta.query("chr1", 5, 9).unwrap(), b"ACGTA");
    }

    #[test]
    fn the_second_contig_is_found_by_its_own_offset() {
        let (source, index) = fixture();
        let mut fasta = IndexedFasta::new(source, index);
        assert_eq!(fasta.query("chr2", 1, 8).unwrap(), b"TTTTTTGG");
    }

    #[test]
    fn an_empty_query_is_legal_and_a_reversed_one_is_not() {
        let (source, index) = fixture();
        let mut fasta = IndexedFasta::new(source, index);
        assert_eq!(fasta.query("chr1", 5, 4).unwrap(), b"");
        assert_eq!(
            fasta.query("chr1", 6, 4),
            Err(FastaIndexError::MalformedQuery { start: 6, stop: 4 })
        );
    }

    #[test]
    fn past_the_end_is_refused_by_the_index_rather_than_by_the_read() {
        let (source, index) = fixture();
        let mut fasta = IndexedFasta::new(source, index);
        assert_eq!(
            fasta.query("chr1", 1, 15),
            Err(FastaIndexError::PastEndOfContig)
        );
        assert_eq!(
            fasta.query("chrX", 1, 2),
            Err(FastaIndexError::UnknownContig("chrX".to_string()))
        );
    }
}
