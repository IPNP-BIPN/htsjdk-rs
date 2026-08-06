//! The three flag words a CRAM record carries, and the mate chain that fills two of them.
//!
//! Ported from `htsjdk.samtools.cram.structure.CRAMCompressionRecord` at htsjdk 4.2.0: its flag
//! predicates, `restoreMateInfo`, `setNextMate`, `computeInsertSize` and `setToDetachedState`.
//!
//! A BAM record has one flag word. A CRAM record has three: the BAM flags it will become, a CRAM
//! flags byte the format adds, and a mate flags byte that duplicates two bits of the first at
//! different positions.
//!
//! # Two bits are stored twice, at different positions
//!
//! Mate unmapped is `0x2` in the mate flags and `0x8` in the BAM flags; mate reverse strand is
//! `0x1` and `0x20`. Setting either through the record sets both, and nothing keeps them in step if
//! the two words are set independently. htsjdk's own comment calls the duplication redundant and
//! the specification ambiguous about it.
//!
//! # The two narrow words are masked to a byte on the way out
//!
//! `getCRAMFlags` and `getMateFlags` return `0xFF & field`, so a value above 255 comes back
//! truncated rather than refused. Measured: 511 reads back as 255, 256 as 0, -1 as 255.
//!
//! # Restoring mate info walks a ring
//!
//! Each record takes its mate's position, reference and two flags from the next; the last takes
//! them from the first, so a chain of N becomes a cycle. A chain of one is left untouched, because
//! the walk returns before it starts.
//!
//! # The template length is computed once and negated once
//!
//! On the first and the last of the chain. Every record between them keeps whatever it was
//! constructed with, which for a triple measured 250, 0 and -250.

/// `CF_QS_PRESERVED_AS_ARRAY`.
pub const CF_QS_PRESERVED_AS_ARRAY: i32 = 0x1;
/// `CF_DETACHED`: the mate is stored literally rather than as a record offset.
pub const CF_DETACHED: i32 = 0x2;
/// `CF_HAS_MATE_DOWNSTREAM`.
pub const CF_HAS_MATE_DOWNSTREAM: i32 = 0x4;
/// `CF_UNKNOWN_BASES`.
pub const CF_UNKNOWN_BASES: i32 = 0x8;

/// `MF_MATE_NEG_STRAND`, which is `0x20` in the BAM flags and `0x1` here.
pub const MF_MATE_NEG_STRAND: i32 = 0x1;
/// `MF_MATE_UNMAPPED`, which is `0x8` in the BAM flags and `0x2` here.
pub const MF_MATE_UNMAPPED: i32 = 0x2;

/// The BAM flags the record reads for itself, by their SAM values.
pub const READ_PAIRED: i32 = 0x1;
pub const READ_UNMAPPED: i32 = 0x4;
pub const MATE_UNMAPPED: i32 = 0x8;
pub const READ_REVERSE_STRAND: i32 = 0x10;
pub const MATE_REVERSE_STRAND: i32 = 0x20;
pub const FIRST_OF_PAIR: i32 = 0x40;
pub const SECOND_OF_PAIR: i32 = 0x80;
pub const SECONDARY_ALIGNMENT: i32 = 0x100;

/// `SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX`.
pub const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;
/// `SAMRecord.NO_ALIGNMENT_START`.
pub const NO_ALIGNMENT_START: i32 = 0;

/// The three flag words, and the questions a record answers from them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Flags {
    pub bam: i32,
    pub cram: i32,
    pub mate: i32,
}

impl Flags {
    /// `getCRAMFlags`, masked to a byte.
    pub fn cram_flags(&self) -> i32 {
        self.cram & 0xFF
    }

    /// `getMateFlags`, masked to a byte.
    pub fn mate_flags(&self) -> i32 {
        self.mate & 0xFF
    }

    pub fn is_detached(&self) -> bool {
        self.cram & CF_DETACHED != 0
    }

    pub fn has_mate_downstream(&self) -> bool {
        self.cram & CF_HAS_MATE_DOWNSTREAM != 0
    }

    pub fn is_force_preserve_quality_scores(&self) -> bool {
        self.cram & CF_QS_PRESERVED_AS_ARRAY != 0
    }

    pub fn is_unknown_bases(&self) -> bool {
        self.cram & CF_UNKNOWN_BASES != 0
    }

    pub fn is_read_paired(&self) -> bool {
        self.bam & READ_PAIRED != 0
    }

    pub fn is_segment_unmapped(&self) -> bool {
        self.bam & READ_UNMAPPED != 0
    }

    pub fn is_first_segment(&self) -> bool {
        self.bam & FIRST_OF_PAIR != 0
    }

    pub fn is_last_segment(&self) -> bool {
        self.bam & SECOND_OF_PAIR != 0
    }

    pub fn is_secondary_alignment(&self) -> bool {
        self.bam & SECONDARY_ALIGNMENT != 0
    }

    pub fn is_negative_strand(&self) -> bool {
        self.bam & READ_REVERSE_STRAND != 0
    }

    /// Read from the mate flags, which is where the reference reads it even though the same fact
    /// is in the BAM flags.
    pub fn is_mate_unmapped(&self) -> bool {
        self.mate & MF_MATE_UNMAPPED != 0
    }

    pub fn is_mate_negative_strand(&self) -> bool {
        self.mate & MF_MATE_NEG_STRAND != 0
    }

    /// `setMateUnmapped`, which writes both words.
    pub fn set_mate_unmapped(&mut self, value: bool) {
        self.mate = set_bit(self.mate, MF_MATE_UNMAPPED, value);
        self.bam = set_bit(self.bam, MATE_UNMAPPED, value);
    }

    /// `setMateNegativeStrand`, which writes both words.
    pub fn set_mate_negative_strand(&mut self, value: bool) {
        self.mate = set_bit(self.mate, MF_MATE_NEG_STRAND, value);
        self.bam = set_bit(self.bam, MATE_REVERSE_STRAND, value);
    }
}

fn set_bit(word: i32, bit: i32, value: bool) -> i32 {
    if value {
        word | bit
    } else {
        word & !bit
    }
}

/// As much of a record as the mate code touches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MateRecord {
    pub flags: Flags,
    pub reference_index: i32,
    pub alignment_start: i32,
    /// Derived from the read features, not stored. `AlignmentContext.NO_ALIGNMENT_END` where the
    /// record is not placed.
    pub alignment_end: i32,
    pub mate_reference_index: i32,
    pub mate_alignment_start: i32,
    pub template_size: i32,
    pub records_to_next_fragment: i32,
}

impl MateRecord {
    /// `isPlaced`, which consults only the alignment start. The two aberrant combinations are
    /// warnings in the reference and nothing here, because a warning is not an output.
    pub fn is_placed(&self) -> bool {
        self.alignment_start != NO_ALIGNMENT_START
    }

    /// `setToDetachedState`: detached, no mate downstream, and no distance to one.
    pub fn set_to_detached_state(&mut self) {
        self.flags.cram = set_bit(self.flags.cram, CF_DETACHED, true);
        self.flags.cram = set_bit(self.flags.cram, CF_HAS_MATE_DOWNSTREAM, false);
        self.records_to_next_fragment = -1;
    }

    /// `setNextMate`: take the mate's position, reference and two flags from `next`.
    ///
    /// A mate on no reference loses its position: the start is forced to zero after being taken,
    /// whatever the mate's own start was.
    fn set_next_mate(&mut self, next: &MateRecord) {
        self.mate_alignment_start = next.alignment_start;
        self.flags
            .set_mate_unmapped(next.flags.is_segment_unmapped());
        self.flags
            .set_mate_negative_strand(next.flags.is_negative_strand());
        self.mate_reference_index = next.reference_index;
        if self.mate_reference_index == NO_ALIGNMENT_REFERENCE_INDEX {
            self.mate_alignment_start = NO_ALIGNMENT_START;
        }
    }
}

/// `restoreMateInfo`, over a chain given in order.
///
/// Each record takes its mate from the next and the last from the first, so the chain becomes a
/// ring. The template length is computed between the first and the last, stored on the first and
/// negated on the last; everything between keeps what it had.
///
/// A chain of fewer than two is left exactly as it was, because the reference returns before the
/// walk when there is no next segment.
pub fn restore_mate_info(records: &mut [MateRecord]) {
    if records.len() < 2 {
        return;
    }

    for index in 0..records.len() - 1 {
        let next = records[index + 1].clone();
        records[index].set_next_mate(&next);
    }
    let first = records[0].clone();
    let last = records.len() - 1;
    records[last].set_next_mate(&first);

    let template_length = compute_insert_size(&records[0], &records[last]);
    records[0].template_size = template_length;
    records[last].template_size = -template_length;
}

/// `computeInsertSize`.
///
/// Zero if either end is unmapped or the two are on different references. Otherwise the distance
/// between the two 5' ends plus a sign, so the shortest template is 1 and never 0. A
/// negative-strand record measures from its alignment end.
pub fn compute_insert_size(first: &MateRecord, last: &MateRecord) -> i32 {
    if first.flags.is_segment_unmapped()
        || last.flags.is_segment_unmapped()
        || first.reference_index != last.reference_index
    {
        return 0;
    }

    let first_five_prime = if first.flags.is_negative_strand() {
        first.alignment_end
    } else {
        first.alignment_start
    };
    let last_five_prime = if last.flags.is_negative_strand() {
        last.alignment_end
    } else {
        last.alignment_start
    };

    let adjustment = if last_five_prime >= first_five_prime {
        1
    } else {
        -1
    };
    last_five_prime
        .wrapping_sub(first_five_prime)
        .wrapping_add(adjustment)
}
