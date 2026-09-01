//! `htsjdk.samtools.SBIIndexWriter`: the splitting index, which says where a BAM may be cut.
//!
//! An SBI is a list of virtual file offsets, one every `granularity` records, so a reader can split
//! a BAM into pieces that each begin at a record boundary. It is not a BAI: it answers "where may I
//! start reading" rather than "where do the records for this region live".
//!
//! Two details decide the bytes, and both are counters rather than format:
//!
//! * **The record count and the offset count are different numbers.** `processRecord` writes an
//!   offset when `recordCount++ % granularity == 0`, so the first record always writes one, and
//!   `finish` then writes the **final** offset as well. An index over `n` records at granularity
//!   `g` therefore holds `ceil(n / g) + 1` offsets, and the header's `totalNumberOfRecords` is `n`,
//!   which is neither of the two.
//! * **The offsets must be non-decreasing**, and an out-of-order one is an
//!   `IllegalArgumentException` naming both in hex rather than a silently sorted list.
//!
//! `gatk-rs` ports this inside `CreateHadoopBamSplittingIndex`, which is the tool that drives it.

/// `SBIIndex.SBI_MAGIC`.
pub const SBI_MAGIC: [u8; 4] = *b"SBI\x01";

/// `SBIIndexWriter.DEFAULT_GRANULARITY`.
pub const DEFAULT_GRANULARITY: u64 = 4096;

/// `SBIIndexWriter`, accumulating offsets rather than spilling them to a temporary file.
///
/// htsjdk writes the offsets to a temp file because an index at low granularity can hold 10^8 of
/// them and a `List<Long>` resizes badly. That is a memory strategy, not a format: the bytes it
/// writes are the same ones a vector produces, and the port says so rather than reproducing the
/// temporary file.
pub struct SbiIndexWriter {
    granularity: u64,
    offsets: Vec<u64>,
    previous: Option<u64>,
    record_count: u64,
}

impl SbiIndexWriter {
    pub fn new(granularity: u64) -> Self {
        SbiIndexWriter {
            granularity,
            offsets: Vec::new(),
            previous: None,
            record_count: 0,
        }
    }

    /// `processRecord(virtualOffset)`: an offset is kept every `granularity` records, starting with
    /// the first.
    pub fn process_record(&mut self, virtual_offset: u64) -> Result<(), String> {
        let take = self.record_count.is_multiple_of(self.granularity);
        self.record_count += 1;
        if take {
            self.write_virtual_offset(virtual_offset)?;
        }
        Ok(())
    }

    /// `writeVirtualOffset`, including its refusal.
    fn write_virtual_offset(&mut self, virtual_offset: u64) -> Result<(), String> {
        if let Some(previous) = self.previous {
            if previous > virtual_offset {
                // The message is htsjdk's, `%#x` and all: a caller comparing refusals compares text.
                return Err(format!(
                    "Offsets not in order: {previous:#x} > {virtual_offset:#x}"
                ));
            }
        }
        self.offsets.push(virtual_offset);
        self.previous = Some(virtual_offset);
        Ok(())
    }

    /// `finish(finalVirtualOffset, dataFileLength, md5, uuid)`: the whole file's bytes.
    ///
    /// `md5` and `uuid` are sixteen zero bytes when absent, which is what `EMPTY_MD5` and
    /// `EMPTY_UUID` are, and a wrong length is htsjdk's `IllegalArgumentException`.
    pub fn finish(
        mut self,
        final_virtual_offset: u64,
        data_file_length: u64,
        md5: Option<&[u8]>,
        uuid: Option<&[u8]>,
    ) -> Result<Vec<u8>, String> {
        if let Some(md5) = md5 {
            if md5.len() != 16 {
                return Err(format!("Invalid MD5 length: {}", md5.len()));
            }
        }
        if let Some(uuid) = uuid {
            if uuid.len() != 16 {
                return Err(format!("Invalid UUID length: {}", uuid.len()));
            }
        }
        // The final offset is written whatever the granularity says, so an index always ends on the
        // end of the data.
        self.write_virtual_offset(final_virtual_offset)?;

        let mut out = Vec::with_capacity(4 + 8 + 16 + 16 + 8 + 8 + 8 + self.offsets.len() * 8);
        out.extend_from_slice(&SBI_MAGIC);
        out.extend_from_slice(&data_file_length.to_le_bytes());
        out.extend_from_slice(md5.unwrap_or(&[0u8; 16]));
        out.extend_from_slice(uuid.unwrap_or(&[0u8; 16]));
        out.extend_from_slice(&self.record_count.to_le_bytes());
        out.extend_from_slice(&self.granularity.to_le_bytes());
        out.extend_from_slice(&(self.offsets.len() as u64).to_le_bytes());
        for offset in &self.offsets {
            out.extend_from_slice(&offset.to_le_bytes());
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_first_record_always_writes_an_offset() {
        let mut writer = SbiIndexWriter::new(4096);
        writer.process_record(100).unwrap();
        assert_eq!(writer.offsets, vec![100]);
    }

    #[test]
    fn an_index_holds_one_more_offset_than_the_granularity_implies() {
        let mut writer = SbiIndexWriter::new(2);
        for record in 0..5u64 {
            writer.process_record(record * 10).unwrap();
        }
        let bytes = writer.finish(999, 12345, None, None).unwrap();
        // ceil(5 / 2) = 3 offsets from the records, plus the final one.
        // Layout: magic 0..4, fileLength 4..12, md5 12..28, uuid 28..44, records 44..52,
        // granularity 52..60, offset count 60..68, then the offsets.
        let count = u64::from_le_bytes(bytes[60..68].try_into().unwrap());
        assert_eq!(count, 4);
        // And the header's record count is the records, which is neither 3 nor 4.
        let records = u64::from_le_bytes(bytes[44..52].try_into().unwrap());
        assert_eq!(records, 5);
        assert_eq!(bytes.len(), 68 + 4 * 8);
    }

    #[test]
    fn an_out_of_order_offset_is_refused_in_htsjdks_words() {
        let mut writer = SbiIndexWriter::new(1);
        writer.process_record(0x20).unwrap();
        let error = writer.process_record(0x10).unwrap_err();
        assert_eq!(error, "Offsets not in order: 0x20 > 0x10");
    }

    #[test]
    fn an_empty_index_still_carries_the_final_offset() {
        let writer = SbiIndexWriter::new(4096);
        let bytes = writer.finish(0, 0, None, None).unwrap();
        assert_eq!(&bytes[0..4], &SBI_MAGIC);
        let count = u64::from_le_bytes(bytes[60..68].try_into().unwrap());
        assert_eq!(count, 1, "the final offset is written even with no records");
    }

    #[test]
    fn a_wrong_length_digest_is_refused() {
        let writer = SbiIndexWriter::new(4096);
        let error = writer.finish(0, 0, Some(&[0u8; 8]), None).unwrap_err();
        assert_eq!(error, "Invalid MD5 length: 8");
    }
}
