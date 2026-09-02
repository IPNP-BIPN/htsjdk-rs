//! `SamFiles.findIndex`: the index file a data file's own name implies.
//!
//! A consumer that opens `reads.bam` and looks for `reads.bam.bai` finds nothing when htsjdk wrote
//! `reads.bai`, and answers an interval query with zero records rather than with a refusal. Both
//! names are htsjdk's, and which of the two wins is not the order a reader would guess.
//!
//! # The order, which is not one rule but three
//!
//! For a name ending in `.bam`, the extension is **replaced** first and appended second:
//! `reads.bai`, then `reads.csi`, then `reads.bam.bai`, then `reads.bam.csi`. The consequence is
//! that a `.csi` beside a `.bam` beats a `.bam.bai` beside the same file, because the replaced
//! pair is exhausted before the appended one is tried at all.
//!
//! For a name ending in `.cram` the replacement is `.crai` only, and it is followed by the
//! *appended* `.cram.crai` before the shared fallthrough. There is no `reads.csi` for a CRAM: the
//! replaced `.csi` belongs to the `.bam` branch alone.
//!
//! For every other name, only the fallthrough runs: `reads.txt.bai`, then `reads.txt.csi`.
//!
//! # Existence is a regular file
//!
//! `Files.isRegularFile` follows symbolic links and answers false for a directory, so a directory
//! named `reads.bai` is not an index and the search moves on. And when nothing is found beside the
//! path as given, the path is resolved through its symbolic links and the whole search runs again
//! at the real location: an index that lives next to the target of a link is found, and one that
//! lives next to the link is found first.
//!
//! Ported from `htsjdk.samtools.SamFiles`.

use std::path::{Path, PathBuf};

/// `FileExtensions`, for the six names this search is spelled out of.
const BAM: &str = ".bam";
const BAI_INDEX: &str = ".bai";
const CRAM: &str = ".cram";
const CRAM_INDEX: &str = ".crai";
const CSI: &str = ".csi";

/// `Files.isRegularFile`, which follows links and refuses a directory.
fn is_regular_file(path: &Path) -> bool {
    std::fs::metadata(path)
        .map(|meta| meta.is_file())
        .unwrap_or(false)
}

/// `Path.resolveSibling`: the name in place of the last component.
///
/// A path with no parent resolves to the bare name, which is what Java answers for the same case.
fn sibling(path: &Path, name: &str) -> PathBuf {
    match path.parent() {
        Some(parent) => parent.join(name),
        None => PathBuf::from(name),
    }
}

/// The first of the candidate names that is a regular file.
fn first_existing(path: &Path, names: &[String]) -> Option<PathBuf> {
    names
        .iter()
        .map(|name| sibling(path, name))
        .find(|candidate| is_regular_file(candidate))
}

/// `lookForIndex`, without the symlink retry.
fn look_for_index(sam_path: &Path) -> Option<PathBuf> {
    let file_name = sam_path.file_name()?.to_str()?.to_string();

    // The replaced pair, which is tried before anything is appended.
    if let Some(stem) = file_name.strip_suffix(BAM) {
        let found = first_existing(
            sam_path,
            &[format!("{stem}{BAI_INDEX}"), format!("{stem}{CSI}")],
        );
        if found.is_some() {
            return found;
        }
    } else if let Some(stem) = file_name.strip_suffix(CRAM) {
        // A CRAM's replaced name is `.crai` alone, and its appended `.crai` comes before the
        // fallthrough rather than after it.
        let found = first_existing(
            sam_path,
            &[
                format!("{stem}{CRAM_INDEX}"),
                format!("{file_name}{CRAM_INDEX}"),
            ],
        );
        if found.is_some() {
            return found;
        }
    }

    // The appended pair, which every name reaches.
    first_existing(
        sam_path,
        &[
            format!("{file_name}{BAI_INDEX}"),
            format!("{file_name}{CSI}"),
        ],
    )
}

/// `SamFiles.findIndex`: the index beside the path, or the one beside what the path really is.
///
/// `None` where the reference returns null, which includes the case where the search found nothing
/// and the case where resolving the links failed.
pub fn find_index(sam_path: &Path) -> Option<PathBuf> {
    look_for_index(sam_path).or_else(|| {
        // `toRealPath` throws on a path that does not exist, and the reference catches it and
        // answers null rather than letting it out.
        let canonical = std::fs::canonicalize(sam_path).ok()?;
        look_for_index(&canonical)
    })
}
