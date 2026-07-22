//! `SAMUtils.calculateReadGroupRecordChecksum`.
//!
//! Ported from `htsjdk.samtools.SAMUtils.calculateReadGroupRecordChecksum` and the
//! `SAMHeaderRecordComparator` it sorts with, tag 4.2.0. This is the whole of what Picard's
//! `CalculateReadGroupChecksum` computes: an MD5 over the read groups' attributes, which changes
//! whenever a read group is added, removed, or edited.
//!
//! Two ordering rules decide the bytes fed to the digest, and both are reproduced exactly.
//!
//! **The read groups are sorted** by a fixed tag priority: `PU, LB, DT, SM, CN, PL, DS, ID`. A null
//! (absent) tag sorts before any present value, and values compare as strings. The `ID` tag is on
//! the list but is inert: htsjdk stores a read group's ID as a separate field, not among the
//! attributes `getAttribute` reads, so `getAttribute("ID")` is null for every record and the tie is
//! never broken by it. The comment in htsjdk calls it there "just in case". Reproduced, including
//! its inertness, since the sort is stable and equal records keep their header order.
//!
//! **Within each read group the attributes are sorted by key** (a `TreeMap`), and for each the key
//! bytes then the value bytes (UTF-8) are fed to the digest. The `ID` is excluded, which is
//! automatic here because it is not one of the stored attributes.
//!
//! The result is `new BigInteger(1, digest).toString(16)` left-padded to 32 characters, which is
//! exactly the lowercase hex of the 16-byte digest with its leading zeros kept.

use md5::{Digest, Md5};

use crate::header::ReadGroup;

/// The tag priority `SAMHeaderRecordComparator` is constructed with for read groups.
const SORT_TAGS: [&str; 8] = ["PU", "LB", "DT", "SM", "CN", "PL", "DS", "ID"];

/// `SAMHeaderRecordComparator.compare` over the read-group sort tags.
///
/// A null value sorts first; two nulls tie and move to the next tag; otherwise the string order of
/// the values decides. `ID` is read from the attributes (where it never is), so it is inert.
fn compare(left: &ReadGroup, right: &ReadGroup) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    for tag in SORT_TAGS {
        let l = left.attributes.get(tag);
        let r = right.attributes.get(tag);
        match (l, r) {
            (None, None) => continue,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some(lv), Some(rv)) => match lv.cmp(rv) {
                Ordering::Equal => continue,
                other => return other,
            },
        }
    }
    Ordering::Equal
}

/// `SAMUtils.calculateReadGroupRecordChecksum`: the 32-character lowercase MD5 hex of the read
/// groups' attributes, in the tool's sort order.
pub fn calculate_read_group_record_checksum(read_groups: &[ReadGroup]) -> String {
    // Collections.sort is stable, and Rust's sort_by is stable too, so equal records keep their
    // header order.
    let mut sorted: Vec<&ReadGroup> = read_groups.iter().collect();
    sorted.sort_by(|a, b| compare(a, b));

    let mut digest = Md5::new();
    for rg in &sorted {
        // TreeMap<String, String>: attributes visited in ascending key order.
        let mut attrs: Vec<(&str, &str)> = rg.attributes.iter().collect();
        attrs.sort_by(|a, b| a.0.cmp(b.0));
        for (key, value) in attrs {
            // Redundant here (ID is not a stored attribute), kept for parity with htsjdk.
            if key != "ID" {
                digest.update(key.as_bytes());
                digest.update(value.as_bytes());
            }
        }
    }

    let bytes = digest.finalize();
    let mut hex = String::with_capacity(32);
    for b in bytes {
        hex.push_str(&format!("{b:02x}"));
    }
    hex
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rg(id: &str, attrs: &[(&str, &str)]) -> ReadGroup {
        let mut r = ReadGroup::new(id);
        for (k, v) in attrs {
            r.attributes.set(k, v);
        }
        r
    }

    #[test]
    fn a_known_single_read_group_hashes_to_its_md5() {
        // One read group with SM:sample1 LB:lib1. The digest is over the key/value bytes in key
        // order: "LB"+"lib1" then "SM"+"sample1". Verified against Java's SAMUtils.
        let rgs = [rg("A", &[("SM", "sample1"), ("LB", "lib1")])];
        let mut d = Md5::new();
        d.update(b"LB");
        d.update(b"lib1");
        d.update(b"SM");
        d.update(b"sample1");
        let expect: String = d.finalize().iter().map(|b| format!("{b:02x}")).collect();
        assert_eq!(calculate_read_group_record_checksum(&rgs), expect);
        assert_eq!(expect.len(), 32);
    }

    #[test]
    fn the_hash_is_independent_of_read_group_input_order() {
        // Two read groups differing in LB sort by LB, so header order does not matter.
        let a = rg("A", &[("LB", "libB"), ("SM", "s")]);
        let b = rg("B", &[("LB", "libA"), ("SM", "s")]);
        let forward = calculate_read_group_record_checksum(&[a.clone(), b.clone()]);
        let reverse = calculate_read_group_record_checksum(&[b, a]);
        assert_eq!(forward, reverse);
    }

    #[test]
    fn a_missing_tag_sorts_before_a_present_one() {
        use std::cmp::Ordering;
        let with_pu = rg("A", &[("PU", "unit")]);
        let without = rg("B", &[("LB", "lib")]);
        // `without` has no PU, so it sorts before `with_pu`.
        assert_eq!(compare(&without, &with_pu), Ordering::Less);
    }
}
