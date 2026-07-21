//! The SAM text header.
//!
//! Ported from `htsjdk.samtools.SAMTextHeaderCodec.encode` and the record classes it walks:
//! `AbstractSAMHeaderRecord`, `SAMFileHeader`, `SAMSequenceRecord`, `SAMReadGroupRecord`,
//! `SAMProgramRecord`.
//!
//! The header is text, so a divergence here is visible rather than hidden. It is still easy to
//! get wrong, because the one thing that decides the byte order of every line is a detail of
//! the Java collection used to hold the attributes:
//!
//! `AbstractSAMHeaderRecord` stores them in a **`LinkedHashMap`**. That means insertion order,
//! not sorted order, and it means that overwriting an existing key leaves it in its
//! **original** position. Both a `BTreeMap` (sorted) and a naive remove-then-append (moves to
//! the end) produce a header that is correct SAM and different bytes.

use std::fmt::Write as _;

/// `SAMFileHeader.CURRENT_VERSION`.
pub const CURRENT_VERSION: &str = "1.6";

/// `SAMFileHeader.VERSION_TAG`.
pub const VERSION_TAG: &str = "VN";

/// `SAMTextHeaderCodec.HEADER_LINE_START`.
const HEADER_LINE_START: char = '@';

/// An attribute list with `LinkedHashMap` semantics.
///
/// Insertion order is preserved, and re-setting an existing key updates its value **in place**
/// rather than moving it to the end. Setting a key to `None` removes it, matching
/// `AbstractSAMHeaderRecord.setAttribute(key, null)`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Attributes {
    entries: Vec<(String, String)>,
}

impl Attributes {
    pub fn new() -> Self {
        Self::default()
    }

    /// `AbstractSAMHeaderRecord.setAttribute`.
    pub fn set(&mut self, key: &str, value: &str) {
        match self.entries.iter_mut().find(|(k, _)| k == key) {
            // LinkedHashMap.put on an existing key does not reposition it.
            Some(slot) => slot.1 = value.to_string(),
            None => self.entries.push((key.to_string(), value.to_string())),
        }
    }

    /// `setAttribute(key, null)`.
    pub fn remove(&mut self, key: &str) {
        self.entries.retain(|(k, _)| k != key);
    }

    pub fn get(&self, key: &str) -> Option<&str> {
        self.entries
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    pub fn iter(&self) -> impl Iterator<Item = (&str, &str)> {
        self.entries.iter().map(|(k, v)| (k.as_str(), v.as_str()))
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// `SAMSequenceRecord`. `SN` and `LN` are separate fields, written before the attributes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    pub name: String,
    pub length: i32,
    pub attributes: Attributes,
}

impl SequenceRecord {
    pub fn new(name: &str, length: i32) -> Self {
        SequenceRecord {
            name: name.to_string(),
            length,
            attributes: Attributes::new(),
        }
    }
}

/// `SAMReadGroupRecord`. `ID` is a separate field, written first.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGroup {
    pub id: String,
    pub attributes: Attributes,
}

impl ReadGroup {
    pub fn new(id: &str) -> Self {
        ReadGroup {
            id: id.to_string(),
            attributes: Attributes::new(),
        }
    }
}

/// `SAMProgramRecord`. `ID` is a separate field, written first.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProgramRecord {
    pub id: String,
    pub attributes: Attributes,
}

impl ProgramRecord {
    pub fn new(id: &str) -> Self {
        ProgramRecord {
            id: id.to_string(),
            attributes: Attributes::new(),
        }
    }
}

/// `SAMFileHeader`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SamHeader {
    /// The `@HD` attributes. `SAMFileHeader()` sets `VN` first, so `VN` leads every header.
    pub attributes: Attributes,
    pub sequences: Vec<SequenceRecord>,
    pub read_groups: Vec<ReadGroup>,
    pub programs: Vec<ProgramRecord>,
    /// Comments, each already carrying its `@CO\t` prefix as `addComment` leaves them.
    pub comments: Vec<String>,
}

impl Default for SamHeader {
    /// `new SAMFileHeader()`, which sets `VN` to [`CURRENT_VERSION`] in its constructor.
    fn default() -> Self {
        let mut attributes = Attributes::new();
        attributes.set(VERSION_TAG, CURRENT_VERSION);
        SamHeader {
            attributes,
            sequences: Vec::new(),
            read_groups: Vec::new(),
            programs: Vec::new(),
            comments: Vec::new(),
        }
    }
}

impl SamHeader {
    pub fn new() -> Self {
        Self::default()
    }

    /// `SAMFileHeader.setSortOrder`.
    pub fn set_sort_order(&mut self, order: &str) {
        self.attributes.set("SO", order);
    }

    /// `SAMFileHeader.setGroupOrder`.
    pub fn set_group_order(&mut self, order: &str) {
        self.attributes.set("GO", order);
    }

    /// `SAMFileHeader.addComment`, which prefixes `@CO\t`.
    pub fn add_comment(&mut self, comment: &str) {
        self.comments.push(format!("@CO\t{comment}"));
    }

    /// `SAMTextHeaderCodec.encode` with `keepExistingVersionNumber = true`, which is what
    /// `BAMFileWriter.writeHeader` passes.
    ///
    /// Line order is fixed by the Java: `@HD`, then every `@SQ`, then every `@RG`, then every
    /// `@PG`, then the comments. Within a line the fixed fields come first and the attributes
    /// follow in insertion order.
    pub fn encode(&self) -> String {
        let mut out = String::new();

        // @HD, even when it carries nothing but VN.
        let mut fields = vec![format!("{HEADER_LINE_START}HD")];
        push_attributes(&mut fields, &self.attributes);
        let _ = writeln!(out, "{}", fields.join("\t"));

        for seq in &self.sequences {
            let mut fields = vec![
                format!("{HEADER_LINE_START}SQ"),
                format!("SN:{}", seq.name),
                format!("LN:{}", seq.length),
            ];
            push_attributes(&mut fields, &seq.attributes);
            let _ = writeln!(out, "{}", fields.join("\t"));
        }

        for rg in &self.read_groups {
            let mut fields = vec![format!("{HEADER_LINE_START}RG"), format!("ID:{}", rg.id)];
            push_attributes(&mut fields, &rg.attributes);
            let _ = writeln!(out, "{}", fields.join("\t"));
        }

        for pg in &self.programs {
            let mut fields = vec![format!("{HEADER_LINE_START}PG"), format!("ID:{}", pg.id)];
            push_attributes(&mut fields, &pg.attributes);
            let _ = writeln!(out, "{}", fields.join("\t"));
        }

        for comment in &self.comments {
            let _ = writeln!(out, "{comment}");
        }

        out
    }
}

/// `SAMTextHeaderCodec.encodeTags`, via `TagValueCodec.encodeUntypedTag`.
fn push_attributes(fields: &mut Vec<String>, attributes: &Attributes) {
    for (k, v) in attributes.iter() {
        fields.push(format!("{k}:{v}"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The load-bearing property of the whole module.
    #[test]
    fn overwriting_an_attribute_keeps_its_original_position() {
        let mut a = Attributes::new();
        a.set("SM", "first");
        a.set("LB", "lib");
        a.set("SM", "second");
        a.set("PL", "ILLUMINA");
        let order: Vec<(&str, &str)> = a.iter().collect();
        assert_eq!(
            order,
            vec![("SM", "second"), ("LB", "lib"), ("PL", "ILLUMINA")],
            "LinkedHashMap.put on an existing key updates in place; it does not move the key"
        );
    }

    #[test]
    fn attributes_are_not_sorted() {
        let mut a = Attributes::new();
        a.set("ZZ", "1");
        a.set("AA", "2");
        let keys: Vec<&str> = a.iter().map(|(k, _)| k).collect();
        assert_eq!(keys, vec!["ZZ", "AA"], "insertion order, not sorted order");
    }

    #[test]
    fn setting_to_none_removes() {
        let mut a = Attributes::new();
        a.set("SM", "x");
        a.remove("SM");
        assert!(a.is_empty());
    }

    /// A fresh header carries VN, and VN comes first because the constructor set it first.
    #[test]
    fn a_new_header_leads_with_the_version() {
        assert_eq!(SamHeader::new().encode(), "@HD\tVN:1.6\n");
    }

    #[test]
    fn the_sort_order_follows_the_version() {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        assert_eq!(h.encode(), "@HD\tVN:1.6\tSO:coordinate\n");
    }

    #[test]
    fn sq_writes_sn_and_ln_before_its_attributes() {
        let mut h = SamHeader::new();
        let mut s = SequenceRecord::new("chr1", 250_000_000);
        s.attributes.set("AS", "GRCh38");
        s.attributes.set("M5", "d41d8cd98f00b204e9800998ecf8427e");
        h.sequences.push(s);
        assert_eq!(
            h.encode(),
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:250000000\tAS:GRCh38\tM5:d41d8cd98f00b204e9800998ecf8427e\n"
        );
    }

    /// The section order is fixed: HD, SQ, RG, PG, CO.
    #[test]
    fn sections_come_in_htsjdks_order() {
        let mut h = SamHeader::new();
        h.add_comment("a comment");
        h.programs.push(ProgramRecord::new("prog1"));
        h.read_groups.push(ReadGroup::new("rg1"));
        h.sequences.push(SequenceRecord::new("chr1", 100));
        let text = h.encode();
        let starts: Vec<&str> = text.lines().map(|l| &l[..3]).collect();
        assert_eq!(starts, vec!["@HD", "@SQ", "@RG", "@PG", "@CO"]);
    }

    #[test]
    fn every_line_ends_with_a_newline_including_the_last() {
        let mut h = SamHeader::new();
        h.add_comment("last");
        assert!(h.encode().ends_with("@CO\tlast\n"));
    }
}
