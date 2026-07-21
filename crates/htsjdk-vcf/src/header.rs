//! The VCF header.
//!
//! Ported from `htsjdk.variant.vcf.VCFHeader`, `VCFHeaderLine`, `VCFCompoundHeaderLine` and
//! `htsjdk.variant.variantcontext.writer.VCFWriter.writeHeader`.
//!
//! ## The counterpart to decision 0009, with the opposite answer
//!
//! Decision 0009 recorded that a SAM header's attribute order is **insertion order**, because
//! `AbstractSAMHeaderRecord` holds them in a `LinkedHashMap`. The same library orders VCF
//! header lines the other way round:
//!
//! ```java
//! public Set<VCFHeaderLine> getMetaDataInSortedOrder() {
//!     return makeGetMetaDataSet(new TreeSet<VCFHeaderLine>(mMetaData));
//! }
//!
//! public int compareTo(Object other) {
//!     return toString().compareTo(other.toString());
//! }
//! ```
//!
//! A `TreeSet` over the **rendered string of the whole line**, payload included. So the order is
//! plain ASCII lexicographic over complete lines: `##FILTER` before `##FORMAT` before `##INFO`
//! before `##alsoUnstructured`, because `I` (0x49) sorts before `a` (0x61).
//!
//! Two header formats in one library, with opposite ordering rules, and neither states its rule
//! anywhere except in the choice of collection. A port that guessed "sorted" for SAM or
//! "insertion" for VCF would be wrong in both directions and produce valid files either way.
//!
//! ## And the VCF rule is not a total order
//!
//! `VCFContigHeaderLine` overrides the comparison:
//!
//! ```java
//! public int compareTo(final Object other) {
//!     if (other instanceof VCFContigHeaderLine)
//!         return contigIndex.compareTo(((VCFContigHeaderLine) other).contigIndex);
//!     else
//!         return super.compareTo(other);   // by rendered string
//! }
//! ```
//!
//! Contigs compare to each other **by index** and to everything else **by string**, which can
//! form a cycle. Measured, not reasoned: three lines whose pairwise comparisons give
//! `aaa < mmm < zzz < aaa` produce **two different headers** depending on the order they were
//! inserted. See decision 0016.
//!
//! Real headers avoid this because no line other than a contig renders a string starting with
//! `contig=`, so the cross comparisons never interleave. That is a property of the key
//! namespace, not of the comparator, and nothing enforces it.

/// `VCFHeader.METADATA_INDICATOR`.
pub const METADATA_INDICATOR: &str = "##";
/// `VCFHeader.HEADER_INDICATOR`.
pub const HEADER_INDICATOR: &str = "#";
/// `VCFConstants.FIELD_SEPARATOR`.
pub const FIELD_SEPARATOR: char = '\t';

/// The eight mandatory columns, `VCFHeader.HEADER_FIELDS`.
pub const HEADER_FIELDS: [&str; 8] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

/// `VCFHeaderLineCount`, the `Number` field of a compound line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cardinality {
    /// A fixed count, rendered as the integer.
    Fixed(i32),
    /// One value per alternate allele.
    A,
    /// One value per genotype.
    G,
    /// One value per allele, reference included.
    R,
    /// Unbounded, rendered as `.`.
    Unbounded,
}

impl Cardinality {
    fn render(self) -> String {
        match self {
            Cardinality::Fixed(n) => n.to_string(),
            Cardinality::A => "A".into(),
            Cardinality::G => "G".into(),
            Cardinality::R => "R".into(),
            Cardinality::Unbounded => ".".into(),
        }
    }
}

/// `VCFHeaderLineType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineType {
    Integer,
    Float,
    String,
    Character,
    Flag,
}

impl LineType {
    fn render(self) -> &'static str {
        match self {
            LineType::Integer => "Integer",
            LineType::Float => "Float",
            LineType::String => "String",
            LineType::Character => "Character",
            LineType::Flag => "Flag",
        }
    }
}

/// One header line.
///
/// `Contig` is a distinct variant rather than a `Structured` line with `key = "contig"`,
/// because htsjdk gives it a different comparison rule and the two are not interchangeable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HeaderLine {
    /// `##key=value`, with the value written raw.
    Unstructured { key: String, value: String },
    /// `##INFO=<...>` or `##FORMAT=<...>`.
    Compound {
        /// `INFO` or `FORMAT`.
        key: String,
        id: String,
        number: Cardinality,
        line_type: LineType,
        description: String,
        /// Extra key/value pairs appended after `Description`, in the order given.
        extra: Vec<(String, String)>,
    },
    /// `##FILTER=<ID=...,Description=...>`.
    Filter { id: String, description: String },
    /// `##contig=<...>`, which sorts by its **index** against other contigs.
    Contig {
        index: i32,
        fields: Vec<(String, String)>,
    },
    /// Any other structured line, rendered from its pairs in order.
    Structured {
        key: String,
        fields: Vec<(String, String)>,
    },
}

/// `VCFHeaderLine.escapeQuotes`.
fn escape_quotes(s: &str) -> String {
    s.replace('"', "\\\"")
}

/// `VCFHeaderLine.toStringEncoding(Map)`.
///
/// A value is quoted when it **contains** a comma or a space, **or** when its key is one of
/// `Description`, `Source` or `Version`. That is an OR of a content rule and a key rule, so
/// `Description=nospace` is quoted while `ID=nospace` is not, for the same value.
fn render_fields(fields: &[(String, String)]) -> String {
    let mut out = String::from("<");
    for (i, (k, v)) in fields.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(k);
        out.push('=');
        let quote = v.contains(',')
            || v.contains(' ')
            || k == "Description"
            || k == "Source"
            || k == "Version";
        if quote {
            out.push('"');
            out.push_str(&escape_quotes(v));
            out.push('"');
        } else {
            out.push_str(v);
        }
    }
    out.push('>');
    out
}

impl HeaderLine {
    /// The key, which is what the writer skips the `fileformat` line by.
    pub fn key(&self) -> &str {
        match self {
            HeaderLine::Unstructured { key, .. } => key,
            HeaderLine::Compound { key, .. } => key,
            HeaderLine::Filter { .. } => "FILTER",
            HeaderLine::Contig { .. } => "contig",
            HeaderLine::Structured { key, .. } => key,
        }
    }

    /// `VCFHeaderLine.toString`, which is also the sort key.
    pub fn render(&self) -> String {
        match self {
            HeaderLine::Unstructured { key, value } => format!("{key}={value}"),
            HeaderLine::Compound {
                key,
                id,
                number,
                line_type,
                description,
                extra,
            } => {
                let mut fields = vec![
                    ("ID".to_string(), id.clone()),
                    ("Number".to_string(), number.render()),
                    ("Type".to_string(), line_type.render().to_string()),
                    ("Description".to_string(), description.clone()),
                ];
                fields.extend(extra.iter().cloned());
                format!("{key}={}", render_fields(&fields))
            }
            HeaderLine::Filter { id, description } => {
                let fields = vec![
                    ("ID".to_string(), id.clone()),
                    ("Description".to_string(), description.clone()),
                ];
                format!("FILTER={}", render_fields(&fields))
            }
            HeaderLine::Contig { fields, .. } => format!("contig={}", render_fields(fields)),
            HeaderLine::Structured { key, fields } => {
                format!("{key}={}", render_fields(fields))
            }
        }
    }

    /// An `##INFO` line.
    pub fn info(id: &str, number: Cardinality, line_type: LineType, description: &str) -> Self {
        HeaderLine::Compound {
            key: "INFO".into(),
            id: id.into(),
            number,
            line_type,
            description: description.into(),
            extra: Vec::new(),
        }
    }

    /// A `##FORMAT` line.
    pub fn format(id: &str, number: Cardinality, line_type: LineType, description: &str) -> Self {
        HeaderLine::Compound {
            key: "FORMAT".into(),
            id: id.into(),
            number,
            line_type,
            description: description.into(),
            extra: Vec::new(),
        }
    }

    /// A `##FILTER` line.
    pub fn filter(id: &str, description: &str) -> Self {
        HeaderLine::Filter {
            id: id.into(),
            description: description.into(),
        }
    }

    /// A `##contig` line at a dictionary index.
    ///
    /// The index is what contigs sort by, so it is required rather than optional: a contig
    /// without one cannot be placed.
    pub fn contig(id: &str, length: i64, index: i32) -> Self {
        HeaderLine::Contig {
            index,
            fields: vec![
                ("ID".to_string(), id.to_string()),
                ("length".to_string(), length.to_string()),
            ],
        }
    }
}

/// `VCFHeader`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VcfHeader {
    pub lines: Vec<HeaderLine>,
    pub samples: Vec<String>,
}

/// The version line `VCFWriter` writes first, before anything else.
pub const VERSION_LINE: &str = "##fileformat=VCFv4.2";

impl VcfHeader {
    pub fn new() -> Self {
        Self::default()
    }

    /// `VCFWriter.writeHeader`.
    ///
    /// The version line is written first and unconditionally. `getMetaDataInSortedOrder` also
    /// *prepends* a `fileformat` line to the sorted set, which the writer then skips with
    /// `if (VCFHeaderVersion.isFormatString(line.getKey())) continue`. So the line is added and
    /// filtered out again; reproducing the skip is what keeps a user-supplied `fileformat` line
    /// from appearing twice.
    pub fn write(&self) -> String {
        let mut out = String::new();
        out.push_str(VERSION_LINE);
        out.push('\n');

        // `VCFHeaderLine.compareTo` is by rendered string, except that contigs compare to each
        // other by index. Reproduced as a sort key: contigs carry their index, everything else
        // carries only its string, so contigs order among themselves by index and against
        // others by string. Where those two disagree the Java comparator has a cycle and its
        // result depends on TreeSet insertion order; see decision 0016 and the module doc.
        let mut sorted: Vec<(SortKey, String)> = self
            .lines
            .iter()
            .filter(|l| !is_format_string(l.key()))
            .map(|l| (l.sort_key(), l.render()))
            .collect();
        sorted.sort_by(|a, b| a.0.cmp(&b.0));
        sorted.dedup_by(|a, b| a.1 == b.1);
        for (_, line) in sorted {
            out.push_str(METADATA_INDICATOR);
            out.push_str(&line);
            out.push('\n');
        }

        out.push_str(HEADER_INDICATOR);
        out.push_str(&HEADER_FIELDS.join(&FIELD_SEPARATOR.to_string()));
        if !self.samples.is_empty() {
            out.push(FIELD_SEPARATOR);
            out.push_str("FORMAT");
            for s in &self.samples {
                out.push(FIELD_SEPARATOR);
                out.push_str(s);
            }
        }
        out.push('\n');
        out
    }
}

/// The key a line sorts by.
///
/// Contigs sort among themselves by index and against everything else by their rendered
/// string, which is what `VCFContigHeaderLine.compareTo` does. Encoding that as a key rather
/// than as a comparator makes the resulting order a genuine total order, which the Java one is
/// not. The difference is only observable on headers that trigger the cycle, and those are
/// exactly the headers on which htsjdk itself is order-dependent.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct SortKey {
    string: String,
    /// `Some` for contigs, so two contigs compare by index once their strings have been
    /// compared. Only reached when two contigs render identically, which the dedup handles.
    index: Option<i32>,
}

impl HeaderLine {
    fn sort_key(&self) -> SortKey {
        match self {
            HeaderLine::Contig { index, .. } => SortKey {
                // Contigs share a prefix so that they group together, then order by index.
                string: format!("contig={:012}", index),
                index: Some(*index),
            },
            other => SortKey {
                string: other.render(),
                index: None,
            },
        }
    }
}

/// `VCFHeaderVersion.isFormatString`.
fn is_format_string(key: &str) -> bool {
    key == "fileformat" || key == "format"
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The ordering rule, stated on its own because it is the opposite of the SAM one.
    #[test]
    fn lines_are_sorted_by_their_whole_rendered_string() {
        let mut h = VcfHeader::new();
        h.lines = vec![
            HeaderLine::info("ZZ", Cardinality::Fixed(1), LineType::Integer, "last by id"),
            HeaderLine::filter("zFilter", "z"),
            HeaderLine::info("AA", Cardinality::Fixed(1), LineType::String, "first by id"),
            HeaderLine::format("ZQ", Cardinality::Fixed(1), LineType::Float, "z format"),
            HeaderLine::filter("aFilter", "a"),
            HeaderLine::format("AQ", Cardinality::Fixed(1), LineType::Float, "a format"),
        ];
        let text = h.write();
        let keys: Vec<&str> = text
            .lines()
            .filter(|l| l.starts_with("##") && !l.contains("fileformat"))
            .map(|l| l.split('=').next().unwrap())
            .collect();
        assert_eq!(
            keys,
            vec!["##FILTER", "##FILTER", "##FORMAT", "##FORMAT", "##INFO", "##INFO"]
        );
    }

    /// Plain ASCII, so an upper-case key sorts before a lower-case one whatever it means.
    #[test]
    fn upper_case_sorts_before_lower_case() {
        let mut h = VcfHeader::new();
        h.lines = vec![
            HeaderLine::Unstructured {
                key: "alsoUnstructured".into(),
                value: "v".into(),
            },
            HeaderLine::info("DP", Cardinality::Fixed(1), LineType::Integer, "d"),
        ];
        let text = h.write();
        let lines: Vec<String> = text
            .lines()
            .filter(|l| l.starts_with("##") && !l.contains("fileformat"))
            .map(str::to_string)
            .collect();
        assert!(lines[0].starts_with("##INFO"), "got {lines:?}");
        assert!(lines[1].starts_with("##alsoUnstructured"));
    }

    /// The quoting rule is an OR of a content test and a key test, so the same value is quoted
    /// under `Description` and bare under `ID`.
    #[test]
    fn quoting_depends_on_the_key_as_well_as_the_value() {
        let line = HeaderLine::info(
            "nospace",
            Cardinality::Fixed(1),
            LineType::String,
            "nospace",
        );
        let rendered = line.render();
        assert!(rendered.contains("ID=nospace,"), "ID is bare: {rendered}");
        assert!(
            rendered.contains(r#"Description="nospace""#),
            "Description is always quoted: {rendered}"
        );
    }

    #[test]
    fn a_value_with_a_comma_or_a_space_is_quoted_whatever_its_key() {
        let line = HeaderLine::Structured {
            key: "contig".into(),
            fields: vec![
                ("ID".into(), "chr 1".into()),
                ("assembly".into(), "a,b".into()),
                ("length".into(), "100".into()),
            ],
        };
        let r = line.render();
        assert!(r.contains(r#"ID="chr 1""#), "{r}");
        assert!(r.contains(r#"assembly="a,b""#), "{r}");
        assert!(
            r.contains("length=100"),
            "an ordinary value stays bare: {r}"
        );
    }

    #[test]
    fn a_quote_inside_a_value_is_backslash_escaped() {
        let line = HeaderLine::info("Q", Cardinality::Fixed(1), LineType::String, "with\"quote");
        assert!(line.render().contains(r#"Description="with\"quote""#));
    }

    #[test]
    fn cardinalities_render_as_their_letters() {
        for (c, expected) in [
            (Cardinality::Fixed(3), "Number=3"),
            (Cardinality::A, "Number=A"),
            (Cardinality::G, "Number=G"),
            (Cardinality::R, "Number=R"),
            (Cardinality::Unbounded, "Number=."),
        ] {
            let line = HeaderLine::info("X", c, LineType::Float, "d");
            assert!(line.render().contains(expected), "{}", line.render());
        }
    }

    #[test]
    fn a_minimal_header_is_the_version_line_and_the_column_line() {
        assert_eq!(
            VcfHeader::new().write(),
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        );
    }

    #[test]
    fn samples_add_a_format_column_and_their_names() {
        let mut h = VcfHeader::new();
        h.samples = vec!["NA12878".into(), "NA12891".into()];
        let last = h.write().lines().last().unwrap().to_string();
        assert!(last.ends_with("\tINFO\tFORMAT\tNA12878\tNA12891"), "{last}");
    }

    /// A user-supplied `fileformat` line is skipped, so it cannot appear twice.
    #[test]
    fn a_supplied_fileformat_line_is_not_written_again() {
        let mut h = VcfHeader::new();
        h.lines = vec![HeaderLine::Unstructured {
            key: "fileformat".into(),
            value: "VCFv4.1".into(),
        }];
        let text = h.write();
        assert_eq!(text.matches("fileformat").count(), 1);
        assert!(text.starts_with("##fileformat=VCFv4.2\n"));
    }

    /// Sorting by the whole line means a prefix ID sorts before a longer one, because `,` is
    /// 0x2C and any letter is above it.
    #[test]
    fn a_prefix_id_sorts_before_a_longer_one() {
        let mut h = VcfHeader::new();
        h.lines = vec![
            HeaderLine::info("AB", Cardinality::Fixed(1), LineType::Integer, "d"),
            HeaderLine::info("A", Cardinality::Fixed(1), LineType::Integer, "d"),
        ];
        let text = h.write();
        let ids: Vec<String> = text
            .lines()
            .filter(|l| l.starts_with("##INFO"))
            .map(|l| l.split(&['=', ','][..]).nth(2).unwrap().to_string())
            .collect();
        assert_eq!(ids, vec!["A", "AB"]);
    }
}
