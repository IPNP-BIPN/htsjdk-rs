//! The `MetricsFile` layout.
//!
//! Ported from `htsjdk.samtools.metrics.MetricsFile.write`, `printHeaders`, `printBeanMetrics`
//! and `printHistogram`.
//!
//! Two properties are worth naming before the code.
//!
//! **The column order is unspecified by the Java language.** `printBeanMetrics` walks
//! `getBeanType().getFields()`, and `Class.getFields()` is documented as returning fields in
//! *no particular order*. On HotSpot it returns declaration order, and every Picard metrics
//! file ever written depends on that. A port must reproduce declaration order and can only
//! claim to match the reference implementation, not the language.
//!
//! **The trailing blank lines are structural.** `write` emits a newline after the headers,
//! after the bean rows, and after the histogram — including when the histogram is empty and
//! prints nothing at all. So a metrics file with no histogram ends with two consecutive blank
//! lines, and a writer that "tidied up" the trailing whitespace would produce a file every
//! parser still reads and no byte comparison accepts.

use crate::format::{format_bool, format_double, format_long};

/// `MetricsFile.MAJOR_HEADER_PREFIX`.
pub const MAJOR_HEADER_PREFIX: &str = "## ";
/// `MetricsFile.MINOR_HEADER_PREFIX`.
pub const MINOR_HEADER_PREFIX: &str = "# ";
/// `MetricsFile.SEPARATOR`.
pub const SEPARATOR: &str = "\t";
/// `MetricsFile.METRIC_HEADER`.
pub const METRIC_HEADER: &str = "## METRICS CLASS\t";
/// `MetricsFile.HISTO_HEADER`.
pub const HISTO_HEADER: &str = "## HISTOGRAM\t";

/// `htsjdk.samtools.metrics.StringHeader`, the only header type Picard tools use.
pub const STRING_HEADER_CLASS: &str = "htsjdk.samtools.metrics.StringHeader";

/// One cell of a metrics row, in the Java types `FormatUtil.format(Object)` dispatches on.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Long(i64),
    Double(f64),
    Bool(bool),
    /// A string, an enum name, or a date already rendered.
    Str(String),
    /// A null field, which `FormatUtil.format` renders as the empty string.
    Null,
}

impl Value {
    /// `FormatUtil.format(Object)`.
    pub fn format(&self) -> String {
        match self {
            Value::Long(v) => format_long(*v),
            Value::Double(v) => format_double(*v),
            Value::Bool(v) => format_bool(*v).to_string(),
            Value::Str(s) => s.clone(),
            // `if (value == null) return "";`
            Value::Null => String::new(),
        }
    }
}

/// A metric bean: its Java class name, its column names, and one row of values.
///
/// Columns and values must be in **declaration order**, which is what HotSpot's
/// `Class.getFields()` returns.
pub trait MetricBean {
    fn class_name(&self) -> &str;
    fn columns(&self) -> &[&'static str];
    fn values(&self) -> Vec<Value>;
}

/// A `Histogram` as `printHistogram` writes it: a bin label, a value label, and sorted bins.
#[derive(Debug, Clone)]
pub struct Histogram {
    pub bin_label: String,
    pub value_label: String,
    /// Java class name of the key type, written into the `## HISTOGRAM` line.
    pub key_class: String,
    /// Bins, which the writer keeps in the order given. htsjdk sorts them through a `TreeSet`
    /// on the histogram's comparator, so the caller supplies them already sorted.
    pub bins: Vec<(String, f64)>,
}

impl Histogram {
    pub fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }
}

/// `MetricsFile`.
#[derive(Debug, Default)]
pub struct MetricsFile {
    /// Header strings, each written as a class line followed by a value line.
    pub headers: Vec<String>,
    class_name: Option<String>,
    columns: Vec<&'static str>,
    rows: Vec<Vec<Value>>,
    pub histograms: Vec<Histogram>,
}

impl MetricsFile {
    pub fn new() -> Self {
        Self::default()
    }

    /// `MetricsFile.addHeader`, for the `StringHeader` case Picard uses.
    pub fn add_header(&mut self, text: &str) {
        self.headers.push(text.to_string());
    }

    /// `MetricsFile.addMetric`.
    pub fn add_metric<B: MetricBean>(&mut self, bean: &B) {
        self.class_name = Some(bean.class_name().to_string());
        self.columns = bean.columns().to_vec();
        self.rows.push(bean.values());
    }

    /// `MetricsFile.write(Writer)`.
    pub fn write(&self) -> String {
        let mut out = String::new();

        // printHeaders: two lines per header, the class name then the value.
        for h in &self.headers {
            out.push_str(MAJOR_HEADER_PREFIX);
            out.push_str(STRING_HEADER_CLASS);
            out.push('\n');
            out.push_str(MINOR_HEADER_PREFIX);
            out.push_str(h);
            out.push('\n');
        }
        out.push('\n');

        // printBeanMetrics: returns immediately when there are no metrics, so the blank line
        // above and the one below end up adjacent.
        if !self.rows.is_empty() {
            out.push_str(METRIC_HEADER);
            out.push_str(self.class_name.as_deref().unwrap_or(""));
            out.push('\n');
            out.push_str(&self.columns.join(SEPARATOR));
            out.push('\n');
            for row in &self.rows {
                let cells: Vec<String> = row.iter().map(|v| v.format()).collect();
                out.push_str(&cells.join(SEPARATOR));
                out.push('\n');
            }
        }
        out.push('\n');

        // printHistogram: empty histograms are dropped, and if none survive it writes nothing.
        let non_empty: Vec<&Histogram> = self.histograms.iter().filter(|h| !h.is_empty()).collect();
        if !non_empty.is_empty() {
            out.push_str(HISTO_HEADER);
            out.push_str(&non_empty[0].key_class);
            out.push('\n');
            out.push_str(&non_empty[0].bin_label);
            for h in &non_empty {
                out.push_str(SEPARATOR);
                out.push_str(&h.value_label);
            }
            out.push('\n');
            // The combined key set, in the order of the first histogram's bins followed by any
            // key only later histograms have. htsjdk builds a TreeSet, so callers supply
            // sorted bins and this preserves them.
            let mut keys: Vec<&str> = Vec::new();
            for h in &non_empty {
                for (k, _) in &h.bins {
                    if !keys.contains(&k.as_str()) {
                        keys.push(k);
                    }
                }
            }
            for key in keys {
                out.push_str(key);
                for h in &non_empty {
                    out.push_str(SEPARATOR);
                    // A key absent from this histogram contributes 0, not a blank.
                    let v = h
                        .bins
                        .iter()
                        .find(|(k, _)| k == key)
                        .map(|(_, v)| *v)
                        .unwrap_or(0.0);
                    out.push_str(&format_double(v));
                }
                out.push('\n');
            }
        }
        out.push('\n');

        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Simple;
    impl MetricBean for Simple {
        fn class_name(&self) -> &str {
            "test.Simple"
        }
        fn columns(&self) -> &[&'static str] {
            &["A", "B"]
        }
        fn values(&self) -> Vec<Value> {
            vec![Value::Long(1), Value::Double(0.5)]
        }
    }

    #[test]
    fn a_header_is_two_lines_the_class_then_the_value() {
        let mut f = MetricsFile::new();
        f.add_header("CollectSomething INPUT=x");
        let text = f.write();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "## htsjdk.samtools.metrics.StringHeader");
        assert_eq!(lines[1], "# CollectSomething INPUT=x");
    }

    /// The blank lines are structural: `write` emits one after the headers, one after the
    /// beans, and one after the histogram even when the histogram printed nothing.
    #[test]
    fn a_file_without_a_histogram_ends_with_two_blank_lines() {
        let mut f = MetricsFile::new();
        f.add_header("h");
        f.add_metric(&Simple);
        let text = f.write();
        assert!(
            text.ends_with("1\t0.5\n\n\n"),
            "trailing blank lines are part of the format, got {:?}",
            &text[text.len().saturating_sub(20)..]
        );
    }

    #[test]
    fn an_empty_file_is_still_three_newlines() {
        assert_eq!(MetricsFile::new().write(), "\n\n\n");
    }

    #[test]
    fn columns_come_out_in_the_order_given() {
        let mut f = MetricsFile::new();
        f.add_metric(&Simple);
        let text = f.write();
        assert!(text.contains("## METRICS CLASS\ttest.Simple\nA\tB\n1\t0.5\n"));
    }

    #[test]
    fn a_null_field_is_the_empty_string() {
        assert_eq!(Value::Null.format(), "");
        assert_eq!(Value::Bool(true).format(), "Y");
        assert_eq!(Value::Bool(false).format(), "N");
    }

    #[test]
    fn an_empty_histogram_prints_nothing_at_all() {
        let mut f = MetricsFile::new();
        f.add_metric(&Simple);
        f.histograms.push(Histogram {
            bin_label: "bin".into(),
            value_label: "count".into(),
            key_class: "java.lang.Integer".into(),
            bins: Vec::new(),
        });
        assert!(!f.write().contains("HISTOGRAM"));
    }

    #[test]
    fn a_histogram_writes_its_key_class_and_labels() {
        let mut f = MetricsFile::new();
        f.add_metric(&Simple);
        f.histograms.push(Histogram {
            bin_label: "quality".into(),
            value_label: "count".into(),
            key_class: "java.lang.Integer".into(),
            bins: vec![("1".into(), 10.0), ("2".into(), 20.5)],
        });
        let text = f.write();
        assert!(text.contains("## HISTOGRAM\tjava.lang.Integer\nquality\tcount\n1\t10\n2\t20.5\n"));
    }

    /// A key present in one histogram and absent from another contributes 0, not a blank.
    #[test]
    fn a_missing_bin_contributes_zero() {
        let mut f = MetricsFile::new();
        f.add_metric(&Simple);
        f.histograms.push(Histogram {
            bin_label: "k".into(),
            value_label: "a".into(),
            key_class: "java.lang.Integer".into(),
            bins: vec![("1".into(), 5.0)],
        });
        f.histograms.push(Histogram {
            bin_label: "k".into(),
            value_label: "b".into(),
            key_class: "java.lang.Integer".into(),
            bins: vec![("2".into(), 7.0)],
        });
        let text = f.write();
        assert!(text.contains("1\t5\t0\n2\t0\t7\n"), "got {text}");
    }
}
