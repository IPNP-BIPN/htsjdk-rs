//! Port of htsjdk's metrics output: number formatting and the `MetricsFile` layout.
//!
//! Ported from htsjdk 4.2.0 `htsjdk.samtools.util.FormatUtil` and
//! `htsjdk.samtools.metrics.MetricsFile`.

pub mod file;
pub mod format;

pub use file::{MetricBean, MetricsFile, Value};
pub use format::{format_bool, format_double, format_long};
