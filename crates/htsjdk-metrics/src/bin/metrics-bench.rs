//! Writing a metrics file, on the beans and histograms `tools/benchmark/MetricsBench.java` writes
//! with htsjdk's `MetricsFile`. The output protocol is `bam-ops-bench`'s.
//!
//! The rows are many more than a Picard tool writes, so the formatting is what the time measures;
//! the two histograms are the size of an insert-size histogram.
//!
//! Usage: metrics-bench dir reps

/// The allocator picard-rs links, so the rows measure the library as it ships there.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static ALLOCATOR: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::path::PathBuf;
use std::time::Instant;

use htsjdk_metrics::file::Histogram;
use htsjdk_metrics::{MetricBean, MetricsFile, Value};

const ROWS: usize = 200_000;
const BINS: usize = 1_000;

struct Bench {
    i: usize,
}

impl MetricBean for Bench {
    fn class_name(&self) -> &str {
        "MetricsBench$BenchMetric"
    }

    fn columns(&self) -> &[&'static str] {
        &["READS", "MEAN", "NAME", "PASS"]
    }

    fn values(&self) -> Vec<Value> {
        let i = self.i;
        vec![
            Value::Long(i as i64 * 37),
            Value::Double(i as f64 / 7.0),
            Value::Str(format!("sample_{i}")),
            Value::Bool(i.is_multiple_of(3)),
        ]
    }
}

fn histogram(value_label: &str, scale: f64) -> Histogram {
    Histogram {
        bin_label: "insert_size".to_string(),
        value_label: value_label.to_string(),
        key_class: "java.lang.Integer".to_string(),
        bins: (0..BINS)
            .map(|k| (k.to_string(), (k as f64 * scale).floor()))
            .collect(),
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let dir = PathBuf::from(args.next().expect("dir"));
    let reps: usize = args.next().map(|a| a.parse().expect("reps")).unwrap_or(3);
    let out = dir.join("port_metrics_write.txt");
    let write = || {
        let mut file = MetricsFile::new();
        file.add_header("benchmark");
        for i in 0..ROWS {
            file.add_metric(&Bench { i });
        }
        file.histograms.push(histogram("count", 1.5));
        file.histograms.push(histogram("count_dups", 0.25));
        let text = file.write();
        std::fs::write(&out, &text).expect("writes");
        text.len()
    };
    let size = write();
    println!("op name=metrics_write bytes={size} records={ROWS} check=bytes files=.txt");
    for run in 0..reps {
        let start = Instant::now();
        write();
        println!(
            "rust_metrics_write_run{run}_secs={:.6}",
            start.elapsed().as_secs_f64()
        );
    }
}
