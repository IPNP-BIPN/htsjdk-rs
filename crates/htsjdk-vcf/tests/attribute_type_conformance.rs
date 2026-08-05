//! Conformance for the typed attribute accessors, against `CommonInfo` and
//! `VCFUtils.parseVcfDouble`.
//!
//! Goldens from `tools/vcf-conformance/VcfAttributeTypeDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the three spellings of missing, which are identical in every
//! rendering and are not identical to the accessors:
//!
//! ```text
//! as  missing-dot    I1  asString  .        as  missing-dot    I1  asInt  E:...NumberFormatException
//! as  missing-empty  I1  asString  .        as  missing-empty  I1  asInt  -1
//! as  missing-bare   I1  asString  .        as  missing-bare   I1  asInt  -1
//! ```
//!
//! And the two accessors over the same conversion that disagree about what missing is:
//!
//! ```text
//! as  missing-empty  I1  asInt         -1
//! as  missing-empty  I1  asDouble      E:java.lang.NumberFormatException:For input string: "."
//! as  missing-empty  I1  asDoubleList  [-1.0]
//! ```

use std::io::Read;

use htsjdk_vcf::attributes::{
    as_boolean, as_double, as_double_list, as_int, as_int_list, as_list, as_string, as_string_list,
    parse_vcf_double, AttributeError,
};
use htsjdk_vcf::header_lines::parse_meta_line;
use htsjdk_vcf::header_parse::read_header_frame;
use htsjdk_vcf::record_parse::decode_line;
use htsjdk_vcf::variant::{Value, VariantContext};
use htsjdk_vcf::VcfHeader;

/// `VcfAttributeTypeDump.HEADER`, one declaration per type the format has.
const HEADER: &str = "##fileformat=VCFv4.2\n\
    ##INFO=<ID=I1,Number=1,Type=Integer,Description=\"one integer\">\n\
    ##INFO=<ID=IA,Number=A,Type=Integer,Description=\"per alt integer\">\n\
    ##INFO=<ID=F1,Number=1,Type=Float,Description=\"one float\">\n\
    ##INFO=<ID=FR,Number=R,Type=Float,Description=\"per allele float\">\n\
    ##INFO=<ID=S1,Number=1,Type=String,Description=\"one string\">\n\
    ##INFO=<ID=C1,Number=1,Type=Character,Description=\"one character\">\n\
    ##INFO=<ID=CU,Number=.,Type=Character,Description=\"characters\">\n\
    ##INFO=<ID=B1,Number=0,Type=Flag,Description=\"a flag\">\n\
    ##contig=<ID=chr1,length=100000>\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Label and INFO field, in the dump's order.
const CASES: &[(&str, &str)] = &[
    ("integer", "I1=42"),
    ("integer-list", "IA=1,2,3"),
    ("float", "F1=0.5"),
    ("float-list", "FR=0.5,0.25"),
    ("string", "S1=hello"),
    ("character", "C1=x"),
    ("character-list", "CU=a,b,c"),
    ("flag", "B1"),
    ("missing-dot", "I1=."),
    ("missing-empty", "I1="),
    ("missing-bare", "I1"),
    ("missing-dot-float", "F1=."),
    ("missing-in-a-list", "IA=1,.,3"),
    ("integer-holding-a-float", "I1=1.5"),
    ("integer-holding-text", "I1=abc"),
    ("float-holding-an-integer", "F1=7"),
    ("float-holding-text", "F1=abc"),
    ("character-holding-a-word", "C1=word"),
    ("string-holding-true", "S1=true"),
    ("string-holding-TRUE", "S1=TRUE"),
    ("string-holding-one", "S1=1"),
    ("float-inf", "F1=inf"),
    ("float-minus-inf", "F1=-inf"),
    ("float-nan", "F1=nan"),
    ("absent", "S1=present"),
];

/// The `dbl` cases, in the dump's order.
const DOUBLES: &[&str] = &[
    "1",
    "1.5",
    "-1.5",
    "1e3",
    "1E3",
    "+1",
    " 1",
    "1 ",
    "1f",
    "1d",
    "0x1p3",
    "Infinity",
    "-Infinity",
    "+Infinity",
    "inf",
    "-inf",
    "+inf",
    "INF",
    "Inf",
    "infinity",
    "nan",
    "NaN",
    "NAN",
    "-nan",
    "+nan",
    ".",
    "",
    "1,5",
];

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/vcf_attribute_type.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn header() -> (VcfHeader, usize) {
    let frame = read_header_frame(HEADER).expect("the fixture header parses");
    let mut header = VcfHeader::new();
    header.samples = frame.samples.clone();
    let mut contigs = 0;
    for line in &frame.meta_lines {
        if let Ok(parsed) = parse_meta_line(line, frame.version, contigs) {
            if matches!(parsed, htsjdk_vcf::HeaderLine::Contig { .. }) {
                contigs += 1;
            }
            header.lines.push(parsed);
        }
    }
    (header, frame.meta_lines.len() + 1)
}

fn record(info: &str) -> VariantContext {
    let (header, line_no) = header();
    decode_line(
        &format!("chr1\t100\t.\tA\tT\t50\tPASS\t{info}"),
        &header,
        line_no,
        htsjdk_vcf::header_parse::VcfVersion::Vcf4_2,
    )
    .expect("the fixture line decodes")
    .expect("a data line")
    .variant
}

/// The dump's escaping, restricted to what these values contain.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            c if (c as u32) < 0x20 || (c as u32) > 0x7e => {
                out.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out
}

/// `String.valueOf` over what an accessor returned, or the exception's class and message.
fn outcome<T>(result: Result<T, AttributeError>, render: impl Fn(T) -> String) -> String {
    match result {
        Ok(value) => escape(&render(value)),
        Err(error) => format!("E:{}:{}", error.class(), escape(&error.message())),
    }
}

/// `Double.toString`, for the shapes these accessors return.
fn java_double(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else if value.is_infinite() {
        if value > 0.0 {
            "Infinity".to_string()
        } else {
            "-Infinity".to_string()
        }
    } else if value == value.trunc() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

fn render_list<T>(values: Vec<T>, render: impl Fn(&T) -> String) -> String {
    format!(
        "[{}]",
        values.iter().map(render).collect::<Vec<_>>().join(",")
    )
}

/// `String.valueOf(Object)` over a stored value, which is how the `attr` row renders it.
fn stored_text(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        Value::Missing => ".".to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Int(n) => n.to_string(),
        Value::Double(d) => htsjdk_vcf::variant::format_vcf_double(*d),
        Value::List(items) => format!(
            "[{}]",
            items.iter().map(stored_text).collect::<Vec<_>>().join(", ")
        ),
    }
}

/// The Java class a stored value reports, which is the row that says the declared Type converted
/// nothing.
fn stored_class(value: Option<&Value>) -> &'static str {
    match value {
        None => "null",
        Some(Value::List(_)) => "ArrayList",
        Some(Value::Bool(_)) => "Boolean",
        _ => "String",
    }
}

#[test]
fn every_accessor_answers_as_the_reference_answers() {
    let corpus = corpus();
    let golden: Vec<&str> = corpus
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    let mut produced: Vec<String> = Vec::new();

    for (label, info) in CASES {
        let variant = record(info);
        let key = if *label == "absent" {
            "MISSING"
        } else {
            info.split(['=', ';']).next().expect("a key")
        };
        let value = variant
            .attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value);

        produced.push(format!(
            "attr\t{label}\t{key}\t{}\t{}",
            stored_class(value),
            escape(&value.map_or_else(|| "null".to_string(), stored_text))
        ));

        let mut row = |accessor: &str, text: String| {
            produced.push(format!("as\t{label}\t{key}\t{accessor}\t{text}"));
        };
        row("asString", escape(&as_string(value, "DEFAULT")));
        row("asInt", outcome(as_int(value, -1), |n| n.to_string()));
        row("asDouble", outcome(as_double(value, -1.0), java_double));
        row(
            "asBoolean",
            outcome(as_boolean(value, false), |b| b.to_string()),
        );
        row("asList", escape(&render_list(as_list(value), stored_text)));
        row(
            "asStringList",
            escape(&render_list(as_string_list(value, "D"), |s| s.clone())),
        );
        row(
            "asIntList",
            outcome(as_int_list(value, -1), |list| {
                render_list(list, |n| n.to_string())
            }),
        );
        row(
            "asDoubleList",
            outcome(as_double_list(value, -1.0), |list| {
                render_list(list, |d| java_double(*d))
            }),
        );
    }

    for raw in DOUBLES {
        produced.push(format!(
            "dbl\t{}\t{}",
            escape(raw),
            outcome(parse_vcf_double(raw), java_double)
        ));
    }

    assert_eq!(
        produced.len(),
        golden.len(),
        "the port produced {} rows and the golden has {}",
        produced.len(),
        golden.len()
    );
    for (index, (mine, theirs)) in produced.iter().zip(golden.iter()).enumerate() {
        assert_eq!(mine, theirs, "row {index}");
    }
}

/// Stated on its own, because it is the reason [`Value::Missing`] exists as a variant distinct from
/// `Value::Str(".")` and a port that collapsed them would pass every other suite.
#[test]
fn the_three_spellings_of_missing_are_two_behaviours() {
    let dot = record("I1=.");
    let empty = record("I1=");
    let bare = record("I1");

    let value = |vc: &VariantContext| vc.attributes.first().map(|(_, value)| value.clone());
    assert!(as_int(value(&dot).as_ref(), -1).is_err());
    assert_eq!(as_int(value(&empty).as_ref(), -1), Ok(-1));
    assert_eq!(as_int(value(&bare).as_ref(), -1), Ok(-1));

    // And all three render identically, which is why nothing downstream can tell them apart.
    for vc in [&dot, &empty, &bare] {
        assert_eq!(as_string(value(vc).as_ref(), "D"), ".");
    }
}
