//! The header lines htsjdk rewrites on the way in.
//!
//! Ported from `htsjdk.variant.vcf.VCFStandardHeaderLines` at htsjdk 4.2.0.
//!
//! `AbstractVCFCodec.doOnTheFlyModifications` defaults to **true**, so every VCF read through the
//! codec goes through `repairStandardHeaderLines`, and the header a reader hands back is not the
//! header the file declared. For eighteen IDs htsjdk holds its own definition, and when the file's
//! count or type disagrees the whole line is **replaced**, description included.
//!
//! Three details decide what survives:
//!
//!  * **the description alone never triggers a repair.** `REPAIR_BAD_DESCRIPTIONS` is a private
//!    `false`, so a file's own wording for a standard ID is kept, and only a count or type
//!    disagreement loses it. That makes the rewrite hard to notice: the line that comes back
//!    usually looks like the one that went in;
//!  * **when a repair does happen the description goes too**, because the replacement is the whole
//!    standard line rather than a patched copy of the file's;
//!  * **the count comparison is in two parts.** A different *kind* of count (`A` against a fixed
//!    number) is a disagreement on its own; a different fixed *value* is only compared when both
//!    are fixed. So `Number=A` against a standard `Number=1` is caught by the first test and never
//!    reaches the second.
//!
//! `DP` appears in both tables with **different descriptions**, so the repaired line depends on
//! whether it was an INFO or a FORMAT line.

use crate::header::{Cardinality, HeaderLine, LineType};

/// One entry of `formatStandards` or `infoStandards`.
struct Standard {
    id: &'static str,
    number: Cardinality,
    line_type: LineType,
    description: &'static str,
}

/// `VCFStandardHeaderLines`' static initialiser, FORMAT half.
const FORMAT_STANDARDS: &[Standard] = &[
    Standard {
        id: "GT",
        number: Cardinality::Fixed(1),
        line_type: LineType::String,
        description: "Genotype",
    },
    Standard {
        id: "GQ",
        number: Cardinality::Fixed(1),
        line_type: LineType::Integer,
        description: "Genotype Quality",
    },
    Standard {
        id: "DP",
        number: Cardinality::Fixed(1),
        line_type: LineType::Integer,
        description: "Approximate read depth (reads with MQ=255 or with bad mates are filtered)",
    },
    Standard {
        id: "PL",
        number: Cardinality::G,
        line_type: LineType::Integer,
        description: "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF \
                      specification",
    },
    Standard {
        id: "AD",
        number: Cardinality::R,
        line_type: LineType::Integer,
        description: "Allelic depths for the ref and alt alleles in the order listed",
    },
    Standard {
        id: "FT",
        number: Cardinality::Unbounded,
        line_type: LineType::String,
        description: "Genotype-level filter",
    },
    Standard {
        id: "PS",
        number: Cardinality::Fixed(1),
        line_type: LineType::Integer,
        description: "Phasing set (typically the position of the first variant in the set)",
    },
    Standard {
        id: "PQ",
        number: Cardinality::Fixed(1),
        line_type: LineType::Float,
        description: "Read-backed phasing quality",
    },
];

/// `VCFStandardHeaderLines`' static initialiser, INFO half.
const INFO_STANDARDS: &[Standard] = &[
    Standard {
        id: "END",
        number: Cardinality::Fixed(1),
        line_type: LineType::Integer,
        description: "Stop position of the interval",
    },
    Standard {
        id: "DB",
        number: Cardinality::Fixed(0),
        line_type: LineType::Flag,
        description: "dbSNP Membership",
    },
    Standard {
        id: "DP",
        number: Cardinality::Fixed(1),
        line_type: LineType::Integer,
        description: "Approximate read depth; some reads may have been filtered",
    },
    Standard {
        id: "SB",
        number: Cardinality::Fixed(1),
        line_type: LineType::Float,
        description: "Strand Bias",
    },
    Standard {
        id: "AF",
        number: Cardinality::A,
        line_type: LineType::Float,
        description: "Allele Frequency, for each ALT allele, in the same order as listed",
    },
    Standard {
        id: "AC",
        number: Cardinality::A,
        line_type: LineType::Integer,
        description: "Allele count in genotypes, for each ALT allele, in the same order as listed",
    },
    Standard {
        id: "AN",
        number: Cardinality::Fixed(1),
        line_type: LineType::Integer,
        description: "Total number of alleles in called genotypes",
    },
    Standard {
        id: "MQ0",
        number: Cardinality::Fixed(1),
        line_type: LineType::Integer,
        description: "Total Mapping Quality Zero Reads",
    },
    Standard {
        id: "MQ",
        number: Cardinality::Fixed(1),
        line_type: LineType::Float,
        description: "RMS Mapping Quality",
    },
    Standard {
        id: "SOMATIC",
        number: Cardinality::Fixed(0),
        line_type: LineType::Flag,
        description: "Somatic event",
    },
];

/// `Standards.repair`, for one line. Anything that is not a compound line is returned untouched.
pub fn repair(line: &HeaderLine) -> HeaderLine {
    let HeaderLine::Compound {
        key,
        id,
        number,
        line_type,
        ..
    } = line
    else {
        return line.clone();
    };

    let table = match key.as_str() {
        "INFO" => INFO_STANDARDS,
        "FORMAT" => FORMAT_STANDARDS,
        // `repairStandardHeaderLines` dispatches on the Java type of the line, and only
        // `VCFInfoHeaderLine` and `VCFFormatHeaderLine` have standards.
        _ => return line.clone(),
    };
    let Some(standard) = table.iter().find(|s| s.id == *id) else {
        return line.clone();
    };

    // The count is compared in two steps, and the second is only reached when both are fixed.
    let bad_count_type = count_kind(*number) != count_kind(standard.number);
    let bad_count = match (number, standard.number) {
        (Cardinality::Fixed(a), Cardinality::Fixed(b)) => !bad_count_type && *a != b,
        _ => false,
    };
    let bad_type = *line_type != standard.line_type;

    // `REPAIR_BAD_DESCRIPTIONS` is false, so a description that disagrees on its own is kept.
    if bad_count_type || bad_count || bad_type {
        HeaderLine::Compound {
            key: key.clone(),
            id: standard.id.to_string(),
            number: standard.number,
            line_type: standard.line_type,
            description: standard.description.to_string(),
            // The replacement is the standard line itself, so any extra tags the file carried on
            // that line are dropped with it.
            extra: Vec::new(),
        }
    } else {
        line.clone()
    }
}

/// `VCFCompoundHeaderLine.getCountType`: which *kind* of count, ignoring a fixed one's value.
fn count_kind(number: Cardinality) -> u8 {
    match number {
        Cardinality::Fixed(_) => 0,
        Cardinality::A => 1,
        Cardinality::G => 2,
        Cardinality::R => 3,
        Cardinality::Unbounded => 4,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn info(number: Cardinality, line_type: LineType, description: &str) -> HeaderLine {
        HeaderLine::info("DP", number, line_type, description)
    }

    #[test]
    fn a_wrong_type_replaces_the_whole_line_including_the_description() {
        let repaired = repair(&info(Cardinality::Fixed(1), LineType::Float, "my depth"));
        let HeaderLine::Compound {
            line_type,
            description,
            ..
        } = &repaired
        else {
            panic!("still a compound line");
        };
        assert_eq!(*line_type, LineType::Integer);
        assert_eq!(
            description,
            "Approximate read depth; some reads may have been filtered"
        );
    }

    #[test]
    fn a_wrong_description_alone_is_kept() {
        let line = info(Cardinality::Fixed(1), LineType::Integer, "my depth");
        assert_eq!(repair(&line), line);
    }

    /// `Number=A` against a fixed standard is caught by the count *kind*, before any value is
    /// compared.
    #[test]
    fn a_different_count_kind_is_a_repair() {
        let repaired = repair(&info(Cardinality::A, LineType::Integer, "d"));
        let HeaderLine::Compound { number, .. } = &repaired else {
            panic!("still a compound line");
        };
        assert_eq!(*number, Cardinality::Fixed(1));
    }

    /// The same ID under the other key gets the other table's description.
    #[test]
    fn dp_is_in_both_tables_with_different_descriptions() {
        let repaired = repair(&HeaderLine::format(
            "DP",
            Cardinality::Fixed(2),
            LineType::Integer,
            "d",
        ));
        let HeaderLine::Compound { description, .. } = &repaired else {
            panic!("still a compound line");
        };
        assert_eq!(
            description,
            "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"
        );
    }

    #[test]
    fn an_id_with_no_standard_is_untouched() {
        let line = HeaderLine::info("XX", Cardinality::Fixed(2), LineType::Float, "Mine");
        assert_eq!(repair(&line), line);
    }
}
