//! Canonicalization: the fields a comparison is allowed to ignore, and the record of which
//! ones it actually ignored.
//!
//! Some fields legitimately vary between two correct runs: timestamps, version strings, the
//! command line recorded in a header. Comparing them raw would fail every time and prove
//! nothing.
//!
//! The hazard is the opposite one. Each rule added here weakens the claim by exactly the
//! ground it covers, and a rule added to make a stubborn test pass is indistinguishable from a
//! rule added for a good reason. So:
//!
//! - Rules are **named and declared**, never applied inline at a comparison site.
//! - Every rule carries the justification for its own existence, in `why`.
//! - A comparison reports which rules fired, so "identical" is never reported without saying
//!   what was excused to get there.

/// One field a comparison may normalise away.
#[derive(Debug, Clone)]
pub struct Rule {
    /// Short identifier reported in results.
    pub name: &'static str,
    /// Why this field is allowed to vary between two correct runs. Not decoration: a rule
    /// whose justification cannot be written is a rule that should not exist.
    pub why: &'static str,
    /// Lines matching this prefix (after leading whitespace) are replaced wholesale.
    pub line_prefix: &'static str,
}

/// Rules for Picard metrics files.
pub const PICARD_METRICS: &[Rule] = &[
    Rule {
        name: "picard.started-on",
        why: "wall-clock time of the run; changes on every invocation by construction",
        line_prefix: "# Started on:",
    },
    Rule {
        name: "picard.command-line",
        why: "records absolute input paths and the tool version, neither of which is part of \
               the computed result",
        line_prefix: "# ",
    },
];

/// Rules for GATK report files.
pub const GATK_REPORT: &[Rule] = &[Rule {
    name: "gatk.argument-table-command-line",
    why: "the argument table echoes the invocation, including paths that differ between the \
           reference run and ours",
    line_prefix: "#:GATKReport",
}];

/// Rules for VCF headers.
pub const VCF_HEADER: &[Rule] = &[
    Rule {
        name: "vcf.file-date",
        why: "date of the run, not of the data",
        line_prefix: "##fileDate",
    },
    Rule {
        name: "vcf.source",
        why: "names the producing tool and its version; the port is deliberately a different \
               tool with a different version string",
        line_prefix: "##source",
    },
];

/// Applies rules to text, returning the normalised bytes and the names of the rules that
/// actually fired.
///
/// A rule that is declared but never matches does **not** appear in the returned list. That
/// distinction matters: it means a result saying "identical, no rules applied" is a genuine
/// byte-identity claim even though rules were available.
pub fn canonicalize(input: &[u8], rules: &[Rule]) -> (Vec<u8>, Vec<String>) {
    if rules.is_empty() {
        return (input.to_vec(), Vec::new());
    }
    let text = String::from_utf8_lossy(input);
    let mut fired: Vec<String> = Vec::new();
    let mut out = String::with_capacity(text.len());

    for line in text.split_inclusive('\n') {
        let trimmed = line.trim_start();
        let matched = rules.iter().find(|r| trimmed.starts_with(r.line_prefix));
        match matched {
            Some(rule) => {
                if !fired.iter().any(|f| f == rule.name) {
                    fired.push(rule.name.to_string());
                }
                out.push_str(&format!("<canonicalized:{}>\n", rule.name));
            }
            None => out.push_str(line),
        }
    }
    (out.into_bytes(), fired)
}

/// Compares two texts under a rule set, preferring the strongest true result.
///
/// Raw equality is tried first, so a pair that happens to be byte-identical is reported as
/// such rather than being downgraded merely because rules were available.
pub fn compare_with_rules(left: &[u8], right: &[u8], rules: &[Rule]) -> crate::Comparison {
    let raw = crate::compare_text(left, right);
    if raw.is_equal() {
        return raw;
    }
    let (cl, fired_l) = canonicalize(left, rules);
    let (cr, fired_r) = canonicalize(right, rules);
    let canon = crate::compare_text(&cl, &cr);
    if canon.is_equal() {
        let mut rules_applied = fired_l;
        for r in fired_r {
            if !rules_applied.contains(&r) {
                rules_applied.push(r);
            }
        }
        rules_applied.sort();
        return crate::Comparison::IdenticalAfterCanonicalization { rules_applied };
    }
    // Report the divergence from the canonicalized text: the raw one would point at a
    // timestamp and hide the real difference further down.
    canon
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Comparison;

    #[test]
    fn byte_identical_is_not_downgraded_when_rules_exist() {
        let a = b"# Started on: Mon\nDATA\t1\n";
        assert_eq!(
            compare_with_rules(a, a, PICARD_METRICS),
            Comparison::ByteIdentical,
            "identical inputs must report the strongest result, not a canonicalized one"
        );
    }

    #[test]
    fn a_varying_timestamp_is_excused_and_reported() {
        let a = b"# Started on: Mon Jan 1\nDATA\t1\n";
        let b = b"# Started on: Tue Feb 2\nDATA\t1\n";
        match compare_with_rules(a, b, PICARD_METRICS) {
            Comparison::IdenticalAfterCanonicalization { rules_applied } => {
                assert!(rules_applied.iter().any(|r| r == "picard.started-on"));
            }
            other => panic!("expected canonicalized equality, got {other}"),
        }
    }

    /// The point of the whole module: a real difference must survive canonicalization.
    #[test]
    fn canonicalization_does_not_hide_a_real_difference() {
        let a = b"# Started on: Mon\nDATA\t1\n";
        let b = b"# Started on: Tue\nDATA\t2\n";
        let c = compare_with_rules(a, b, PICARD_METRICS);
        assert!(
            !c.is_equal(),
            "a differing data row must not be excused by a timestamp rule; got {c}"
        );
    }

    #[test]
    fn no_rules_means_raw_comparison() {
        let a = b"x\n";
        let b = b"y\n";
        assert!(!compare_with_rules(a, b, &[]).is_equal());
        let (out, fired) = canonicalize(a, &[]);
        assert_eq!(out, a);
        assert!(fired.is_empty());
    }

    #[test]
    fn every_rule_states_why_it_exists() {
        for set in [PICARD_METRICS, GATK_REPORT, VCF_HEADER] {
            for rule in set {
                assert!(
                    rule.why.len() > 20,
                    "rule `{}` has no real justification; a rule that cannot be justified \
                     should not exist",
                    rule.name
                );
            }
        }
    }
}
