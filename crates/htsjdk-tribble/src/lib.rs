//! Tribble feature codecs, ported from `htsjdk.tribble` (htsjdk 4.2.0).
//!
//! A Tribble codec turns a line of a feature file into a feature, and GATK's `-L` accepts three
//! such files (`.bed`, `.interval_list`, VCF). The interval list has its own reader in
//! [`htsjdk_bam`]; this crate is the rest.

pub mod bed;
