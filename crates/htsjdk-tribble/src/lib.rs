//! Tribble feature codecs, ported from `htsjdk.tribble` (htsjdk 4.2.0).
//!
//! A Tribble codec turns a line of a feature file into a feature, and GATK's `-L` accepts three
//! such files (`.bed`, `.interval_list`, VCF).
//!
//! The interval list has **two** parsers in htsjdk: the reader in [`htsjdk_bam::interval`], which
//! `IntervalList.fromReader` uses, and the codec in [`interval_list`], which Tribble uses. They
//! disagree on which files are valid, and `-L` goes through the codec.

pub mod bed;
pub mod interval_list;
