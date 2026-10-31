//! A deflate compressor written to reproduce another implementation's bytes, not to be fast.
//!
//! ## Why this exists
//!
//! GATK writes BAM files through Intel's GKL, and GKL does not produce the bytes the JDK's
//! deflater produces. htsjdk-rs decision 0029 read the level branch out of `libgkl_compression.so`
//! and found two backends behind one `Deflater`:
//!
//! | level | backend |
//! |---|---|
//! | 1, 2 | ISA-L igzip |
//! | 3 to 9 | zlib 1.2.13 carrying Intel's `deflate_medium` patch |
//!
//! htsjdk's BGZF default is **level 5**, so the second one decides the bytes of every BAM GATK
//! writes without `--use-jdk-deflater`.
//!
//! ## What is here
//!
//! Both flavours, verified two different ways.
//!
//! [`Flavour::Jdk`] is stock zlib, which is what `java.util.zip.Deflater` reaches. It is checked
//! against the C zlib the workspace already links, byte for byte, at every level.
//!
//! [`Flavour::Gkl`] is Intel's fork: `deflate_medium` at levels 4 to 6, and a CRC-32C positional
//! hash in place of zlib's multiplicative rolling one. It is checked against the hashes GKL itself
//! produced in the pinned container, so the assertion is against the real library rather than
//! against a reading of its source.
//!
//! **GKL's output depends on the CPU, and this is where that becomes visible.** The fork selects
//! its hash on `x86_cpu_has_sse42` at load time, so a host without SSE4.2 fills the chains
//! differently and emits different bytes at every level from 3 up. `Flavour::Gkl { sse42 }` makes
//! that a parameter rather than a hidden assumption; the default of `true` is the column every
//! oracle run has been measured in, because every host so far reports SSE4.2.
//!
//! ## What is not here
//!
//! **igzip, for levels 1 and 2.** Those are the only levels GKL does not route through zlib, and
//! nothing in htsjdk asks for them by default. [`deflate_gkl`] refuses them rather than quietly
//! answering with zlib's bytes, which would be a wrong answer that looks like a right one.

mod deflate;
#[cfg(feature = "isal")]
mod igzip;
#[cfg(feature = "isal")]
mod igzip_canary;
mod trees;

pub use deflate::Flavour;

/// Whether `deflate_gkl` can answer levels 1 and 2 on this build and this host.
///
/// False when the crate is built without the `isal` feature, and false when ISA-L is present but
/// fell back to its readable C, which happens without an assembler, without SSE4.2, and on any
/// non-x86 architecture. In that state the library still returns valid deflate; it is simply not
/// GKL's, so those levels refuse. Callers that would rather degrade than fail can ask first.
pub fn igzip_available() -> bool {
    #[cfg(feature = "isal")]
    {
        igzip::usable()
    }
    #[cfg(not(feature = "isal"))]
    {
        false
    }
}

/// Compress `data` as a raw deflate stream the way `java.util.zip.Deflater` does, which is the
/// `nowrap` mode BGZF blocks use.
///
/// `level` is 1 to 9. Level 0 (store only) is not implemented: htsjdk never asks for it, and a
/// path nothing exercises is a path nothing checks.
pub fn deflate(data: &[u8], level: usize) -> Vec<u8> {
    assert!(
        (1..=9).contains(&level),
        "level must be 1 to 9, got {level}"
    );
    deflate_flavour(data, level, Flavour::Jdk)
}

/// [`deflate_flavour`] over and over without allocating: a thread's deflater, reset between
/// blocks. The bytes are those of a new one, which `a_reused_deflater_answers_as_a_new_one` checks
/// for both flavours.
pub struct ReusableDeflater {
    state: deflate::Deflater,
    fresh: bool,
}

impl ReusableDeflater {
    /// `level` is 1 to 9 for [`Flavour::Jdk`] and 3 to 9 for [`Flavour::Gkl`], as in
    /// [`deflate_flavour`].
    pub fn new(level: usize, flavour: Flavour) -> Self {
        match flavour {
            Flavour::Jdk => assert!(
                (1..=9).contains(&level),
                "level must be 1 to 9, got {level}"
            ),
            Flavour::Gkl { .. } => {
                assert!(
                    (3..=9).contains(&level),
                    "GKL levels 1 and 2 are igzip's; got {level}"
                )
            }
        }
        ReusableDeflater {
            state: deflate::Deflater::new(level, flavour),
            fresh: true,
        }
    }

    /// Replaces `out`'s contents with `data`'s deflate stream. `out`'s allocation joins the state's
    /// and comes back, so a caller that keeps passing the same vector allocates nothing.
    pub fn deflate_into(&mut self, data: &[u8], out: &mut Vec<u8>) {
        if !std::mem::take(&mut self.fresh) {
            self.state.reset();
        }
        std::mem::swap(&mut self.state.out, out);
        self.state.out.clear();
        self.state.run(data);
        std::mem::swap(&mut self.state.out, out);
    }
}

/// Compress `data` with an explicit flavour, for callers that need the non-default CPU branch.
pub fn deflate_flavour(data: &[u8], level: usize, flavour: Flavour) -> Vec<u8> {
    match flavour {
        Flavour::Jdk => assert!(
            (1..=9).contains(&level),
            "level must be 1 to 9, got {level}"
        ),
        Flavour::Gkl { .. } => assert!(
            (3..=9).contains(&level),
            "GKL routes levels 1 and 2 through igzip, which is not implemented; got {level}"
        ),
    }
    let mut deflater = deflate::Deflater::new(level, flavour);
    deflater.run(data);
    std::mem::take(&mut deflater.out)
}

/// Level 0 is store only, through GKL's zlib: what `deflate_stored` writes for one `Z_FINISH` over
/// a BGZF block, a single final stored block (header byte, length, its complement, the bytes).
/// Measured against Picard 3.5.0 at `COMPRESSION_LEVEL=0` with GKL and with libdeflate, which
/// write the same bytes. Anything longer than one stored block's 65,535 bytes is refused: BGZF
/// never asks for it, and how this zlib would split one is not measured.
pub fn stored(data: &[u8]) -> Vec<u8> {
    let length = u16::try_from(data.len()).expect("a stored block holds at most 65,535 bytes");
    let mut out = Vec::with_capacity(data.len() + 5);
    out.push(0x01);
    out.extend_from_slice(&length.to_le_bytes());
    out.extend_from_slice(&(!length).to_le_bytes());
    out.extend_from_slice(data);
    out
}

/// Compress `data` the way GKL's `IntelDeflater` does, on a CPU reporting SSE4.2.
///
/// `level` is 0 to 9; 0 is [`stored`]. Levels 1 and 2 go through ISA-L itself rather than through a port of it, for
/// the reason decision 0034 gives: ISA-L's own readable C finds different matches from the
/// assembly kernels it ships, so translating the C would answer confidently and wrongly.
///
/// Without the `isal` feature those two levels are refused rather than approximated.
pub fn deflate_gkl(data: &[u8], level: usize) -> Vec<u8> {
    if level == 0 {
        return stored(data);
    }
    if level == 1 || level == 2 {
        #[cfg(feature = "isal")]
        {
            // GKL passes the level through to ISA-L (decision 0043, correcting 0031, whose four
            // fixtures happen to compress identically at ISA-L levels 1 and 2). The call refuses
            // if this build cannot reproduce GKL; see `igzip::usable`.
            return igzip::deflate(data, level);
        }
        #[cfg(not(feature = "isal"))]
        panic!(
            "GKL routes levels 1 and 2 through igzip; build with the `isal` feature, or use \
             level 3 or above. Answering with zlib's bytes would be a wrong answer wearing a \
             right one's shape."
        );
    }
    assert!(
        (3..=9).contains(&level),
        "level must be 1 to 9, got {level}"
    );
    deflate_flavour(data, level, Flavour::Gkl { sse42: true })
}

#[cfg(test)]
mod reuse_tests {
    use super::*;

    /// A state that has already deflated something answers as a new one does: the window's bytes
    /// past the data, the chains and the bit buffer all start over.
    #[test]
    fn a_reused_deflater_answers_as_a_new_one() {
        let mut seed = 7u64;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 33) as u8
        };
        let bases: Vec<u8> = (0..65_498)
            .map(|_| b"ACGT"[(next() & 3) as usize])
            .collect();
        let noise: Vec<u8> = (0..65_498).map(|_| next()).collect();
        let runs = vec![b'A'; 65_498];
        let short = b"the quick brown fox".to_vec();
        let inputs = [&bases, &noise, &runs, &short, &bases, &short, &noise];
        for flavour in [
            Flavour::Jdk,
            Flavour::Gkl { sse42: true },
            Flavour::Gkl { sse42: false },
        ] {
            let levels = if flavour == Flavour::Jdk {
                1..=9
            } else {
                3..=9
            };
            for level in levels {
                let mut reused = ReusableDeflater::new(level, flavour);
                let mut out = Vec::new();
                for input in inputs {
                    reused.deflate_into(input, &mut out);
                    assert_eq!(
                        out,
                        deflate_flavour(input, level, flavour),
                        "{flavour:?} {level}"
                    );
                }
            }
        }
    }
}
