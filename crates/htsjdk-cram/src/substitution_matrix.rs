//! The CRAM substitution matrix: five bytes inside the preservation map that decide which
//! substitution gets the shortest code.
//!
//! Ported from `htsjdk.samtools.cram.structure.SubstitutionMatrix` and `SubstitutionBase` at
//! htsjdk 4.2.0, and from `java.util.TimSort`'s small-array path, which the ranking depends on.
//!
//! [`crate::preservation_map`] measured the `SM` key as five bytes and left them opaque. Each byte
//! is a packed vector of four two-bit codes, one per possible substitute of that reference base, in
//! the order `A, C, G, T, N` minus the base itself. The codes are **ranks by observed frequency**,
//! so the commonest substitution gets code 0 and therefore the shortest ITF8 in the read features.
//!
//! # The sort is run twice, and the second run is a sort by ordinal in disguise
//!
//! `substitutionCodeVector` sorts by frequency, writes each entry's rank, then **zeroes every
//! frequency** and sorts again with the same comparator, which now falls through to the ordinal
//! tie-break on every pair. A port that sorts once and packs in sorted order puts the right codes
//! in the wrong slots.
//!
//! # The comparator subtracts two longs and casts to int
//!
//! `(int) (o2.freq - o1.freq)`. Two frequencies whose difference is a non-zero multiple of 2^32
//! compare **equal**, and the ordinal tie-break decides instead. Measured, and it is not subtle:
//! with `C` substituted 4294967296 times and nothing else substituted at all, reference base `G`
//! gives the code vector **27**, in which `C` ranks *second* behind a substitution never observed.
//! One more occurrence, 4294967297, gives **75**, in which `C` ranks first. The commonest
//! substitution in the file loses the shortest code to one that never happened.
//!
//! # A lower case reference base can be decoded but not encoded
//!
//! The reading constructor copies each row to its lower case twin, so [`SubstitutionMatrix::base`]
//! accepts `a`; [`SubstitutionMatrix::code`] refuses it by name.
//!
//! # The failure message blames the wrong argument
//!
//! When `base` lands on `NO_BASE` it formats the **reference** base into
//! `Attempt to retrieve a substitution base for invalid base`, even though the reference base is
//! the one thing that was valid enough to index with. The code is what was wrong.
//!
//! # The default matrix is not zeroes
//!
//! With no substitutions observed every frequency ties, the ordinal order wins, and each byte
//! becomes **0x1b**: the ranks 0, 1, 2, 3 packed two bits apiece. That is the `1b1b1b1b1b` the
//! preservation map suite measured in every file.

/// `SubstitutionBase.values()`, whose ordinal order is the tie-break and the packing order.
pub const BASES: [u8; 5] = *b"ACGTN";
/// `SubstitutionMatrix.BASES_SIZE`.
pub const BASES_SIZE: usize = 5;
/// `SubstitutionMatrix.CODES_PER_BASE`, which is one fewer: a base is not its own substitute.
pub const CODES_PER_BASE: usize = BASES_SIZE - 1;
/// `SubstitutionMatrix.SYMBOL_SPACE_SIZE`. Both lookup tables are this square so a base can index
/// them directly, which is why an invalid base is caught by a sign test and not by a bounds check.
pub const SYMBOL_SPACE_SIZE: usize = 128;
/// `SubstitutionMatrix.NO_BASE`.
pub const NO_BASE: u8 = 0;

/// What the two accessors refuse with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixError {
    /// `code` on a reference base that is not positive, or is lower case.
    InvalidOrLowerCaseReferenceBase(u8),
    /// `code` on a read base that is not positive.
    InvalidReadBase(u8),
    /// `base` on a reference base that is not positive.
    InvalidReferenceBase(u8),
    /// `base` where the lookup landed on `NO_BASE`. The message names the **reference** base, not
    /// the code, which is the argument that was actually wrong.
    NoSubstitutionBase(u8),
}

/// Java's `(char)` cast from a `byte`: sign-extend to `int`, then truncate to sixteen bits.
///
/// `%c` on a byte of -1 formats U+FFFF, not U+00FF, and on a byte of 0 formats a NUL. Both appear
/// in these messages, and both are what a port has to reproduce to carry the message verbatim.
fn java_char(value: u8) -> char {
    let widened = (value as i8) as i32 as u32 & 0xFFFF;
    char::from_u32(widened).unwrap_or(char::REPLACEMENT_CHARACTER)
}

impl MatrixError {
    pub fn message(&self) -> String {
        match self {
            MatrixError::InvalidOrLowerCaseReferenceBase(base) => format!(
                "CRAM: Attempt to generate a substitution code for invalid or lower case reference base '{}'",
                java_char(*base)
            ),
            MatrixError::InvalidReadBase(base) => format!(
                "CRAM: Attempt to generate a substitution code for an invalid read base value '{}'",
                java_char(*base)
            ),
            MatrixError::InvalidReferenceBase(base) => format!(
                "CRAM: Attempt to generate a substitution code for invalid reference base '{}'",
                java_char(*base)
            ),
            MatrixError::NoSubstitutionBase(base) => format!(
                "CRAM: Attempt to retrieve a substitution base for invalid base '{}'",
                java_char(*base)
            ),
        }
    }
}

/// The matrix, in both the form it is serialised in and the two lookups it is used through.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubstitutionMatrix {
    encoded: [u8; BASES_SIZE],
    /// `(refBase, readBase) -> code`.
    code_by_base: Vec<[u8; SYMBOL_SPACE_SIZE]>,
    /// `(refBase, code) -> base`. Filled for lower case reference bases too.
    base_by_code: Vec<[u8; SYMBOL_SPACE_SIZE]>,
}

impl SubstitutionMatrix {
    fn empty() -> Self {
        Self {
            encoded: [0u8; BASES_SIZE],
            code_by_base: vec![[0u8; SYMBOL_SPACE_SIZE]; SYMBOL_SPACE_SIZE],
            base_by_code: vec![[0u8; SYMBOL_SPACE_SIZE]; SYMBOL_SPACE_SIZE],
        }
    }

    /// `SubstitutionMatrix(byte[])`: the reading constructor.
    ///
    /// It unpacks each byte into four `(code, base)` pairs, copies each row to its lower case twin,
    /// and only then derives the forward lookup from the reverse one.
    pub fn from_encoded(encoded: [u8; BASES_SIZE]) -> Self {
        let mut out = Self::empty();
        out.encoded = encoded;

        // The substitutes of each base, in the packing order: the base list without the base.
        for (index, reference) in BASES.iter().enumerate() {
            let byte = encoded[index];
            let mut shift = 6i32;
            for substitute in BASES.iter().filter(|b| *b != reference) {
                let code = ((byte >> shift) & 3) as usize;
                out.base_by_code[*reference as usize][code] = *substitute;
                shift -= 2;
            }
            let lower = reference.to_ascii_lowercase() as usize;
            // `System.arraycopy(baseByCode[R], 0, baseByCode[r], 0, 4)`: four entries, not the row.
            // `N` has no lower case twin in htsjdk's unrolled constructor, so it gets none here.
            if *reference != b'N' {
                for code in 0..CODES_PER_BASE {
                    out.base_by_code[lower][code] = out.base_by_code[*reference as usize][code];
                }
            }
        }

        for reference in BASES {
            for code in 0..CODES_PER_BASE {
                let base = out.base_by_code[reference as usize][code] as usize;
                out.code_by_base[reference as usize][base] = code as u8;
            }
        }
        out
    }

    /// `SubstitutionMatrix(List<CRAMCompressionRecord>)`, given the frequencies it would build.
    ///
    /// `frequencies[refBase][readBase]` is the count of that substitution.
    pub fn from_frequencies(frequencies: &[[i64; SYMBOL_SPACE_SIZE]]) -> Self {
        let mut out = Self::empty();
        for (index, reference) in BASES.iter().enumerate() {
            out.encoded[index] =
                out.substitution_code_vector(*reference, &frequencies[*reference as usize]);
        }
        for reference in BASES {
            for substitute in BASES {
                if reference != substitute {
                    let code = out.code_by_base[reference as usize][substitute as usize] as usize;
                    out.base_by_code[reference as usize][code] = substitute;
                    out.base_by_code[reference.to_ascii_lowercase() as usize][code] = substitute;
                }
            }
        }
        out
    }

    /// `substitutionCodeVector`, side effect included: it fills this base's row of the forward
    /// lookup as well as returning the packed byte.
    pub fn substitution_code_vector(
        &mut self,
        reference: u8,
        frequencies: &[i64; SYMBOL_SPACE_SIZE],
    ) -> u8 {
        // (base, ordinal, frequency, rank)
        let mut codes: Vec<(u8, usize, i64, u8)> = BASES
            .iter()
            .enumerate()
            .filter(|(_, base)| **base != reference)
            .map(|(ordinal, base)| (*base, ordinal, frequencies[*base as usize], 0u8))
            .collect();

        java_sort(&mut codes, compare);
        for (rank, entry) in codes.iter_mut().enumerate() {
            entry.3 = rank as u8;
        }
        // Every frequency is zeroed, so the second sort is the ordinal tie-break on every pair.
        for entry in codes.iter_mut() {
            entry.2 = 0;
        }
        java_sort(&mut codes, compare);

        let mut vector = 0u8;
        for entry in &codes {
            vector = (vector << 2) | entry.3;
        }
        for entry in &codes {
            self.code_by_base[reference as usize][entry.0 as usize] = entry.3;
        }
        vector
    }

    /// The five bytes as they are serialised into the preservation map's `SM` key.
    pub fn encoded(&self) -> [u8; BASES_SIZE] {
        self.encoded
    }

    /// `SubstitutionMatrix.code`.
    pub fn code(&self, reference: u8, read: u8) -> Result<u8, MatrixError> {
        if reference == 0 || reference > 127 || reference.is_ascii_lowercase() {
            return Err(MatrixError::InvalidOrLowerCaseReferenceBase(reference));
        }
        if read == 0 || read > 127 {
            return Err(MatrixError::InvalidReadBase(read));
        }
        Ok(self.code_by_base[reference as usize][read as usize])
    }

    /// `SubstitutionMatrix.base`, whose failure names the reference base rather than the code.
    pub fn base(&self, reference: u8, code: u8) -> Result<u8, MatrixError> {
        if reference == 0 || reference > 127 {
            return Err(MatrixError::InvalidReferenceBase(reference));
        }
        let base = self.base_by_code[reference as usize][code as usize];
        if base == NO_BASE {
            return Err(MatrixError::NoSubstitutionBase(reference));
        }
        Ok(base)
    }

    /// `SubstitutionMatrix.toString`: the upper case rows, then the lower case ones.
    pub fn display(&self) -> String {
        let mut out = String::new();
        for reference in BASES {
            out.push(reference as char);
            out.push(':');
            for code in 0..CODES_PER_BASE {
                out.push(self.base_by_code[reference as usize][code] as char);
            }
            out.push('\t');
        }
        for reference in BASES {
            let lower = reference.to_ascii_lowercase();
            out.push(lower as char);
            out.push(':');
            for code in 0..CODES_PER_BASE {
                out.push(self.base_by_code[lower as usize][code] as char);
            }
            out.push('\t');
        }
        out
    }
}

/// The comparator, overflow included: a `long` difference narrowed to an `int`.
fn compare(left: &(u8, usize, i64, u8), right: &(u8, usize, i64, u8)) -> std::cmp::Ordering {
    if left.2 != right.2 {
        // `(int) (o2.freq - o1.freq)`: the subtraction is 64-bit and wraps, then the cast keeps the
        // low 32 bits. A difference that is a non-zero multiple of 2^32 becomes 0, which the sort
        // reads as "equal" even though the branch was taken because they are not.
        let difference = right.2.wrapping_sub(left.2) as i32;
        return difference.cmp(&0);
    }
    // The spec's base order.
    (left.1 as i32).cmp(&(right.1 as i32))
}

/// `java.util.Arrays.sort(T[], Comparator)` for an array below `MIN_MERGE`, which is the only path
/// a four-element array can take.
///
/// Ported rather than delegated to Rust's sort because the comparator above is **not** a total
/// order: with an inconsistent comparator two sorts can disagree, so the one that decides the bytes
/// has to be the one htsjdk runs. Rust's `sort_by` is also entitled to panic on an inconsistent
/// comparator, and htsjdk's does not.
fn java_sort<T, F>(values: &mut [T], compare: F)
where
    F: Fn(&T, &T) -> std::cmp::Ordering,
{
    use std::cmp::Ordering;
    let length = values.len();
    if length < 2 {
        return;
    }
    // `countRunAndMakeAscending`.
    let mut run_hi = 1usize;
    if compare(&values[run_hi], &values[0]) == Ordering::Less {
        run_hi += 1;
        while run_hi < length && compare(&values[run_hi], &values[run_hi - 1]) == Ordering::Less {
            run_hi += 1;
        }
        values[..run_hi].reverse();
    } else {
        run_hi += 1;
        while run_hi < length && compare(&values[run_hi], &values[run_hi - 1]) != Ordering::Less {
            run_hi += 1;
        }
    }

    // `binarySort`, starting past the run that is already ordered.
    for start in run_hi..length {
        let (mut left, mut right) = (0usize, start);
        while left < right {
            let mid = (left + right) >> 1;
            if compare(&values[start], &values[mid]) == Ordering::Less {
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        values[left..=start].rotate_right(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat(counts: [(u8, i64); 5]) -> [i64; SYMBOL_SPACE_SIZE] {
        let mut out = [0i64; SYMBOL_SPACE_SIZE];
        for (base, count) in counts {
            out[base as usize] = count;
        }
        out
    }

    /// No substitutions observed: every frequency ties, the ordinal order wins, and the byte is the
    /// `1b` the preservation map suite found in every file.
    #[test]
    fn the_default_matrix_is_one_b_and_not_zero() {
        let mut matrix = SubstitutionMatrix::empty();
        let frequencies = [0i64; SYMBOL_SPACE_SIZE];
        for reference in BASES {
            assert_eq!(
                matrix.substitution_code_vector(reference, &frequencies),
                0x1b,
                "reference {}",
                reference as char
            );
        }
    }

    /// The commonest substitution in the file loses the shortest code to one that never happened,
    /// because the comparator narrows a long difference to an int.
    #[test]
    fn a_frequency_difference_of_two_to_the_thirty_two_compares_equal() {
        let mut matrix = SubstitutionMatrix::empty();
        let overflowing = flat([
            (b'A', 0),
            (b'C', 4294967296),
            (b'G', 0),
            (b'T', 0),
            (b'N', 0),
        ]);
        assert_eq!(matrix.substitution_code_vector(b'G', &overflowing), 27);

        let mut matrix = SubstitutionMatrix::empty();
        let one_more = flat([
            (b'A', 0),
            (b'C', 4294967297),
            (b'G', 0),
            (b'T', 0),
            (b'N', 0),
        ]);
        assert_eq!(matrix.substitution_code_vector(b'G', &one_more), 75);
    }

    /// A lower case reference base decodes and does not encode.
    #[test]
    fn lower_case_decodes_but_does_not_encode() {
        let matrix = SubstitutionMatrix::from_encoded([0x1b; BASES_SIZE]);
        assert_eq!(matrix.base(b'a', 0), Ok(b'C'));
        assert_eq!(
            matrix.code(b'a', b'C'),
            Err(MatrixError::InvalidOrLowerCaseReferenceBase(b'a'))
        );
    }

    /// The message names the reference base, which is the argument that was fine.
    #[test]
    fn the_failure_blames_the_reference_base_rather_than_the_code() {
        let matrix = SubstitutionMatrix::from_encoded([0x1b; BASES_SIZE]);
        let error = matrix.base(b'A', 100).expect_err("no such code");
        assert_eq!(error, MatrixError::NoSubstitutionBase(b'A'));
        assert!(error.message().contains("invalid base 'A'"));
    }
}
