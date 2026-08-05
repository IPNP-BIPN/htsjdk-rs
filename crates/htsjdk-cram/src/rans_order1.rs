//! rANS 4x8 order 1: the same codec with a context, and a different arithmetic underneath.
//!
//! Ported from `htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Encode.compressOrder1Way4`
//! and `RANS4x8Decode.uncompressOrder1Way4` at htsjdk 4.2.0.
//!
//! [`crate::rans`] pinned order 0. Order 1 is not order 0 with more tables: five things change,
//! and three of them are not in the specification.
//!
//! # The frequency table counts three bytes that are not bigrams
//!
//! Each of the four lanes starts with context 0, and only the first lane's first byte actually
//! follows nothing. So `calcFrequenciesOrder1` adds `F[0][x] += 1` for the byte at each quarter
//! boundary, three counts the input does not contain, and adds 3 to `T[0]`.
//!
//! Measured on the four bytes `ACGT`: context 0 ends up holding `A`, `C`, `G` and `T` at one
//! apiece, so the table says that after nothing, `C` is exactly as likely as `A`, which the input
//! never shows. A port that counts only the bigrams it can see produces a different table on every
//! input.
//!
//! # The normalisation is floating point, where order 0's is fixed point
//!
//! `p = 4096.0 / T[i]`, then `F[i][j] *= p`, which in Java is a compound assignment on an `int`
//! and so truncates toward zero. Two normalisations in the same class, two arithmetics. Order 0's
//! `(F[j] * tr) >> 31` is entirely integer; this one is not.
//!
//! # The frequency table has two levels of run length
//!
//! One over contexts, one over symbols within a context, each with order 0's rule: the run byte
//! appears only once the previous entry was also present, and the reader infers it by peeking at
//! whether the next byte is the current index plus one.
//!
//! # A frequency byte of zero means 4096 on the way in
//!
//! `readStatsOrder1` carries `if (D[i].frequencies[j] == 0) D[i].frequencies[j] = TOTAL_FREQ;` and
//! `readStatsOrder0` has no such line. The writer never emits a zero, so this is a **reader-only**
//! rule. Measured on a table built by hand, because no input can produce one: the five bytes
//! `00 41 00 00 00` give context 0, symbol 65 a frequency of **4096**.
//!
//! # The four lanes read four quarters
//!
//! Order 0's lanes take every fourth byte; order 1's take a contiguous quarter each. That is why
//! the contexts have to be seeded at the quarter boundaries and why the remainder past four
//! quarters belongs to the last lane alone.

use crate::rans::{
    EncodingSymbol, Order, LOWER_BOUND, NUMBER_OF_SYMBOLS, PREFIX_LENGTH, TOTAL_FREQ,
    TOTAL_FREQ_SHIFT,
};

/// One frequency table: 256 contexts of 256 symbols. Heap-allocated because the array is 256 KiB.
pub type Table = Vec<[i32; NUMBER_OF_SYMBOLS]>;

fn empty_table() -> Table {
    vec![[0i32; NUMBER_OF_SYMBOLS]; NUMBER_OF_SYMBOLS]
}

/// `RANS4x8Encode.calcFrequenciesOrder1`.
///
/// The three extra counts at the quarter boundaries are not a rounding detail: they are what makes
/// the four lanes decodable, since three of them begin in the middle of the input with no context
/// to inherit.
pub fn calc_frequencies_order1(input: &[u8]) -> Table {
    let mut frequencies = empty_table();
    let mut totals = [0i32; NUMBER_OF_SYMBOLS];

    let mut last = 0usize;
    for byte in input {
        frequencies[last][*byte as usize] += 1;
        totals[last] += 1;
        last = *byte as usize;
    }

    // The bytes three of the four lanes start on, counted as if they followed nothing.
    let quarter = input.len() >> 2;
    for multiple in 1..=3 {
        frequencies[0][input[quarter * multiple] as usize] += 1;
    }
    totals[0] += 3;

    for context in 0..NUMBER_OF_SYMBOLS {
        if totals[context] == 0 {
            continue;
        }
        // Floating point, and the compound assignment truncates toward zero.
        let p = f64::from(TOTAL_FREQ) / f64::from(totals[context]);
        let mut sum = 0i32;
        let mut max_count = 0i32;
        let mut max_symbol = 0usize;
        for (symbol, frequency) in frequencies[context].iter_mut().enumerate() {
            if *frequency == 0 {
                continue;
            }
            // The maximum is taken on the raw count, in the same pass and before the scaling, so a
            // tie goes to the lowest symbol.
            if max_count < *frequency {
                max_count = *frequency;
                max_symbol = symbol;
            }
            *frequency = (f64::from(*frequency) * p) as i32;
            if *frequency == 0 {
                *frequency = 1;
            }
            sum += *frequency;
        }
        // htsjdk's `else` subtracts zero when the sum is already exact, which is why this is not
        // written as three cases.
        if sum < TOTAL_FREQ {
            frequencies[context][max_symbol] += TOTAL_FREQ - sum;
        } else {
            frequencies[context][max_symbol] -= sum - TOTAL_FREQ;
        }
    }
    frequencies
}

/// `RANS4x8Encode.writeFrequenciesOrder1`. Returns how many bytes were appended.
///
/// Two levels of run length, and the outer one is keyed on the context's **total**, which after
/// normalisation is 4096 for every context that exists and 0 for every context that does not.
pub fn write_frequencies_order1(frequencies: &Table, out: &mut Vec<u8>) -> usize {
    let start = out.len();
    let mut totals = [0i32; NUMBER_OF_SYMBOLS];
    for (context, total) in totals.iter_mut().enumerate() {
        *total = frequencies[context].iter().sum();
    }

    let mut run_contexts = 0i32;
    for context in 0..NUMBER_OF_SYMBOLS {
        if totals[context] == 0 {
            continue;
        }
        if run_contexts != 0 {
            run_contexts -= 1;
        } else {
            out.push(context as u8);
            if context != 0 && totals[context - 1] != 0 {
                let mut end = context + 1;
                while end < NUMBER_OF_SYMBOLS && totals[end] != 0 {
                    end += 1;
                }
                run_contexts = (end - context - 1) as i32;
                out.push(run_contexts as u8);
            }
        }

        let mut run_symbols = 0i32;
        for symbol in 0..NUMBER_OF_SYMBOLS {
            let frequency = frequencies[context][symbol];
            if frequency == 0 {
                continue;
            }
            if run_symbols != 0 {
                run_symbols -= 1;
            } else {
                out.push(symbol as u8);
                if symbol != 0 && frequencies[context][symbol - 1] != 0 {
                    let mut end = symbol + 1;
                    while end < NUMBER_OF_SYMBOLS && frequencies[context][end] != 0 {
                        end += 1;
                    }
                    run_symbols = (end - symbol - 1) as i32;
                    out.push(run_symbols as u8);
                }
            }
            if frequency < 128 {
                out.push(frequency as u8);
            } else {
                out.push((128 | (frequency >> 8)) as u8);
                out.push((frequency & 0xFF) as u8);
            }
        }
        // Every context's symbol list is terminated on its own.
        out.push(0);
    }
    // Then the context list.
    out.push(0);
    out.len() - start
}

/// One encoding symbol per (context, symbol), indexed `[context][symbol]`.
///
/// `RANSEncode.buildSymsOrder1` runs order 0's cumulative walk once per context, so a context's
/// symbols are laid out exactly as order 0 lays out a whole table.
pub fn build_symbols_order1(frequencies: &Table) -> Vec<[EncodingSymbol; NUMBER_OF_SYMBOLS]> {
    let mut symbols = vec![[EncodingSymbol::default(); NUMBER_OF_SYMBOLS]; NUMBER_OF_SYMBOLS];
    for context in 0..NUMBER_OF_SYMBOLS {
        let mut cumulative = 0i32;
        for symbol in 0..NUMBER_OF_SYMBOLS {
            if frequencies[context][symbol] != 0 {
                symbols[context][symbol] =
                    EncodingSymbol::set(cumulative, frequencies[context][symbol], TOTAL_FREQ_SHIFT);
                cumulative += frequencies[context][symbol];
            }
        }
    }
    symbols
}

/// `RANS4x8Encode.compressOrder1Way4`.
///
/// The caller must have checked the length: htsjdk refuses order 1 below four bytes by silently
/// using order 0, which [`crate::rans::order_used`] is.
pub fn compress_order1(input: &[u8]) -> Vec<u8> {
    let frequencies = calc_frequencies_order1(input);
    let symbols = build_symbols_order1(&frequencies);

    let mut out = vec![0u8; PREFIX_LENGTH];
    let frequency_table_size = write_frequencies_order1(&frequencies, &mut out);

    let mut blob: Vec<u8> = Vec::with_capacity(input.len());
    let in_size = input.len();
    let quarter = (in_size >> 2) as isize;
    let mut rans = [LOWER_BOUND; 4];

    // The lane cursors, running backwards, two before the end of each quarter. On a four-byte
    // input `quarter` is 1 and every cursor starts negative, so the main loop never runs and the
    // whole input is encoded by the remainder and the four final symbols.
    let mut cursors: [isize; 4] = [
        quarter - 2,
        2 * quarter - 2,
        3 * quarter - 2,
        4 * quarter - 2,
    ];

    // The symbol each lane will encode against, seeded from one past its cursor.
    let mut last = [0u8; 4];
    for lane in 0..3 {
        if cursors[lane] + 1 >= 0 {
            last[lane] = input[(cursors[lane] + 1) as usize];
        }
    }
    last[3] = input[in_size - 1];

    // `symbols[context][symbol]`: going backwards, the byte at the cursor is the CONTEXT and the
    // byte after it is the symbol being encoded. The naming in htsjdk reads the other way round.
    let mut i3 = in_size as isize - 2;
    while i3 > 4 * quarter - 2 && i3 >= 0 {
        let context = input[i3 as usize];
        rans[3] = symbols[context as usize][last[3] as usize].put(rans[3], &mut blob);
        last[3] = context;
        i3 -= 1;
    }
    cursors[3] = i3;

    while cursors[0] >= 0 {
        let context: [u8; 4] = [
            input[cursors[0] as usize],
            input[cursors[1] as usize],
            input[cursors[2] as usize],
            input[cursors[3] as usize],
        ];
        for lane in (0..4).rev() {
            rans[lane] =
                symbols[context[lane] as usize][last[lane] as usize].put(rans[lane], &mut blob);
        }
        last = context;
        for cursor in cursors.iter_mut() {
            *cursor -= 1;
        }
    }

    // The first symbol of each lane, whose context is 0 because there is nothing before it.
    for lane in (0..4).rev() {
        rans[lane] = symbols[0][last[lane] as usize].put(rans[lane], &mut blob);
    }

    for lane in (0..4).rev() {
        blob.extend_from_slice(&(rans[lane] as u32).to_be_bytes());
    }
    blob.reverse();
    out.extend_from_slice(&blob);

    out[0] = Order::One as u8;
    let compressed_size = (frequency_table_size + blob.len()) as i32;
    out[1..5].copy_from_slice(&compressed_size.to_le_bytes());
    out[5..9].copy_from_slice(&(in_size as i32).to_le_bytes());
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The table holds three counts the input does not contain, one per lane that starts in the
    /// middle of it.
    /// The table holds three counts the input does not contain, one per lane that starts in the
    /// middle of it.
    #[test]
    fn the_table_counts_three_bytes_that_are_not_bigrams() {
        let frequencies = calc_frequencies_order1(b"ACGT");
        // A follows nothing for real; C, G and T are the quarter boundaries at 1, 2 and 3.
        for symbol in *b"ACGT" {
            assert_eq!(
                frequencies[0][symbol as usize], 1024,
                "context 0, symbol {symbol}"
            );
        }
        // The bigrams the input does contain take their whole context.
        assert_eq!(frequencies[b'A' as usize][b'C' as usize], TOTAL_FREQ);
        assert_eq!(frequencies[b'C' as usize][b'G' as usize], TOTAL_FREQ);
        assert_eq!(frequencies[b'G' as usize][b'T' as usize], TOTAL_FREQ);
        // T is last, so it is a context of nothing.
        assert_eq!(frequencies[b'T' as usize].iter().sum::<i32>(), 0);
    }
}
