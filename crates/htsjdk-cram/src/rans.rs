//! rANS 4x8 order 0: the entropy codec an ordinary CRAM actually uses.
//!
//! Ported from `htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Encode`,
//! `RANS4x8Decode`, `RANSEncodingSymbol`, `RANSDecodingSymbol` and `Constants` at htsjdk 4.2.0.
//!
//! [`crate::block`] measured the compression methods present in a four-read CRAM and found RAW,
//! GZIP and rANS. This is the rANS, at the order every writer reaches for first.
//!
//! A stream is a nine-byte prefix, a frequency table, and a blob. What follows is what measurement
//! says rather than what the layout says.
//!
//! # An empty input produces zero bytes
//!
//! Not a nine-byte prefix over an empty table: `compress` returns an empty buffer before it does
//! anything at all. A decoder that reads the prefix unconditionally fails on the one input it
//! should find easiest, so `uncompress` returns empty for empty in the same way.
//!
//! # The requested order is not always the written order
//!
//! Below `MINIMUM_ORDER_1_SIZE` = 4 bytes, `compress` uses order 0 whatever the parameters say, and
//! the order byte records the order it **used**. Measured: asking for order 1 on 0, 1, 2 or 3 bytes
//! writes a 0; asking on 4 writes a 1. [`order_used`] is that rule alone.
//!
//! # The four final states are written big-endian and the blob is then reversed
//!
//! Two reversals that cancel: the states arrive at the head of the blob little-endian, and in the
//! order `rans0, rans1, rans2, rans3`, which is the opposite of the order they were written in.
//!
//! # The normalisation is fixed point, and one symbol absorbs the whole rounding error
//!
//! `tr = (4096 << 31) / T + (1 << 30) / T`, then `F[j] = (F[j] * tr) >> 31`, with a floor of 1 for
//! any symbol that rounds to zero. The residue is then dumped on a single symbol: the one with the
//! largest **raw** count, ties going to the **lowest** index, adjusted up or down so the total is
//! exactly 4096.
//!
//! It is not a small correction. Measured on an input holding symbol `i` exactly `i` times, symbol
//! 254 normalises to 31 and symbol 255, whose fair share is the same 31, is written as **152**.
//!
//! # The frequency table's run lengths start at the second consecutive symbol
//!
//! The run byte is written only when `F[j-1]` is also non-zero, so a run of exactly two symbols
//! writes a run byte of **0**, and the decoder infers that a run byte follows by peeking at whether
//! the next symbol byte is the current symbol plus one. The marker is inferred from the data and
//! never signalled.
//!
//! # A uniform input's blob is four unchanged states
//!
//! Measured: 1000 identical bytes compress to 29. A single symbol normalises to the whole 4096, so
//! its complement frequency is zero and its bias is zero, the state never moves, and the blob is
//! the four initial lower bounds. The frequency table is four bytes and the input length is carried
//! by the prefix.

/// `Constants.TOTAL_FREQ_SHIFT`.
pub const TOTAL_FREQ_SHIFT: u32 = 12;
/// `Constants.TOTAL_FREQ`.
pub const TOTAL_FREQ: i32 = 1 << TOTAL_FREQ_SHIFT;
/// `Constants.NUMBER_OF_SYMBOLS`.
pub const NUMBER_OF_SYMBOLS: usize = 256;
/// `Constants.RANS_4x8_LOWER_BOUND`.
pub const LOWER_BOUND: u64 = 1 << 23;
/// `Constants.RANS_4x8_PREFIX_BYTE_LENGTH`: the order byte, then two int32 lengths.
pub const PREFIX_LENGTH: usize = 1 + 4 + 4;
/// `RANS4x8Encode.MINIMUM_ORDER_1_SIZE`.
pub const MINIMUM_ORDER_1_SIZE: usize = 4;

/// `RANSParams.ORDER`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Order {
    Zero = 0,
    One = 1,
}

impl Order {
    /// `RANSParams.ORDER.fromInt`, which indexes `values()` and turns the bounds failure into an
    /// `IllegalArgumentException`.
    pub fn from_int(value: i32) -> Option<Self> {
        match value {
            0 => Some(Order::Zero),
            1 => Some(Order::One),
            _ => None,
        }
    }
}

/// What a rANS stream is refused with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RansError {
    /// `RANSParams.ORDER.fromInt` on a byte that is neither 0 nor 1.
    UnknownOrder(i32),
    /// The declared compressed length disagrees with what is there. htsjdk's message, verbatim.
    BadLength,
    /// The stream ended inside the frequency table or the blob.
    Truncated,
    /// Order 1 is a separate slice of the port and is not here yet. htsjdk encodes and decodes it;
    /// this says so rather than answering wrongly.
    OrderOneNotPorted,
}

impl RansError {
    pub fn message(&self) -> String {
        match self {
            RansError::UnknownOrder(value) => format!("Unknown rANS order: {value}"),
            RansError::BadLength => {
                "Invalid input length detected in a CRAM rans 4x8 input stream.".to_string()
            }
            RansError::Truncated => "the rANS stream ends inside its own data".to_string(),
            RansError::OrderOneNotPorted => "rANS 4x8 order 1 is not ported".to_string(),
        }
    }
}

/// The order the writer will use, given the order it was asked for.
///
/// Below [`MINIMUM_ORDER_1_SIZE`] the answer is not the question: there is not enough symbol
/// context for order 1, so order 0 is used and the order byte says so.
pub fn order_used(requested: Order, input_length: usize) -> Order {
    if input_length < MINIMUM_ORDER_1_SIZE {
        Order::Zero
    } else {
        requested
    }
}

/// `RANSEncodingSymbol`: the fixed-point reciprocal that lets the hot loop divide by multiplying.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EncodingSymbol {
    /// Exclusive upper bound of the pre-normalisation interval.
    pub x_max: u64,
    /// Fixed-point reciprocal frequency. Held unsigned: Java holds it in an `int` and masks it to
    /// 32 bits at every use, so the sign it carries there is not the number it stands for.
    pub rcp_freq: u32,
    pub bias: u32,
    /// `(1 << scaleBits) - freq`.
    pub cmpl_freq: u32,
    /// Already carries the `+ 32` that saves a shift in the hot loop.
    pub rcp_shift: u32,
}

impl EncodingSymbol {
    /// `RANSEncodingSymbol.set`.
    ///
    /// `shift = ceil(log2(freq))` is computed by counting rather than by a logarithm, and the
    /// reciprocal is Alverson's `((1 << (shift + 31)) + freq - 1) / freq`. A frequency below 2 takes
    /// a different branch entirely: the reciprocal is all ones and the bias carries the rounding.
    pub fn set(start: i32, freq: i32, scale_bits: u32) -> Self {
        let x_max = (1u64 << (31 - scale_bits)) * freq as u64;
        let cmpl_freq = ((1i64 << scale_bits) - freq as i64) as u32;
        if freq < 2 {
            Self {
                x_max,
                rcp_freq: u32::MAX,
                bias: (start + (1 << scale_bits) - 1) as u32,
                cmpl_freq,
                // `rcpShift = 0`, then the unconditional `+= 32`.
                rcp_shift: 32,
            }
        } else {
            let mut shift = 0u32;
            while freq as u64 > (1u64 << shift) {
                shift += 1;
            }
            Self {
                x_max,
                rcp_freq: (((1i64 << (shift + 31)) + freq as i64 - 1) / freq as i64) as u32,
                bias: start as u32,
                cmpl_freq,
                rcp_shift: shift - 1 + 32,
            }
        }
    }

    /// `RANSEncodingSymbol.putSymbol4x8`: renormalise by emitting at most two bytes, then advance.
    ///
    /// At most two, not a loop: the state is bounded so a third byte is never needed, and htsjdk
    /// writes the two tests out by hand rather than looping.
    pub fn put(&self, r: u64, out: &mut Vec<u8>) -> u64 {
        let mut r = r;
        if r >= self.x_max {
            out.push((r & 0xFF) as u8);
            r >>= 8;
            if r >= self.x_max {
                out.push((r & 0xFF) as u8);
                r >>= 8;
            }
        }
        let q = (r * u64::from(self.rcp_freq)) >> self.rcp_shift;
        r + u64::from(self.bias) + q * u64::from(self.cmpl_freq)
    }
}

/// `RANS4x8Encode.calcFrequenciesOrder0`: count, scale by a fixed-point reciprocal, then hand the
/// entire residue to one symbol.
pub fn calc_frequencies_order0(input: &[u8]) -> [i32; NUMBER_OF_SYMBOLS] {
    let total = input.len() as i64;
    let mut frequencies = [0i32; NUMBER_OF_SYMBOLS];
    for byte in input {
        frequencies[*byte as usize] += 1;
    }
    if total == 0 {
        return frequencies;
    }

    // The symbol with the largest *raw* count, ties to the lowest index, chosen before anything is
    // scaled. It is the one that will absorb the rounding error.
    let mut max_count = 0;
    let mut max_symbol = 0usize;
    for (symbol, count) in frequencies.iter().enumerate() {
        if max_count < *count {
            max_count = *count;
            max_symbol = symbol;
        }
    }

    let tr = ((i64::from(TOTAL_FREQ) << 31) / total) + ((1i64 << 30) / total);
    let mut sum = 0i32;
    for frequency in frequencies.iter_mut() {
        if *frequency == 0 {
            continue;
        }
        *frequency = ((i64::from(*frequency) * tr) >> 31) as i32;
        if *frequency == 0 {
            // A symbol that occurs must not scale to nothing, or the decoder cannot spell it.
            *frequency = 1;
        }
        sum += *frequency;
    }

    // htsjdk's `else` runs when the sum is already exact and subtracts zero, which is why this is
    // not written as three cases.
    if sum < TOTAL_FREQ {
        frequencies[max_symbol] += TOTAL_FREQ - sum;
    } else {
        frequencies[max_symbol] -= sum - TOTAL_FREQ;
    }
    frequencies
}

/// `RANS4x8Encode.writeFrequenciesOrder0`. Returns how many bytes were appended.
pub fn write_frequencies_order0(
    frequencies: &[i32; NUMBER_OF_SYMBOLS],
    out: &mut Vec<u8>,
) -> usize {
    let start = out.len();
    let mut rle = 0i32;
    for symbol in 0..NUMBER_OF_SYMBOLS {
        let frequency = frequencies[symbol];
        if frequency == 0 {
            continue;
        }
        if rle != 0 {
            // Inside a run: the symbol byte is implied by the run that announced it.
            rle -= 1;
        } else {
            out.push(symbol as u8);
            // A run byte only once the *previous* symbol was also present, so a run of two writes a
            // run byte of zero. htsjdk's `rle == 0` here is already known and is left out.
            if symbol != 0 && frequencies[symbol - 1] != 0 {
                let mut end = symbol + 1;
                while end < NUMBER_OF_SYMBOLS && frequencies[end] != 0 {
                    end += 1;
                }
                rle = (end - symbol - 1) as i32;
                out.push(rle as u8);
            }
        }
        if frequency < 128 {
            out.push(frequency as u8);
        } else {
            out.push((128 | (frequency >> 8)) as u8);
            out.push((frequency & 0xFF) as u8);
        }
    }
    // The zero that ends the table, which is also why a symbol-zero entry cannot end it.
    out.push(0);
    out.len() - start
}

/// The encoding symbols of a normalised frequency table, in cumulative order.
///
/// `RANSEncode.updateEncodingSymbols`: only present symbols get one, and the cumulative frequency
/// advances only over them.
pub fn build_symbols_order0(
    frequencies: &[i32; NUMBER_OF_SYMBOLS],
) -> [EncodingSymbol; NUMBER_OF_SYMBOLS] {
    let mut symbols = [EncodingSymbol::default(); NUMBER_OF_SYMBOLS];
    let mut cumulative = 0i32;
    for symbol in 0..NUMBER_OF_SYMBOLS {
        if frequencies[symbol] != 0 {
            symbols[symbol] =
                EncodingSymbol::set(cumulative, frequencies[symbol], TOTAL_FREQ_SHIFT);
            cumulative += frequencies[symbol];
        }
    }
    symbols
}

/// `RANS4x8Encode.compressOrder0Way4`.
pub fn compress_order0(input: &[u8]) -> Vec<u8> {
    if input.is_empty() {
        return Vec::new();
    }

    let frequencies = calc_frequencies_order0(input);
    let symbols = build_symbols_order0(&frequencies);

    let mut out = Vec::with_capacity(input.len() + 256 * 3 + PREFIX_LENGTH + 16);
    out.resize(PREFIX_LENGTH, 0);
    let frequency_table_size = write_frequencies_order0(&frequencies, &mut out);

    // The blob is built forwards and reversed whole at the end, so everything below writes in the
    // opposite order to the one it will be read in.
    let mut blob: Vec<u8> = Vec::with_capacity(input.len());
    let in_size = input.len();
    let mut rans = [LOWER_BOUND; 4];

    // The tail, whose four cases fall through in htsjdk and so are three tests here. The indices
    // are the last `in_size & 3` bytes, taken from the highest downwards but assigned to states 2,
    // 1 and 0 in that order.
    let tail = in_size & 3;
    if tail == 3 {
        rans[2] = symbols[input[in_size - (tail - 2)] as usize].put(rans[2], &mut blob);
    }
    if tail >= 2 {
        rans[1] = symbols[input[in_size - (tail - 1)] as usize].put(rans[1], &mut blob);
    }
    if tail >= 1 {
        rans[0] = symbols[input[in_size - tail] as usize].put(rans[0], &mut blob);
    }

    let mut i = in_size & !3;
    while i > 0 {
        for lane in (0..4).rev() {
            rans[lane] = symbols[input[i - 4 + lane] as usize].put(rans[lane], &mut blob);
        }
        i -= 4;
    }

    // Big-endian, and reversed a moment later into little-endian at the head of the blob.
    for lane in (0..4).rev() {
        blob.extend_from_slice(&(rans[lane] as u32).to_be_bytes());
    }
    blob.reverse();

    out.extend_from_slice(&blob);

    // The prefix, written last because its lengths are only known now.
    out[0] = Order::Zero as u8;
    let compressed_size = (frequency_table_size + blob.len()) as i32;
    out[1..5].copy_from_slice(&compressed_size.to_le_bytes());
    out[5..9].copy_from_slice(&(in_size as i32).to_le_bytes());
    out
}

/// `RANS4x8Encode.compress`, including the order it silently substitutes on short input.
pub fn compress(input: &[u8], requested: Order) -> Result<Vec<u8>, RansError> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    Ok(match order_used(requested, input.len()) {
        Order::Zero => compress_order0(input),
        Order::One => crate::rans_order1::compress_order1(input),
    })
}

/// One symbol's range, as the decoder holds it.
#[derive(Debug, Clone, Copy, Default)]
struct DecodingSymbol {
    start: i32,
    freq: i32,
}

impl DecodingSymbol {
    /// `RANSDecodingSymbol.advanceSymbolStep`.
    ///
    /// The subtraction cannot go below zero on a well-formed stream, because the cumulative
    /// frequency it was looked up by lies in `[start, start + freq)`. It wraps rather than panics
    /// so a malformed stream produces a wrong answer where htsjdk produces a wrong answer, instead
    /// of an abort where htsjdk has none.
    fn advance_step(&self, r: u64) -> u64 {
        let mask = (1u64 << TOTAL_FREQ_SHIFT) - 1;
        ((self.freq as u64) * (r >> TOTAL_FREQ_SHIFT) + (r & mask)).wrapping_sub(self.start as u64)
    }
}

/// `Utils.RANSDecodeRenormalize4x8`: pull bytes until the state is back above the lower bound.
fn renormalise(r: u64, blob: &[u8], at: &mut usize) -> Result<u64, RansError> {
    let mut r = r;
    while r < LOWER_BOUND {
        let byte = *blob.get(*at).ok_or(RansError::Truncated)?;
        *at += 1;
        r = (r << 8) | u64::from(byte);
    }
    Ok(r)
}

/// The frequency table, and the reverse lookup that maps a cumulative frequency back to a symbol.
struct Stats {
    symbols: [DecodingSymbol; NUMBER_OF_SYMBOLS],
    reverse_lookup: Vec<u8>,
}

/// `RANS4x8Decode.readStatsOrder0`.
///
/// The run marker is a peek: when the next symbol byte is the current symbol plus one, the byte
/// after it is a run length. Nothing in the stream says so.
fn read_stats_order0(bytes: &[u8], at: &mut usize) -> Result<Stats, RansError> {
    let next = |at: &mut usize| -> Result<u8, RansError> {
        let byte = *bytes.get(*at).ok_or(RansError::Truncated)?;
        *at += 1;
        Ok(byte)
    };

    let mut stats = Stats {
        symbols: [DecodingSymbol::default(); NUMBER_OF_SYMBOLS],
        reverse_lookup: vec![0u8; TOTAL_FREQ as usize],
    };
    let mut rle = 0i32;
    let mut cumulative = 0i32;
    let mut symbol = next(at)? as usize;
    loop {
        let first = next(at)?;
        let frequency = if first >= 0x80 {
            (i32::from(first & 0x7F) << 8) | i32::from(next(at)?)
        } else {
            i32::from(first)
        };

        stats.symbols[symbol] = DecodingSymbol {
            start: cumulative,
            freq: frequency,
        };
        let from = cumulative.max(0) as usize;
        let to = (cumulative + frequency).max(0) as usize;
        if to > stats.reverse_lookup.len() {
            return Err(RansError::Truncated);
        }
        stats.reverse_lookup[from..to].fill(symbol as u8);
        cumulative += frequency;

        let peek = usize::from(*bytes.get(*at).ok_or(RansError::Truncated)?);
        if rle == 0 && symbol + 1 == peek {
            symbol = next(at)? as usize;
            rle = i32::from(next(at)?);
        } else if rle != 0 {
            rle -= 1;
            symbol += 1;
        } else {
            symbol = next(at)? as usize;
        }
        if symbol == 0 {
            break;
        }
    }
    Ok(stats)
}

/// `RANS4x8Decode.uncompress`, order 0.
///
/// An empty input decodes to an empty output, which is the counterpart of an empty input
/// compressing to nothing.
pub fn uncompress(input: &[u8]) -> Result<Vec<u8>, RansError> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    if input.len() < PREFIX_LENGTH {
        return Err(RansError::Truncated);
    }

    let order =
        Order::from_int(i32::from(input[0])).ok_or(RansError::UnknownOrder(i32::from(input[0])))?;
    let compressed_size = i32::from_le_bytes(input[1..5].try_into().expect("four bytes"));
    // htsjdk compares against `remaining() - 4` after reading the order and the compressed length,
    // which is the whole input less the nine-byte prefix.
    if compressed_size as i64 != (input.len() - PREFIX_LENGTH) as i64 {
        return Err(RansError::BadLength);
    }
    let out_size = i32::from_le_bytes(input[5..9].try_into().expect("four bytes")).max(0) as usize;
    if order == Order::One {
        return Err(RansError::OrderOneNotPorted);
    }

    let mut at = PREFIX_LENGTH;
    let stats = read_stats_order0(input, &mut at)?;

    let mut rans = [0u64; 4];
    for state in rans.iter_mut() {
        let bytes: [u8; 4] = input
            .get(at..at + 4)
            .ok_or(RansError::Truncated)?
            .try_into()
            .expect("four bytes");
        *state = u64::from(u32::from_le_bytes(bytes));
        at += 4;
    }

    let mut out = vec![0u8; out_size];
    let end = out_size & !3;
    let mut i = 0usize;
    while i < end {
        for lane in 0..4 {
            let symbol = stats.reverse_lookup[(rans[lane] & (TOTAL_FREQ as u64 - 1)) as usize];
            out[i + lane] = symbol;
        }
        for lane in 0..4 {
            let symbol = out[i + lane];
            rans[lane] = stats.symbols[symbol as usize].advance_step(rans[lane]);
            rans[lane] = renormalise(rans[lane], input, &mut at)?;
        }
        i += 4;
    }

    // The tail, decoded from states 0, 1 and 2 in that order, which is where the encoder put it.
    for lane in 0..(out_size & 3) {
        let symbol = stats.reverse_lookup[(rans[lane] & (TOTAL_FREQ as u64 - 1)) as usize];
        out[end + lane] = symbol;
        let advanced = stats.symbols[symbol as usize].advance_step(rans[lane]);
        rans[lane] = renormalise(advanced, input, &mut at)?;
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The order a caller asks for is not the order that gets written.
    #[test]
    fn order_one_is_refused_below_four_bytes() {
        for length in 0..4 {
            assert_eq!(
                order_used(Order::One, length),
                Order::Zero,
                "{length} bytes"
            );
        }
        assert_eq!(order_used(Order::One, 4), Order::One);
        assert_eq!(order_used(Order::Zero, 4), Order::Zero);
    }

    #[test]
    fn an_empty_input_produces_no_bytes_at_all() {
        assert!(compress_order0(&[]).is_empty());
        assert_eq!(uncompress(&[]), Ok(Vec::new()));
    }

    /// One symbol takes the whole table, so its complement frequency is zero, the state never
    /// moves, and the blob is the four initial lower bounds.
    #[test]
    fn a_uniform_input_never_moves_its_states() {
        let frequencies = calc_frequencies_order0(&[b'A'; 1000]);
        assert_eq!(frequencies[b'A' as usize], TOTAL_FREQ);
        let symbol = EncodingSymbol::set(0, TOTAL_FREQ, TOTAL_FREQ_SHIFT);
        assert_eq!(symbol.cmpl_freq, 0);
        assert_eq!(symbol.bias, 0);
        let compressed = compress_order0(&[b'A'; 1000]);
        assert_eq!(compressed.len(), PREFIX_LENGTH + 4 + 16);
    }

    /// The residue lands on one symbol, and it is not a rounding nudge.
    #[test]
    fn one_symbol_absorbs_the_whole_rounding_error() {
        let mut input = Vec::new();
        for symbol in 0..256usize {
            input.extend(std::iter::repeat_n(symbol as u8, symbol));
        }
        let frequencies = calc_frequencies_order0(&input);
        assert_eq!(frequencies[254], 31);
        assert_eq!(
            frequencies[255], 152,
            "the maximum symbol carries the residue"
        );
        assert_eq!(frequencies.iter().sum::<i32>(), TOTAL_FREQ);
    }

    /// A run of exactly two symbols writes a run byte of zero.
    #[test]
    fn a_run_of_two_writes_a_run_length_of_zero() {
        let mut frequencies = [0i32; NUMBER_OF_SYMBOLS];
        frequencies[5] = 2048;
        frequencies[6] = 2048;
        let mut out = Vec::new();
        write_frequencies_order0(&frequencies, &mut out);
        // symbol 5, freq 2048 in two bytes, symbol 6, run length 0, freq 2048, terminator.
        assert_eq!(out, vec![5, 0x88, 0x00, 6, 0, 0x88, 0x00, 0]);
    }

    #[test]
    fn every_shape_round_trips() {
        let mut inputs: Vec<Vec<u8>> = vec![
            b"A".to_vec(),
            b"AB".to_vec(),
            b"ACG".to_vec(),
            b"ACGT".to_vec(),
            vec![b'A'; 1000],
            (0..=255u8).collect(),
        ];
        for length in [1000, 1001, 1002, 1003] {
            inputs.push(
                b"ACGT"
                    .iter()
                    .copied()
                    .cycle()
                    .take(length)
                    .collect::<Vec<u8>>(),
            );
        }
        let mut seed = 0x1234_5678u64;
        let noise: Vec<u8> = (0..10_000)
            .map(|_| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (seed >> 33) as u8
            })
            .collect();
        inputs.push(noise);

        for input in inputs {
            let compressed = compress_order0(&input);
            assert_eq!(
                uncompress(&compressed).map_err(|e| e.message()),
                Ok(input.clone()),
                "{} bytes",
                input.len()
            );
        }
    }

    #[test]
    fn a_length_that_disagrees_with_the_stream_is_refused() {
        let mut compressed = compress_order0(b"ACGTACGT");
        compressed[1] = compressed[1].wrapping_add(1);
        assert_eq!(uncompress(&compressed), Err(RansError::BadLength));
    }
}
