//! The CRAM 3.1 fqzcomp quality codec, block compression method 7.
//!
//! Ported from `htsjdk.samtools.cram.compression.fqzcomp` at htsjdk 5.0.0. Each quality is coded
//! with an adaptive byte model chosen by a 16-bit context built from the previous qualities, the
//! position in the read and a running count of changes; the read lengths that delimit records
//! are coded inline, because the codec needs them and nothing else in the block carries them.

use crate::range::{ByteModel, RangeCoder};
use crate::rans_nx16::{Input, NxError};

const SUPPORTED_VERSION: u8 = 5;

const MULTI_PARAM_FLAG_MASK: u8 = 0x01;
const SELECTOR_TABLE_FLAG_MASK: u8 = 0x02;
const DO_REVERSE_FLAG_MASK: u8 = 0x04;

const DEDUP_FLAG_MASK: u8 = 0x02;
const FIXED_LEN_FLAG_MASK: u8 = 0x04;
const SEL_FLAG_MASK: u8 = 0x08;
const QMAP_FLAG_MASK: u8 = 0x10;
const PTAB_FLAG_MASK: u8 = 0x20;
const DTAB_FLAG_MASK: u8 = 0x40;
const QTAB_FLAG_MASK: u8 = 0x80;

/// `FQZParam`: one parameter block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Param {
    pub(crate) context: u32,
    pub(crate) flags: u8,
    pub(crate) do_dedup: bool,
    /// Positive until the first record sets it to the negated length; zero when lengths vary.
    pub(crate) fixed_len: i32,
    pub(crate) do_sel: bool,
    pub(crate) do_qmap: bool,
    pub(crate) do_pos: bool,
    pub(crate) do_delta: bool,
    pub(crate) do_qtab: bool,
    pub(crate) max_symbols: u32,
    pub(crate) quality_context_bits: u32,
    pub(crate) quality_context_shift: u32,
    pub(crate) quality_context_location: u32,
    pub(crate) selector_context_location: u32,
    pub(crate) position_context_location: u32,
    pub(crate) delta_context_location: u32,
    pub(crate) quality_map: Vec<u32>,
    pub(crate) quality_context_table: Vec<u32>,
    pub(crate) position_context_table: Vec<u32>,
    pub(crate) delta_context_table: Vec<u32>,
}

impl Param {
    fn read(input: &mut Input<'_>) -> Result<Self, NxError> {
        let context = u32::from(input.byte()?) | (u32::from(input.byte()?) << 8);
        let flags = input.byte()?;
        let max_symbols = u32::from(input.byte()?);
        let x = input.byte()?;
        let y = input.byte()?;
        let z = input.byte()?;
        let mut param = Param {
            context,
            flags,
            do_dedup: flags & DEDUP_FLAG_MASK != 0,
            fixed_len: i32::from(flags & FIXED_LEN_FLAG_MASK),
            do_sel: flags & SEL_FLAG_MASK != 0,
            do_qmap: flags & QMAP_FLAG_MASK != 0,
            do_pos: flags & PTAB_FLAG_MASK != 0,
            do_delta: flags & DTAB_FLAG_MASK != 0,
            do_qtab: flags & QTAB_FLAG_MASK != 0,
            max_symbols,
            quality_context_bits: u32::from(x >> 4),
            quality_context_shift: u32::from(x & 0x0F),
            quality_context_location: u32::from(y >> 4),
            selector_context_location: u32::from(y & 0x0F),
            position_context_location: u32::from(z >> 4),
            delta_context_location: u32::from(z & 0x0F),
            quality_map: (0..256).collect(),
            quality_context_table: vec![0; 1024],
            position_context_table: Vec::new(),
            delta_context_table: Vec::new(),
        };
        if param.do_qmap {
            for i in 0..max_symbols as usize {
                param.quality_map[i] = u32::from(input.byte()?);
            }
        }
        if param.quality_context_bits > 0 && param.do_qtab {
            read_array(input, &mut param.quality_context_table, 256)?;
        } else {
            for i in 0..256 {
                param.quality_context_table[i] = i as u32;
            }
        }
        if param.do_pos {
            param.position_context_table = vec![0; 1024];
            read_array(input, &mut param.position_context_table, 1024)?;
        }
        if param.do_delta {
            param.delta_context_table = vec![0; 256];
            read_array(input, &mut param.delta_context_table, 256)?;
        }
        Ok(param)
    }
}

/// `FQZParams`: the global flags, the selector table and every parameter block.
pub(crate) struct Params {
    pub(crate) flags: u8,
    pub(crate) params: Vec<Param>,
    pub(crate) max_selector: u32,
    pub(crate) max_symbol: u32,
    pub(crate) selector_table: Vec<u32>,
}

impl Params {
    fn read(input: &mut Input<'_>) -> Result<Self, NxError> {
        let flags = input.byte()?;
        // A signed byte in the reference.
        let count = if flags & MULTI_PARAM_FLAG_MASK != 0 {
            i32::from(input.byte()? as i8)
        } else {
            1
        };
        if count < 0 {
            return Err(NxError::Cram(format!("Illegal Capacity: {count}")));
        }
        let mut max_selector = if count > 1 { count as u32 } else { 0 };
        let mut selector_table = vec![0u32; 256];
        if flags & SELECTOR_TABLE_FLAG_MASK != 0 {
            max_selector = u32::from(input.byte()?);
            read_array(input, &mut selector_table, 256)?;
        } else {
            for (i, slot) in selector_table.iter_mut().enumerate() {
                *slot = if (i as i32) < count {
                    i as u32
                } else {
                    (count - 1).max(0) as u32
                };
            }
        }
        let mut params = Vec::with_capacity(count as usize);
        let mut max_symbol = 0;
        for _ in 0..count {
            let param = Param::read(input)?;
            max_symbol = max_symbol.max(param.max_symbols);
            params.push(param);
        }
        Ok(Self {
            flags,
            params,
            max_selector,
            max_symbol,
            selector_table,
        })
    }
}

/// `FQZUtils.readArray`: a doubly run-length coded table.
fn read_array(input: &mut Input<'_>, table: &mut [u32], size: usize) -> Result<(), NxError> {
    let mut runs: Vec<u32> = Vec::with_capacity(1024);
    let mut z = 0usize;
    let mut last: i32 = -1;
    while z < size {
        let run = u32::from(input.byte()?);
        runs.push(run);
        z += run as usize;
        if run as i32 == last {
            let copy = u32::from(input.byte()?);
            z += (run * copy) as usize;
            for _ in 0..copy {
                runs.push(run);
            }
        }
        last = run as i32;
    }
    let mut i = 0u32;
    let mut j = 0usize;
    z = 0;
    while z < size {
        let mut length = 0u32;
        loop {
            let part = *runs.get(j).ok_or(NxError::Truncated)?;
            j += 1;
            length += part;
            if part != 255 {
                break;
            }
        }
        for _ in 0..length {
            *table.get_mut(z).ok_or(NxError::Truncated)? = i;
            z += 1;
        }
        i += 1;
    }
    Ok(())
}

/// The adaptive models, the 65536 quality ones built when first used.
pub(crate) struct Models {
    quality: Vec<Option<ByteModel>>,
    quality_symbols: usize,
    pub(crate) length: Vec<ByteModel>,
    pub(crate) reverse: ByteModel,
    pub(crate) duplicate: ByteModel,
    pub(crate) selector: Option<ByteModel>,
}

impl Models {
    pub(crate) fn new(max_symbol: u32, max_selector: u32) -> Self {
        Self {
            quality: (0..1 << 16).map(|_| None).collect(),
            quality_symbols: max_symbol as usize + 1,
            length: vec![ByteModel::new(256); 4],
            reverse: ByteModel::new(2),
            duplicate: ByteModel::new(2),
            selector: (max_selector > 0).then(|| ByteModel::new(max_selector as usize + 1)),
        }
    }

    pub(crate) fn quality(&mut self, context: usize) -> &mut ByteModel {
        let symbols = self.quality_symbols;
        self.quality[context].get_or_insert_with(|| ByteModel::new(symbols))
    }
}

/// `FQZState`.
#[derive(Default)]
pub(crate) struct State {
    pub(crate) quality_context: u32,
    pub(crate) previous_quality: u32,
    pub(crate) delta: u32,
    pub(crate) bases: i32,
    pub(crate) selector: u32,
    pub(crate) selector_table: u32,
    pub(crate) record_length: i32,
    pub(crate) is_duplicate: bool,
    pub(crate) read_ordinal: usize,
    pub(crate) quality_lengths: Vec<i32>,
    pub(crate) reverse: Vec<bool>,
}

/// `fqzUpdateContext`.
pub(crate) fn update_context(param: &Param, state: &mut State, quality: u32) -> usize {
    let mut last = param.context;
    state.quality_context = (state.quality_context << param.quality_context_shift)
        .wrapping_add(param.quality_context_table[quality as usize]);
    last = last.wrapping_add(
        (state.quality_context & ((1u32 << param.quality_context_bits) - 1))
            << param.quality_context_location,
    );
    if param.do_pos {
        let at = state.bases.clamp(0, 1023) as usize;
        last =
            last.wrapping_add(param.position_context_table[at] << param.position_context_location);
    }
    if param.do_delta {
        let at = state.delta.min(255) as usize;
        last = last.wrapping_add(param.delta_context_table[at] << param.delta_context_location);
        if state.previous_quality != quality {
            state.delta += 1;
        }
        state.previous_quality = quality;
    }
    if param.do_sel {
        last = last.wrapping_add(state.selector << param.selector_context_location);
    }
    state.bases -= 1;
    (last & 0xffff) as usize
}

/// `FQZCompDecode.uncompress`.
pub fn uncompress(bytes: &[u8]) -> Result<Vec<u8>, NxError> {
    let mut input = Input::new(bytes);
    let length = input.uint7()? as usize;
    let version = input.byte()?;
    if version != SUPPORTED_VERSION {
        return Err(NxError::Cram(format!(
            "Invalid FQZComp format version number: {version}"
        )));
    }
    let mut params = Params::read(&mut input)?;
    let mut models = Models::new(params.max_symbol, params.max_selector);
    let mut state = State::default();
    let mut coder = RangeCoder::new();
    coder.decode_start(&mut input)?;
    let mut out = vec![0u8; length];
    let do_reverse = params.flags & DO_REVERSE_FLAG_MASK != 0;

    let mut current = 0usize;
    let mut last = 0usize;
    let mut i = 0usize;
    while i < length {
        if state.bases == 0 {
            new_record(
                &mut input,
                &mut coder,
                &mut models,
                &mut params,
                &mut state,
                do_reverse,
            )?;
            if state.is_duplicate {
                let n = state.record_length.max(0) as usize;
                if n > i || i + n > length {
                    return Err(NxError::Truncated);
                }
                out.copy_within(i - n..i, i);
                i += n;
                last = 0;
                state.bases = 0;
                state.is_duplicate = false;
                state.quality_lengths.push(state.record_length);
                continue;
            }
            state.quality_lengths.push(state.record_length);
            current = state.selector_table as usize;
            // `last = fqzState.getContext()`, which was just reset to zero: the first quality of
            // a record is coded in context 0, not in the parameter block's.
            last = 0;
        }
        let param = params.params.get(current).ok_or(NxError::Truncated)?;
        let quality = models.quality(last).decode(&mut input, &mut coder)?;
        out[i] = if param.do_qmap {
            param.quality_map[quality as usize] as u8
        } else {
            quality as u8
        };
        i += 1;
        last = update_context(param, &mut state, quality);
    }

    if do_reverse {
        let mut index = 0usize;
        for record in 0..state.read_ordinal {
            if index >= length {
                break;
            }
            let n = state.quality_lengths[record].max(0) as usize;
            if state.reverse.get(record).copied().unwrap_or(false) {
                let end = (index + n).min(length);
                out[index..end].reverse();
            }
            index += n;
        }
    }
    Ok(out)
}

/// `decodeFQZNewRecord`.
fn new_record(
    input: &mut Input<'_>,
    coder: &mut RangeCoder,
    models: &mut Models,
    params: &mut Params,
    state: &mut State,
    do_reverse: bool,
) -> Result<(), NxError> {
    state.selector = match models.selector.as_mut() {
        Some(model) if params.max_selector > 0 => model.decode(input, coder)?,
        _ => 0,
    };
    state.selector_table = *params
        .selector_table
        .get(state.selector as usize)
        .ok_or(NxError::Truncated)?;
    let param = params
        .params
        .get_mut(state.selector_table as usize)
        .ok_or(NxError::Truncated)?;
    let length = if param.fixed_len >= 0 {
        let mut length = models.length[0].decode(input, coder)? as i32;
        length |= (models.length[1].decode(input, coder)? as i32) << 8;
        length |= (models.length[2].decode(input, coder)? as i32) << 16;
        length |= (models.length[3].decode(input, coder)? as i32) << 24;
        if param.fixed_len > 0 {
            param.fixed_len = -length;
        }
        length
    } else {
        -param.fixed_len
    };
    state.record_length = length;
    if do_reverse {
        let reverse = models.reverse.decode(input, coder)? != 0;
        state.reverse.push(reverse);
    }
    state.is_duplicate = false;
    if param.do_dedup && models.duplicate.decode(input, coder)? != 0 {
        state.is_duplicate = true;
    }
    state.bases = length;
    state.delta = 0;
    state.quality_context = 0;
    state.previous_quality = 0;
    state.read_ordinal += 1;
    Ok(())
}
