//! bzip2 as Apache commons-compress 1.26.2's `BZip2CompressorOutputStream` writes it at block
//! size 9, which is what htsjdk 5.0.0's `BZIP2ExternalCompressor` hands a CRAM block to.
//!
//! The format is bzip2's, and libbzip2 decodes it, but libbzip2's encoder makes other bytes from
//! the same input. commons-compress is a port of libbzip2 0.9.x/1.0 that drifted in a few places,
//! each of which changes the stream:
//!
//! * a block holds up to `900000 - 20` run-length-coded bytes before the next run, not
//!   `900000 - 19`;
//! * the coding tables are built with a longest code of 20 bits, not 17, and the rescaling pass
//!   that shortens a table too long skips the last symbol;
//! * each group of fifty symbols takes the cheapest table counting down from the last, so a tie
//!   goes to the highest table rather than the lowest;
//! * `mainQSort3` pushes its three partitions in a fixed order, and a block sort that exceeds
//!   `30 * last` comparisons starts over with `fallbackSort`.
//!
//! Every routine below is the Java one, statement for statement, with its array layouts: `block`
//! holds byte `i` at `i + 1` and the last byte again at 0, and the quadrant array is `sfmap`.

const BASE_BLOCK_SIZE: i32 = 100_000;
const BLOCK_SIZE_100K: i32 = 9;
const MAX_ALPHA_SIZE: usize = 258;
const RUNA: u16 = 0;
const RUNB: u16 = 1;
const N_GROUPS: usize = 6;
const G_SIZE: i32 = 50;
const N_ITERS: usize = 4;
const MAX_SELECTORS: usize = 18002;
const NUM_OVERSHOOT_BYTES: i32 = 20;
const GREATER_ICOST: u8 = 15;
const LESSER_ICOST: u8 = 0;

const FALLBACK_QSORT_SMALL_THRESH: i32 = 10;
const SMALL_THRESH: i32 = 20;
const DEPTH_THRESH: i32 = 10;
const WORK_FACTOR: i32 = 30;
const SETMASK: i32 = 1 << 21;
const CLEARMASK: i32 = !SETMASK;
const QSORT_STACK_SIZE: usize = 1000;
const INCS: [i32; 14] = [
    1, 4, 13, 40, 121, 364, 1093, 3280, 9841, 29524, 88573, 265720, 797161, 2391484,
];

/// `BZip2CompressorOutputStream(out)` over `data`, then `finish()`: the whole stream.
pub fn compress(data: &[u8]) -> Vec<u8> {
    let mut stream = Stream::new();
    for &byte in data {
        stream.write0(i32::from(byte));
    }
    stream.finish()
}

/// bzip2's CRC-32: MSB first, polynomial 0x04c11db7, as `CRC.update(int, int)` applies it.
struct Crc {
    table: [u32; 256],
    value: u32,
}

impl Crc {
    fn new() -> Self {
        let mut table = [0u32; 256];
        for (i, entry) in table.iter_mut().enumerate() {
            let mut c = (i as u32) << 24;
            for _ in 0..8 {
                c = if c & 0x8000_0000 != 0 {
                    (c << 1) ^ 0x04c1_1db7
                } else {
                    c << 1
                };
            }
            *entry = c;
        }
        Self {
            table,
            value: 0xffff_ffff,
        }
    }

    fn update(&mut self, byte: i32, repeat: i32) {
        for _ in 0..repeat {
            let index = ((self.value >> 24) ^ byte as u32) & 0xff;
            self.value = (self.value << 8) ^ self.table[index as usize];
        }
    }

    fn get(&self) -> u32 {
        !self.value
    }
}

/// `BZip2CompressorOutputStream` with its `Data` and its `BlockSort`.
struct Stream {
    out: Vec<u8>,
    bs_buff: u32,
    bs_live: u32,
    last: i32,
    allowable_block_size: i32,
    crc: Crc,
    combined_crc: u32,
    current_char: i32,
    run_length: i32,
    n_in_use: usize,
    n_mtf: i32,
    // Data
    in_use: [bool; 256],
    unseq_to_seq: [u8; 256],
    mtf_freq: [i32; MAX_ALPHA_SIZE],
    selector: Vec<u8>,
    selector_mtf: Vec<u8>,
    len: [[u8; MAX_ALPHA_SIZE]; N_GROUPS],
    rfreq: [[i32; MAX_ALPHA_SIZE]; N_GROUPS],
    code: [[i32; MAX_ALPHA_SIZE]; N_GROUPS],
    block: Vec<u8>,
    fmap: Vec<i32>,
    /// `sfmap`, which `BlockSort` also uses as its `quadrant`.
    sfmap: Vec<u16>,
    orig_ptr: i32,
    // BlockSort
    work_done: i32,
    work_limit: i32,
    first_attempt: bool,
    ftab: Vec<i32>,
    eclass: Vec<i32>,
}

impl Stream {
    fn new() -> Self {
        let n = (BLOCK_SIZE_100K * BASE_BLOCK_SIZE) as usize;
        let mut stream = Self {
            out: Vec::new(),
            bs_buff: 0,
            bs_live: 0,
            last: -1,
            allowable_block_size: BLOCK_SIZE_100K * BASE_BLOCK_SIZE - 20,
            crc: Crc::new(),
            combined_crc: 0,
            current_char: -1,
            run_length: 0,
            n_in_use: 0,
            n_mtf: 0,
            in_use: [false; 256],
            unseq_to_seq: [0; 256],
            mtf_freq: [0; MAX_ALPHA_SIZE],
            selector: vec![0; MAX_SELECTORS],
            selector_mtf: vec![0; MAX_SELECTORS],
            len: [[0; MAX_ALPHA_SIZE]; N_GROUPS],
            rfreq: [[0; MAX_ALPHA_SIZE]; N_GROUPS],
            code: [[0; MAX_ALPHA_SIZE]; N_GROUPS],
            block: vec![0; n + 1 + NUM_OVERSHOOT_BYTES as usize],
            fmap: vec![0; n],
            sfmap: vec![0; 2 * n],
            orig_ptr: 0,
            work_done: 0,
            work_limit: 0,
            first_attempt: false,
            ftab: vec![0; 65537],
            eclass: Vec::new(),
        };
        stream.out.extend_from_slice(b"BZh");
        stream.out.push(b'0' + BLOCK_SIZE_100K as u8);
        stream.init_block();
        stream
    }

    fn init_block(&mut self) {
        self.crc = Crc::new();
        self.last = -1;
        self.in_use = [false; 256];
    }

    fn write0(&mut self, b: i32) {
        if self.current_char != -1 {
            let b = b & 0xff;
            if self.current_char == b {
                self.run_length += 1;
                if self.run_length > 254 {
                    self.write_run();
                    self.current_char = -1;
                    self.run_length = 0;
                }
            } else {
                self.write_run();
                self.run_length = 1;
                self.current_char = b;
            }
        } else {
            self.current_char = b & 0xff;
            self.run_length += 1;
        }
    }

    fn write_run(&mut self) {
        let last = self.last;
        if last < self.allowable_block_size {
            let current = self.current_char;
            self.in_use[current as usize] = true;
            let ch = current as u8;
            let mut run_length = self.run_length;
            self.crc.update(current, run_length);
            let at = (last + 2) as usize;
            match run_length {
                1 => {
                    self.block[at] = ch;
                    self.last = last + 1;
                }
                2 => {
                    self.block[at] = ch;
                    self.block[at + 1] = ch;
                    self.last = last + 2;
                }
                3 => {
                    self.block[at] = ch;
                    self.block[at + 1] = ch;
                    self.block[at + 2] = ch;
                    self.last = last + 3;
                }
                _ => {
                    run_length -= 4;
                    self.in_use[run_length as usize] = true;
                    self.block[at..at + 4].fill(ch);
                    self.block[at + 4] = run_length as u8;
                    self.last = last + 5;
                }
            }
        } else {
            self.end_block();
            self.init_block();
            self.write_run();
        }
    }

    fn finish(mut self) -> Vec<u8> {
        if self.run_length > 0 {
            self.write_run();
        }
        self.current_char = -1;
        self.end_block();
        // endCompression
        for byte in [0x17, 0x72, 0x45, 0x38, 0x50, 0x90] {
            self.bs_w(8, byte);
        }
        self.bs_put_int(self.combined_crc);
        while self.bs_live > 0 {
            self.out.push((self.bs_buff >> 24) as u8);
            self.bs_buff <<= 8;
            self.bs_live = self.bs_live.saturating_sub(8);
        }
        self.out
    }

    fn end_block(&mut self) {
        let block_crc = self.crc.get();
        self.combined_crc = self.combined_crc.rotate_left(1) ^ block_crc;
        if self.last == -1 {
            return;
        }
        self.block_sort();
        for byte in [0x31, 0x41, 0x59, 0x26, 0x53, 0x59] {
            self.bs_w(8, byte);
        }
        self.bs_put_int(block_crc);
        self.bs_w(1, 0);
        // moveToFrontCodeAndSend
        self.bs_w(24, self.orig_ptr as u32);
        self.generate_mtf_values();
        self.send_mtf_values();
    }

    /// `bsPutInt`: four bytes, most significant first.
    fn bs_put_int(&mut self, u: u32) {
        for shift in [24, 16, 8, 0] {
            self.bs_w(8, (u >> shift) & 0xff);
        }
    }

    /// `bsW(n, v)`: `n` bits of `v`, most significant first. `n` is at most 24.
    fn bs_w(&mut self, n: u32, v: u32) {
        while self.bs_live >= 8 {
            self.out.push((self.bs_buff >> 24) as u8);
            self.bs_buff <<= 8;
            self.bs_live -= 8;
        }
        self.bs_buff |= v << (32 - self.bs_live - n);
        self.bs_live += n;
    }

    fn generate_mtf_values(&mut self) {
        let mut n_in_use = 0usize;
        for i in 0..256 {
            if self.in_use[i] {
                self.unseq_to_seq[i] = n_in_use as u8;
                n_in_use += 1;
            }
        }
        self.n_in_use = n_in_use;
        let eob = n_in_use + 1;
        self.mtf_freq[..=eob].fill(0);
        let mut yy = [0u8; 256];
        for (i, y) in yy.iter_mut().enumerate().take(n_in_use) {
            *y = i as u8;
        }
        let mut wr = 0usize;
        let mut z_pend = 0i32;
        for i in 0..=self.last as usize {
            let ll_i = self.unseq_to_seq[self.block[self.fmap[i] as usize] as usize];
            let mut tmp = yy[0];
            let mut j = 0usize;
            while ll_i != tmp {
                j += 1;
                let tmp2 = tmp;
                tmp = yy[j];
                yy[j] = tmp2;
            }
            yy[0] = tmp;
            if j == 0 {
                z_pend += 1;
            } else {
                if z_pend > 0 {
                    self.flush_runs(&mut wr, z_pend - 1);
                    z_pend = 0;
                }
                self.sfmap[wr] = (j + 1) as u16;
                wr += 1;
                self.mtf_freq[j + 1] += 1;
            }
        }
        if z_pend > 0 {
            self.flush_runs(&mut wr, z_pend - 1);
        }
        self.sfmap[wr] = eob as u16;
        self.mtf_freq[eob] += 1;
        self.n_mtf = wr as i32 + 1;
    }

    /// A run of `z_pend + 1` zeros as RUNA and RUNB digits.
    fn flush_runs(&mut self, wr: &mut usize, mut z_pend: i32) {
        loop {
            let symbol = if z_pend & 1 == 0 { RUNA } else { RUNB };
            self.sfmap[*wr] = symbol;
            *wr += 1;
            self.mtf_freq[symbol as usize] += 1;
            if z_pend >= 2 {
                z_pend = (z_pend - 2) >> 1;
            } else {
                break;
            }
        }
    }

    fn send_mtf_values(&mut self) {
        let alpha_size = self.n_in_use + 2;
        for t in 0..N_GROUPS {
            self.len[t][..alpha_size].fill(GREATER_ICOST);
        }
        let n_groups = match self.n_mtf {
            n if n < 200 => 2,
            n if n < 600 => 3,
            n if n < 1200 => 4,
            n if n < 2400 => 5,
            _ => 6,
        };
        self.send_mtf_values0(n_groups, alpha_size);
        let n_selectors = self.send_mtf_values1(n_groups, alpha_size);
        self.send_mtf_values2(n_groups, n_selectors);
        self.send_mtf_values3(n_groups, alpha_size);
        self.send_mtf_values4();
        self.send_mtf_values5(n_groups, n_selectors);
        self.send_mtf_values6(n_groups, alpha_size);
        self.send_mtf_values7();
    }

    /// The initial tables: each covers a run of symbols holding its share of the frequencies.
    fn send_mtf_values0(&mut self, n_groups: usize, alpha_size: usize) {
        let mut rem_f = self.n_mtf;
        let mut gs = 0i32;
        let mut n_part = n_groups;
        while n_part > 0 {
            let t_freq = rem_f / n_part as i32;
            let mut ge = gs - 1;
            let mut a_freq = 0;
            let ax = alpha_size as i32 - 1;
            while a_freq < t_freq && ge < ax {
                ge += 1;
                a_freq += self.mtf_freq[ge as usize];
            }
            if ge > gs && n_part != n_groups && n_part != 1 && (n_groups - n_part) & 1 != 0 {
                a_freq -= self.mtf_freq[ge as usize];
                ge -= 1;
            }
            let len = &mut self.len[n_part - 1];
            for v in (0..alpha_size as i32).rev() {
                len[v as usize] = if v >= gs && v <= ge {
                    LESSER_ICOST
                } else {
                    GREATER_ICOST
                };
            }
            gs = ge + 1;
            rem_f -= a_freq;
            n_part -= 1;
        }
    }

    /// Four rounds of assigning each group of fifty symbols to its cheapest table and rebuilding
    /// the tables from what they were given.
    fn send_mtf_values1(&mut self, n_groups: usize, alpha_size: usize) -> usize {
        let n_mtf = self.n_mtf;
        let mut n_selectors = 0usize;
        for _ in 0..N_ITERS {
            let mut fave = [0i32; N_GROUPS];
            for t in (0..n_groups).rev() {
                fave[t] = 0;
                self.rfreq[t][..alpha_size].fill(0);
            }
            n_selectors = 0;
            let mut gs = 0i32;
            while gs < n_mtf {
                let ge = (gs + G_SIZE - 1).min(n_mtf - 1);
                let mut cost = [0i16; N_GROUPS];
                if n_groups == N_GROUPS {
                    for i in gs..=ge {
                        let icv = self.sfmap[i as usize] as usize;
                        for (t, c) in cost.iter_mut().enumerate() {
                            *c = c.wrapping_add(i16::from(self.len[t][icv] as i8));
                        }
                    }
                } else {
                    for i in gs..=ge {
                        let icv = self.sfmap[i as usize] as usize;
                        for t in (0..n_groups).rev() {
                            cost[t] = cost[t].wrapping_add(i16::from(self.len[t][icv] as i8));
                        }
                    }
                }
                let mut bt = -1i32;
                let mut bc = 999_999_999i32;
                for t in (0..n_groups).rev() {
                    let c = i32::from(cost[t]);
                    if c < bc {
                        bc = c;
                        bt = t as i32;
                    }
                }
                let bt = bt as usize;
                fave[bt] += 1;
                self.selector[n_selectors] = bt as u8;
                n_selectors += 1;
                for i in gs..=ge {
                    self.rfreq[bt][self.sfmap[i as usize] as usize] += 1;
                }
                gs = ge + 1;
            }
            for t in 0..n_groups {
                let freq = self.rfreq[t];
                hb_make_code_lengths(&mut self.len[t], &freq, alpha_size, 20);
            }
        }
        n_selectors
    }

    fn send_mtf_values2(&mut self, n_groups: usize, n_selectors: usize) {
        let mut pos = [0u8; N_GROUPS];
        for (i, p) in pos.iter_mut().enumerate().take(n_groups) {
            *p = i as u8;
        }
        for i in 0..n_selectors {
            let ll_i = self.selector[i];
            let mut tmp = pos[0];
            let mut j = 0usize;
            while ll_i != tmp {
                j += 1;
                let tmp2 = tmp;
                tmp = pos[j];
                pos[j] = tmp2;
            }
            pos[0] = tmp;
            self.selector_mtf[i] = j as u8;
        }
    }

    fn send_mtf_values3(&mut self, n_groups: usize, alpha_size: usize) {
        for t in 0..n_groups {
            let mut min_len = 32;
            let mut max_len = 0;
            for &l in &self.len[t][..alpha_size] {
                max_len = max_len.max(l);
                min_len = min_len.min(l);
            }
            let mut vec = 0i32;
            for n in min_len..=max_len {
                for i in 0..alpha_size {
                    if self.len[t][i] == n {
                        self.code[t][i] = vec;
                        vec += 1;
                    }
                }
                vec <<= 1;
            }
        }
    }

    fn send_mtf_values4(&mut self) {
        let mut in_use16 = [false; 16];
        for (i, used) in in_use16.iter_mut().enumerate() {
            *used = self.in_use[i * 16..i * 16 + 16].iter().any(|&u| u);
        }
        for &used in &in_use16 {
            self.bs_w(1, u32::from(used));
        }
        for (i, &used) in in_use16.iter().enumerate() {
            if used {
                for j in 0..16 {
                    self.bs_w(1, u32::from(self.in_use[i * 16 + j]));
                }
            }
        }
    }

    fn send_mtf_values5(&mut self, n_groups: usize, n_selectors: usize) {
        self.bs_w(3, n_groups as u32);
        self.bs_w(15, n_selectors as u32);
        for i in 0..n_selectors {
            for _ in 0..self.selector_mtf[i] {
                self.bs_w(1, 1);
            }
            self.bs_w(1, 0);
        }
    }

    fn send_mtf_values6(&mut self, n_groups: usize, alpha_size: usize) {
        for t in 0..n_groups {
            let mut curr = i32::from(self.len[t][0]);
            self.bs_w(5, curr as u32);
            for i in 0..alpha_size {
                let lti = i32::from(self.len[t][i]);
                while curr < lti {
                    self.bs_w(2, 2);
                    curr += 1;
                }
                while curr > lti {
                    self.bs_w(2, 3);
                    curr -= 1;
                }
                self.bs_w(1, 0);
            }
        }
    }

    fn send_mtf_values7(&mut self) {
        let n_mtf = self.n_mtf;
        let mut sel_ctr = 0usize;
        let mut gs = 0i32;
        while gs < n_mtf {
            let ge = (gs + G_SIZE - 1).min(n_mtf - 1);
            let table = self.selector[sel_ctr] as usize;
            while gs <= ge {
                let symbol = self.sfmap[gs as usize] as usize;
                let (len, code) = (self.len[table][symbol], self.code[table][symbol]);
                self.bs_w(u32::from(len), code as u32);
                gs += 1;
            }
            gs = ge + 1;
            sel_ctr += 1;
        }
    }

    // ---- BlockSort ----

    fn block_sort(&mut self) {
        let last = self.last;
        self.work_limit = WORK_FACTOR * last;
        self.work_done = 0;
        self.first_attempt = true;
        if last + 1 < 10000 {
            self.fallback_sort_block(last);
        } else {
            self.main_sort(last);
            if self.first_attempt && self.work_done > self.work_limit {
                self.fallback_sort_block(last);
            }
        }
        self.orig_ptr = -1;
        for i in 0..=last as usize {
            if self.fmap[i] == 0 {
                self.orig_ptr = i as i32;
                break;
            }
        }
    }

    /// `fallbackSort(Data, int)`: sorts the rotations of `block[0..=last]`, which starts with the
    /// last byte, then shifts the result back by one.
    fn fallback_sort_block(&mut self, last: i32) {
        self.block[0] = self.block[(last + 1) as usize];
        let nblock = (last + 1) as usize;
        if self.eclass.len() < nblock {
            self.eclass = vec![0; self.sfmap.len() / 2];
        }
        fallback_sort(&mut self.fmap, &self.block, &mut self.eclass, nblock);
        for f in &mut self.fmap[..nblock] {
            *f -= 1;
        }
        if let Some(f) = self.fmap[..nblock].iter_mut().find(|f| **f == -1) {
            *f = last;
        }
    }

    fn main_sort(&mut self, last: i32) {
        let mut running_order = [0i32; 256];
        let mut copy = [0i32; 256];
        let mut big_done = [false; 256];
        let work_limit = self.work_limit;
        let first_attempt = self.first_attempt;
        self.ftab.fill(0);
        for i in 0..NUM_OVERSHOOT_BYTES {
            self.block[(last + i + 2) as usize] = self.block[((i % (last + 1)) + 1) as usize];
        }
        for i in (0..(last + NUM_OVERSHOOT_BYTES + 1) as usize).rev() {
            self.sfmap[i] = 0;
        }
        self.block[0] = self.block[(last + 1) as usize];

        let mut c1 = i32::from(self.block[0]);
        for i in 0..=last as usize {
            let c2 = i32::from(self.block[i + 1]);
            self.ftab[((c1 << 8) + c2) as usize] += 1;
            c1 = c2;
        }
        for i in 1..=65536 {
            self.ftab[i] += self.ftab[i - 1];
        }
        c1 = i32::from(self.block[1]);
        for i in 0..last {
            let c2 = i32::from(self.block[(i + 2) as usize]);
            let slot = ((c1 << 8) + c2) as usize;
            self.ftab[slot] -= 1;
            self.fmap[self.ftab[slot] as usize] = i;
            c1 = c2;
        }
        let slot =
            ((i32::from(self.block[(last + 1) as usize]) << 8) + i32::from(self.block[1])) as usize;
        self.ftab[slot] -= 1;
        self.fmap[self.ftab[slot] as usize] = last;

        for i in (0..256).rev() {
            big_done[i] = false;
            running_order[i] = i as i32;
        }
        let ftab = &self.ftab;
        let size = |ro: i32| ftab[((ro + 1) << 8) as usize] - ftab[(ro << 8) as usize];
        let mut h = 364i32;
        while h != 1 {
            h /= 3;
            for i in h..=255 {
                let vv = running_order[i as usize];
                let a = size(vv);
                let b = h - 1;
                let mut j = i;
                let mut ro = running_order[(j - h) as usize];
                while size(ro) > a {
                    running_order[j as usize] = ro;
                    j -= h;
                    if j <= b {
                        break;
                    }
                    ro = running_order[(j - h) as usize];
                }
                running_order[j as usize] = vv;
            }
        }

        for i in 0..=255usize {
            let ss = running_order[i];
            // Step 1: sort the small buckets [ss, j] no earlier scan has finished.
            for j in 0..=255 {
                let sb = ((ss << 8) + j) as usize;
                let ftab_sb = self.ftab[sb];
                if ftab_sb & SETMASK != SETMASK {
                    let lo = ftab_sb & CLEARMASK;
                    let hi = (self.ftab[sb + 1] & CLEARMASK) - 1;
                    if hi > lo {
                        self.main_qsort3(lo, hi, 2, last);
                        if first_attempt && self.work_done > work_limit {
                            return;
                        }
                    }
                    self.ftab[sb] = ftab_sb | SETMASK;
                }
            }
            // Step 2: the sorted bucket orders the small buckets [t, ss].
            for j in 0..=255 {
                copy[j as usize] = self.ftab[((j << 8) + ss) as usize] & CLEARMASK;
            }
            let start = self.ftab[(ss << 8) as usize] & CLEARMASK;
            let end = self.ftab[((ss + 1) << 8) as usize] & CLEARMASK;
            for j in start..end {
                let fmap_j = self.fmap[j as usize];
                let c = self.block[fmap_j as usize] as usize;
                if !big_done[c] {
                    self.fmap[copy[c] as usize] = if fmap_j == 0 { last } else { fmap_j - 1 };
                    copy[c] += 1;
                }
            }
            for j in (0..256).rev() {
                self.ftab[((j << 8) + ss) as usize] |= SETMASK;
            }
            // Step 3: the quadrant descriptors of the finished bucket.
            big_done[ss as usize] = true;
            if i < 255 {
                let bb_start = self.ftab[(ss << 8) as usize] & CLEARMASK;
                let bb_size = (self.ftab[((ss + 1) << 8) as usize] & CLEARMASK) - bb_start;
                let mut shifts = 0;
                while (bb_size >> shifts) > 65534 {
                    shifts += 1;
                }
                for j in 0..bb_size {
                    let a2update = self.fmap[(bb_start + j) as usize];
                    let q = (j >> shifts) as u16;
                    self.sfmap[a2update as usize] = q;
                    if a2update < NUM_OVERSHOOT_BYTES {
                        self.sfmap[(a2update + last + 1) as usize] = q;
                    }
                }
            }
        }
    }

    fn main_qsort3(&mut self, lo_st: i32, hi_st: i32, d_st: i32, last: i32) {
        let mut stack_ll = [0i32; QSORT_STACK_SIZE];
        let mut stack_hh = [0i32; QSORT_STACK_SIZE];
        let mut stack_dd = [0i32; QSORT_STACK_SIZE];
        stack_ll[0] = lo_st;
        stack_hh[0] = hi_st;
        stack_dd[0] = d_st;
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let (lo, hi, d) = (stack_ll[sp], stack_hh[sp], stack_dd[sp]);
            if hi - lo < SMALL_THRESH || d > DEPTH_THRESH {
                if self.main_simple_sort(lo, hi, d, last) {
                    return;
                }
                continue;
            }
            let d1 = d + 1;
            let at = |fmap: &[i32], block: &[u8], x: i32| {
                i32::from(block[(fmap[x as usize] + d1) as usize])
            };
            let med = med3(
                at(&self.fmap, &self.block, lo),
                at(&self.fmap, &self.block, hi),
                at(&self.fmap, &self.block, ((lo + hi) as u32 >> 1) as i32),
            );
            let (mut un_lo, mut un_hi, mut lt_lo, mut gt_hi) = (lo, hi, lo, hi);
            loop {
                while un_lo <= un_hi {
                    let n = at(&self.fmap, &self.block, un_lo) - med;
                    if n == 0 {
                        self.fmap.swap(un_lo as usize, lt_lo as usize);
                        un_lo += 1;
                        lt_lo += 1;
                    } else if n < 0 {
                        un_lo += 1;
                    } else {
                        break;
                    }
                }
                while un_lo <= un_hi {
                    let n = at(&self.fmap, &self.block, un_hi) - med;
                    if n == 0 {
                        self.fmap.swap(un_hi as usize, gt_hi as usize);
                        un_hi -= 1;
                        gt_hi -= 1;
                    } else if n > 0 {
                        un_hi -= 1;
                    } else {
                        break;
                    }
                }
                if un_lo > un_hi {
                    break;
                }
                self.fmap.swap(un_lo as usize, un_hi as usize);
                un_lo += 1;
                un_hi -= 1;
            }
            if gt_hi < lt_lo {
                stack_ll[sp] = lo;
                stack_hh[sp] = hi;
                stack_dd[sp] = d1;
                sp += 1;
            } else {
                let n = (lt_lo - lo).min(un_lo - lt_lo);
                vswap(&mut self.fmap, lo, un_lo - n, n);
                let m = (hi - gt_hi).min(gt_hi - un_hi);
                vswap(&mut self.fmap, un_lo, hi - m + 1, m);
                let n = lo + un_lo - lt_lo - 1;
                let m = hi - (gt_hi - un_hi) + 1;
                for (l, h, dd) in [(lo, n, d), (n + 1, m - 1, d1), (m, hi, d)] {
                    stack_ll[sp] = l;
                    stack_hh[sp] = h;
                    stack_dd[sp] = dd;
                    sp += 1;
                }
            }
        }
    }

    /// Shell sort of `fmap[lo..=hi]` by the rotations from depth `d`, with `mainGtU` inline.
    /// Returns whether the work limit is spent.
    fn main_simple_sort(&mut self, lo: i32, hi: i32, d: i32, last: i32) -> bool {
        let big_n = hi - lo + 1;
        if big_n < 2 {
            return self.first_attempt && self.work_done > self.work_limit;
        }
        let mut hp = 0usize;
        while INCS[hp] < big_n {
            hp += 1;
        }
        let first_attempt = self.first_attempt;
        let work_limit = self.work_limit;
        let mut work_done = self.work_done;
        let fmap = &mut self.fmap;
        let block = &self.block;
        let quadrant = &self.sfmap;

        'hp: while hp > 0 {
            hp -= 1;
            let h = INCS[hp];
            let mj = lo + h - 1;
            let mut i = lo + h;
            while i <= hi {
                let mut k = 3;
                while i <= hi && {
                    k -= 1;
                    k >= 0
                } {
                    let v = fmap[i as usize];
                    let vd = v + d;
                    let mut j = i;
                    let mut once_run = false;
                    let mut a = 0;
                    loop {
                        if once_run {
                            fmap[j as usize] = a;
                            j -= h;
                            if j <= mj {
                                break;
                            }
                        } else {
                            once_run = true;
                        }
                        a = fmap[(j - h) as usize];
                        if !greater(a + d, vd, block, quadrant, last, &mut work_done) {
                            break;
                        }
                    }
                    fmap[j as usize] = v;
                    i += 1;
                }
                if first_attempt && i <= hi && work_done > work_limit {
                    break 'hp;
                }
            }
        }
        self.work_done = work_done;
        first_attempt && work_done > work_limit
    }
}

/// The inline `mainGtU` of `mainSimpleSort`: whether the rotation at `i1` sorts after the one at
/// `i2`, counting each four-byte step of the long comparison as one unit of work.
fn greater(
    mut i1: i32,
    mut i2: i32,
    block: &[u8],
    quadrant: &[u16],
    last: i32,
    work_done: &mut i32,
) -> bool {
    let b = |x: i32| block[x as usize];
    let q = |x: i32| quadrant[x as usize];
    for k in 1..=5 {
        if b(i1 + k) != b(i2 + k) {
            return b(i1 + k) > b(i2 + k);
        }
    }
    i1 += 6;
    i2 += 6;
    if b(i1) != b(i2) {
        return b(i1) > b(i2);
    }
    let last_plus1 = last + 1;
    let mut x = last;
    while x > 0 {
        x -= 4;
        for k in 0..4 {
            if b(i1 + k + 1) != b(i2 + k + 1) {
                return b(i1 + k + 1) > b(i2 + k + 1);
            }
            if q(i1 + k) != q(i2 + k) {
                return q(i1 + k) > q(i2 + k);
            }
        }
        i1 += 4;
        if i1 >= last_plus1 {
            i1 -= last_plus1;
        }
        i2 += 4;
        if i2 >= last_plus1 {
            i2 -= last_plus1;
        }
        *work_done += 1;
    }
    false
}

fn med3(a: i32, b: i32, c: i32) -> i32 {
    if a < b {
        if b < c {
            b
        } else if a < c {
            c
        } else {
            a
        }
    } else if b > c {
        b
    } else if a > c {
        c
    } else {
        a
    }
}

fn vswap(fmap: &mut [i32], mut p1: i32, mut p2: i32, n: i32) {
    let end = p1 + n;
    while p1 < end {
        fmap.swap(p1 as usize, p2 as usize);
        p1 += 1;
        p2 += 1;
    }
}

/// `fallbackSort(int[], byte[], int)`: libbzip2's prefix-doubling sort of the rotations of
/// `block[..nblock]`, with a `BitSet` for the bucket heads.
fn fallback_sort(fmap: &mut [i32], block: &[u8], eclass: &mut [i32], nblock: usize) {
    let mut ftab = [0i32; 257];
    eclass[..nblock].fill(0);
    for &byte in &block[..nblock] {
        ftab[byte as usize] += 1;
    }
    for i in 1..257 {
        ftab[i] += ftab[i - 1];
    }
    for (i, &byte) in block[..nblock].iter().enumerate() {
        let j = byte as usize;
        ftab[j] -= 1;
        fmap[ftab[j] as usize] = i as i32;
    }
    let mut bhtab = vec![false; nblock + 64 + 2];
    for &f in &ftab[..256] {
        bhtab[f as usize] = true;
    }
    for i in 0..32 {
        bhtab[nblock + 2 * i] = true;
        bhtab[nblock + 2 * i + 1] = false;
    }
    let next_clear =
        |bh: &[bool], from: usize| (from..).find(|&k| k >= bh.len() || !bh[k]).unwrap();
    let next_set = |bh: &[bool], from: usize| (from..bh.len()).find(|&k| bh[k]);
    let mut h = 1usize;
    loop {
        let mut j = 0i32;
        for i in 0..nblock {
            if bhtab[i] {
                j = i as i32;
            }
            let mut k = fmap[i] - h as i32;
            if k < 0 {
                k += nblock as i32;
            }
            eclass[k as usize] = j;
        }
        let mut n_not_done = 0usize;
        let mut r: i64 = -1;
        loop {
            let k = next_clear(&bhtab, (r + 1) as usize);
            let l = k as i64 - 1;
            if l >= nblock as i64 {
                break;
            }
            let Some(k) = next_set(&bhtab, k + 1) else {
                break;
            };
            r = k as i64 - 1;
            if r >= nblock as i64 {
                break;
            }
            if r > l {
                let (l, r) = (l as i32, r as i32);
                n_not_done += (r - l + 1) as usize;
                fallback_qsort3(fmap, eclass, l, r);
                let mut cc = -1;
                for i in l..=r {
                    let cc1 = eclass[fmap[i as usize] as usize];
                    if cc != cc1 {
                        bhtab[i as usize] = true;
                        cc = cc1;
                    }
                }
            }
        }
        h *= 2;
        if h > nblock || n_not_done == 0 {
            break;
        }
    }
}

fn fallback_qsort3(fmap: &mut [i32], eclass: &[i32], lo_st: i32, hi_st: i32) {
    let mut stack = vec![(lo_st, hi_st)];
    let mut r: i64 = 0;
    while let Some((lo, hi)) = stack.pop() {
        if hi - lo < FALLBACK_QSORT_SMALL_THRESH {
            fallback_simple_sort(fmap, eclass, lo, hi);
            continue;
        }
        r = ((r * 7621) + 1) % 32768;
        let ec = |x: i32| eclass[fmap[x as usize] as usize];
        let med = match r % 3 {
            0 => ec(lo),
            1 => ec(((lo + hi) as u32 >> 1) as i32),
            _ => ec(hi),
        };
        let (mut un_lo, mut lt_lo, mut un_hi, mut gt_hi) = (lo, lo, hi, hi);
        loop {
            while un_lo <= un_hi {
                let n = eclass[fmap[un_lo as usize] as usize] - med;
                if n == 0 {
                    fmap.swap(un_lo as usize, lt_lo as usize);
                    lt_lo += 1;
                    un_lo += 1;
                    continue;
                }
                if n > 0 {
                    break;
                }
                un_lo += 1;
            }
            while un_lo <= un_hi {
                let n = eclass[fmap[un_hi as usize] as usize] - med;
                if n == 0 {
                    fmap.swap(un_hi as usize, gt_hi as usize);
                    gt_hi -= 1;
                    un_hi -= 1;
                    continue;
                }
                if n < 0 {
                    break;
                }
                un_hi -= 1;
            }
            if un_lo > un_hi {
                break;
            }
            fmap.swap(un_lo as usize, un_hi as usize);
            un_lo += 1;
            un_hi -= 1;
        }
        if gt_hi < lt_lo {
            continue;
        }
        let n = (lt_lo - lo).min(un_lo - lt_lo);
        vswap(fmap, lo, un_lo - n, n);
        let m = (hi - gt_hi).min(gt_hi - un_hi);
        vswap(fmap, un_hi + 1, hi - m + 1, m);
        let n = lo + un_lo - lt_lo - 1;
        let m = hi - (gt_hi - un_hi) + 1;
        if n - lo > hi - m {
            stack.push((lo, n));
            stack.push((m, hi));
        } else {
            stack.push((m, hi));
            stack.push((lo, n));
        }
    }
}

fn fallback_simple_sort(fmap: &mut [i32], eclass: &[i32], lo: i32, hi: i32) {
    if lo == hi {
        return;
    }
    if hi - lo > 3 {
        let mut i = hi - 4;
        while i >= lo {
            let tmp = fmap[i as usize];
            let ec_tmp = eclass[tmp as usize];
            let mut j = i + 4;
            while j <= hi && ec_tmp > eclass[fmap[j as usize] as usize] {
                fmap[(j - 4) as usize] = fmap[j as usize];
                j += 4;
            }
            fmap[(j - 4) as usize] = tmp;
            i -= 1;
        }
    }
    let mut i = hi - 1;
    while i >= lo {
        let tmp = fmap[i as usize];
        let ec_tmp = eclass[tmp as usize];
        let mut j = i + 1;
        while j <= hi && ec_tmp > eclass[fmap[j as usize] as usize] {
            fmap[(j - 1) as usize] = fmap[j as usize];
            j += 1;
        }
        fmap[(j - 1) as usize] = tmp;
        i -= 1;
    }
}

/// `hbMakeCodeLengths`: Huffman code lengths for `freq`, rescaled until none exceeds `max_len`.
/// The rescaling pass stops before the last symbol, as commons-compress's does.
fn hb_make_code_lengths(len: &mut [u8], freq: &[i32], alpha_size: usize, max_len: i32) {
    let mut heap = [0i32; MAX_ALPHA_SIZE + 2];
    let mut weight = [0i32; MAX_ALPHA_SIZE * 2];
    let mut parent = [0i32; MAX_ALPHA_SIZE * 2];
    for i in (0..alpha_size).rev() {
        weight[i + 1] = (if freq[i] == 0 { 1 } else { freq[i] }) << 8;
    }
    let mut too_long = true;
    while too_long {
        too_long = false;
        let mut n_nodes = alpha_size as i32;
        let mut n_heap = 0usize;
        heap[0] = 0;
        weight[0] = 0;
        parent[0] = -2;
        for i in 1..=alpha_size {
            parent[i] = -1;
            n_heap += 1;
            heap[n_heap] = i as i32;
            let mut zz = n_heap;
            let tmp = heap[zz];
            while weight[tmp as usize] < weight[heap[zz >> 1] as usize] {
                heap[zz] = heap[zz >> 1];
                zz >>= 1;
            }
            heap[zz] = tmp;
        }
        while n_heap > 1 {
            let n1 = heap[1];
            heap[1] = heap[n_heap];
            n_heap -= 1;
            down_heap(&mut heap, &weight, n_heap);
            let n2 = heap[1];
            heap[1] = heap[n_heap];
            n_heap -= 1;
            down_heap(&mut heap, &weight, n_heap);
            n_nodes += 1;
            parent[n1 as usize] = n_nodes;
            parent[n2 as usize] = n_nodes;
            let (w1, w2) = (weight[n1 as usize], weight[n2 as usize]);
            weight[n_nodes as usize] =
                ((w1 & !0xff) + (w2 & !0xff)) | (1 + (w1 & 0xff).max(w2 & 0xff));
            parent[n_nodes as usize] = -1;
            n_heap += 1;
            heap[n_heap] = n_nodes;
            let mut zz = n_heap;
            let tmp = heap[zz];
            while weight[tmp as usize] < weight[heap[zz >> 1] as usize] {
                heap[zz] = heap[zz >> 1];
                zz >>= 1;
            }
            heap[zz] = tmp;
        }
        for i in 1..=alpha_size {
            let mut j = 0;
            let mut k = i as i32;
            while parent[k as usize] >= 0 {
                k = parent[k as usize];
                j += 1;
            }
            len[i - 1] = j as u8;
            if j > max_len {
                too_long = true;
            }
        }
        if too_long {
            for w in &mut weight[1..alpha_size] {
                let j = *w >> 8;
                *w = (1 + (j >> 1)) << 8;
            }
        }
    }
}

fn down_heap(heap: &mut [i32], weight: &[i32], n_heap: usize) {
    let mut zz = 1usize;
    let tmp = heap[zz];
    loop {
        let mut yy = zz << 1;
        if yy > n_heap {
            break;
        }
        if yy < n_heap && weight[heap[yy + 1] as usize] < weight[heap[yy] as usize] {
            yy += 1;
        }
        if weight[tmp as usize] < weight[heap[yy] as usize] {
            break;
        }
        heap[zz] = heap[yy];
        zz = yy;
    }
    heap[zz] = tmp;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// The inputs the expectations below were made from: SA-tag-like text on both sides of the
    /// 10000-byte switch to `mainSort` and past one 900k block, periodic data that exhausts
    /// `mainSort`'s work limit, every run-length edge, and random and low-entropy bytes.
    fn corpus() -> Vec<(&'static str, Vec<u8>)> {
        let mut seed = 0x5eed_u64;
        let mut next = move |bound: u64| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) % bound
        };
        let mut sa = |n: usize| -> Vec<u8> {
            let mut out = Vec::new();
            for _ in 0..n {
                let strand = if next(2) == 0 { '+' } else { '-' };
                out.extend(
                    format!(
                        "{},{},{strand},{}M{}S,{},{};",
                        1 + next(22),
                        1 + next(60_000_000),
                        20 + next(130),
                        next(100),
                        next(61),
                        next(10)
                    )
                    .bytes(),
                );
            }
            out.push(0);
            out
        };
        let sa_small = sa(5);
        let sa_mid = sa(400);
        let sa_switch = sa(2000);
        let sa_big = sa(80_000);
        let mut next = {
            let mut seed = 0xb2_u64;
            move |bound: u64| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (seed >> 33) % bound
            }
        };
        let random: Vec<u8> = (0..50_000).map(|_| next(256) as u8).collect();
        let low: Vec<u8> = (0..200_000)
            .map(|_| b"AAAAAAAAAAAACGTN"[next(16) as usize])
            .collect();
        vec![
            ("empty", vec![]),
            ("one", b"A".to_vec()),
            ("two", b"AB".to_vec()),
            ("run4", b"AAAA".to_vec()),
            ("run5", b"AAAAA".to_vec()),
            ("run255", vec![b'C'; 255]),
            ("run256", vec![b'C'; 256]),
            ("run600", vec![b'C'; 600]),
            ("sa_small", sa_small),
            ("sa_mid", sa_mid),
            ("sa_9999", sa_switch[..9999].to_vec()),
            ("sa_10000", sa_switch[..10000].to_vec()),
            ("sa_big", sa_big),
            ("period3", b"ABC".repeat(7000)),
            ("period7", b"ACGTTGA".repeat(3001)[..20000].to_vec()),
            ("zeros", vec![0; 300_000]),
            ("all256", (0..=255u8).collect::<Vec<_>>().repeat(40)),
            ("random", random),
            ("low", low),
        ]
    }

    /// commons-compress 1.26.2's stream for each input, as its length and CRC-32: the oracle
    /// image's `BZip2CompressorOutputStream` over the files `dump_corpus` writes.
    const JAVA: &[(&str, usize, u32)] = &[
        ("empty", 14, 0x039e0d3a),
        ("one", 37, 0x929ce026),
        ("two", 37, 0x446b2eaa),
        ("run4", 40, 0x03814163),
        ("run5", 40, 0x2ca563bb),
        ("run255", 40, 0x1edc4dd3),
        ("run256", 40, 0xa8ceb1e1),
        ("run600", 44, 0xb949b556),
        ("sa_small", 117, 0xa503e731),
        ("sa_mid", 4005, 0xce58ba7b),
        ("sa_9999", 3913, 0xc371d64b),
        ("sa_10000", 3919, 0x8b8c383d),
        ("sa_big", 625838, 0x854dc5c4),
        ("period3", 47, 0x68e3f1b5),
        ("period7", 58, 0x27d3e563),
        ("zeros", 47, 0xb9427245),
        ("all256", 753, 0xc340954f),
        ("random", 50579, 0xc0a25e77),
        ("low", 38513, 0xdc0a7aea),
    ];

    #[test]
    fn every_stream_is_the_one_commons_compress_writes() {
        for (name, data) in corpus() {
            let ours = compress(&data);
            let &(_, length, crc) = JAVA.iter().find(|e| e.0 == name).expect("an expectation");
            assert_eq!(ours.len(), length, "{name}");
            assert_eq!(crate::compression_header::crc32(&ours), crc, "{name}");
            let mut back = Vec::new();
            bzip2::read::BzDecoder::new(ours.as_slice())
                .read_to_end(&mut back)
                .expect("a stream libbzip2 reads");
            assert_eq!(back, data, "{name}");
        }
    }

    /// Writes the corpus to `$BZ_DIR`, for the Java side of [`JAVA`].
    #[test]
    #[ignore]
    fn dump_corpus() {
        let dir = std::env::var("BZ_DIR").expect("BZ_DIR");
        for (name, data) in corpus() {
            std::fs::write(format!("{dir}/{name}.in"), data).unwrap();
        }
    }
}
