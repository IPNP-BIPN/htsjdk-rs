//! The Huffman side of deflate: the tables, the two trees, and the bit writer.
//!
//! Ported from zlib's `trees.c` (zlib licence, Jean-loup Gailly and Mark Adler). The structure is
//! kept rather than improved, because the output has to match bit for bit and every choice here is
//! observable: which tree wins a block, how ties inside the heap are broken, the order the bit
//! length codes are sent in.
//!
//! The static tables are **generated at load time by the same algorithm zlib uses**, not copied
//! from `trees.h`. Copying the numbers would hide the rule that produced them, and the rule is
//! short: see [`Tables::new`].

/// Deflate's alphabet sizes, named as zlib names them.
pub const LITERALS: usize = 256;
pub const LENGTH_CODES: usize = 29;
pub const L_CODES: usize = LITERALS + 1 + LENGTH_CODES;
pub const D_CODES: usize = 30;
pub const BL_CODES: usize = 19;
pub const HEAP_SIZE: usize = 2 * L_CODES + 1;
pub const MAX_BITS: usize = 15;
pub const END_BLOCK: usize = 256;

const STORED_BLOCK: u32 = 0;
const STATIC_TREES: u32 = 1;
const DYN_TREES: u32 = 2;

const EXTRA_LBITS: [u8; LENGTH_CODES] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const EXTRA_DBITS: [u8; D_CODES] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];
const EXTRA_BLBITS: [u8; BL_CODES] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 7];

/// The bit length codes are sent in order of decreasing probability, so that the trailing unused
/// ones can be dropped from the header.
const BL_ORDER: [u8; BL_CODES] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

/// Repeat codes, as they appear in the bit length alphabet.
const REP_3_6: usize = 16;
const REPZ_3_10: usize = 17;
const REPZ_11_138: usize = 18;

/// One node: a frequency (or, once the tree is built, a code) and a length (or a parent).
///
/// zlib overlays these two meanings in a union. Keeping both fields costs a few kilobytes and
/// removes a class of bug that a union invites, so they are kept apart.
#[derive(Clone, Copy, Default)]
pub struct Node {
    pub freq: u16,
    pub code: u16,
    pub len: u16,
    pub dad: u16,
}

/// The tables `tr_static_init` builds, and the static trees a block may be sent with.
pub struct Tables {
    /// `length - 3` (0..=255) to its length code (0..=28).
    pub length_code: [u8; 256],
    /// A distance to its distance code, the first 256 directly and the rest by `dist >> 7`.
    pub dist_code: [u8; 512],
    pub base_length: [u16; LENGTH_CODES],
    pub base_dist: [u16; D_CODES],
    pub static_ltree: [Node; L_CODES + 2],
    pub static_dtree: [Node; D_CODES],
}

impl Tables {
    /// Ported from zlib's `tr_static_init`.
    fn new() -> Self {
        let mut tables = Tables {
            length_code: [0; 256],
            dist_code: [0; 512],
            base_length: [0; LENGTH_CODES],
            base_dist: [0; D_CODES],
            static_ltree: [Node::default(); L_CODES + 2],
            static_dtree: [Node::default(); D_CODES],
        };

        // length (0..255, meaning a match of 3..258) -> length code (0..28).
        let mut length = 0usize;
        for (code, &extra) in EXTRA_LBITS.iter().enumerate().take(LENGTH_CODES - 1) {
            tables.base_length[code] = length as u16;
            for _ in 0..(1u32 << extra) {
                tables.length_code[length] = code as u8;
                length += 1;
            }
        }
        // Length 255 (a match of 258) has two encodings, code 284 plus five extra bits or code
        // 285. The last write wins and picks the shorter one.
        tables.length_code[length - 1] = (LENGTH_CODES - 1) as u8;

        // distance (0..32K) -> distance code (0..29), the top half indexed by dist >> 7.
        let mut dist = 0usize;
        for (code, &extra) in EXTRA_DBITS.iter().enumerate().take(16) {
            tables.base_dist[code] = dist as u16;
            for _ in 0..(1u32 << extra) {
                tables.dist_code[dist] = code as u8;
                dist += 1;
            }
        }
        dist >>= 7;
        for (code, &extra) in EXTRA_DBITS.iter().enumerate().skip(16) {
            tables.base_dist[code] = (dist << 7) as u16;
            for _ in 0..(1u32 << (extra - 7)) {
                tables.dist_code[256 + dist] = code as u8;
                dist += 1;
            }
        }

        // The static literal tree: four runs of fixed lengths, defined by the deflate spec.
        let mut bl_count = [0u16; MAX_BITS + 1];
        for n in 0..=143 {
            tables.static_ltree[n].len = 8;
            bl_count[8] += 1;
        }
        for n in 144..=255 {
            tables.static_ltree[n].len = 9;
            bl_count[9] += 1;
        }
        for n in 256..=279 {
            tables.static_ltree[n].len = 7;
            bl_count[7] += 1;
        }
        for n in 280..=287 {
            tables.static_ltree[n].len = 8;
            bl_count[8] += 1;
        }
        // The static tree is built with a length of 9 rather than MAX_BITS, so that the codes come
        // out in the order the spec fixes rather than in canonical order.
        gen_codes(&mut tables.static_ltree, L_CODES + 1, &mut bl_count);

        for n in 0..D_CODES {
            tables.static_dtree[n].len = 5;
            tables.static_dtree[n].code = bit_reverse(n as u32, 5) as u16;
        }
        tables
    }

    #[inline]
    fn d_code(&self, dist: usize) -> usize {
        if dist < 256 {
            self.dist_code[dist] as usize
        } else {
            self.dist_code[256 + (dist >> 7)] as usize
        }
    }
}

/// Built once. `tr_static_init` in zlib is guarded by a flag and does the same thing.
pub static TABLES: std::sync::LazyLock<Tables> = std::sync::LazyLock::new(Tables::new);

/// Reverse the low `len` bits of `code`. Huffman codes travel most-significant-bit first, and
/// deflate writes bits least-significant first, so every code is stored reversed.
fn bit_reverse(mut code: u32, len: u32) -> u32 {
    let mut result = 0u32;
    for _ in 0..len {
        result |= code & 1;
        code >>= 1;
        result <<= 1;
    }
    result >> 1
}

/// Ported from zlib's `gen_codes`: turn a table of code lengths into canonical codes.
fn gen_codes(tree: &mut [Node], max_code: usize, bl_count: &mut [u16; MAX_BITS + 1]) {
    let mut next_code = [0u16; MAX_BITS + 1];
    let mut code = 0u16;
    for bits in 1..=MAX_BITS {
        code = (code + bl_count[bits - 1]) << 1;
        next_code[bits] = code;
    }
    for node in tree[..=max_code].iter_mut() {
        let len = node.len as usize;
        if len == 0 {
            continue;
        }
        node.code = bit_reverse(next_code[len] as u32, len as u32) as u16;
        next_code[len] += 1;
    }
}

/// Which alphabet a tree is for. The three differ in their extra-bit tables and in how many
/// elements they have, and `build_tree` needs all three facts.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Literal,
    Distance,
    BitLength,
}

impl Kind {
    fn extra_bits(self, n: usize) -> u32 {
        match self {
            // The literal alphabet's extra bits start at the first length code.
            Kind::Literal => {
                if n > LITERALS {
                    EXTRA_LBITS[n - LITERALS - 1] as u32
                } else {
                    0
                }
            }
            Kind::Distance => EXTRA_DBITS[n] as u32,
            Kind::BitLength => EXTRA_BLBITS[n] as u32,
        }
    }

    fn extra_base(self) -> usize {
        match self {
            Kind::Literal => LITERALS + 1,
            _ => 0,
        }
    }

    fn max_length(self) -> usize {
        match self {
            Kind::BitLength => 7,
            _ => MAX_BITS,
        }
    }
}

/// Everything a block needs while it is being accumulated.
pub struct Trees {
    pub dyn_ltree: [Node; HEAP_SIZE],
    pub dyn_dtree: [Node; 2 * D_CODES + 1],
    pub bl_tree: [Node; 2 * BL_CODES + 1],

    l_max_code: usize,
    d_max_code: usize,
    bl_max_code: usize,

    heap: [usize; HEAP_SIZE],
    heap_len: usize,
    heap_max: usize,
    depth: [u8; HEAP_SIZE],

    /// zlib holds these in `ulg` and lets them wrap: `build_tree` decrements them before
    /// `gen_bitlen` adds anything, so on an empty tree they briefly go below zero. Signed here,
    /// which reaches the same answer without relying on wrapping.
    opt_len: i64,
    static_len: i64,

    /// The symbols of the current block, three bytes each: distance low, distance high, then
    /// either the literal or `length - 3`. A distance of zero marks a literal.
    pub sym_buf: Vec<u8>,
    pub sym_next: usize,
    pub sym_end: usize,

    bi_buf: u32,
    /// Signed, because zlib's `bi_valid += length - Buf_size` is a negative step in the flushing
    /// branch. Unsigned here would underflow where the C simply subtracts.
    bi_valid: i32,

    /// zlib keeps this inside the tree descriptor. Here it is one field, written by `gen_bitlen`
    /// and read by `gen_codes` immediately afterwards.
    bl_count: [u16; MAX_BITS + 1],
}

impl Trees {
    pub fn new(lit_bufsize: usize) -> Self {
        let mut trees = Trees {
            dyn_ltree: [Node::default(); HEAP_SIZE],
            dyn_dtree: [Node::default(); 2 * D_CODES + 1],
            bl_tree: [Node::default(); 2 * BL_CODES + 1],
            l_max_code: 0,
            d_max_code: 0,
            bl_max_code: 0,
            heap: [0; HEAP_SIZE],
            heap_len: 0,
            heap_max: 0,
            depth: [0; HEAP_SIZE],
            opt_len: 0,
            static_len: 0,
            sym_buf: vec![0; lit_bufsize * 3],
            sym_next: 0,
            sym_end: (lit_bufsize - 1) * 3,
            bi_buf: 0,
            bi_valid: 0,
            bl_count: [0; MAX_BITS + 1],
        };
        trees.init_block();
        trees
    }

    fn init_block(&mut self) {
        for node in self.dyn_ltree.iter_mut() {
            node.freq = 0;
        }
        for node in self.dyn_dtree.iter_mut() {
            node.freq = 0;
        }
        for node in self.bl_tree.iter_mut() {
            node.freq = 0;
        }
        self.dyn_ltree[END_BLOCK].freq = 1;
        self.opt_len = 0;
        self.static_len = 0;
        self.sym_next = 0;
    }

    /// Record a literal. Returns true when the symbol buffer is full and the block must be closed.
    pub fn tally_lit(&mut self, c: u8) -> bool {
        self.sym_buf[self.sym_next] = 0;
        self.sym_buf[self.sym_next + 1] = 0;
        self.sym_buf[self.sym_next + 2] = c;
        self.sym_next += 3;
        self.dyn_ltree[c as usize].freq += 1;
        self.sym_next == self.sym_end
    }

    /// Record a match of `len` bytes at `dist` back. `len` arrives already reduced by MIN_MATCH,
    /// as zlib's `_tr_tally_dist` expects.
    pub fn tally_dist(&mut self, dist: usize, len: usize) -> bool {
        self.sym_buf[self.sym_next] = dist as u8;
        self.sym_buf[self.sym_next + 1] = (dist >> 8) as u8;
        self.sym_buf[self.sym_next + 2] = len as u8;
        self.sym_next += 3;
        let dist = dist - 1;
        self.dyn_ltree[TABLES.length_code[len] as usize + LITERALS + 1].freq += 1;
        self.dyn_dtree[TABLES.d_code(dist)].freq += 1;
        self.sym_next == self.sym_end
    }

    /// Ported from zlib's `pqdownheap`. Ties are broken by `depth`, which is what keeps two trees
    /// with the same frequencies from coming out differently.
    fn pqdownheap(&mut self, which: Which, mut k: usize) {
        let v = self.heap[k];
        let mut j = k << 1;
        while j <= self.heap_len {
            if j < self.heap_len && self.smaller(which, self.heap[j + 1], self.heap[j]) {
                j += 1;
            }
            if self.smaller(which, v, self.heap[j]) {
                break;
            }
            self.heap[k] = self.heap[j];
            k = j;
            j <<= 1;
        }
        self.heap[k] = v;
    }

    fn smaller(&self, which: Which, n: usize, m: usize) -> bool {
        let tree = self.tree(which);
        tree[n].freq < tree[m].freq
            || (tree[n].freq == tree[m].freq && self.depth[n] <= self.depth[m])
    }

    fn tree(&self, which: Which) -> &[Node] {
        match which {
            Which::Literal => &self.dyn_ltree,
            Which::Distance => &self.dyn_dtree,
            Which::BitLength => &self.bl_tree,
        }
    }

    fn tree_mut(&mut self, which: Which) -> &mut [Node] {
        match which {
            Which::Literal => &mut self.dyn_ltree,
            Which::Distance => &mut self.dyn_dtree,
            Which::BitLength => &mut self.bl_tree,
        }
    }

    /// Ported from zlib's `build_tree`: build the Huffman tree, then its code lengths, then its
    /// codes, and leave `opt_len` and `static_len` updated so the caller can pick a block type.
    fn build_tree(&mut self, which: Which) {
        let kind = which.kind();
        let elems = which.elems();
        let mut max_code: isize = -1;

        self.heap_len = 0;
        self.heap_max = HEAP_SIZE;

        for n in 0..elems {
            if self.tree(which)[n].freq != 0 {
                self.heap_len += 1;
                self.heap[self.heap_len] = n;
                max_code = n as isize;
                self.depth[n] = 0;
            } else {
                self.tree_mut(which)[n].len = 0;
            }
        }

        // A tree needs at least two codes. If the block used fewer, invent them: the extra code
        // costs a bit in the header and saves the encoder from a special case.
        while self.heap_len < 2 {
            self.heap_len += 1;
            let node = if max_code < 2 {
                max_code += 1;
                max_code as usize
            } else {
                0
            };
            self.heap[self.heap_len] = node;
            self.tree_mut(which)[node].freq = 1;
            self.depth[node] = 0;
            self.opt_len -= 1;
            if let Some(static_tree) = which.static_tree() {
                self.static_len -= static_tree[node].len as i64;
            }
        }
        let max_code = max_code as usize;
        which.set_max_code(self, max_code);

        for n in (1..=self.heap_len / 2).rev() {
            self.pqdownheap(which, n);
        }

        // Repeatedly join the two rarest nodes.
        let mut node = elems;
        loop {
            let n = self.heap[1];
            self.heap[1] = self.heap[self.heap_len];
            self.heap_len -= 1;
            self.pqdownheap(which, 1);
            let m = self.heap[1];

            self.heap_max -= 1;
            self.heap[self.heap_max] = n;
            self.heap_max -= 1;
            self.heap[self.heap_max] = m;

            let freq = self.tree(which)[n].freq + self.tree(which)[m].freq;
            self.tree_mut(which)[node].freq = freq;
            self.depth[node] = self.depth[n].max(self.depth[m]) + 1;
            self.tree_mut(which)[n].dad = node as u16;
            self.tree_mut(which)[m].dad = node as u16;
            self.heap[1] = node;
            node += 1;
            self.pqdownheap(which, 1);
            if self.heap_len < 2 {
                break;
            }
        }
        self.heap_max -= 1;
        self.heap[self.heap_max] = self.heap[1];

        self.gen_bitlen(which, kind, max_code);
        let mut bl_count = self.bl_count;
        let tree = self.tree_mut(which);
        gen_codes(tree, max_code, &mut bl_count);
    }

    /// Ported from zlib's `gen_bitlen`. The overflow loop at the end is the part worth reading: a
    /// Huffman tree can be deeper than `max_length`, and deflate cannot express that, so the tree
    /// is flattened by moving leaves up one level at a time.
    fn gen_bitlen(&mut self, which: Which, kind: Kind, max_code: usize) {
        let mut bl_count = [0u16; MAX_BITS + 1];
        let max_length = kind.max_length();
        let extra_base = kind.extra_base();
        let mut overflow: isize = 0;

        let root = self.heap[self.heap_max];
        self.tree_mut(which)[root].len = 0;

        for h in self.heap_max + 1..HEAP_SIZE {
            let n = self.heap[h];
            let dad = self.tree(which)[n].dad as usize;
            let mut bits = self.tree(which)[dad].len as usize + 1;
            if bits > max_length {
                bits = max_length;
                overflow += 1;
            }
            self.tree_mut(which)[n].len = bits as u16;
            if n > max_code {
                continue;
            }
            bl_count[bits] += 1;
            let extra = if n >= extra_base {
                kind.extra_bits(n)
            } else {
                0
            } as usize;
            let freq = self.tree(which)[n].freq as usize;
            self.opt_len += (freq * (bits + extra)) as i64;
            if let Some(static_tree) = which.static_tree() {
                self.static_len += (freq * (static_tree[n].len as usize + extra)) as i64;
            }
        }

        if overflow != 0 {
            loop {
                let mut bits = max_length - 1;
                while bl_count[bits] == 0 {
                    bits -= 1;
                }
                bl_count[bits] -= 1;
                bl_count[bits + 1] += 2;
                bl_count[max_length] -= 1;
                overflow -= 2;
                if overflow <= 0 {
                    break;
                }
            }

            // Re-issue the lengths in decreasing order, only touching nodes that are too deep.
            let mut h = HEAP_SIZE;
            for bits in (1..=max_length).rev() {
                let mut n = bl_count[bits];
                while n != 0 {
                    h -= 1;
                    let m = self.heap[h];
                    if m > max_code {
                        continue;
                    }
                    if self.tree(which)[m].len as usize != bits {
                        self.opt_len += ((bits - self.tree(which)[m].len as usize)
                            * self.tree(which)[m].freq as usize)
                            as i64;
                        self.tree_mut(which)[m].len = bits as u16;
                    }
                    n -= 1;
                }
            }
        }
        self.bl_count = bl_count;
    }

    /// Ported from zlib's `scan_tree`: count the runs the bit length alphabet will encode.
    fn scan_tree(&mut self, which: Which, max_code: usize) {
        let mut tree: Vec<Node> = self.tree(which)[..max_code + 2].to_vec();
        // zlib's guard: a length no real code can have, so the last run always closes.
        tree[max_code + 1].len = 0xffff;
        let mut prevlen: isize = -1;
        let mut nextlen = tree[0].len as usize;
        let mut count = 0;
        let mut max_count = 7;
        let mut min_count = 4;
        if nextlen == 0 {
            max_count = 138;
            min_count = 3;
        }
        // A sentinel past the end, so the last run is closed like any other.
        let mut curlen;

        for n in 0..=max_code {
            curlen = nextlen;
            nextlen = tree[n + 1].len as usize;
            count += 1;
            if count < max_count && curlen == nextlen {
                continue;
            } else if count < min_count {
                self.bl_tree[curlen].freq += count as u16;
            } else if curlen != 0 {
                if curlen as isize != prevlen {
                    self.bl_tree[curlen].freq += 1;
                }
                self.bl_tree[REP_3_6].freq += 1;
            } else if count <= 10 {
                self.bl_tree[REPZ_3_10].freq += 1;
            } else {
                self.bl_tree[REPZ_11_138].freq += 1;
            }
            count = 0;
            prevlen = curlen as isize;
            if nextlen == 0 {
                max_count = 138;
                min_count = 3;
            } else if curlen == nextlen {
                max_count = 6;
                min_count = 3;
            } else {
                max_count = 7;
                min_count = 4;
            }
        }
    }

    /// Ported from zlib's `send_tree`, the emitting twin of `scan_tree`. The two must walk the
    /// tree identically or the header describes a different tree from the one that follows.
    fn send_tree(&mut self, out: &mut Vec<u8>, which: Which, max_code: usize) {
        let mut tree: Vec<Node> = self.tree(which)[..max_code + 2].to_vec();
        tree[max_code + 1].len = 0xffff;
        let bl_tree = self.bl_tree;
        let mut prevlen: isize = -1;
        let mut nextlen = tree[0].len as usize;
        let mut count = 0;
        let mut max_count = 7;
        let mut min_count = 4;
        if nextlen == 0 {
            max_count = 138;
            min_count = 3;
        }

        for n in 0..=max_code {
            let curlen = nextlen;
            nextlen = tree[n + 1].len as usize;
            count += 1;
            if count < max_count && curlen == nextlen {
                continue;
            } else if count < min_count {
                for _ in 0..count {
                    self.send_code(out, curlen, &bl_tree);
                }
            } else if curlen != 0 {
                if curlen as isize != prevlen {
                    self.send_code(out, curlen, &bl_tree);
                    count -= 1;
                }
                self.send_code(out, REP_3_6, &bl_tree);
                self.send_bits(out, count as u32 - 3, 2);
            } else if count <= 10 {
                self.send_code(out, REPZ_3_10, &bl_tree);
                self.send_bits(out, count as u32 - 3, 3);
            } else {
                self.send_code(out, REPZ_11_138, &bl_tree);
                self.send_bits(out, count as u32 - 11, 7);
            }
            count = 0;
            prevlen = curlen as isize;
            if nextlen == 0 {
                max_count = 138;
                min_count = 3;
            } else if curlen == nextlen {
                max_count = 6;
                min_count = 3;
            } else {
                max_count = 7;
                min_count = 4;
            }
        }
    }

    /// Ported from zlib's `build_bl_tree`. Returns the index of the last bit length code that has
    /// to be sent; the ones after it in `BL_ORDER` are all zero and are dropped.
    fn build_bl_tree(&mut self) -> usize {
        let (l_max, d_max) = (self.l_max_code, self.d_max_code);
        self.scan_tree(Which::Literal, l_max);
        self.scan_tree(Which::Distance, d_max);
        self.build_tree(Which::BitLength);

        let mut max_blindex = BL_CODES - 1;
        loop {
            if self.bl_tree[BL_ORDER[max_blindex] as usize].len != 0 || max_blindex < 3 {
                break;
            }
            max_blindex -= 1;
        }
        self.opt_len += (3 * (max_blindex + 1) + 5 + 5 + 4) as i64;
        max_blindex
    }

    fn send_all_trees(&mut self, out: &mut Vec<u8>, lcodes: usize, dcodes: usize, blcodes: usize) {
        self.send_bits(out, lcodes as u32 - 257, 5);
        self.send_bits(out, dcodes as u32 - 1, 5);
        self.send_bits(out, blcodes as u32 - 4, 4);
        for &order in BL_ORDER.iter().take(blcodes) {
            let len = self.bl_tree[order as usize].len as u32;
            self.send_bits(out, len, 3);
        }
        let (l_max, d_max) = (self.l_max_code, self.d_max_code);
        self.send_tree(out, Which::Literal, l_max);
        self.send_tree(out, Which::Distance, d_max);
    }

    #[inline]
    fn send_code(&mut self, out: &mut Vec<u8>, c: usize, tree: &[Node]) {
        self.send_bits(out, tree[c].code as u32, tree[c].len as u32);
    }

    /// Ported from zlib's `send_bits`. Bits accumulate low-end first in a 16-bit window, and are
    /// flushed a whole `u16` at a time, little-endian.
    #[inline]
    fn send_bits(&mut self, out: &mut Vec<u8>, value: u32, length: u32) {
        const BUF_SIZE: i32 = 16;
        let length = length as i32;
        if self.bi_valid > BUF_SIZE - length {
            self.bi_buf |= (value << self.bi_valid) & 0xffff;
            out.push(self.bi_buf as u8);
            out.push((self.bi_buf >> 8) as u8);
            self.bi_buf = (value & 0xffff) >> (BUF_SIZE - self.bi_valid);
            self.bi_valid += length - BUF_SIZE;
        } else {
            self.bi_buf |= (value << self.bi_valid) & 0xffff;
            self.bi_valid += length;
        }
    }

    /// Flush whatever bits remain and align to a byte boundary.
    fn bi_windup(&mut self, out: &mut Vec<u8>) {
        if self.bi_valid > 8 {
            out.push(self.bi_buf as u8);
            out.push((self.bi_buf >> 8) as u8);
        } else if self.bi_valid > 0 {
            out.push(self.bi_buf as u8);
        }
        self.bi_buf = 0;
        self.bi_valid = 0;
    }

    /// Ported from zlib's `compress_block`: replay the symbol buffer through a pair of trees.
    fn compress_block(&mut self, out: &mut Vec<u8>, ltree: &[Node], dtree: &[Node]) {
        if self.sym_next != 0 {
            let mut sx = 0;
            while sx < self.sym_next {
                let mut dist = self.sym_buf[sx] as usize;
                dist += (self.sym_buf[sx + 1] as usize) << 8;
                let lc = self.sym_buf[sx + 2] as usize;
                sx += 3;
                if dist == 0 {
                    self.send_code(out, lc, ltree);
                } else {
                    let code = TABLES.length_code[lc] as usize;
                    self.send_code(out, code + LITERALS + 1, ltree);
                    let extra = EXTRA_LBITS[code] as u32;
                    if extra != 0 {
                        self.send_bits(out, (lc - TABLES.base_length[code] as usize) as u32, extra);
                    }
                    let dist = dist - 1;
                    let code = TABLES.d_code(dist);
                    self.send_code(out, code, dtree);
                    let extra = EXTRA_DBITS[code] as u32;
                    if extra != 0 {
                        self.send_bits(out, (dist - TABLES.base_dist[code] as usize) as u32, extra);
                    }
                }
            }
        }
        self.send_code(out, END_BLOCK, ltree);
    }

    /// Ported from zlib's `_tr_stored_block`.
    fn stored_block(&mut self, out: &mut Vec<u8>, buf: &[u8], last: bool) {
        self.send_bits(out, (STORED_BLOCK << 1) + last as u32, 3);
        self.bi_windup(out);
        let len = buf.len() as u16;
        out.push(len as u8);
        out.push((len >> 8) as u8);
        out.push(!len as u8);
        out.push(((!len) >> 8) as u8);
        out.extend_from_slice(buf);
    }

    /// Ported from zlib's `_tr_flush_block`: build both trees, then pick the cheapest of stored,
    /// static and dynamic. The three-way comparison is where most of a deflate stream's identity
    /// comes from, and it is decided in whole bytes, so a tie goes to the static tree.
    pub fn flush_block(
        &mut self,
        out: &mut Vec<u8>,
        buf: Option<&[u8]>,
        stored_len: usize,
        last: bool,
    ) {
        self.build_tree(Which::Literal);
        self.build_tree(Which::Distance);
        let max_blindex = self.build_bl_tree();
        let opt = ((self.opt_len + 3 + 7) >> 3) as usize;
        let static_lenb = ((self.static_len + 3 + 7) >> 3) as usize;
        let opt_lenb = if static_lenb <= opt { static_lenb } else { opt };

        if let (Some(buf), true) = (buf, stored_len + 4 <= opt_lenb) {
            self.stored_block(out, buf, last);
        } else if static_lenb == opt_lenb {
            self.send_bits(out, (STATIC_TREES << 1) + last as u32, 3);
            let (ltree, dtree) = (TABLES.static_ltree, TABLES.static_dtree);
            self.compress_block(out, &ltree, &dtree);
        } else {
            self.send_bits(out, (DYN_TREES << 1) + last as u32, 3);
            let (l, d, b) = (self.l_max_code + 1, self.d_max_code + 1, max_blindex + 1);
            self.send_all_trees(out, l, d, b);
            let (ltree, dtree) = (self.dyn_ltree, self.dyn_dtree);
            self.compress_block(out, &ltree, &dtree);
        }
        self.init_block();
        if last {
            self.bi_windup(out);
        }
    }
}

/// Which of the three trees a call is about.
#[derive(Clone, Copy, PartialEq)]
enum Which {
    Literal,
    Distance,
    BitLength,
}

impl Which {
    fn kind(self) -> Kind {
        match self {
            Which::Literal => Kind::Literal,
            Which::Distance => Kind::Distance,
            Which::BitLength => Kind::BitLength,
        }
    }

    fn elems(self) -> usize {
        match self {
            Which::Literal => L_CODES,
            Which::Distance => D_CODES,
            Which::BitLength => BL_CODES,
        }
    }

    fn static_tree(self) -> Option<&'static [Node]> {
        match self {
            Which::Literal => Some(&TABLES.static_ltree),
            Which::Distance => Some(&TABLES.static_dtree),
            Which::BitLength => None,
        }
    }

    fn set_max_code(self, trees: &mut Trees, max_code: usize) {
        match self {
            Which::Literal => trees.l_max_code = max_code,
            Which::Distance => trees.d_max_code = max_code,
            Which::BitLength => trees.bl_max_code = max_code,
        }
    }
}
