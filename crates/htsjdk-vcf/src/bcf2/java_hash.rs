//! The iteration order of a `java.util.HashMap`, which is the order BCF2 writes its INFO fields.
//!
//! `BCF2Writer.buildInfo` walks `vc.getAttributes().entrySet()` and writes each field in turn, so
//! the INFO block of a BCF record is in **hash order**, not in the file's order and not sorted.
//! Every map an attribute passes through on its way there is a `HashMap`, built either by
//! `put` one key at a time or by `putAll` into an empty map, and the two size their tables
//! differently: `new HashMap<>()` starts at 16 and doubles past three quarters full, while
//! `putAll` into an empty map pre-sizes to `tableSizeFor((int) (n / 0.75f + 1))`, which for one
//! entry is a table of two. A key lands in bucket `spread(hashCode) & (capacity - 1)` and a
//! bucket keeps insertion order, so the same keys come out in a different order from each
//! table size, and the order out of one map is the insertion order into the next.
//!
//! [`JavaHashMap`] replays those operations on the keys alone. Tree bins are not modelled: a
//! bucket becomes a tree only at eight colliding keys in a table of at least 64, which an INFO
//! column does not reach.

/// `String.hashCode`: over UTF-16 code units, wrapping.
pub fn string_hash(s: &str) -> i32 {
    s.encode_utf16()
        .fold(0i32, |h, unit| h.wrapping_mul(31).wrapping_add(unit as i32))
}

/// `HashMap.hash`: the high half folded into the low half.
fn spread(h: i32) -> u32 {
    let h = h as u32;
    h ^ (h >> 16)
}

/// `HashMap.tableSizeFor`: the next power of two, at least one.
fn table_size_for(cap: usize) -> usize {
    cap.max(1).next_power_of_two().min(1 << 30)
}

/// The bucket table of a `HashMap`, holding indexes into a key list the caller owns.
#[derive(Debug, Clone)]
pub struct JavaHashMap {
    table: Vec<Vec<usize>>,
    threshold: usize,
    size: usize,
}

impl JavaHashMap {
    /// `new HashMap<>()`.
    pub fn new() -> Self {
        Self {
            table: Vec::new(),
            threshold: 0,
            size: 0,
        }
    }

    /// `new HashMap<>(initialCapacity)`.
    pub fn with_capacity(initial_capacity: usize) -> Self {
        Self {
            table: Vec::new(),
            threshold: table_size_for(initial_capacity),
            size: 0,
        }
    }

    /// `HashMap.resize`, which keeps each bucket's order when it splits it.
    fn resize(&mut self, hashes: &[u32]) {
        let old_cap = self.table.len();
        let old_thr = self.threshold;
        let mut new_thr = 0;
        let new_cap;
        if old_cap > 0 {
            new_cap = old_cap << 1;
            if old_cap >= 16 {
                new_thr = old_thr << 1;
            }
        } else if old_thr > 0 {
            new_cap = old_thr;
        } else {
            new_cap = 16;
            new_thr = 12;
        }
        if new_thr == 0 {
            new_thr = (new_cap as f32 * 0.75) as usize;
        }
        self.threshold = new_thr;
        let mut table = vec![Vec::new(); new_cap];
        for bucket in &self.table {
            for &index in bucket {
                table[hashes[index] as usize & (new_cap - 1)].push(index);
            }
        }
        self.table = table;
    }

    /// `put(key, value)` for the key at `index`; a key already present keeps its place.
    pub fn put(&mut self, index: usize, keys: &[String], hashes: &[u32]) {
        if self.table.is_empty() {
            self.resize(hashes);
        }
        let bucket = hashes[index] as usize & (self.table.len() - 1);
        if self.table[bucket].iter().any(|&i| keys[i] == keys[index]) {
            return;
        }
        self.table[bucket].push(index);
        self.size += 1;
        if self.size > self.threshold {
            self.resize(hashes);
        }
    }

    /// `putAll(source)` into this map, with the source's iteration order given.
    pub fn put_all(&mut self, order: &[usize], keys: &[String], hashes: &[u32]) {
        let s = order.len();
        if s > 0 && self.table.is_empty() {
            let ft = s as f32 / 0.75 + 1.0;
            let t = if ft < (1u32 << 30) as f32 {
                ft as usize
            } else {
                1 << 30
            };
            if t > self.threshold {
                self.threshold = table_size_for(t);
            }
        }
        for &index in order {
            self.put(index, keys, hashes);
        }
    }

    /// The iteration order: bucket by bucket, each in insertion order.
    pub fn order(&self) -> Vec<usize> {
        self.table.iter().flatten().copied().collect()
    }
}

impl Default for JavaHashMap {
    fn default() -> Self {
        Self::new()
    }
}

/// One way a map is filled from the order of the one before it.
#[derive(Debug, Clone, Copy)]
pub enum Fill {
    /// `new HashMap<>()` then `put` per key.
    Put,
    /// `new HashMap<>(capacity)` then `put` per key.
    PutSized(usize),
    /// `new HashMap<>()` then `putAll`.
    PutAll,
}

/// The order `keys` come out in after passing through `fills`, one map after another, starting
/// from the order they are given in. Keys are distinct.
pub fn order_through(keys: &[String], fills: &[Fill]) -> Vec<usize> {
    let hashes: Vec<u32> = keys.iter().map(|k| spread(string_hash(k))).collect();
    let mut order: Vec<usize> = (0..keys.len()).collect();
    for fill in fills {
        let mut map = match fill {
            Fill::Put | Fill::PutAll => JavaHashMap::new(),
            Fill::PutSized(capacity) => JavaHashMap::with_capacity(*capacity),
        };
        match fill {
            Fill::PutAll => map.put_all(&order, keys, &hashes),
            _ => order.iter().for_each(|&i| map.put(i, keys, &hashes)),
        }
        order = map.order();
    }
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(keys: &[&str]) -> Vec<String> {
        keys.iter().map(|k| k.to_string()).collect()
    }

    #[test]
    fn string_hash_is_javas() {
        assert_eq!(string_hash(""), 0);
        assert_eq!(string_hash("DP"), 2188);
        assert_eq!(string_hash("hello"), 99162322);
        // Overflows wrap.
        assert_eq!(string_hash("polygenelubricants"), i32::MIN);
    }

    /// `{DP, AC, AF, AN}` in a table of 16: bucket of each key's spread hash, ties in order.
    #[test]
    fn a_default_map_orders_by_bucket() {
        let keys = strings(&["DP", "AC", "AF", "AN"]);
        let order = order_through(&keys, &[Fill::Put]);
        let names: Vec<&str> = order.iter().map(|&i| keys[i].as_str()).collect();
        // DP=2188 -> 12, AC=2082 -> 2, AF=2085 -> 5, AN=2093 -> 13.
        assert_eq!(names, ["AC", "AF", "DP", "AN"]);
        let mut sorted = keys.clone();
        sorted.sort();
        assert_eq!(order_through(&sorted, &[Fill::Put]).len(), 4);
    }

    /// `putAll` of one key makes a table of two, so a second map can reorder what the first kept.
    #[test]
    fn put_all_presizes_small_tables() {
        let mut map = JavaHashMap::new();
        let keys = strings(&["A"]);
        let hashes = [spread(string_hash("A"))];
        map.put_all(&[0], &keys, &hashes);
        assert_eq!(map.table.len(), 2);
        let keys = strings(&["a", "b", "c", "d", "e", "f"]);
        let hashes: Vec<u32> = keys.iter().map(|k| spread(string_hash(k))).collect();
        let mut map = JavaHashMap::new();
        map.put_all(&[0, 1, 2, 3, 4, 5], &keys, &hashes);
        assert_eq!(map.table.len(), 16);
    }

    #[test]
    fn growth_keeps_every_key() {
        let keys: Vec<String> = (0..40).map(|i| format!("K{i}")).collect();
        let order = order_through(&keys, &[Fill::Put, Fill::PutSized(10), Fill::PutAll]);
        let mut seen = order.clone();
        seen.sort();
        assert_eq!(seen, (0..40).collect::<Vec<_>>());
    }
}
