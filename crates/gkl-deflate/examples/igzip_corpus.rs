//! Writes `<file>.r1` and `<file>.r2` next to every file named: this crate's levels 1 and 2. The
//! oracle in `tools/gkl-probe` writes `.l1` and `.l2` beside them from ISA-L itself.
fn main() {
    for path in std::env::args().skip(1) {
        let data = std::fs::read(&path).unwrap();
        for level in 1..=2 {
            std::fs::write(
                format!("{path}.r{level}"),
                gkl_deflate::deflate_gkl(&data, level),
            )
            .unwrap();
        }
    }
}
