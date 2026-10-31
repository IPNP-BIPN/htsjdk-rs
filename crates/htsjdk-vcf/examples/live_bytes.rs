//! Live heap per BAM record and per VCF record, counted by the allocator rather than reasoned out.
//! `live_bytes bam file.bam` or `live_bytes vcf file.vcf`.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

struct Counting;
static LIVE: AtomicUsize = AtomicUsize::new(0);
static COUNT: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        LIVE.fetch_add(layout.size(), Relaxed);
        COUNT.fetch_add(1, Relaxed);
        System.alloc(layout)
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        LIVE.fetch_sub(layout.size(), Relaxed);
        COUNT.fetch_sub(1, Relaxed);
        System.dealloc(ptr, layout)
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args[0] == "bam" {
        let data = htsjdk_bgzf::decompress_all(&std::fs::read(&args[1]).unwrap()).unwrap();
        let reader = htsjdk_bam::BamReader::new(&data).unwrap();
        let before = (LIVE.load(Relaxed), COUNT.load(Relaxed));
        let records: Vec<_> = reader.map(|r| r.unwrap()).collect();
        report(
            "BAM",
            records.len(),
            before,
            std::mem::size_of_val(&records[0]),
        );
        println!("on disk per record {}", data.len() / records.len());
        println!(
            "TagValue {} entry {}",
            std::mem::size_of::<htsjdk_bam::tag::TagValue>(),
            std::mem::size_of::<(htsjdk_bam::tag::Tag, htsjdk_bam::tag::TagValue)>()
        );
    } else {
        let text = std::fs::read_to_string(&args[1]).unwrap();
        let before = (LIVE.load(Relaxed), COUNT.load(Relaxed));
        let file = htsjdk_vcf::reader::read_vcf(&text).unwrap();
        report(
            "VCF",
            file.records.len(),
            before,
            std::mem::size_of_val(&file.records[0]),
        );
        println!("text per record {}", text.len() / file.records.len());
        println!(
            "Genotype {} Allele {} Value {}",
            std::mem::size_of::<htsjdk_vcf::variant::Genotype>(),
            std::mem::size_of::<htsjdk_vcf::allele::Allele>(),
            std::mem::size_of::<htsjdk_vcf::variant::Value>()
        );
    }
}

fn report(what: &str, n: usize, before: (usize, usize), inline: usize) {
    println!(
        "{what}: {n} records, {} live bytes each, {} allocations each, {inline} inline",
        (LIVE.load(Relaxed) - before.0) / n,
        (COUNT.load(Relaxed) - before.1) / n
    );
}
