//! `bam_to_cram <in.bam> <out.cram> <ref.fasta> <file id>`: SamFormatConverter's CRAM output.
use htsjdk_cram::reader::FastaReference;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let raw = std::fs::read(&args[1]).unwrap();
    let plain = htsjdk_bgzf::decompress_all(&raw).unwrap();
    let reader = htsjdk_bam::reader::BamReader::new(&plain).unwrap();
    let header = reader.header.text.clone();
    let records: Vec<_> = reader.map(|r| r.unwrap()).collect();
    let mut fasta = FastaReference::open(std::path::Path::new(&args[3]));
    match htsjdk_cram::writer::write_cram(&header, &records, &args[4], &mut fasta) {
        Ok(bytes) => std::fs::write(&args[2], bytes).unwrap(),
        Err(e) => {
            eprintln!("Exception in thread \"main\" {}", e.java_line());
            std::process::exit(1);
        }
    }
}
