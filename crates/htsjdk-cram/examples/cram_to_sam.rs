//! `cram_to_sam <in.cram> [ref.fasta]`: the records as SAM text on stdout, for comparing a read
//! against the reference's.
use htsjdk_cram::reader::{CramReader, FastaReference, NoReference, ReferenceSource};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let bytes = std::fs::read(&args[1]).expect("the CRAM reads");
    let mut reader = CramReader::new(&bytes).unwrap_or_else(|e| panic!("{}", e.java_line()));
    let mut fasta;
    let mut none = NoReference;
    let source: &mut dyn ReferenceSource = match args.get(2) {
        Some(path) => {
            fasta = FastaReference::open(std::path::Path::new(path));
            &mut fasta
        }
        None => &mut none,
    };
    match reader.read_all(source) {
        Ok(records) => print!(
            "{}",
            htsjdk_bam::sam_file::write_sam(&reader.header, &records).unwrap()
        ),
        Err(e) => {
            eprintln!("Exception in thread \"main\" {}", e.java_line());
            std::process::exit(1);
        }
    }
}
