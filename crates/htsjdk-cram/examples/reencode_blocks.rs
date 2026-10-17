//! `reencode_blocks <in.cram>...`: decompress every external block and compress it again with the
//! encoder for its method, reporting any block whose bytes differ from the file's.
use htsjdk_cram::block::{BlockHeader, CompressionMethod, ContentType};
use htsjdk_cram::compression_header::CompressionHeader;
use htsjdk_cram::container::{ContainerHeader, FileDefinition, FILE_DEFINITION_LENGTH};
use htsjdk_cram::file::{decompress, read_container};
use htsjdk_cram::record_read::{read_record, RecordReaders, SliceContext};

fn main() {
    let (mut same, mut differ) = (0usize, 0usize);
    for path in std::env::args().skip(1) {
        let cram = std::fs::read(&path).unwrap();
        let definition = FileDefinition::read(&cram).unwrap();
        let major = definition.major;
        // The header container's one block.
        let header = ContainerHeader::read(&cram[FILE_DEFINITION_LENGTH..], major).unwrap();
        let at = FILE_DEFINITION_LENGTH + header.byte_length;
        let block = BlockHeader::read(&cram[at..]).unwrap();
        let raw =
            &cram[at + block.byte_length..at + block.byte_length + block.compressed_size as usize];
        let content = decompress(block.method().unwrap(), raw, block.uncompressed_size).unwrap();
        let ok = htsjdk_cram::gzip::compress(&content, 5) == raw;
        println!(
            "{path}: header block {}",
            if ok { "same" } else { "DIFFERS" }
        );

        let mut at = at + header.blocks_byte_size as usize;
        loop {
            let container_header = ContainerHeader::read(&cram[at..], major).unwrap();
            if container_header.is_eof() {
                break;
            }
            let blocks_start = at + container_header.byte_length;
            let end = blocks_start + container_header.blocks_byte_size as usize;
            let compression =
                CompressionHeader::read_block(&cram[blocks_start..end], major).unwrap();
            let (walk, _) = read_container(&cram, at, major).unwrap();
            // Walk the raw blocks again for their compressed bytes, slice by slice.
            let mut cursor = blocks_start
                + BlockHeader::read(&cram[blocks_start..])
                    .unwrap()
                    .total_length(major);
            for slice in &walk.slices {
                let readers = RecordReaders::new(&compression.encodings).unwrap();
                let blocks = slice.block_bytes();
                let mut streams = htsjdk_cram::external_codecs::SliceReadStreams::new(&blocks);
                let context = SliceContext {
                    reference_context: slice.header.reference_context_id,
                    alignment_start: slice.header.alignment_start,
                };
                let mut lengths = Vec::new();
                let mut flags = Vec::new();
                let mut previous = slice.header.alignment_start;
                for _ in 0..slice.header.record_count {
                    let r = read_record(&readers, &compression, &context, &mut streams, previous)
                        .unwrap();
                    previous = r.alignment_start;
                    if r.flags.is_force_preserve_quality_scores() && r.read_length != 0 {
                        lengths.push(r.read_length);
                        flags.push(r.flags.bam);
                    }
                }
                // The slice header block, then the slice's blocks.
                cursor += BlockHeader::read(&cram[cursor..])
                    .unwrap()
                    .total_length(major);
                for _ in 0..slice.header.block_count {
                    let h = BlockHeader::read(&cram[cursor..]).unwrap();
                    let start = cursor + h.byte_length;
                    let raw = &cram[start..start + h.compressed_size as usize];
                    cursor += h.total_length(major);
                    if h.content_type() != Some(ContentType::External) {
                        continue;
                    }
                    let method = h.method().unwrap();
                    let content = decompress(method, raw, h.uncompressed_size).unwrap();
                    let candidates: Vec<Vec<u8>> = match method {
                        CompressionMethod::Raw => vec![content.clone()],
                        CompressionMethod::Gzip => vec![htsjdk_cram::gzip::compress(&content, 5)],
                        CompressionMethod::RansNx16 => [0u8, 1]
                            .iter()
                            .map(|f| htsjdk_cram::rans_nx16::compress(&content, *f))
                            .collect(),
                        CompressionMethod::Fqzcomp => {
                            let sum: i32 = lengths.iter().sum();
                            let l = if sum as usize == content.len() {
                                lengths.clone()
                            } else {
                                vec![content.len() as i32]
                            };
                            let f = if sum as usize == content.len() {
                                flags.clone()
                            } else {
                                vec![0]
                            };
                            vec![htsjdk_cram::fqzcomp::compress(&content, &l, Some(&f))]
                        }
                        CompressionMethod::NameTokeniser => {
                            vec![htsjdk_cram::name_tokeniser::compress(&content, 0).unwrap()]
                        }
                        _ => vec![],
                    };
                    if candidates.iter().any(|c| c == raw) {
                        same += 1;
                    } else {
                        differ += 1;
                        println!(
                            "{path}: container@{at} block id={} {} DIFFERS (len {} vs {:?})",
                            h.content_id,
                            method.name(),
                            raw.len(),
                            candidates.iter().map(Vec::len).collect::<Vec<_>>()
                        );
                    }
                }
            }
            at = end;
        }
    }
    println!("same={same} differ={differ}");
}
