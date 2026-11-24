# The bug fixes of htsjdk 5.0.1, 5.0.2 and 6.0.0, triaged

htsjdk-rs ports htsjdk 5.0.0, the version GATK 4.7.0.0 and Picard 3.5.0 run and the one the port is
checked against byte for byte. Upstream has since released 5.0.1, 5.0.2 and 6.0.0. A fix from them
comes in only if it can't move a byte GATK or Picard writes on an input they process. Each fix
listed in the three releases' notes is sorted into one of three classes:

- **(a)** Changes no output on a file GATK or Picard accept: a crash, a wrong refusal, a leak, a
  race, or speed. These are ported where htsjdk-rs has the bug.
- **(b)** Changes output GATK 4.7.0.0 or Picard 3.5.0 would write, or the values they decode. These are
  not ported. Each one is an upstream divergence that htsjdk 6 fixes. Adopt it when GATK and Picard
  move to htsjdk 6. A fix marked **⚑** repairs a crash or data loss on input GATK accepts, so it is
  raised with the maintainer rather than decided here.
- **(c)** API-only (Java `Path`, unmodifiable collections, Java serialization, the class path).
  There is nothing to port.

"PR" is the upstream pull request in samtools/htsjdk. A 5.0.x fix names its backport too. Features
that fix nothing (VCF 4.4/4.5, CSI, CRAI, bgzipped SAM, new comparators and profiles) aren't listed.
Upstream's own list of what they change is in 6.0.0's "Output that differs from 5.x".

The **Port** column gives the htsjdk-rs result for class (a): **ported** (with its htsjdk-rs PR),
**not present** (the Rust code doesn't have the bug), or **no counterpart** (htsjdk-rs doesn't port
that code). For a (b) row it notes what 5.0.0 does, or that the port already behaves as 6.0.0
does, which is a divergence from the oracle in its own right (see the last section).

## 5.0.1 and 5.0.2

| Fix | PR | Class | Port |
|---|---|---|---|
| BCF decoding reads record bytes directly (speed) | #1856 | a | no counterpart (Java stream locking) |
| A truncated or corrupt BCF record throws instead of reading missing bytes as 0xFF | #1856 | b | 5.0.0 finishes and writes the `-1`s; the port reproduces it (`bcf2/types.rs` `Cursor::byte`) |
| VCF writing with many samples is faster (speed) | #1857, #1858 | a | no counterpart (Java `Formatter`/lock costs) |
| A record that fails to encode no longer leaves its columns at the start of the next line | #1857 | a | not present: `write_vcf` builds the record whole and returns `Err` |
| A CRAM `ReferenceSource` shared between threads no longer returns another contig's bases | #1792, #1868 | a | not present: `sequence_bases` takes `&mut self` |

## 6.0.0: SAM, BAM and CRAM

| Fix | PR | Class | Port |
|---|---|---|---|
| A mapped CRAM read without bases keeps its CIGAR, NM and MD (5.x wrote it back as `0M`) | #1837 | b ⚑ | ported as a divergence, #298 |
| CRAMs written readable by htsjdk 4.x (empty blocks RAW, not GZIP) | #1838 | b | |
| A CRAM labelled 3.0 holds no 3.1 codecs (FAST profile); valid 3.0 from a default strategy | #1838, #1936 | b | |
| Mate template lengths derived and stored as htslib does | #1843 | b | |
| Unsigned `B` arrays keep their type through CRAM; bytes 0x80–0xFF in read names and tags kept | #1843, #1840 | b ⚑ | ported as a divergence, #298 |
| CRAM NM counts every inserted base; bases past the reference end kept with NM and MD | #1850 | b | |
| CRAM reading matches samtools for MD past the reference end, `cF`, unpaired PNEXT | #1850, #1887 | b | |
| The CRAM writer keeps an unpaired read's pair/mate flags and PNEXT | #1887 | b ⚑ | ported as a divergence, #298 (not the STRICT check) |
| The CRAM writer keeps an `RG` tag that names no `@RG` line | #1891 | b ⚑ | ported as a divergence, #298 (not the warning) |
| The CRAM writer takes an empty quality array as missing qualities | #1850 | b | |
| The CRAM writer no longer crashes on a slice whose placed reads are all unmapped | #1855 | a | ported, #290: only the spans 5.0.0 threw on (the port wrote them negative); upstream's `max(end, start)` would move spans 5.0.0 writes |
| `@SQ` lines get an `M5` when a CRAM is written | #1855 | b | |
| A CRAM with an EOF container partway through is read to its end | #1876 | b ⚑ | ported as a divergence, #298 |
| CRAM quality scores with a byte-array HUFFMAN encoding can be read | #1854 | a | ported, #290 |
| CRAM iterators over one seekable stream no longer cut each other short | #1908 | a | no counterpart: each reader owns its bytes |
| A CRAM opened by `Path` is read through a buffer (speed on NFS) | #1822 | a | not present: the stream reader fills 1 MiB at a time |
| `ReferenceSource` keyed on the contig's name: one source serving CRAMs with `@SQ` in different orders | #1875 | b ⚑ | already 6.0.0's: the port's cache is keyed by name |
| `CRAMFileReader` constructors without a stringency use `SILENT` | #1909 | c | |
| The beta CRAM encoders write the version they name | #1906 | c | |
| A CRAM reader on a stream that can't seek throws `IllegalStateException`, not NPE | #1791 | a | no counterpart: no query on a stream |
| A failed CRAM reference lookup names the sequence and the files it needs | #1790, #1793 | b | message text |
| CRAM's rANS tables no longer fill the heap (one cache per writer) | #1815 | a | no counterpart (Java heap) |
| Multi-interval queries return every record when references are not in dictionary order | #1877 | b ⚑ | ported as a divergence, #294; gatk-rs's own copy in `reads.rs` follows in its bump |
| `queryUnmapped()` on a CSI-indexed BAM no longer seeks into the header | #1809 | b | |
| A read past 2^29 on a longer reference is refused by the BAI indexer, not a panic (5.0.0's `SAMException`, a 0-byte `.bai`); a read below it on such a reference binned by position | — | a | ported, #297: the port panicked in the writer, `build_bam_index` and the CRAM BAI |
| A BAM writer that can't build or write its index refuses later writes | #1809 | a | not present: 5.0.0 already fails the write, and the port has no later write |
| A BAI/TBI metadata bin with the wrong chunk count suggests a CSI | #1821 | b | message text |
| A sort-order check no longer fails because the caller changed a record it was given | #1844 | a | no counterpart (no asserting iterator) |
| `makeDefault()` uses the inflater factory in force when called | #1844 | c | |
| A plain SAM file can be iterated again after its first iterator closes | #1915 | a | not present: `read_sam_with` returns the records |
| A reader whose header or iterator fails no longer leaks its file handle | #1913 | a | not present (Rust drops it) |
| `SAMFileWriterFactory` copies keep the deflater factory and flag format | — | c | |
| `isUnsignedArrayAttribute` no longer throws NPE | #1885 | a | no counterpart: the `TagValue` variant carries the sign |
| Overlap clipping clips part of a `=` or `X` element; mates on different references aren't clipped | #1802 | b | |
| An out-of-order record error names the record before the offending one | #1797 | b | message text |
| Bgzipped SAM read with the factory's `InflaterFactory` | #1814 | c | |
| Dates with a space before the time, or a lower-case `t` or `z`, are parsed | #1918 | b | 5.0.0 keeps such a DT at LENIENT/SILENT and would rewrite a parsed one |
| `@RG DT` kept as written, not rewritten in the JVM's time zone | #1918 | b | |
| Header text UTF-8; read names and `Z`/`A` tags byte for byte | #1840 | b | |
| A CIGAR too long for BAM goes in `CG:B:I`, not `CG:B:i` | #1840 | b | |
| A U in SEQ is read and written as T | #1890 | b | |
| `x.sam.gz` written as bgzipped SAM; names in any case | caa921cc | b | |
| Records the spec allows on unmapped/unpaired reads are accepted by validation | #1917 | b | |
| `@SQ` names outside the RNAME pattern reported by stringency | #1920 | b | |
| A SAM record with a trailing tab is accepted; line numbers from the file's start | #1929 | a | ported, #292: 5.0.0 itself reads one trailing tab, the port refused it |
| A SAM record with RNAME or RNEXT `*` gets the unplaced index | #1933 | b | |
| `SamLocusIterator` and `EdgeReadIterator` no longer slow down as reads get longer | #1912 | a | not present: picard-rs's walkers queue in a `VecDeque` or `BTreeMap` |
| An indexed CRAM is written with a `.crai`, not a `.cram.bai` | #1798 | b | |
| BAI/TBI/CSI content as htslib's; TBI record counts | #1817, #1834 | b | |

## 6.0.0: VCF and BCF

| Fix | PR | Class | Port |
|---|---|---|---|
| The BCF writer keeps interior missing values, sizes FORMAT vectors and UTF-8 strings, writes String INFO whatever its comma count, accepts `Flag` with `Number`≠0, masks n_sample to 24 bits | #1828, #1873, #1892 | b ⚑ | interior missing values ported as a divergence, #295 (an Integer list with a null threw NPE); the rest unported |
| A missing value inside a BCF vector no longer shifts the values after it; one at a 2.2 vector's end is missing, not padding | #1827, #1889 | b ⚑ | ported for BCF 2.1 as a divergence, #295; `0x81` as END_OF_VECTOR not ported (5.0.0 writes -127 so) |
| A BCF GT starting with END_OF_VECTOR reads as a sample without GT | #1894 | b | |
| INFO values typed by their INFO line, FORMAT values by their FORMAT line | #1879 | b ⚑ | not ported: not in the approved list |
| `getAttributeAs{Int,Double,IntList,DoubleList}` return the default for a missing value | #1849, #1899 | a | ported, #291: 5.0.0 threw `NumberFormatException` |
| `getAttributeAsBoolean` returns the default for a missing value | #1899 | b | 5.0.0 returns `false` |
| Spaces in INFO values are accepted | #1825 | a | ported, #291: a tab is still refused (a sites-only line's ninth column) |
| A malformed genotype is a `TribbleException` naming the sample and key | #1849 | b | message text |
| `VCFFileReader` opens VCF 3.2/3.3; a 3.x `FILTER=PASS` passes | #1826 | b | |
| `smartMergeHeaders` promotes Integer/Float either order; keeps every contig when dictionaries differ | #1878 | b ⚑ | contigs ported as a divergence, #296 (the Integer/Float promotion is a 5.0.0 warning, not ported); gatk-rs's `multi_variant.rs` follows in its bump |
| `sameGenotype` ignoring phase tells polyploid genotypes apart | #1881 | b | |
| `getStructuralVariantType()` is `DEL` for `SVTYPE=DEL:ME` | #1830 | b | |
| `hasID()` false for a deserialized record without an ID | #1882 | c | |
| An empty VCF is reported as empty | #1849 | b | message text |
| `toString()` no longer decodes lazy genotypes | #1849 | a | not present: `Debug` doesn't decode |
| Closing a VCF/BCF writer or `BlockCompressedOutputStream` twice does nothing | #1845 | a | not present: BGZF `finish` is idempotent |
| VCF and BCF work without commons-jexl on the class path | #1880 | c | |
| `VCFFileReader` opens a path in a zip file system | — | c | |
| `VariantContextWriterBuilder` keeps a caller's `IndexCreator` for bgzipped VCF | #1807 | c | |
| Writing out of order with on-the-fly indexing throws | #1845 | b | |
| A header keeps the first INFO/FORMAT/FILTER/structured line per ID; numbers its contigs | #1916 | b | |
| Repeated unstructured lines kept; contig md5, URL, species carried | #1846 | b | |
| A standard line declared with another type kept as declared; `PQ` is Integer | #1816 | b | |
| A header keeps its declared version | #1816 | b | |
| A record's end from `END`, `SVLEN` and `LEN` as htslib does | #1823 | b | |
| A genotype keeps each allele's phase | #1820 | b | |
| Negative doubles written like positive ones; missing GT as `.`; all-missing sample as `.` | #1849, #1888, #1895 | b | |
| VCF, GFF3 and Tribble text read and written as UTF-8 | #1824 | b | |
| Sample and header collections unmodifiable | #1934 | c | |
| VCF/BCF recognised by one naming rule | #1923 | b | |
| The VCF codecs are thread-safe once the header is read | #1812 | a | not present (no shared codec state) |

## 6.0.0: indexes, tabix, FASTA, I/O

| Fix | PR | Class | Port |
|---|---|---|---|
| A tabix index with a sequence that has no records is written correctly | #1806 | b | the port's creator never writes one |
| `TabixReader` returns a last line without newline, strips CR, rejects plain gzip, names the indexes it tried | #1806 | b | no counterpart (no `TabixReader`) |
| `TabixReader` region strings accept contig names with `:` or `-` | #1852 | a | no counterpart |
| A Tribble index with a byte ≥ 0x80 in a name is read | #1831 | a | not present: names are read lossily |
| An empty or inverted query on a Tribble `.idx` returns nothing | #1911 | b | |
| A Tribble index can't be built over BGZF BCF | #1907 | b | |
| BAI/CSI read lazily, sparse binning index (speed, memory) | #1810, #1806 | a | no counterpart |
| `FastaSequenceIndexCreator` indexes blank lines and a last line without newline | #1848 | a | no counterpart (no index creator) |
| A sequence longer than 2^30 bases without a dictionary no longer fails | #1848 | a | not present: `read_fasta` grows a `Vec` |
| A gzipped (not bgzipped) FASTA with a `.fai` read sequentially | #1848 | a | no counterpart |
| A FASTA whose first lines are `;` comments can be read and indexed | #1922 | a | ported, #292 |
| An empty indexed FASTA with an empty `.fai` opens | #1871 | a | not present: the bug is in 6.0.0's own new check |
| Opening an indexed FASTA checks it matches its `.fai` | #1745 | b | new refusal |
| Concurrent indexed FASTA lookups; bgzipped lookups don't re-inflate | #1842 | a | not present |
| `.gzi` written so htsjdk can read it; bgzipped `FastaReferenceWriter` output indexable | #1842 | b | |
| `GZIIndex.buildIndex` handles more than 2 GiB | #1886 | a | no counterpart (no GZI code) |
| `FastaReferenceWriter` with no sequences writes an empty FASTA | #1925 | b | |
| A short sequence's own length as its `.fai` line width | #1926 | b | |
| Gzipped and bgzipped input from pipes, sockets and URLs no longer truncated | #1691 | a | not present: BGZF is read a block at a time |
| A BGZF file with an empty block not at the end is read to its end | 24abe9ee | b ⚑ | already 6.0.0's: `BgzfReader::fill` skips empty blocks |
| `IOUtil.isGZIPInputStream` checks only the gzip header | #1691 | a | no counterpart |
| FTP fixes; format of an FTP file by its name | #1898 | c | |
| `SeekableMemoryStream.read(buf, off, 0)` returns 0; `SeekablePathStream` retries an empty read | #1810, #1832 | c | |
| `IOUtil.deletePaths` given one `Path` | #1872 | c | |
| `assertFileIsWritable` on NFS | #1902 | c | |
| `ProcessExecutor` no longer deadlocks on large output | #1805 | c | |
| `SortingCollection` buffers temp files; `SortingLongCollection` clamps its buffer | #1853, #1883 | a | no counterpart |
| `TagFilter` accepts any `List` | #1853 | c | |
| `Histogram.getMin`/`getMax`/`getMode` throw `IllegalStateException`, not NPE, on an empty histogram | #1900 | a | not present: they return `Option`, and Picard never asks an empty one |
| Malformed BED/GFF3 lines reported with file, column and text | #1904 | b | message text |
| `.md5` of a gzipped file is of the compressed bytes; FASTQ ends lines with LF everywhere | #1851 | b | |
| Warnings through `Log`, not `System.err` | #1910 | b | stderr text |
| The POM lists all three licences | #1847 | c | |

## Divergences to adopt with htsjdk 6

Every class (b) row is one. When GATK and Picard pin htsjdk 6, each becomes a port to make, and the
goldens it moves get refrozen against the new oracle. Error and warning text counts as output here,
because the oracle prints it.

The ⚑ rows repair data loss or a crash on input GATK accepts. The maintainer chose to port the CRAM,
BCF interior-missing, multi-interval query and contig-merge ones as **deliberate divergences from
GATK 4.7.0.0 and Picard 3.5.0**: #294, #295, #296 and #298. Each suite row they move is checked as a
listed divergence with its count, never refrozen. `#1879` (INFO typed by a FORMAT line sharing its
ID) and the BCF String-INFO and `Flag` writer fixes stay unported.

| Divergence | Inputs it changes | Upstream | htsjdk-rs | Rows lifted |
|---|---|---|---|---|
| A mapped CRAM read without bases keeps CIGAR, NM, MD | CRAM with SEQ `*` on mapped reads | #1837 | #298 | `cram_read_features` `no-sequence` (1) |
| Unsigned `B` arrays and bytes 0x80–0xFF kept in CRAM | CRAM with `B:C/S/I`, non-ASCII names or tags | #1843, #1840 | #298 | — |
| Unpaired mate flags, PNEXT and orphan `RG` kept in CRAM | CRAM of htslib-style unpaired reads | #1887, #1891 | #298 | — |
| A CRAM is read past a mid-file EOF container | CRAMs joined by `samtools cat` < 1.13 | #1876 | #298 | — |
| An interior missing BCF value stays in place | BCF 2.1 with `1,.,3` or AD `5,.,7` | #1827, #1828, #1889 | #295 | `bcf2_conformance` `refuse-nullint` (1); picard-rs `bcf_tools` `refuse-nullint/converter-bcf` |
| Multi-interval queries keep a cursor per reference | Indexed BAM/CRAM with references out of `@SQ` order | #1877 | #294 | `query_conformance` `filter` rows 4, 6, 8 × 2 (6) |
| Merged headers keep every contig | VCF inputs whose dictionaries differ | #1878 | #296 | gatk-rs `multi_variant_walker` header rows |

The bump PRs in gatk-rs and picard-rs list the rows they lift there.

## Where htsjdk-rs already differs from 5.0.0

The checks behind this table found four places where the port doesn't do what 5.0.0 does. Two of
them are 6.0.0's behaviour already. None of them came from this triage's ports.

- **An empty BGZF block that isn't the last** (⚑ above). 5.0.0 stops there; the port reads on.
- **`ReferenceSource` keyed on index** (⚑ above). 5.0.0 can return the wrong contig; the port can't.
- **`@RG DT`**. 5.0.0 rewrites every DT it can parse as `yyyy-MM-dd'T'HH:mm:ssZ` in the JVM's
  time zone (`2016-01-01` becomes `2016-01-01T00:00:00+0000` under UTC), so a re-encoded header differs.
  The port copies DT as written, which is what 6.0.0 does (#1918).
- **An empty interior SAM field**. 5.0.0 reports it by stringency; the port refuses it at every
  stringency (`text_parse.rs`, `EmptyField`).

Each would be its own fix against the oracle, not part of this upstream triage.
