# goose-seqair

Experimental BAM/SAM/CRAM reading and writing, and raw pileups in [Goose](https://github.com/aardappel/goose),
inspired by [Seqair](https://github.com/Softleif/seqair). The parser and pileup
engine are Goose; [goose-bgzf](https://github.com/sstadick/goose-bgzf)
supplies shared BGZF reading/writing, backed by
[goose-deflate](https://github.com/sstadick/goose-deflate). Small native kernels
handle base packing/unpacking and packed rANS decoding. No htslib or Rust runtime dependency.

**[API and executable examples](API.md)** · **[Performance baseline](bench/README.md)**

Reads BAM sequentially or loads regions through **BAI/CSI indexes**, then
links retained mates and produces raw pileup columns. Indexed queries reuse
the open reader, parsed index and region slabs, including backward queries.
The API is experimental.

## Try it locally

Run from this checkout. The generated fixture needs only Python's standard
library; substitute your own coordinate-sorted BAM and contig/coordinates.

```sh
python3 tests/fixtures.py tests/data/example.bam
goose run examples/records.goose -- tests/data/example.bam
goose run examples/pileup.goose -- tests/data/example.bam chrA 5 18
goose run examples/pileup.goose -- tests/data/example.bam chrA 5 18 tests/data/example.bam.bai
```

`examples/records.goose` prints name, reference ID, start, exclusive end, and
sequence length. `examples/pileup.goose` prints zero-based position and raw
depth. Depth includes deletion and reference-skip entries.

The [complete pileup example](examples/pileup.goose) declares all transitive
imports required by Goose. Its main API calls are:

```goose
let reader = bam::Reader { .. };
defer reader.close();
assert(reader.open("tests/data/example.bam") == bam::Status.Record);
let index = bam_index::Index { .. };
assert(index.open("tests/data/example.bam.bai") == bam::Status.Record);
let store = pileup::Store { .. };
assert(store.load_indexed(reader, index, reader.reference_id("chrA"), 5, 18) == bam::Status.Record);
assert(store.mate_stats.ambiguous_qnames == 0);
store.columns() { pos, entries =>
    for &entry in entries {
        let record = store.record(entry.record_index);
        // record.name(), record.fields.flags, entry.base, entry.quality,
        // entry.qpos, entry.op, entry.indel and record.fields.mate_index are available.
        assert(record.fields.pos <= pos);
    }
};
```

Requires Goose, CMake 3.20+, and a native C compiler, including for `goose run`
(libdeflate is compiled natively). Tested with Goose
[`74c641c5`](https://github.com/aardappel/goose/commit/74c641c5ac6d762e57e3177678fe9dea414d042a).
Runtime dependencies are pinned in the source and example imports.
Root programs must also declare `import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";`
and select goose-io v0.3.0 for descriptor-owned seeking.
The former internal BGZF implementation is now that standalone package.

CRAM imports goose-simd v0.3.0. Its two shared native-only headers are vendored
unchanged at the same revision, with provenance and license in `vendor/goose-simd`:
Goose's separately compiled C sources cannot resolve another package's headers
by package name. No sibling SIMD checkout is required.

## Ownership, coordinates, and errors

- `bam::Reader` is a root owner with reusable BGZF, header, reference-dictionary,
  and cross-block scratch buffers. Like Seqair, records wholly within one
  decoded BGZF block borrow it without copying; only cross-block records use
  scratch. Close explicitly, normally with `defer`. Reopening closes the previous
  file; closing twice is safe.
- `reader.record()` borrows bytes until the next read/seek/open/close. Store a local
  record view before calling accessors. Copying the view does not own its bytes.
- `pileup::Store` follows Seqair's slab model: compact metadata with offsets
  into separate names, decoded IUPAC bases, native CIGAR words, qualities, and
  raw auxiliary bytes. Only overlapping, accepted records are retained; no
  packed BAM copy is kept. Its `pileup::Record` views survive closing the reader,
  but not the next load or slab mutation. Treat implementation fields as read-only.
- `store.columns()` borrows its entry slice for one callback. Cursors advance
  through each read's CIGAR monotonically and cache slab offsets. Scratch buffers
  are reused across columns; the hot loop does not reconstruct record views.
  Entry order follows input order; uncovered positions are omitted.
- Coordinates are zero-based and half-open. Reference IDs follow the BAM
  header. Reader views decode packed sequence on access; region stores expand
  it once at load time. Both retain BAM orientation: reverse-strand alignments
  are **not** reverse-complemented a second time. All 16 IUPAC codes are preserved.
- `bam::Status.Record` means success for open/load and one available alignment
  for read. `End` means clean EOF; `Invalid`, `Truncated`, `IoError`, `Closed`,
  and `TooLarge` distinguish failures. Reader EOF/errors remain sticky until
  reopen or seek. A failed region load clears the store rather than returning a
  partial result. `Unsorted` detects decreasing mapped coordinates in visited
  data only; indexed queries do not validate skipped blocks.
- Default allocation limits are 16 MiB per BAM record/header text/reference-name
  buffer, 1,048,576 references, and 4 GiB minus one byte of variable slab payload
  per region. Slabs grow on demand; the limit does not preallocate memory.
  The latter counts name + decoded bases + CIGAR + qualities + auxiliary bytes,
  not packed record size or total memory: metadata and buffer capacity add overhead.
  Set `Store.max_bytes` lower before loading to impose a smaller budget. Region
  payload cannot exceed 2^32-1 bytes, keeping slab offsets representable.

`bam_index::Index.open` takes an explicit BAI or CSI path, raw or BGZF-compressed.
It bounds decompressed index input to 256 MiB by default (`max_bytes`), validates
counts/chunks, excludes pseudo-bin statistics, and merges overlapping query
chunks. BAI supports coordinates below 2^29; use CSI for larger references.
The index must match the BAM: reference counts are checked, but these formats
have no file-identity checksum, so stale indexes cannot always be detected.
Index creation is not included; use `samtools index` or `samtools index -c`.
The sequential `store.load` remains available without an index; reopen the
reader before each independent sequential query.

Indexed loading follows Seqair's staged decoding: inspect core fields and the
wire CIGAR span to reject non-overlapping/filtered reads before full decoding.
BGZF integrity and record framing are always checked, but names, auxiliary tags
and query lengths of discarded records are **not** fully validated. Retained
records still receive the full validation used by `reader.read()`. A long-CIGAR
escape is resolved before rejecting it by span, rather than trusting its
placeholder. Use sequential reading to validate every record, not a region
query. Loading invalidates the reader's current record view; call `read()`
before borrowing another one with `record()`.

Both loaders automatically link eligible primary mapped mates already in the
store. `record.fields.mate_index` is the partner's store index, or -1. Linking
matches name bytes, same-reference mate fields and reciprocal positions, like
Seqair; READ1/READ2 bits alone do not determine a pair. Secondary/supplementary,
unmapped/mate-unmapped, empty-name, and cross-reference mates do not link.
`store.mate_stats` reports pairs and ambiguous same-name arrivals; it does not
silently collapse extra primary alignments. Links are replaced on load and can
be rebuilt with `store.link_mates()`. Mates outside the region are not fetched,
and linking does not change column depth or remove overlapping observations.

BGZF checks include its `BC` subfield, compressed/uncompressed sizes, CRC, and
exact decompressor consumption. Records can cross any number of blocks. Empty
blocks are skipped, including internal EOF markers; a missing final marker is
accepted at a complete block boundary. BAM validation checks layout, lengths,
auxiliary types, and mapped CIGAR/query lengths, not every SAM semantic rule.
Auxiliary bytes remain raw, including `CG:B:I` long-CIGAR data; the CIGAR
accessor transparently uses that array when the BAM long-CIGAR escape applies.

## Raw pileups, not samtools mpileup defaults

Default read filtering excludes unmapped, secondary, QC-failed, and duplicate
reads (`0x704`); supplementary reads remain. Set `Filter.min_mapq` and
`Filter.exclude_flags` for other policies. Unmapped and zero-reference-span
records never create columns, even with an exclusion mask of zero.

`Entry.op` distinguishes a base, deletion, and reference skip. `indel` is the
signed length anchored after the current column: positive insertion, negative
deletion. Insertions/clips do not create reference columns. `qpos` is -1 for
deletions/skips or absent sequence; unavailable base/quality are `N`/255.
For absent sequence this deliberately differs from htslib's hypothetical query
offset. Quality 255 means missing, not a high-quality observation.

There is no base-quality filtering, BAQ, mate-overlap adjustment, orphan
filtering, or depth cap. Callers can count only `Op.Base` and apply their own
base-quality policy. There is no typed auxiliary lookup
or parallel BGZF decompression yet. Seqair's persistent distant-bin
cache, per-record extras and realignment are not implemented. The storage,
indexed loading and forward-cursor model do not imply full feature parity.

## SAM reading

Import `seqair.sam` and use `sam::Reader` with the same `open/read/record/close`
pattern as BAM. Plain text and BGZF SAM are detected by their contents.
Header-only and headerless unmapped files, CRLF, final unterminated lines,
standard auxiliary types, long CIGARs and CG promotion are supported. Ordinary
gzip SAM is not implemented yet. `line_number` includes header lines and helps
locate a parsing error.

Like Seqair's `SamRecordEncoder`, the parser reuses a packed BAM payload; its
`bam::Record` views expire on the next reader operation. `store.load(reader, ...)`
accepts either reader and feeds the same decoded slabs, mate links and forward
CIGAR cursors. SAM region loading is sequential; there is no text index yet.

Parsing follows Seqair's decimal-number checks, reference/position fixups,
unsigned-first auxiliary integers and SAM-to-IUPAC mapping. It rejects junk
instead of adopting htslib's permissive numeric defaults. Unlike upstream's
recoverable malformed-line errors, errors here stay sticky, matching the Goose
BAM API. Empty/duplicate reference names, nonpositive reference lengths, lengths
beyond BAM's signed-32-bit limit, and deprecated CIGAR `B` are rejected. Header
and line allocations default to 16 MiB; record limits match BAM.

## BAM writing

Import `seqair.bam_writer` and, for newly constructed alignments,
`seqair.bam_encode`. The writer follows Seqair's eager-header, reusable-buffer,
single-threaded BGZF path, with libdeflate levels 0–12 (default 6). It preserves
caller order and keeps ordinary records in one block; oversized records span
blocks. It does not sort, build an index, or compress in parallel.

```goose
let writer = bam_writer::Writer { .. };
defer writer.close();
assert(writer.open("copy.bam", reader.header, reader.references, reader.names) == bam::Status.Record);
while reader.read() == bam::Status.Record {
    let record = reader.record();
    assert(writer.write(record) == bam::Status.Record);
}
assert(reader.status == bam::Status.End);
assert(writer.close() == bam::Status.Record);
```

Use the same reference dictionary/order as the input; IDs are not remapped.
Raw reader records are copied byte-for-byte. Slab records are repacked directly
from their decoded bases and native CIGAR words. `bam_encode::Alignment` accepts
borrowed fields for new records: ASCII IUPAC sequence in stored orientation,
raw Phred qualities (empty means missing), native CIGAR words and raw BAM aux.
No input views are retained after `write`. Long CIGARs use the BAM CG escape.

**Check `close()`**, not only `write()`: it flushes the last block and writes
the EOF member, and can fail. Deferred cleanup alone cannot report that failure.
A failed write poisons the writer and closes without a success EOF marker.
Output creation truncates the destination; do not use the input file as output.

## CRAM reading

`cram::Reader` reads **CRAM 3.0 and 3.1** into the same borrowed `bam::Record`
interface. The implementation follows Seqair's container → compression header
→ slice → decoded slabs → mate resolution model. It reconstructs all records
in a slice before exposing any, including mate chains, lossy names, missing
sequence/quality, auxiliary tags, read groups and MD/NM synthesis. It is not an
htslib wrapper. CRAM framing, data-series decoding, reconstruction and entropy
tables are Goose; native kernels implement Seqair's packed scalar, AVX2,
first-pass AVX-512 and NEON rANS paths with checked tails and runtime CPU
dispatch. Shared goose-simd operations supply masked expansion and the
BAM/CRAM nibble decoder. Arithmetic coding uses Seqair's four-entry scalar
prefix and depth-selected paired-vector search; FQZ/run models stay scalar.

Supply a reference provider with `fetch(name, begin, end, output) -> bool`.
It must replace `output` with the exact zero-based half-open window requested.
For goose-bio's indexed FASTA reader, the adapter is:

```goose
struct Reference { reader: faidx::Reader& }
fn fetch(source: Reference&, name: const u8[:], begin: i64, end: i64, output: u8[>..]&) -> bool {
    source.reader.fetch(name, begin, end, output) == faidx::Status.Ok
}
```

After opening the FASTA and creating `source`, stream with
`reader.read(source)` or load a region with
`cram::load(store, reader, source, tid, begin, end, filter)`. Then use the same
`store.columns()` and `store.record()` as for BAM. Region loading is sequential,
requires coordinate-sorted input, and links retained mates. Reopen before an
independent query. The [CRAM benchmark](bench/cram.goose) has complete imports
and a working faidx adapter using goose-bio v0.4.1; no sibling Bio checkout is needed.

For indexed access, open a matching plain or gzip/BGZF CRAI with
`cram_index::Index.open`, select a region with
`cursor.fetch(index, reader, tid, begin, end)`, then call
`cursor.read(reader, source)` until `End`. Records use `reader.record()`.
Each worker owns its reader and cursor; backward and repeated queries reuse
container and codec buffers. Only selected slices are decoded. Nested spans
are retained, duplicate slice addresses removed, and records stay in file order.
CRAI has no identity checksum; skipped slices and stale indexes cannot be
fully validated. Index creation is not included.

Single-reference slices fetch their declared window. Multi-reference slices
pre-scan data series and fetch one union window per contig rather than whole
contigs. External reference bases are uppercased; nonzero single-reference
slice MD5s are checked. Embedded and no-reference slices never call `fetch`;
a provider returning false is sufficient for those files. Embedded consensus
is not checked against an unrelated FASTA. No reference downloading occurs.

Supported block codecs: raw, gzip, bzip2, XZ, rANS4x8, rANS Nx16, arithmetic,
FQZComp and Tok3. **Bzip2/XZ require system libbz2/liblzma**, loaded on first use;
these are standard codec libraries, not htslib. Missing backends return `Invalid`;
direct bzip2/XZ blocks also set `unsupported_method`. Container and decoded-block CRCs are checked;
like Seqair, EOF is recognized by its header without decoding its dummy blocks.

Streaming and CRAI-indexed access are experimental.
CRAM 2.x and ordinary-gzip SAM remain unsupported.
Limits bound individual containers, decoded buffers and records, not aggregate
process memory. Linux x86-64 is tested; the optional native loader is POSIX-only.

## SAM and CRAM writing

Import `seqair.sam_writer` or `seqair.cram_writer`. Both use the same
`open(path, header, references, names)`, `write(record)` and `close()` pattern
shown for BAM, accepting borrowed BAM records, slab records, or
`bam_encode::Alignment`. Input views are never retained. Check the final close
status; output can fail during flushing. Records must use the supplied reference
dictionary. Neither writer sorts, builds indexes, or writes in parallel.

**SAM:** plain text with Seqair's header, CIGAR, sequence, quality and auxiliary
formatting. Missing SQ lines are synthesized from the binary dictionary. Long
CIGAR escapes expand into the CIGAR column; integer tags normalize to `i`.
Scalar and array floats match Seqair/htslib's six-significant-digit formatting,
so SAM is not a bit-exact interchange for arbitrary floating-point tags.
The native float formatter uses a thread-local C numeric locale (POSIX).
Defaults are 16 MiB input payloads and 64 MiB output lines, including newline;
set `max_record_bytes`/`max_line_bytes` before opening to change them.

**CRAM:** reference-free **3.0/3.1** output; no FASTA is needed to read it.
`open(..., level=6, minor_version=1)` selects libdeflate level 0–12 and version.
Level 0 emits raw streams; other levels use gzip when smaller. Encoding uses
the inverse of Seqair's container/slice/data-series model: independent reusable
streams, explicit bases, ordered tag dictionaries, detached mate fields, CRCs
and the canonical EOF container. Seqair has no CRAM writer at the pinned revision;
this is a new encoder, not a port of an upstream encoder or an htslib wrapper.

The writer stages up to `records_per_slice` (default 10,000), one multi-reference
slice per container. `flush()` emits a pending container and flushes buffered
I/O; it is not fsync. `max_slice_bytes` defaults to 64 MiB with conservative
record-size admission, and `max_record_bytes` to 16 MiB. These limit payloads,
not total process memory. Reuse retains allocated buffers.

CRAM decoders normalize CIGAR `=`/`X` to `M`, remove zero-length operations and
merge adjacent equal operations. Bases, qualities, ordinary tags, names and
detached mate fields are preserved. Inputs that cannot be preserved by this
encoder return `Invalid`: mapped records without CIGAR/reference/position,
unmapped records with CIGAR or nonzero MAPQ, empty names, reserved `cF` tags,
nonstandard double (`d`) tags, unpaired records with a mate target/reverse flag,
and zero-based positions ≥ 2,147,483,647. Header SQ entries must agree with the
binary dictionary. Other failures poison writes until close/reopen.

Reference-based encoding, embedded references, entropy-codec selection,
mate-chain/name compression and CRAI creation are not implemented. This first
encoder does **not** claim htslib compression-ratio or encoding-speed parity.
The [writer benchmarks](bench/README.md) report measured speed and output size.

## Development

```sh
python3 -m venv .venv
.venv/bin/pip install -r tests/requirements.txt
.venv/bin/python tests/run.py
python3 tests/native_simd.py --sanitize
# On the corresponding hardware; these fail on an unsupported host:
python3 tests/native_simd.py --require-isa neon
python3 tests/native_simd.py --require-isa avx512
scripts/docs.sh
scripts/docs.sh --check
scripts/doctest.sh
```

`GOOSE` selects the compiler. `GOOSE_DOC` names goose-doc's `main.goose`, default
`../goose-doc/main.goose`; documentation was checked with goose-doc
[`bb399143`](https://github.com/sstadick/goose-doc/commit/bb399143866c0c42268f23c8dfda7a447173f93d).
Source comments generate `API.md`; `api.txt` selects documented public names,
not language-enforced visibility. Doctests generate their own small BAM fixture.

Tests independently encode BAM/BGZF with Python `struct`/`zlib`, then compare
per-read pileups against htslib through pinned pysam. They exercise cross-block
records, malformed lengths, truncation, CRCs, empty input, auxiliary data,
65,536-operation CIGARs, clipping, padding, adjacent insertions, D/N, missing
sequence/quality, flags/MAPQ, multi-reference input, unsorted input, memory
limits, and region boundaries through sequential, BAI and CSI paths. Tests also
cover seek/EOF recovery, index truncation, overlapping chunks, pseudo-bins,
coarse bins, CSI beyond 512 Mbp, skipping unrelated BGZF blocks, and mate-link
eligibility/ambiguity/reuse. htslib BAQ, overlap correction, and base-quality
filtering are disabled to match the raw-column contract.

Local validation covers Linux x86-64, TinyCC at both Goose inlining levels,
Clang debug and release, and AddressSanitizer/UndefinedBehaviorSanitizer.
CRAM tests compare fields and per-read pileup entries with htslib for external,
embedded, no-reference, multi-reference and multi-slice inputs. Additional
validation uses Seqair's fidelity/conformance corpus and 189 official codec
vectors, with byte comparisons and sanitizer runs. These do not establish
complete CRAM conformance or cross-platform performance parity. Goose now
matches or beats Seqair on the five measured CRAM scans; see the benchmark
checkpoints for timings and remaining implementation differences.
The shared SIMD tests, 189 codec vectors and 2,048 native rANS cases also pass
in a Clang AArch64 build under QEMU. This is functional NEON coverage, not
native hardware performance. AVX-512 compiles but has not executed here;
`--require-isa` prevents a CI host from silently substituting another path.
Native ARM/AVX-512 runs and large real-world datasets remain pending. Benchmarks
explicitly compile native code with **Clang `-O3 -march=native`**, including I/O and libdeflate;
Goose's `-O2` controls frontend inlining only.

## Provenance and license

MIT; see [LICENSE](LICENSE). BAM/BGZF parsing implements the
[format specification](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf).
CRAM reconstruction/codecs and region storage follow Seqair, but this is not
a complete Seqair port. Reference source was inspected at
[`8401961d`](https://github.com/Softleif/seqair/tree/8401961de45883d4772226d91da95c9965104d59).
The test fixture generator is original and does not redistribute upstream BAMs.
