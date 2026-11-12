# Native benchmarks

## BAM writer checkpoint

Add `--writer --seqair` to the command below to benchmark raw-record writing
at libdeflate levels 0, 1, 6 and 12. Both implementations preload validated
records outside the timer, then time creation, header, writes and checked
finish to a regular file (no fsync). Each iteration reopens/truncates the file;
iteration zero warms up. Output is independently decompressed and compared
byte-for-byte with the original BAM after every sample. This isolates writing,
not BAM-to-BAM transcoding or serialization from decoded slab records.

Measured on the same local machine, CPU 4, five alternating rounds of ten
iterations, using upstream `tests/data/test.bam` (29,358 records):

| Level | Goose ms | Seqair ms | Output bytes (both) |
| ---: | ---: | ---: | ---: |
| 0 | 12.861 | 13.193 | 7,310,092 |
| 1 | 26.195 | 26.680 | 1,742,856 |
| 6 | 71.599 | 72.447 | 1,585,426 |
| 12 | 1,721.581 | 1,812.566 | 1,478,436 |

This is write parity on this small warm-cache dataset, not a universal speed
claim. The same run measured Goose reads at 6.654 ms, sequential region load
at 0.316 ms and columns at 0.100 ms, versus the pre-writer checkpoint's
6.54/0.313/0.098 ms. No material read/pileup regression was observed.

## SAM and CRAM writer checkpoints

`--format-writers --seqair` measures SAM formatting against Seqair and the new
reference-free CRAM encoder at levels 0, 1, 6 and 12. It requires pysam for
untimed, independent decoding. Records are preloaded outside timing; each
iteration includes a fresh writer, header, encoding, file writes and checked
close, without fsync. The final runner alternates SAM executable order and
checks output after every sample. SAM bytes must agree exactly. CRAM fields
must match the input after documented CIGAR normalization (`=/X` → `M`, zero
operations removed, adjacent equal operations merged).

AMD Ryzen 9 3950X, CPU 4, three rounds × three timed iterations with one warmup
per sample. Clang **`-O3 -march=native`**, including native adapters and
libdeflate; Rust opt-level 3, target-cpu native, linker Clang. Input is the same
29,358-record Seqair `tests/data/test.bam` (SHA-256
`90d52f3cf8a119ddd0b3cabfbc4e98e231e98847fb12cc12bdf130e5fa0a60be`).
Milliseconds per complete output:

| Writer | First implementation | Bulk quality copy | Final accounting | Output bytes |
| --- | ---: | ---: | ---: | ---: |
| Goose SAM | 15.626 | 15.242 | 15.799 | 9,019,872 |
| Seqair SAM | 18.768 | 18.713 | 18.615 | 9,019,872 |
| Goose CRAM, level 0 | 18.214 | 14.988 | 13.610 | 7,831,571 |
| Goose CRAM, level 1 | 27.424 | 25.431 | 23.961 | 1,278,551 |
| Goose CRAM, level 6 | 98.037 | 95.894 | 94.005 | 1,067,394 |
| Goose CRAM, level 12 | 2,327.393 | 2,578.468 | 2,595.636 | 949,742 |

The second checkpoint replaces per-quality pushes with contiguous copying
for CRAM and a pre-sized ASCII conversion for SAM. The final checkpoint reuses
the conservative record-size bound for slice admission instead of scanning
26 stream lengths after every record; absent SEQ uses the CIGAR-derived length
before allocating dummy bases/qualities. Output sizes are unchanged on this
fixture. These are separate runs, not interleaved before/after binaries; small
differences include machine variation. Level 12 became slower across runs,
and its size saving costs substantially more time than levels 1 or 6.

**Seqair has no CRAM writer at this revision.** There is no matched upstream
CRAM encoding-time claim, nor an htslib compression-ratio/speed comparison.
This encoder uses raw/gzip streams, no external reference and no entropy-codec
selection. These warm-cache numbers do not establish large-file, cold-I/O,
multi-threaded or cross-platform performance. Independent tests decode both
CRAM versions with htslib, Goose and Seqair; they also cover long CIGARs, missing
SEQ/QUAL, mates, tags (including exact NaN payloads), error/close paths and limits.

```sh
taskset -c 4 python3 bench/run.py input.bam chr19 6105700 6106700 \
  --seqair --format-writers --rounds 3 --iterations 3 --output writers.json
```

## SAM checkpoints

`--sam input.sam --seqair` adds a plain-text SAM scan. Both readers count every
record and its sequence length, including header/open costs and excluding the
warm-up iteration. Use a SAM generated with `samtools view -h --no-PG` from the
same BAM; the JSON records its SHA-256. Correctness tests separately compare all
record fields, auxiliary bytes and individual pileup entries with htslib.

On the dataset/affinity/rounds above (29,358 records, 2,310,244 bases):

| Checkpoint | Goose ms | Seqair ms |
| --- | ---: | ---: |
| Initial reusable SAM encoder | 19.204 | 15.540 |
| Table lookup, avoid repeated ordinary-record validation | 18.087 | 15.330 |
| Direct native pair packing, no intermediate decoded sequence | 15.689 | 15.401 |

The final checkpoint is within 2% on this workload. Parsing, record metadata,
slab storage and pileups remain Goose; the small base-packing loop is C, like
the existing base-unpacking loop, compiled with Clang `-O3 -march=native`.

## Initial CRAM parity checkpoint

Both consumers use Seqair-style streaming reconstruction and produce validated
packed BAM records. Timings include opening the CRAM/FASTA/index, decoding all
containers/slices, resolving mates and serializing records. They exclude
output-file writing. These are fresh-reader, single-threaded, warm-cache scans,
not indexed CRAM queries or pileup timings. CPU 4; five alternating rounds of
ten iterations, excluding one warm-up iteration per sample.

| Input | Records | Goose before this pass ms | Goose after this pass ms | Seqair ms |
| --- | ---: | ---: | ---: | ---: |
| `tests/data/test_v30.cram` | 29,358 | 34.982 | 29.636 | 42.645 |
| `hts-specs/cram/3.1/passed/level-1.cram` | 20,000 | 28.949 | 20.471 | 25.019 |
| `hts-specs/cram/3.1/passed/level-2.cram` | 20,000 | 37.543 | 24.122 | 26.828 |
| `hts-specs/cram/3.1/passed/level-3.cram` | 20,000 | 64.521 | 58.929 | 60.542 |
| `hts-specs/cram/3.1/passed/level-4.cram` | 20,000 | 67.591 | 63.271 | 66.774 |

Paths are relative to pinned Seqair's `tests/data`, except the first row,
which includes that prefix. All three executables were measured together;
both Goose versions include reference MD5 validation and bounded reference
windows. The before executable is the reader at the start of this parity pass,
not the earlier exploratory implementation.

This pass carried over Seqair's performance paths, not just its decoding results:

- Resolve external data-series descriptors once per slice; shared content IDs
  still share one cursor, and bindings reset for every slice.
- Pack each rANS lookup into one symbol/frequency/bias word. Use the actual
  precision as row stride, with reusable scratch and the unpacked fallback for
  non-exact tables or absent initial context.
- Keep four-state scalar loops specialized. Use runtime-dispatched AVX2 for
  32 states, eight-lane shuffle renormalization, 32×32 order-1 output tiles,
  and checked scalar tails. The native kernel does not call htscodecs.
- Expand short RLE runs/literal stretches with fixed 16-byte stores into
  allocated slack; expand PACK bytes through 256-entry tables above Seqair's
  1,024-byte threshold. Framing, tables and transforms remain Goose.

**That checkpoint matched or beat Seqair on these five scans:** about 0.69× its
time on 3.0 and 0.82× / 0.90× / 0.97× / 0.95× on the four 3.1 files. This is
not full implementation or cross-platform parity. Other CPUs, cold I/O and
large real-world CRAMs remain unmeasured. At that checkpoint NEON/AVX-512 and
Seqair's adaptive range-coder SIMD search were not ported, and 4-bit PACK used
byte-pair tables. The follow-up below addresses those implementation gaps.

Untimed validation checked byte-identical length-prefixed BAM payloads against
Seqair for 15 CRAM files, plus 189 official codec vectors. ASan/UBSan/leak checks
passed, including 15,120 mutated/truncated vectors and 1,454 targeted lane,
precision, renormalization, tile and PACK boundary cases. Mutation verdicts and
accepted bytes matched the before decoder exactly; forced scalar dispatch also
matched AVX2 on the official vectors and boundary cases. The timed consumers
check record/base counts, not byte identity.

## Shared SIMD follow-up

The sibling `../goose-simd` checkout now supplies the unsigned prefix-sum,
comparison/reduction and masked-expansion operations, plus a shared nibble
decoder for BAM (high-first) and CRAM PACK (low-first). Arithmetic coding
follows Seqair's four-entry scalar prefix and depth threshold of 16 × 32;
FQZ and run models retain scalar searches. A range-search-only incremental
check measured level 4 at 62.764 ms versus Seqair's 65.999 ms.

First-pass rANS kernels follow Seqair's ISA-specific structure: AVX-512 uses
16-lane widening/expansion renormalization; NEON uses four-lane expansion in
16-state batches, vector narrowing/next-index calculation and 32×32 output
tiles. Both keep the existing normalized-state/headroom guards and checked
scalar tails. Four-state decoding stays scalar. Nibble decoding includes
NEON and AVX-512 shuffles/interleaves and overlapping final vectors.

Final rerun on the same Ryzen 9 3950X/CPU 4, five alternating rounds × ten
iterations, fresh readers and warm cache, with no concurrent validation:

| Input | Goose ms | Seqair ms |
| --- | ---: | ---: |
| `test_v30.cram` | 30.305 | 42.515 |
| 3.1 level 1 | 20.778 | 24.802 |
| 3.1 level 2 | 23.644 | 26.663 |
| 3.1 level 3 | 58.129 | 59.793 |
| 3.1 level 4 | 63.610 | 66.373 |

This retains measured parity (0.71–0.97× Seqair's time); it is not a claim
that every new path speeds up this AVX2 machine. A BAM/pileup rerun of
`test.bam`, chr19:6105700–6106700, measured Goose read/load/columns at
6.708/0.309/0.098 ms versus Seqair sequential 8.342/0.846/0.153 ms.
Indexed load/columns were 0.293/0.098 versus 0.307/0.153 ms. Counts and
checksums matched in every round.

Validation after integration: byte-identical BAM payloads against Seqair on
15 CRAM files, 189 official codec vectors under ASan/UBSan, the public-reader
and pileup suite, and 2,048 independent native rANS decode/consumption/tail
cases in scalar and AVX2 builds. Clang AArch64 builds passed the 189 vectors,
2,048 cases and shared SIMD tests under QEMU. **NEON has functional emulation
coverage, not native timing; AVX-512 is compile-checked but unexecuted.**
Use `python3 tests/native_simd.py --require-isa neon` or `--require-isa avx512`
on supported CI hardware; requesting the wrong ISA fails instead of skipping.
The shared library's CI has an ARM runner and a configurable AVX-512 runner
entry, but neither new CI configuration has run yet.

All current native measurements use **Clang `-O3 -march=native`**, including
libdeflate, bzip2 1.0.8 and XZ 5.2.5. The latter two were built as shared
libraries and selected with `LD_LIBRARY_PATH`; actual loaded paths were checked.
Early exploratory level-3/4 runs used system codecs of unknown build provenance
and are deliberately omitted from this table. Rust uses opt-level 3,
target-cpu=native and the Clang linker; native Cargo dependencies use the same
C flags. Goose's `-O2` below is frontend inlining, not native optimization.

The CRAM benchmark currently needs the sibling `../goose-bio` checkout with
its indexed FASTA reader and the unreleased `../goose-simd` changes.
Build the pinned Seqair consumer, then run:

```sh
CC=clang CFLAGS='-O3 -march=native' \
RUSTFLAGS='-C opt-level=3 -C target-cpu=native -C linker=clang' \
cargo +1.98.1 build --locked --release --manifest-path bench/seqair/Cargo.toml --bin seqair-cram

# For bzip2/XZ inputs, point this at libbz2.so.1.0 and liblzma.so.5 built
# with Clang -O3 -march=native, not unverified system binaries.
export LD_LIBRARY_PATH=/path/to/clang-native-codecs
taskset -c 4 python3 bench/cram.py input.cram reference.fa \
  --seqair-bin bench/seqair/target/release/seqair-cram --output cram.json
```

`bench/cram.py` builds Goose, its I/O/FASTA adapters and libdeflate with the
explicit native flags, then saves input hashes, compiler details, affinity,
library search path, counts and all timing samples. The FASTA needs a matching
`.fai`. `seqair-cram input.cram reference.fa 0 records.raw` additionally dumps
length-prefixed BAM payloads as a byte-fidelity oracle.

## Run the benchmarks

Run from the repository root with Python 3.11+, Goose, CMake, and Clang:

```sh
python3 bench/run.py input.bam chr19 6105700 6106700 --output results.json
```

Add `--index input.bam.bai` (or a CSI path) to measure Goose indexed loading too.
Add `--seqair` to compare both sequential and indexed paths in Goose and the
pinned Rust Seqair. This needs rustup/Rust 1.98.1 and a BAM index beside the input;
the runner finds the index automatically unless `--index` names one explicitly.
It passes that exact index to both implementations and records its hash for
both; Rust does not independently choose a different adjacent CSI/BAI:

```sh
rustup toolchain install 1.98.1 --profile minimal
python3 bench/run.py input.bam chr19 6105700 6106700 --seqair --output comparison.json
```

The Rust consumer is in [seqair/main.rs](seqair/main.rs), with its upstream
revision and transitive dependencies pinned by its Cargo manifest/lockfile.
The runner sets Rust **`-C opt-level=3 -C target-cpu=native -C linker=clang`**
and native **`CC=clang CFLAGS="-O3 -march=native"`** explicitly.
No upstream source is modified. Both use libdeflate 1.26. Cargo's build cache
defaults to `bench/seqair/target`; `CARGO_TARGET_DIR` can move it elsewhere.

The script emits standalone C, compiles it, the seek/base/SAM-format adapters and goose-io's native code with
**Clang `-O3 -march=native`**, and independently builds the pinned libdeflate
with **Clang `-O3 -march=native -DNDEBUG`**. The runner accepts only Clang.
Goose's `-O2` setting is frontend inlining, not the native optimization level.
It does not use `goose build` for timings, because that command currently fixes
the native optimization level at `-O2`. Actual native commands and codec build
flags are printed; the JSON records compiler versions, commands, input hash,
per-round samples, counts, and medians. These benchmark binaries target the
build machine's CPU and are not portable release binaries. The runner
currently targets Linux and macOS Clang command lines, not MSVC.

Each sample performs one untimed warmup followed by ten timed iterations;
seven samples are collected by default. Override `--iterations` and `--rounds`.
All builds finish before timing; program order reverses every other round.
Times exclude compilation and process startup. The three measurements are:

- **Read:** open the BAM and sequentially read all records, validating framing,
  checksums, record layout, and CIGARs; sum sequence lengths without expanding
  every packed base.
- **Load:** reopen and load the requested region into reusable slabs, including
  mate linking. Sequential cases scan from the start. Indexed cases parse the
  index and seek to candidate chunks: Goose uses `load_indexed`, Seqair uses
  `IndexedBamReader::fetch_into_customized`. Index parsing is charged to both,
  even though the Goose API allows callers to reuse a parsed index/open reader.
- **Columns:** traverse the loaded store, generating every entry and summing
  position, record index, query position, base, quality, operation, and indel.
  Parsing and decompression are outside this measurement.

Seqair's sequential adapter uses `BamStreamReader::read_record`, the same
sorted-prefix stopping/filtering rules as Goose, and `RecordStore::push_raw`.
This bridge parses retained records again. Seqair's normal indexed path avoids
that bridge. Both indexed paths borrow records inside a decoded BGZF block and
copy only cross-block records. Both reject non-overlaps using raw core fields
and CIGAR spans before full decoding. Goose still validates retained records,
including auxiliary fields and resolved long CIGARs; indexed loading is not a
whole-file validation pass. Both stores eagerly expand bases into separate slabs alongside
names, native CIGAR words, qualities, and auxiliary bytes. Goose's cursors cache
slab offsets; columns do not reconstruct packed records or unpack bases.
The base-unpacking kernel now uses goose-simd's shared nibble decoder, with
runtime ISA selection and checked buffer lengths. It writes the same
0–15 IUPAC codes as Seqair's SIMD decoder, with no overread or padded tail.
Both link mates using names and reciprocal positions before columns are timed.
Goose uses its standard dictionary over borrowed name slices; Seqair caches a
name hash and additionally builds interval metadata. Its required
`prepare_for_pileup` is charged to load. **Load + columns** is the useful
end-to-end comparison; sharing the storage model does not make the
implementations feature-equivalent.

All cases use one decoding thread, exclude flags `0x704`, apply no MAPQ/base
quality threshold or depth cap, and do no BAQ or overlap correction. The same
half-open region is translated to Seqair's inclusive-end convention. The read
measurement always scans the whole file, even in the indexed-labeled run.

Counts and checksums must match across all implementations and rounds. The
checksum consumes the fields named above, so this measures real per-entry work,
not just depth. The Rust `dump`/`indexed-dump` modes emit per-read results for
untimed exact comparison with `tests/probe.goose`; `mates`/`indexed-mates`
emit every mate index and linking statistics. Indexed modes accept an explicit
index as the final argument. Its base accessor collapses ambiguous IUPAC
codes to N; the measured regions contain only A/C/G/T/N at contributing bases.
Arbitrary inputs with other IUPAC codes may fail comparison for that semantic
difference. No output formatting or terminal I/O is inside timed sections.

## Goose versus Seqair, Linux x86-64

Measured on an AMD Ryzen 9 3950X, pinned to logical CPU 4 with `taskset -c 4`,
warm filesystem cache, nine rounds × thirty iterations. Goose `74c641c5` uses
Clang 18.1.3 at `-O3 -march=native`. Seqair
[`8401961d`](https://github.com/Softleif/seqair/tree/8401961de45883d4772226d91da95c9965104d59)
uses Rust 1.98.1 at `-C opt-level=3 -C target-cpu=native`, with Clang as its
linker and its native codec built by Clang at `-O3 -march=native`.
Input is Seqair's
[`tests/data/test.bam`](https://github.com/Softleif/seqair/blob/8401961de45883d4772226d91da95c9965104d59/tests/data/test.bam)
(1,581,777 compressed bytes, 29,358 records), SHA-256
`90d52f3cf8a119ddd0b3cabfbc4e98e231e98847fb12cc12bdf130e5fa0a60be`.

**Whole-file sequential reading:** 8.678 ms for Rust Seqair and 6.591 ms for
Goose/Clang (sequential-adapter medians from the early-window run).
Goose/Clang takes 0.76× Seqair's time on this input.

Region measurements, milliseconds per iteration. Total is the median of each
sample's load + columns, not the sum of independently computed medians.

| chr19 region | Implementation | Load | Columns | Total |
| --- | --- | ---: | ---: | ---: |
| `[6105700,6106700)` | Seqair sequential | 0.864 | 0.152 | 1.016 |
| 233 reads / 16,506 entries | Seqair indexed | 0.311 | 0.152 | 0.463 |
| | Goose/Clang sequential | 0.316 | 0.099 | 0.415 |
| | Goose/Clang indexed | 0.295 | 0.099 | 0.394 |
| `[6100000,6144000)` | Seqair sequential | 4.394 | 8.527 | 12.912 |
| 12,273 reads / 968,169 entries | Seqair indexed | 3.614 | 8.551 | 12.157 |
| | Goose/Clang sequential | 3.667 | 5.880 | 9.543 |
| | Goose/Clang indexed | 3.720 | 5.885 | 9.606 |
| `[6142000,6143000)` | Seqair sequential | 4.353 | 0.212 | 4.565 |
| 333 reads / 23,135 entries | Seqair indexed | 1.073 | 0.211 | 1.283 |
| | Goose/Clang sequential | 2.839 | 0.142 | 2.983 |
| | Goose/Clang indexed | 1.058 | 0.143 | 1.201 |
| `[6142000,6142001)` | Seqair sequential | 4.251 | 0.002 | 4.253 |
| 27 reads / 27 entries | Seqair indexed | 0.940 | 0.002 | 0.942 |
| | Goose/Clang sequential | 2.701 | <0.001 | 2.702 |
| | Goose/Clang indexed | 0.924 | <0.001 | 0.925 |

Goose indexed load + columns takes **0.85× / 0.79× / 0.94× / 0.98×** Seqair
indexed time on the early/wide/late/one-base windows. The one-base result is
near parity, not a meaningful lead. The late query is **2.48× faster** than
Goose's sequential path in the same run. Indexing does not help the wide
window here: its candidate ranges cover much of the same prefix, and index
parsing has a cost. These results do not establish parity on all workloads or
imply that mate linking is free; it is included in every load.

A separate paired control interleaved the saved pre-optimization binary, the
current binary and Seqair, with the same native flags, BAI, CPU and nine rounds
of thirty iterations. It isolates the change from run-to-run drift:

| Indexed window | Goose before | Goose after | Seqair |
| --- | ---: | ---: | ---: |
| Early | 0.436 | 0.393 | 0.461 |
| Wide | 9.862 | 9.560 | 12.132 |
| Late | 1.403 | 1.197 | 1.279 |
| One base | 1.121 | 0.921 | 0.945 |

Values are load + columns in milliseconds. The late and one-base cases take
15% and 18% less time than before. The changes borrow within-block records,
reduce redundant byte-load bounds checks, expand bases directly into the slab
with a vectorizable native loop, and reject indexed non-overlaps before full
decoding. This preserves Seqair's decoded-slab/forward-cursor model; it does
not bring back packed-record storage. Sequential reading still validates every
record; indexed loading fully validates only retained candidates.

All four paths produced equal retained-read counts, entry counts,
and checksums in every round. Untimed comparisons additionally checked every
per-read position/name/query-position/op/indel/base/quality row for the first
three real regions, 100-base windows on both synthetic inputs below, and cases
drawn from Seqair's indel and mate-linking tests. Both BAI and CSI matched:
**2,167,610 identical rows per executable**, plus every stored mate index and
linking statistic. Real-region pair counts were 93 / 6,053 / 124, with zero
ambiguous arrivals. These comparisons passed in optimized and ASan/UBSan Goose
builds. A malformed auxiliary tag on a discarded non-overlap also confirmed
Seqair's staged-validation behavior; Goose still rejects it when retained.

Reproduce each row with the same BAI and its begin/end:

```sh
taskset -c 4 python3 bench/run.py input.bam chr19 6142000 6143000 \
  --index input.bam.bai --seqair --rounds 9 --iterations 30 --output comparison.json
```

These are **fresh-reader, warm-cache results**, not general performance parity,
a disk-throughput test, or a matched htslib performance comparison. Large real
BAMs, cold I/O, and persistent indexed-reader reuse remain unmeasured. Goose
does not yet have Seqair's persistent distant-bin/block/resume caching; these
timings cannot establish parity for repeated adjacent queries. The source BAM
is not redistributed here.

## Seqair's DNA and spliced benchmark shapes

[`synthetic.py`](synthetic.py) adapts
[`synthetic_input`](https://github.com/Softleif/seqair/blob/8401961de45883d4772226d91da95c9965104d59/crates/seqair/benches/bam.rs)
with the same seed, RNG draw order, 200,000 sorted starts over 1 Mb, CIGARs,
bases, qualities, names and flags. DNA is mostly 150M, with short insertions
and deletions. Spliced reads contain 150 aligned bases and one or two introns
of 100–4,999 bases. Unlike the upstream in-memory benchmark, this encodes a
BAM and BAI/CSI so both implementations read exactly the same input. We query
`chrA:[500000,510000)`, not the upstream benchmark's whole-contig interval.

Same CPU/compiler settings as above; seven rounds × ten iterations, warm
cache, fresh readers and the same BAI. Milliseconds per iteration:

| Shape | Implementation | Load | Columns | Total |
| --- | --- | ---: | ---: | ---: |
| DNA: 2,050 reads / 302,652 entries | Seqair indexed | 2.725 | 2.495 | 5.216 |
| | Goose indexed | 2.622 | 1.767 | 4.389 |
| Spliced: 2,739 reads / 7,026,148 entries | Seqair indexed | 4.018 | 45.216 | 49.234 |
| | Goose indexed | 3.955 | 33.982 | 37.952 |

Goose takes **0.84× / 0.77×** Seqair's load + columns time. All four paths
(including sequential) matched counts/checksums in every round. These inputs
have no paired reads; mate-link correctness also uses the real BAM and crafted
paired cases above. Skips contribute entries, so the spliced case exercises
millions of forward-cursor steps, not just aligned bases.

With `tests/requirements.txt` installed in the active Python environment:

```sh
python3 bench/synthetic.py dna tests/data/dna.bam
python3 bench/synthetic.py spliced tests/data/spliced.bam
taskset -c 4 python3 bench/run.py tests/data/dna.bam chrA 500000 510000 \
  --seqair --index tests/data/dna.bam.bai --output dna.json
taskset -c 4 python3 bench/run.py tests/data/spliced.bam chrA 500000 510000 \
  --seqair --index tests/data/spliced.bam.bai --output spliced.json
```

Measured BAM SHA-256 hashes (pysam 0.24.1 / htslib 1.24):

- DNA, 31,175,720 bytes: `003b878ce7bee6173de41235fc2bcbb142616b9344705a3de57913a3124320bf`.
- Spliced, 33,188,370 bytes: `6107075a28342b0d287533ab0adb13618d0af5bb217d6049e45f113dc03e0c81`.

## Historical packed-store synthetic checks

These results predate the slab conversion above; they are not measurements of
the current store. Additional warm-cache checks used the same CPU, native flags, and seven rounds
of ten iterations, alternating the execution order of the old Goose binary,
the new Goose binary, and both Seqair paths. All counts/checksums matched.
The inputs use the independent encoder in `tests/fixtures.py`, repeating
ACGTN sequence, its position-dependent qualities, `NMi:2` and `RGZ:group` tags,
and unique names `run:1:flowcell:2:1101:{i:08d}`. They deliberately isolate
ordinary matches, high depth, and long reference skips; they are not substitutes
for representative real datasets.

| Workload | Records / start positions | CIGAR | chrA region | Entries |
| --- | --- | --- | --- | ---: |
| Short reads | 60,000 / `i * 3` | `150M` | `[100000,110000)` | 500,000 |
| High depth | 10,000 / `i // 10` | `3S70M2I40M3D35M` | `[400,700)` | 444,000 |
| Spliced reads | 10,000 / `i * 4` | `60M1000N90M` | `[20000,22000)` | 575,000 |

Times below are milliseconds; Seqair columns use its sequential path.

| Workload | Goose read | Seqair read | Goose columns before | Goose columns after | Seqair columns |
| --- | ---: | ---: | ---: | ---: | ---: |
| Short reads | 5.723 | 6.676 | 5.088 | 3.379 | 3.990 |
| High depth | 1.082 | 2.066 | 4.936 | 3.894 | 4.747 |
| Spliced reads | 1.060 | 2.055 | 5.201 | 3.342 | 3.858 |

## Standalone BGZF extraction check

On 2026-11-12, the private BGZF reader/writer moved into
[goose-bgzf v0.1.0](https://github.com/sstadick/goose-bgzf/releases/tag/v0.1.0).
The standalone package records its own
[read/write throughput baseline](https://github.com/sstadick/goose-bgzf/blob/main/bench/README.md).
The native benchmark scripts now link the seek adapter from the exact resolved
BGZF dependency instead of the deleted local copy.

A warm-cache integration smoke measurement on the same Ryzen 9 3950X used
Clang 18.1.3 `-O3 -march=native` for all timed C and libdeflate, three rounds
of three iterations, upstream `test.bam`, and chr19:[0,61431566). Milliseconds:

| BGZF implementation | BAM read | Region load | Columns |
| --- | ---: | ---: | ---: |
| Private copy | 7.154 | 3.787 | 5.828 |
| Published package | 7.347 | 3.918 | 5.995 |

Both produced 29,358 input records, 12,273 retained records, 968,169 entries
and checksum 5,934,549,274,535. These short, unpinned runs were not alternated
and overlapped other validation, so the roughly 3% timing difference does not
establish a regression or equivalence. The extraction adds no new decoding
algorithm. Full debug/release interoperability tests also pass with the released
package; no upstream Seqair performance claim is inferred from this check.
