// Read a BAM index and find file ranges that may contain alignments for a
// genomic region. This module reads the index, not the BAM records themselves.
namespace bam_index;
import std;
import os;
import bam;
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";
import io from "github.com/sstadick/goose-io@v0.3.0";

// A candidate range of records in a BAM file. begin is inclusive and end is
// exclusive. These are BGZF virtual offsets: the high 48 bits locate a
// compressed block in the file, and the low 16 bits locate a byte inside the
// decompressed block. They are not genomic coordinates or byte distances.
//
// Call reader.seek(chunk.begin), then read records while reader.tell() is less
// than chunk.end. Check each record's reference ID and alignment coordinates:
// a chunk may include records outside the genomic region you requested.
// A chunk describes where to read; it does not hold those records in memory.
//
// ```goose
// let chunk = bam_index::Chunk { begin: 100, end: 200 };
// assert(chunk.begin < chunk.end);
// ```
struct Chunk {
    begin: u64,
    end:   u64,
}

struct Bin {
    id:      u32,
    loffset: u64,
    chunks:  i64,
    count:   i64,
}

struct Reference {
    bins:         i64,
    count:        i64,
    linear:       i64,
    linear_count: i64,
}

// An in-memory BAM index, reusable across region queries. open reads a .bai or
// .csi file and closes it after parsing. Both uncompressed and BGZF-compressed
// index files are supported. Pass the index path, not the BAM path.
//
// max_bytes limits the decompressed index input to 256 MiB by default; parsed
// tables and temporary buffers require additional memory. Set this limit
// before open and otherwise treat the fields as read-only.
// BAI supports positions below 2^29. CSI can cover longer references; accepted
// CSI headers have depth 0..10 and min_shift + 3*depth at most 62.
//
// ```goose
// let index = bam_index::Index { .. };
// assert(index.open("tests/data/example.bam.bai") == bam::Status.Ok);
// assert(index.references.len == 2);
// ```
root struct Index {
    references: Reference[>..] = [],
    bins:       Bin[>..]       = [],
    chunks:     Chunk[>..]     = [],
    linear:     u64[>..]       = [],
    min_shift:  i64            = 14,
    depth:      i64            = 5,
    csi:        bool           = false,
    status:     bam::Status    = bam::Status.Closed,
    max_bytes:  i64            = 268435456,
}

struct Parser { bytes: const u8[:], at: i64 = 0, status: bam::Status = bam::Status.Ok }
fn take(parser: Parser&, size: i64) -> const u8[:] {
    if parser.status != bam::Status.Ok { return parser.bytes[..0]; }
    if size < 0 || size > parser.bytes.len - parser.at {
        parser.status = bam::Status.Truncated;
        return parser.bytes[..0];
    }
    let start = parser.at;
    parser.at += size;
    parser.bytes[start..parser.at]
}
fn word(parser: Parser&) -> u32 {
    let bytes = parser.take(4);
    if bytes.len == 4 { bgzf::le32(bytes, 0) } else { 0 }
}
fn wide(parser: Parser&) -> u64 {
    let bytes = parser.take(8);
    if bytes.len == 8 { bgzf::le32(bytes, 0) | ((bgzf::le32(bytes, 4) as u64) << 32) } else { 0 }
}
fn count(parser: Parser&, width: i64) -> i64 {
    let n = parser.word() as! i32;
    if n < 0 { parser.status = bam::Status.Invalid; return 0; }
    if n > (parser.bytes.len - parser.at) / width { parser.status = bam::Status.Truncated; return 0; }
    n
}
fn clear(index: Index&) {
    index.references.clear(); index.bins.clear(); index.chunks.clear(); index.linear.clear();
}

// Replace the index, returning Ok on success. Failure clears any previous data.
// The index must belong to the BAM used for loading; the formats have no file
// identity checksum, so a stale index cannot always be detected.
//
// ```goose
// let index = bam_index::Index { .. };
// assert(index.open("tests/data/example.bam.bai") == bam::Status.Ok);
// assert(index.open("") == bam::Status.IoError);
// assert(index.references.len == 0);
// ```
fn open(index: Index&, path: const u8[:]) -> bam::Status {
    index.clear();
    var bytes: u8[>..] = [];
    index.status = read_bytes(path, bytes, index.max_bytes);
    if index.status == bam::Status.Ok { index.status = parse(index, bytes); }
    if index.status != bam::Status.Ok { index.clear(); }
    index.status
}

// Discover an adjacent index without opening the BAM. Try path + ".csi",
// then replace a trailing ".bam" with ".csi", then the same two BAI names.
// The first existing file wins; parsing or reading errors do not fall back.
// Returns Ok on success; no index is created. Missing indexes return IoError
// and clear old data.
// Use open(index_path) to select an explicit index instead.
//
// ```goose
// let index = bam_index::Index { .. };
// index.open_for("tests/data/example.bam").check(bam::Status.Ok);
// assert(index.references.len == 2);
// ```
fn open_for(index: Index&, bam_path: const u8[:]) -> bam::Status {
    for suffix in [".csi", ".bai"] {
        let appended = str(bam_path, suffix);
        if file_exists(appended) { return index.open(appended); }
        if bam_path.len >= 4 && bam_path[bam_path.len - 4..] == ".bam" {
            let replaced = str(bam_path[..bam_path.len - 4], suffix);
            if file_exists(replaced) { return index.open(replaced); }
        }
    }
    index.clear();
    index.status = bam::Status.IoError;
    index.status
}

fn read_bytes(path: const u8[:], bytes: u8[>..]&, limit: i64) -> bam::Status {
    let file, ok = io::open_read(path);
    if !ok { return bam::Status.IoError; }
    defer io::close(file);
    var buffer: u8[65536] = [0; 65536];
    var n = io::read(file, buffer);
    if n < 0 { return bam::Status.IoError; }
    if n >= 2 && buffer[0] == 31 && buffer[1] == 139 {
        io::close(file);
        var stream = bgzf::Reader { .. };
        defer stream.close();
        let opened = stream.open(path);
        if opened != bgzf::Status.Ready { return bam::from_bgzf(opened); }
        while true {
            let status = stream.block_read();
            if status == bgzf::Status.End { break; }
            if status != bgzf::Status.Ready { return bam::from_bgzf(status); }
            if stream.end > limit - bytes.len { return bam::Status.TooLarge; }
            bytes.append(stream.decoded[..stream.end]);
        }
    } else {
        while n > 0 {
            if n > limit - bytes.len { return bam::Status.TooLarge; }
            bytes.append(buffer[..n]);
            n = io::read(file, buffer);
            if n < 0 { return bam::Status.IoError; }
        }
    }
    bam::Status.Ok
}

fn parse(index: Index&, bytes: const u8[:]) -> bam::Status {
    var parser = Parser { bytes: bytes, .. };
    let magic = parser.take(4);
    if parser.status != bam::Status.Ok { return parser.status; }
    index.csi = magic == "CSI\x01";
    if !index.csi && magic != "BAI\x01" { return bam::Status.Invalid; }
    index.min_shift = 14; index.depth = 5;
    if index.csi {
        index.min_shift = parser.word() as! i32;
        index.depth = parser.word() as! i32;
        if parser.status != bam::Status.Ok { return parser.status; }
        if index.min_shift < 0 || index.depth < 0 || index.depth > 10 ||
           index.min_shift + index.depth * 3 > 62 { return bam::Status.Invalid; }
        let aux = parser.count(1);
        parser.take(aux);
    }
    let refs = parser.count(4);
    if refs > 1048576 { return bam::Status.TooLarge; }
    let bin_limit = ((1 as u64) << (3 * (index.depth + 1))) / 7;
    for r in refs {
        let bins = parser.count(if index.csi { 16 } else { 8 });
        var ref = Reference { bins: index.bins.len, count: 0, linear: index.linear.len, linear_count: 0 };
        for b in bins {
            let id = parser.word();
            let loffset = if index.csi { parser.wide() } else { 0 as u64 };
            let n = parser.count(16);
            if parser.status != bam::Status.Ok { return parser.status; }
            if id == bin_limit + 1 {
                // Pseudo-bin: file span and mapped/unmapped counts, not chunks.
                if n != 2 { return bam::Status.Invalid; }
                parser.take(32);
                continue;
            }
            if id >= bin_limit { return bam::Status.Invalid; }
            let at = index.chunks.len;
            for c in n {
                let start = parser.wide();
                let end = parser.wide();
                if start > end { return bam::Status.Invalid; }
                index.chunks.push(Chunk { begin: start, end: end });
            }
            index.bins.push(Bin { id: id, loffset: loffset, chunks: at, count: n });
        }
        ref.count = index.bins.len - ref.bins;
        sort(index.bins[ref.bins..]) { a, b => a.id < b.id };
        if !index.csi {
            ref.linear_count = parser.count(8);
            for i in ref.linear_count { index.linear.push(parser.wide()); }
        }
        if parser.status != bam::Status.Ok { return parser.status; }
        index.references.push(ref);
    }
    if parser.status != bam::Status.Ok { return parser.status; }
    // The optional final u64 counts unplaced unmapped reads.
    let remaining = bytes.len - parser.at;
    if remaining != 0 && remaining != 8 { return bam::Status.Invalid; }
    bam::Status.Ok
}

fn lower_bound(bins: const Bin[:], id: u64) -> i64 {
    var lo: i64 = 0;
    var hi = bins.len;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if bins[mid].id < id { lo = mid + 1; } else { hi = mid; }
    }
    lo
}

// Find file ranges that may contain records overlapping one genomic region.
// tid is the reference's zero-based index in the BAM header, not a name or
// an offset into the header's name array. begin and end are zero-based genomic
// coordinates: begin is included and end is excluded.
//
// Replaces output with chunks sorted by file position, merging overlapping
// ranges so the same file range is not visited twice within this query.
// Separate queries can still return overlapping chunks. The chunks are only
// candidates: after reading, keep records with the requested tid and with
// record.fields.pos < end and record.fields.end > begin.
//
// Returns Ok on success, including an empty region with no chunks.
// Invalid reference IDs or coordinates outside the index's range fail.
// The index cannot determine whether its BAM file has changed since indexing.
//
// ```goose
// let index = bam_index::Index { .. };
// assert(index.open("tests/data/example.bam.bai") == bam::Status.Ok);
// var chunks: bam_index::Chunk[>..] = [];
// assert(index.query(0, 5, 7, chunks) == bam::Status.Ok && chunks.len > 0);
// assert(index.query(0, 7, 7, chunks) == bam::Status.Ok && chunks.len == 0);
// ```
fn query(index: const Index&, tid: i32, begin: i64, end: i64, output: Chunk[>..]&) -> bam::Status {
    output.clear();
    if index.status != bam::Status.Ok { return index.status; }
    let limit = (1 as u64) << (index.min_shift + index.depth * 3);
    if tid < 0 || tid >= index.references.len || begin < 0 || end < begin || end as u64 > limit {
        return bam::Status.Invalid;
    }
    if begin == end { return bam::Status.Ok; }
    let ref = index.references[tid as i64];
    let bins = index.bins[ref.bins..ref.bins + ref.count];
    var minimum: u64 = 0;
    if index.csi {
        var id = (((1 as u64) << (3 * index.depth)) - 1) / 7 + ((begin as u64) >> index.min_shift);
        while true {
            let at = lower_bound(bins, id);
            if at < bins.len && bins[at].id == id { minimum = bins[at].loffset; break; }
            if id == 0 { break; }
            id = (id - 1) / 8;
        }
    } else {
        let window = begin >> 14;
        if window < ref.linear_count { minimum = index.linear[ref.linear + window]; }
    }
    var offset: u64 = 0;
    for level in index.depth + 1 {
        let shift = index.min_shift + 3 * (index.depth - level);
        let first = offset + ((begin as u64) >> shift);
        let last = offset + (((end - 1) as u64) >> shift);
        var at = lower_bound(bins, first);
        while at < bins.len && bins[at].id <= last {
            let bin = bins[at];
            for chunk in index.chunks[bin.chunks..bin.chunks + bin.count] {
                if chunk.end > minimum && chunk.begin < chunk.end {
                    output.push(Chunk { begin: max(chunk.begin, minimum), end: chunk.end });
                }
            }
            at += 1;
        }
        offset += (1 as u64) << (3 * level);
    }
    sort(output) { a, b => a.begin < b.begin };
    var kept: i64 = 0;
    for chunk in output {
        if kept > 0 && chunk.begin <= output[kept - 1].end {
            output[kept - 1].end = max(output[kept - 1].end, chunk.end);
        } else { output[kept] = chunk; kept += 1; }
    }
    output.resize(kept);
    bam::Status.Ok
}

// Select and visit alignments overlapping [begin, end) on reference tid.
// Coordinates are zero-based and half-open. Reuses the open reader and index,
// including for backward queries or after EOF. Requires coordinate-sorted BAM
// and a matching index. Sorting is checked only in the visited records.
// Returns Ok on completion, including empty results; failures may follow
// earlier callbacks and do not undo their effects. Every visited record is
// fully decoded and validated. Unlike pileup, no flag/MAPQ filters are applied:
// positioned unmapped and zero-span records use BAM's minimum span of one.
// Each record borrows reader storage until the next reader operation. Do not
// read, seek, close or fetch on that reader inside the callback. An inline
// return exits the enclosing function and permits stopping early.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// reader.open("tests/data/example.bam").check(bam::Status.Ok);
// let index = bam_index::Index { .. };
// index.open_for("tests/data/example.bam").check(bam::Status.Ok);
// var count = 0;
// let status = index.fetch(reader, 0 as i32, 5, 7) { record =>
//     assert(record.name() == "match");
//     count += 1;
// };
// status.check(bam::Status.Ok);
// assert(count == 1);
// ```
fn fetch<F: fn(const bam::Record&)>(index: const Index&, reader: bam::Reader&,
        tid: i32, begin: i64, end: i64) -> bam::Status {
    if reader.stream.input.fd < 0 { return bam::Status.Closed; }
    if index.status != bam::Status.Ok { return index.status; }
    if index.references.len != reader.references.len || tid < 0 || tid >= reader.references.len ||
       begin < 0 || end < begin || end > reader.references[tid as i64].length {
        return bam::Status.Invalid;
    }
    var chunks: Chunk[>..] = [];
    let status = index.query(tid, begin, end, chunks);
    if status != bam::Status.Ok { return status; }
    var last_tid: i32 = -1;
    var last_pos: i64 = -1;
    for chunk, i in chunks {
        if i == 0 || reader.tell() < chunk.begin {
            let status = reader.seek(chunk.begin);
            if status != bam::Status.Ok {
                return if status == bam::Status.End { bam::Status.Truncated } else { status };
            }
        }
        while reader.tell() < chunk.end {
            let status = reader.read();
            if status != bam::Status.Record {
                return if status == bam::Status.End { bam::Status.Truncated } else { status };
            }
            let record = reader.record();
            let f = record.fields;
            if f.tid < 0 || f.pos < 0 { continue; }
            if f.tid < last_tid || (f.tid == last_tid && f.pos < last_pos) { return bam::Status.Unsorted; }
            last_tid = f.tid; last_pos = f.pos;
            if f.tid > tid || (f.tid == tid && f.pos >= end) { return bam::Status.Ok; }
            if f.tid == tid && f.end > begin { F(record); }
        }
    }
    bam::Status.Ok
}

// Fetch by reference name instead of ID. Unknown names return Invalid without
// invoking the callback. All other behavior is the same as the tid overload.
// Returns Ok on completion, including empty results.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// reader.open("tests/data/example.bam").check(bam::Status.Ok);
// let index = bam_index::Index { .. };
// index.open_for("tests/data/example.bam").check(bam::Status.Ok);
// let status = index.fetch(reader, "chrA", 5, 7) { record =>
//     assert(record.name() == "match");
// };
// status.check(bam::Status.Ok);
// ```
fn fetch<F: fn(const bam::Record&)>(index: const Index&, reader: bam::Reader&,
        name: const u8[:], begin: i64, end: i64) -> bam::Status {
    index.fetch(reader, reader.reference_id(name), begin, end, F)
}
