// BAI/CSI parsing and half-open region queries. No BAM decoding lives here.
namespace bam_index;
import std;
import bam;
import bgzf;
import io from "github.com/sstadick/goose-io@v0.2.2";

// A half-open file range in BAM virtual offsets, not ordinary byte offsets.
//
// ```goose
// let chunk = bam_index::Chunk { begin: 100, end: 200 };
// assert(chunk.begin < chunk.end);
// ```
struct Chunk { begin: u64, end: u64 }
struct Bin { id: u32, loffset: u64, chunks: i64, count: i64 }
struct Reference { bins: i64, count: i64, linear: i64, linear_count: i64 }

// Owned, reusable BAI or CSI index. open accepts an explicit index path;
// raw and BGZF-compressed indexes are supported. No file remains open.
// max_bytes limits decompressed input (default 256 MiB), not total memory.
// Treat fields other than max_bytes as read-only. CSI depth is 0..10 (u32
// bin IDs), with min_shift + 3*depth <= 62. BAI covers coordinates below 2^29.
//
// ```goose
// let index = bam_index::Index { .. };
// assert(index.open("tests/data/example.bam.bai") == bam::Status.Record);
// assert(index.references.len == 2);
// ```
root struct Index {
    references: Reference[>..] = [], bins: Bin[>..] = [],
    chunks:     Chunk[>..]     = [], linear: u64[>..] = [],
    min_shift:  i64            = 14, depth: i64 = 5, csi: bool = false,
    status:     bam::Status    = bam::Status.Closed,
    max_bytes:  i64            = 268435456,
}

struct Parser { bytes: const u8[:], at: i64 = 0, status: bam::Status = bam::Status.Record }
fn take(parser: Parser&, size: i64) -> const u8[:] {
    if parser.status != bam::Status.Record { return parser.bytes[..0]; }
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

// Replace the index. Failure clears any previous data; errors use bam::Status.
// The index must belong to the BAM used for loading; the formats have no file
// identity checksum, so a stale index cannot always be detected.
//
// ```goose
// let index = bam_index::Index { .. };
// assert(index.open("tests/data/example.bam.bai") == bam::Status.Record);
// assert(index.open("") == bam::Status.IoError);
// assert(index.references.len == 0);
// ```
fn open(index: Index&, path: const u8[:]) -> bam::Status {
    index.clear();
    var bytes: u8[>..] = [];
    index.status = read_bytes(path, bytes, index.max_bytes);
    if index.status == bam::Status.Record { index.status = parse(index, bytes); }
    if index.status != bam::Status.Record { index.clear(); }
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
        var stream = bgzf::Stream { .. };
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
    bam::Status.Record
}

fn parse(index: Index&, bytes: const u8[:]) -> bam::Status {
    var parser = Parser { bytes: bytes, .. };
    let magic = parser.take(4);
    if parser.status != bam::Status.Record { return parser.status; }
    index.csi = magic == "CSI\x01";
    if !index.csi && magic != "BAI\x01" { return bam::Status.Invalid; }
    index.min_shift = 14; index.depth = 5;
    if index.csi {
        index.min_shift = parser.word() as! i32;
        index.depth = parser.word() as! i32;
        if parser.status != bam::Status.Record { return parser.status; }
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
            if parser.status != bam::Status.Record { return parser.status; }
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
        if parser.status != bam::Status.Record { return parser.status; }
        index.references.push(ref);
    }
    if parser.status != bam::Status.Record { return parser.status; }
    // The optional final u64 counts unplaced unmapped reads.
    let remaining = bytes.len - parser.at;
    if remaining != 0 && remaining != 8 { return bam::Status.Invalid; }
    bam::Status.Record
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

// Replace output with sorted, merged candidate chunks for [begin, end).
// Candidates still require BAM overlap/filter checks. Empty regions succeed
// with no chunks; invalid references or out-of-range coordinates fail.
//
// ```goose
// let index = bam_index::Index { .. };
// assert(index.open("tests/data/example.bam.bai") == bam::Status.Record);
// var chunks: bam_index::Chunk[>..] = [];
// assert(index.query(0, 5, 7, chunks) == bam::Status.Record && chunks.len > 0);
// assert(index.query(0, 7, 7, chunks) == bam::Status.Record && chunks.len == 0);
// ```
fn query(index: const Index&, tid: i32, begin: i64, end: i64, output: Chunk[>..]&) -> bam::Status {
    output.clear();
    if index.status != bam::Status.Record { return index.status; }
    let limit = (1 as u64) << (index.min_shift + index.depth * 3);
    if tid < 0 || tid >= index.references.len || begin < 0 || end < begin || end as u64 > limit {
        return bam::Status.Invalid;
    }
    if begin == end { return bam::Status.Record; }
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
    bam::Status.Record
}
