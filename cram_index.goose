// CRAI slice lookup. Coordinates are normalized to zero-based half-open ranges.
namespace cram_index;
import std;
import bam;
import cram;
import pileup;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";

struct Entry { tid: i32, begin: i64, end: i64, container: i64, slice: i64, size: i64, max_end: i64 = 0 }

// Owned CRAI index. Plain text and concatenated gzip/BGZF members are accepted.
// max_bytes bounds compressed and decoded index bytes separately, not total
// allocations. No handle stays open. A CRAI has no file identity checksum.
//
// ```goose
// let index = cram_index::Index { .. };
// assert(index.status == bam::Status.Closed);
// assert(index.open("") == bam::Status.IoError);
// ```
root struct Index { entries: Entry[>..] = [], status: bam::Status = bam::Status.Closed, max_bytes: i64 = 268435456 }

fn read_bytes(path: const u8[:], bytes: u8[>..]&, limit: i64) -> bam::Status {
    let file, ok = io::open_read(path); if !ok { return bam::Status.IoError; }
    defer io::close(file);
    var raw: u8[>..] = []; var buffer: u8[65536] = [0; 65536];
    loop {
        let n = io::read(file, buffer); if n < 0 { return bam::Status.IoError; } if n == 0 { break; }
        if n > limit - raw.len { return bam::Status.TooLarge; } raw.append(buffer[..n]);
    }
    if raw.len < 2 || raw[0] != 31 || raw[1] != 139 { bytes.append(raw); return bam::Status.Record; }
    let decoder, allocated = deflate::decompressor(); if !allocated { return bam::Status.IoError; }
    defer decoder.close(); var at: i64 = 0; var decoded: u8[>..] = [];
    while at < raw.len {
        var size = min(65536, limit - bytes.len);
        loop {
            decoded.resize(size, 0);
            let result = decoder.decompress(deflate::Format.Gzip, raw[at..], decoded);
            if result.status == deflate::Status.Success {
                if result.consumed <= 0 { return bam::Status.Invalid; }
                bytes.append(decoded[..result.written]); at += result.consumed; break;
            }
            if result.status != deflate::Status.OutputTooSmall { return bam::Status.Invalid; }
            if size == limit - bytes.len { return bam::Status.TooLarge; }
            size = min(limit - bytes.len, size * 2);
        }
    }
    bam::Status.Record
}

fn parse(index: Index&, bytes: const u8[:]) -> bam::Status {
    var valid = true;
    bstr::split_str(bytes, "\n") { raw =>
        let line = if raw.len > 0 && raw[raw.len - 1] == '\r' { raw[..raw.len - 1] } else { raw };
        if line.len > 0 && valid {
            var values: i64[6] = [0; 6]; var n: i64 = 0;
            bstr::split_str(line, "\t") { field =>
                let v, ok = parse_int(field);
                if !ok || n >= 6 { valid = false; } else { values[n] = v; } n += 1;
            };
            if n != 6 || values[0] < -1 || values[0] >= 1048576 || values[1] < 0 || values[2] < 0 ||
               (values[0] >= 0 && values[1] == 0) || values[3] < 26 || values[4] < 0 || values[5] <= 0 ||
               values[1] > 9223372036854775807 - values[2] || values[4] > 9223372036854775807 - values[5] { valid = false; }
            if valid {
                index.entries.push(Entry { tid: values[0] as i32, begin: max(0, values[1] - 1),
                    end: max(0, values[1] - 1) + values[2], container: values[3], slice: values[4], size: values[5], .. });
            }
        }
    };
    if !valid { return bam::Status.Invalid; }
    sort(index.entries) { a, b => a.tid < b.tid || (a.tid == b.tid && a.begin < b.begin) };
    var tid: i32 = -2; var end: i64 = 0;
    for &entry in index.entries {
        if entry.tid != tid { tid = entry.tid; end = 0; }
        end = max(end, entry.end); entry.max_end = end;
    }
    bam::Status.Record
}

// Replace an index; errors clear the previous entries. Empty indexes are valid.
//
// ```goose
// let index = cram_index::Index { .. };
// assert(index.open("") == bam::Status.IoError && index.entries.len == 0);
// ```
fn open(index: Index&, path: const u8[:]) -> bam::Status {
    index.entries.clear(); var bytes: u8[>..] = [];
    index.status = if index.max_bytes < 0 { bam::Status.Invalid } else { read_bytes(path, bytes, index.max_bytes) };
    if index.status == bam::Status.Record { index.status = parse(index, bytes); }
    if index.status != bam::Status.Record { index.entries.clear(); }
    index.status
}

// Reusable per-worker region cursor. fetch replaces its selection; read returns
// only overlapping records in physical slice/record order, including placed
// unmapped records. Call reader.record() only after read returns Record. Do not
// interleave cursors on the same reader. End/errors persist until fetch.
//
// ```goose
// let cursor = cram_index::Cursor { .. };
// assert(cursor.status == bam::Status.Closed);
// ```
root struct Cursor {
    slices: Entry[>..] = [], at: i64 = 0, active: bool = false,
    tid: i32 = -1, begin: i64 = 0, end: i64 = 0,
    last_tid: i32 = -1, last_pos: i64 = -1, status: bam::Status = bam::Status.Closed,
}

// Select candidate slices for [begin,end). Nested spans are retained, duplicate
// slice addresses removed, and selected slices sorted in file order. The caller
// must use the index for this CRAM; skipped slices are not validated.
//
// ```goose
// let cursor = cram_index::Cursor { .. }; let index = cram_index::Index { .. };
// let reader = cram::Reader { .. };
// assert(cursor.fetch(index, reader, 0, 0, 1) == bam::Status.Closed);
// ```
fn fetch(cursor: Cursor&, index: const Index&, reader: cram::Reader&, tid: i32, begin: i64, end: i64) -> bam::Status {
    cursor.slices.clear(); cursor.at = 0; cursor.active = false; cursor.last_tid = -1; cursor.last_pos = -1;
    reader.bytes.clear();
    cursor.tid = tid; cursor.begin = begin; cursor.end = end;
    cursor.status = bam::Status.Closed;
    if reader.input.input.fd < 0 { return cursor.status; }
    cursor.status = index.status; if cursor.status != bam::Status.Record { return cursor.status; }
    if tid < 0 || tid >= reader.references.len || begin < 0 || end < begin || end > reader.references[tid as i64].length {
        cursor.status = bam::Status.Invalid; return cursor.status;
    }
    if begin == end { return cursor.status; }
    var lo: i64 = 0; var hi = index.entries.len;
    while lo < hi {
        let mid = lo + (hi - lo) / 2; let entry = index.entries[mid];
        if entry.tid < tid || (entry.tid == tid && entry.begin < end) { lo = mid + 1; } else { hi = mid; }
    }
    while lo > 0 {
        lo -= 1; let entry = index.entries[lo];
        if entry.tid != tid || entry.max_end <= begin { break; }
        if entry.end > begin { cursor.slices.push(entry); }
    }
    sort(cursor.slices) { a, b => a.container < b.container || (a.container == b.container && a.slice < b.slice) };
    var n: i64 = 0;
    for entry in cursor.slices {
        if n == 0 || cursor.slices[n - 1].container != entry.container || cursor.slices[n - 1].slice != entry.slice {
            cursor.slices[n] = entry; n += 1;
        } else if cursor.slices[n - 1].size != entry.size { cursor.status = bam::Status.Invalid; return cursor.status; }
    }
    cursor.slices.resize(n); cursor.status
}

// Decode only the selected slices, retaining the reader's reusable container
// and codec storage. A fetch after EOF/error may seek backward or cross contigs.
// Reference-provider semantics are the same as cram::read.
//
// ```goose
// let cursor = cram_index::Cursor { .. }; let reader = cram::Reader { .. };
// let source = ExampleReference { .. };
// assert(cursor.read(reader, source) == bam::Status.Closed);
// ```
fn read<R>(cursor: Cursor&, reader: cram::Reader&, source: R&) -> bam::Status {
    if cursor.status != bam::Status.Record { return cursor.status; }
    reader.bytes.clear();
    while cursor.at < cursor.slices.len {
        if !cursor.active {
            let entry = cursor.slices[cursor.at];
            cursor.status = reader.seek_slice(entry.container, entry.slice, entry.size, source);
            if cursor.status != bam::Status.Record { return cursor.status; }
            cursor.active = true;
        }
        if reader.next_record >= reader.slice_records.len { cursor.at += 1; cursor.active = false; continue; }
        cursor.status = reader.read(source); if cursor.status != bam::Status.Record { return cursor.status; }
        let f = reader.fields;
        if f.tid < 0 || f.pos < 0 { continue; }
        if f.tid < cursor.last_tid || (f.tid == cursor.last_tid && f.pos < cursor.last_pos) {
            cursor.status = bam::Status.Unsorted; return cursor.status;
        }
        cursor.last_tid = f.tid; cursor.last_pos = f.pos;
        if f.tid == cursor.tid && f.pos < cursor.end && f.end > cursor.begin { return cursor.status; }
    }
    reader.bytes.clear(); cursor.status = bam::Status.End; cursor.status
}
