// Sequential BAM reader. Coordinates are zero-based, with exclusive ends.
namespace bam;
import std;
import bgzf;
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

// Record means a complete validated alignment; End is normal EOF. All other
// statuses are sticky failures until open/seek. TooLarge is a configured size limit.
// Unsorted is returned by pileup loading, not the sequential reader.
//
// ```goose
// let reader = bam::Reader { .. };
// assert(reader.read() == bam::Status.Closed);
// ```
enum Status { Record, End, Invalid, Truncated, IoError, Closed, TooLarge, Unsorted }

// Reference name location in Reader.names and its declared length.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let ref = reader.references[0];
// assert(reader.names[ref.name..ref.name + ref.name_len] == "chrA");
// assert(ref.length == 200000);
// ```
struct Reference { name: i64, name_len: i64, length: i64 }

// Parsed offsets are relative to the record's bytes, not the file. Use the
// Record view's accessors rather than constructing this metadata yourself.
// reference_len is the CIGAR reference span; end uses BAM's minimum span of 1.
// tid/mate_tid and pos/mate_pos use -1 for missing references/positions.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let f = reader.record().fields;
// assert(f.pos == 5 && f.end == 18 && f.reference_len == 13);
// assert(f.seq_len == 13 && f.mate_tid == -1);
// ```
struct Fields {
    tid:      i32, pos: i64, end: i64, reference_len: i64, flags: u16, mapq: u8,
    mate_tid: i32, mate_pos: i64, template_len: i32,
    name_len: i64, cigar: i64, cigar_len: i64, seq: i64, seq_len: i64,
    qual:     i64, aux: i64,
}

// A borrowed packed record. Reader views expire on the next read/seek/open/close.
// Region stores expose a separate slab-backed pileup::Record view.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.name() == "match" && record.fields.tid == 0);
// ```
struct Record { fields: Fields, bytes: const u8[:] }

// Owns stream, header, reference dictionary and a cross-block scratch buffer.
// Records within one decoded BGZF block borrow that block without a copy.
// Construct with Reader { .. }, then open. Do not mutate implementation fields.
// Limits bound allocations from untrusted lengths and may be raised explicitly.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.references.len == 2);
// assert(reader.header[..3] == "@HD");
// ```
root struct Reader {
    stream:           bgzf::Stream,
    status:           Status         = Status.Closed,
    header:           u8[>..]        = [],
    names:            u8[>..]        = [],
    references:       Reference[>..] = [],
    bytes:            u8[>..]        = [],
    block_begin:      i64 = -1,
    fields:           Fields,
    data_start:       u64 = 0,
    max_record_bytes: i64 = 16777216,
    max_header_bytes: i64 = 16777216,
    max_references:   i64 = 1048576,
}

inline fn signed32(bytes: const u8[:], at: i64) -> i32 { bgzf::le32(bytes, at) as! i32 }

fn from_bgzf(status: bgzf::Status) -> Status {
    match status {
        Ready => Status.Record,
        End => Status.End,
        Truncated => Status.Truncated,
        Invalid => Status.Invalid,
        IoError => Status.IoError,
        Closed => Status.Closed,
    }
}

fn required(reader: Reader&, output: u8[:]) -> bool {
    reader.status = from_bgzf(reader.stream.read(output));
    if reader.status == Status.End { reader.status = Status.Truncated; }
    reader.status == Status.Record
}

// Close the file and decompressor. Header storage remains readable until open.
//
// ```goose
// let reader = bam::Reader { .. };
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// reader.close();
// reader.close();
// assert(reader.read() == bam::Status.Closed);
// ```
fn close(reader: Reader&) {
    reader.stream.close();
    reader.status = Status.Closed;
    reader.bytes.clear();
    reader.block_begin = -1;
}

// Open and read the BAM header. Record means ready to read the first alignment;
// it does not mean an alignment is already available. Failure closes the file.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.name() == "match");
// ```
fn open(reader: Reader&, path: const u8[:]) -> Status {
    reader.close();
    reader.header.clear();
    reader.names.clear();
    reader.references.clear();
    reader.data_start = 0;
    reader.status = from_bgzf(reader.stream.open(path));
    if reader.status != Status.Record { return reader.status; }
    if !read_header(reader) {
        reader.stream.close();
        return reader.status;
    }
    reader.data_start = reader.stream.tell();
    Status.Record
}

// Virtual offset of the next alignment, usable with seek on this same BAM.
// The low 16 bits are the uncompressed position inside a compressed BGZF block.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let start = reader.tell();
// assert(reader.read() == bam::Status.Record && reader.tell() > start);
// ```
fn tell(reader: const Reader&) -> u64 { reader.stream.tell() }

// Seek to an alignment boundary obtained from tell or a matching BAM index.
// Invalidates record views and resets EOF/errors, retaining the header. An
// arbitrary byte offset is not an alignment boundary and may fail on read.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let start = reader.tell();
// assert(reader.read() == bam::Status.Record);
// assert(reader.seek(start) == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.name() == "match");
// ```
fn seek(reader: Reader&, offset: u64) -> Status {
    reader.bytes.clear();
    reader.block_begin = -1;
    if reader.stream.input.fd < 0 { reader.status = Status.Closed; }
    else if offset < reader.data_start { reader.status = Status.Invalid; }
    else { reader.status = from_bgzf(reader.stream.seek(offset)); }
    reader.status
}

fn read_header(reader: Reader&) -> bool {
    var word: u8[4] = [0; 4];
    if !reader.required(word) { return false; }
    if word != "BAM\x01" { reader.status = Status.Invalid; return false; }
    if !reader.required(word) { return false; }
    let text_len = signed32(word, 0) as i64;
    if text_len < 0 { reader.status = Status.Invalid; return false; }
    if text_len > reader.max_header_bytes { reader.status = Status.TooLarge; return false; }
    reader.header.resize(text_len, 0);
    if !reader.required(reader.header) || !reader.required(word) { return false; }
    let refs = signed32(word, 0) as i64;
    if refs < 0 { reader.status = Status.Invalid; return false; }
    if refs > reader.max_references { reader.status = Status.TooLarge; return false; }
    for i in refs {
        if !reader.required(word) { return false; }
        let len = signed32(word, 0) as i64;
        if len < 2 { reader.status = Status.Invalid; return false; }
        if len > reader.max_header_bytes - reader.names.len { reader.status = Status.TooLarge; return false; }
        let start = reader.names.len;
        reader.names.resize(start + len, 0);
        if !reader.required(reader.names[start..]) { return false; }
        if reader.names[start + len - 1] != 0 { reader.status = Status.Invalid; return false; }
        for byte in reader.names[start..start + len - 1] {
            if byte == 0 { reader.status = Status.Invalid; return false; }
        }
        if !reader.required(word) { return false; }
        let length = signed32(word, 0) as i64;
        if length < 0 { reader.status = Status.Invalid; return false; }
        reader.references.push(Reference { name: start, name_len: len - 1, length: length });
    }
    true
}

// Look up a reference by its byte name; -1 means absent.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.reference_id("chrB") == 1);
// assert(reader.reference_id("absent") == -1);
// ```
fn reference_id(reader: const Reader&, name: const u8[:]) -> i32 {
    for ref, i in reader.references {
        if reader.names[ref.name..ref.name + ref.name_len] == name { return i as i32; }
    }
    -1
}

// Borrow the last successfully read record. Call only after read returns Record.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let first = reader.record();
// assert(first.name() == "match");
// assert(reader.read() == bam::Status.Record);
// let second = reader.record();
// assert(second.name() == "reverse"); // first's bytes are no longer valid.
// ```
fn record(reader: const Reader&) -> Record {
    let bytes = if reader.block_begin >= 0 {
        reader.stream.decoded[reader.block_begin..reader.stream.begin]
    } else { reader.bytes[..] };
    Record { fields: reader.fields, bytes: bytes }
}

// Read one alignment, reusing bytes. Record bytes are invalidated even when a
// read fails. End and failures are sticky until open/seek. Records can span
// many BGZF blocks.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// var count = 0;
// var status = reader.read();
// while status == bam::Status.Record { count += 1; status = reader.read(); }
// assert(count == 21 && status == bam::Status.End);
// assert(reader.read() == bam::Status.End);
// ```
fn read(reader: Reader&) -> Status {
    let status = reader.read_payload();
    if status != Status.Record { return status; }
    reader.decode_record()
}

// Internal framing step for indexed filtering. fields is not valid until
// decode_record succeeds; only record().bytes may be inspected between them.
fn read_payload(reader: Reader&) -> Status {
    if reader.status != Status.Record { return reader.status; }
    reader.bytes.clear();
    reader.block_begin = -1;
    var len: i64 = 0;
    if reader.stream.end - reader.stream.begin >= 4 {
        len = signed32(reader.stream.decoded, reader.stream.begin);
        reader.stream.begin += 4;
    } else {
        var word: u8[4] = [0; 4];
        reader.status = from_bgzf(reader.stream.read(word));
        if reader.status != Status.Record { return reader.status; }
        len = signed32(word, 0);
    }
    if len < 32 { reader.status = Status.Invalid; return reader.status; }
    if len > reader.max_record_bytes { reader.status = Status.TooLarge; return reader.status; }
    if len <= reader.stream.end - reader.stream.begin {
        reader.block_begin = reader.stream.begin;
        reader.stream.begin += len;
    } else {
        reader.bytes.resize(len, 0);
        if !reader.required(reader.bytes) { return reader.status; }
    }
    Status.Record
}

fn decode_record(reader: Reader&) -> Status {
    let raw = reader.record();
    let fields, ok = decode(raw.bytes, reader.references.len);
    if !ok { reader.status = Status.Invalid; return reader.status; }
    reader.fields = fields;
    reader.status
}

// Conservative prefilter for mapped records before a query's start. Only
// inspect the fixed core and wire CIGAR, as Seqair's indexed reaches() does.
// A malformed CIGAR or long-CIGAR escape falls through to full validation;
// never discard a CG alignment based on its placeholder's declared span.
fn may_reach(bytes: const u8[:], begin: i64) -> bool {
    let pos = signed32(bytes, 4) as i64;
    if pos >= begin { return true; }
    let start = 32 + bytes[8] as i64;
    let count = bgzf::le16(bytes, 12) as i64;
    if start > bytes.len || count > (bytes.len - start) / 4 { return true; }
    let cigar = bytes[start..start + count * 4];
    if count == 2 && (bgzf::le32(cigar, 0) & 15) == 4 &&
       (bgzf::le32(cigar, 4) & 15) == 3 { return true; }
    var span: i64 = 0;
    for i in count {
        let word = bgzf::le32(cigar, i * 4);
        let code = word & 15;
        if code > 8 { return true; }
        if code == 0 || code == 2 || code == 3 || code == 7 || code == 8 { span += word >> 4; }
    }
    pos + max(span, 1) > begin
}

// End of one BAM auxiliary field, or -1 for invalid/truncated framing.
// Writers use this to preserve tag order without interpreting unknown layouts.
fn aux_end(bytes: const u8[:], start: i64) -> i64 {
    if start < 0 || start > bytes.len || bytes.len - start < 3 { return -1; }
    let kind = bytes[start + 2]; var at = start + 3; var size: i64 = 0;
    if kind == 'A' || kind == 'c' || kind == 'C' { size = 1; }
    else if kind == 's' || kind == 'S' { size = 2; }
    else if kind == 'i' || kind == 'I' || kind == 'f' { size = 4; }
    else if kind == 'd' { size = 8; }
    else if kind == 'Z' || kind == 'H' {
        let end = bstr::find_byte(bytes[at..], 0); if end < 0 { return -1; } size = end + 1;
    } else if kind == 'B' {
        if bytes.len - at < 5 { return -1; }
        let sub = bytes[at]; let count = signed32(bytes, at + 1);
        if count < 0 { return -1; }
        let width = if sub == 'c' || sub == 'C' { 1 } else if sub == 's' || sub == 'S' { 2 }
            else if sub == 'i' || sub == 'I' || sub == 'f' { 4 } else { 0 };
        if width == 0 { return -1; }
        size = count as i64 * width; at += 5;
    } else { return -1; }
    if size > bytes.len - at { -1 } else { at + size }
}

// Validate auxiliary fields while locating the long-CIGAR escape, if present.
// Unknown tag types cannot be skipped safely. B arrays have checked lengths.
fn aux_cigar(bytes: const u8[:], start: i64) -> i64, i64, bool {
    var at = start;
    var cg:    i64 = 0;
    var count: i64 = 0;
    while at < bytes.len {
        if bytes.len - at < 3 { return 0, 0, false; }
        let tag = bytes[at..at + 3];
        let is_cg = tag[0] == 'C' && tag[1] == 'G';
        let kind = tag[2];
        at += 3;
        var size: i64 = 0;
        if kind == 'A' || kind == 'c' || kind == 'C' { size = 1; }
        else if kind == 's' || kind == 'S' { size = 2; }
        else if kind == 'i' || kind == 'I' || kind == 'f' { size = 4; }
        else if kind == 'd' { size = 8; }
        else if kind == 'Z' || kind == 'H' {
            let end = bstr::find_byte(bytes[at..], 0);
            if end < 0 { return 0, 0, false; }
            size = end + 1;
        } else if kind == 'B' {
            if bytes.len - at < 5 { return 0, 0, false; }
            let array_header = bytes[at..at + 5];
            let sub = array_header[0];
            let n = signed32(array_header, 1) as i64;
            if n < 0 { return 0, 0, false; }
            var width: i64 = 0;
            if sub == 'c' || sub == 'C' { width = 1; }
            else if sub == 's' || sub == 'S' { width = 2; }
            else if sub == 'i' || sub == 'I' || sub == 'f' { width = 4; }
            else { return 0, 0, false; }
            at += 5;
            size = n * width;
            if is_cg && cg == 0 && sub == 'I' { cg = at; count = n; }
        } else { return 0, 0, false; }
        if size > bytes.len - at { return 0, 0, false; }
        at += size;
    }
    return cg, count, true;
}

fn decode(bytes: const u8[:], references: i64) -> Fields, bool {
    var f = Fields { .. };
    f.tid = signed32(bytes, 0);
    f.pos = signed32(bytes, 4);
    f.mapq = bytes[9];
    f.flags = bgzf::le16(bytes, 14) as u16;
    f.seq_len = signed32(bytes, 16);
    f.mate_tid = signed32(bytes, 20);
    f.mate_pos = signed32(bytes, 24);
    f.template_len = signed32(bytes, 28);
    f.name_len = bytes[8];
    f.cigar = 32 + f.name_len;
    f.cigar_len = bgzf::le16(bytes, 12);
    f.seq = f.cigar + 4 * f.cigar_len;
    f.qual = f.seq + (f.seq_len + 1) / 2;
    f.aux = f.qual + f.seq_len;
    if f.tid < -1 || f.tid >= references || f.mate_tid < -1 || f.mate_tid >= references ||
       f.pos < -1 || f.mate_pos < -1 || f.seq_len < 0 || f.name_len < 1 || f.aux > bytes.len {
        return f, false;
    }
    if bytes[f.cigar - 1] != 0 { return f, false; }
    if bstr::find_byte(bytes[32..f.cigar - 1], 0) >= 0 { return f, false; }
    let cg, count, aux_ok = aux_cigar(bytes, f.aux);
    if !aux_ok { return f, false; }
    if cg > 0 && count > 65535 && f.cigar_len == 2 &&
       bgzf::le32(bytes, f.cigar) == ((f.seq_len as u32) << 4 | 4) &&
       (bgzf::le32(bytes, f.cigar + 4) & 15) == 3 {
        f.cigar = cg;
        f.cigar_len = count;
    }
    var query: i64 = 0;
    var span:  i64 = 0;
    for i in f.cigar_len {
        let op = bgzf::le32(bytes, f.cigar + i * 4);
        let code = op & 15;
        let len = (op >> 4) as i64;
        if code > 8 { return f, false; }
        if code == 0 || code == 1 || code == 4 || code == 7 || code == 8 { query += len; }
        if code == 0 || code == 2 || code == 3 || code == 7 || code == 8 { span += len; }
    }
    if f.flags & 4 == 0 && f.cigar_len > 0 && f.seq_len > 0 && query != f.seq_len { return f, false; }
    f.reference_len = span;
    f.end = f.pos + (if f.flags & 4 != 0 { 1 } else { max(span, 1) });
    return f, true;
}

// Read name without its trailing NUL.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.name() == "match");
// ```
fn name(record: const Record&) -> const u8[:] { record.bytes[32..31 + record.fields.name_len] }

// CIGAR word, with length in its upper 28 bits and SAM operation in its low 4.
// index must be in [0, fields.cigar_len). Long-CIGAR CG arrays are transparent.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.cigar(0) == 13 << 4); // 13M
// ```
fn cigar(record: const Record&, index: i64) -> u32 {
    assert(index >= 0 && index < record.fields.cigar_len);
    bgzf::le32(record.bytes, record.fields.cigar + index * 4)
}

// Decode one packed base to its SAM IUPAC byte. Stored BAM orientation is kept;
// reverse-strand records must not be reverse-complemented again for pileups.
// index must be in [0, fields.seq_len).
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.base(0) == 'A' && record.base(1) == 'C');
// ```
fn base(record: const Record&, index: i64) -> u8 {
    assert(index >= 0 && index < record.fields.seq_len);
    let packed = record.bytes[record.fields.seq + index / 2];
    let code = if index & 1 == 0 { packed >> 4 } else { packed & 15 };
    "=ACMGRSVTWYHKDBN"[code as i64]
}

// Raw Phred quality; 255 means unavailable, not a high-quality base.
// index must be in [0, fields.seq_len).
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.quality(0) == 13 && record.quality(1) == 20);
// ```
fn quality(record: const Record&, index: i64) -> u8 {
    assert(index >= 0 && index < record.fields.seq_len);
    record.bytes[record.fields.qual + index]
}

// Raw BAM auxiliary fields, retained without decoding or re-encoding.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.aux() == "NMi\x03\x00\x00\x00RGZgroup\x00");
// ```
fn aux(record: const Record&) -> const u8[:] { record.bytes[record.fields.aux..] }
