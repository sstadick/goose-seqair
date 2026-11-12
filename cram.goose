// CRAM 3.0/3.1 streaming input, following Seqair's slice reconstruction model.
namespace cram;
import std;
import bam;
import bgzf;
import sam;
import cram_encoding;
import cram_rans;
import cram_nx16;
import cram_arith;
import cram_fqz;
import cram_tok3;
import cram_slice;
import pileup;
import bam_encode;
import deflate from "github.com/sstadick/goose-deflate@8fcf5832a3fa680f4c6a5332e2436f8a6dd5c299";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

extern fn reference_md5(input: const u8[:], output: u8[:]) = "goose_seqair_reference_md5";

struct Container {
    length: i64, tid: i32, start: i64, span: i64, records: i64,
    counter: i64, bases: i64, blocks: i64,
}

struct Window { begin: i64 = -1, end: i64 = 0, offset: i64 = 0 }

// Owns the input, reusable codec scratch and decoded slice slabs. Header and
// reference names live until open; records live until the next read/open/close.
// Limits bound individual containers, decoded blocks/slabs and records, not
// total process memory. End/errors are sticky until open. No CRAI seeking.
//
// ```goose
// let reader = cram::Reader { .. };
// assert(reader.status == bam::Status.Closed);
// assert(reader.max_container_bytes == 268435456);
// ```
root struct Reader {
    input: bgzf::Stream,
    status: bam::Status = bam::Status.Closed,
    header: u8[>..] = [], names: u8[>..] = [], references: bam::Reference[>..] = [],
    container: u8[>..] = [], container_header: u8[>..] = [], landmarks: i64[>..] = [],
    info: Container,
    decoded: u8[>..] = [], blocks: cram_encoding::Block[>..] = [],
    rans_symbols: u8[>..] = [], rans_steps: u32[>..] = [],
    nx_frequencies: u32[>..] = [], nx_cumulative: u32[>..] = [], nx_symbols: u8[>..] = [],
    max_container_bytes: i64 = 268435456, max_decoded_bytes: i64 = 268435456,
    max_header_bytes: i64 = 16777216, max_references: i64 = 1048576,
    unsupported_method: i32 = -1,
    max_record_bytes: i64 = 16777216,
    compression: cram_encoding::Compression,
    encodings: cram_encoding::Encoding[>..] = [], symbols: cram_encoding::Symbol[>..] = [],
    tags: cram_encoding::Tag[>..] = [], tag_uses: cram_encoding::Tag[>..] = [],
    tag_lines: cram_encoding::TagLine[>..] = [],
    groups: bam::Reference[>..] = [], group_names: u8[>..] = [], prefix: u8[>..] = [],
    scratch_name: u8[>..] = [], scratch_bases: u8[>..] = [], scratch_cigar: u32[>..] = [],
    scratch_qualities: u8[>..] = [], scratch_aux: u8[>..] = [],
    scratch_payload: u8[>..] = [], scratch_md: u8[>..] = [],
    slice_records: pileup::SlimRecord[>..] = [], slice_mates: cram_slice::Mate[>..] = [],
    slice_names: u8[>..] = [], slice_bases: u8[>..] = [], slice_cigars: u32[>..] = [],
    slice_qualities: u8[>..] = [], slice_aux: u8[>..] = [],
    reference_bases: u8[>..] = [], reference_tid: i32 = -1,
    windows: Window[>..] = [], window_bases: u8[>..] = [],
    next_slice: i64 = 0, next_record: i64 = 0,
    bytes: u8[>..] = [], fields: bam::Fields,
}

// Close input and invalidate the current record; safe to call repeatedly.
//
// ```goose
// let reader = cram::Reader { .. };
// reader.close(); reader.close();
// assert(reader.status == bam::Status.Closed);
// ```
fn close(reader: Reader&) {
    reader.input.close(); reader.status = bam::Status.Closed;
    reader.container.clear(); reader.decoded.clear(); reader.blocks.clear();
    reader.landmarks.clear(); reader.clear_slice(); reader.bytes.clear();
    reader.next_slice = 0; reader.next_record = 0;
}

fn required(reader: Reader&, out: u8[:]) -> bool {
    let status = bam::from_bgzf(reader.input.input_exact(out));
    if status != bam::Status.Record { reader.status = if status == bam::Status.End { bam::Status.Truncated } else { status }; }
    status == bam::Status.Record
}

fn field(reader: Reader&, long: bool = false) -> i64 {
    if reader.status != bam::Status.Record { return 0; }
    var bytes: u8[9] = [0; 9];
    if !reader.required(bytes[..1]) { return 0; }
    let first = bytes[0];
    var extra = 0;
    if long {
        var mask: u8 = 128;
        while extra < 8 && (first & mask) != 0 { extra += 1; mask >>= 1; }
    } else { extra = if first < 128 { 0 } else if first < 192 { 1 } else if first < 224 { 2 } else if first < 240 { 3 } else { 4 }; }
    if extra > 0 && !reader.required(bytes[1..extra + 1]) { return 0; }
    reader.container_header.append(bytes[..extra + 1]);
    let c = cram_encoding::Cursor { bytes: bytes[..extra + 1] };
    if long { c.ltf8() } else { c.itf8() as i64 }
}

fn next_container(reader: Reader&) -> bam::Status {
    if reader.status != bam::Status.Record { return reader.status; }
    reader.container.clear(); reader.container_header.clear(); reader.landmarks.clear();
    var word: u8[4] = [0; 4];
    reader.status = bam::from_bgzf(reader.input.input_exact(word));
    if reader.status != bam::Status.Record { return reader.status; }
    reader.container_header.append(word);
    let size = (bgzf::le32(word, 0) as! i32) as i64;
    reader.info.length = size;
    reader.info.tid = reader.field() as i32;
    reader.info.start = reader.field(); reader.info.span = reader.field();
    reader.info.records = reader.field(); reader.info.counter = reader.field(true);
    reader.info.bases = reader.field(true); reader.info.blocks = reader.field();
    let n = reader.field();
    if reader.status != bam::Status.Record { return reader.status; }
    if size < 0 || reader.info.tid < -2 || reader.info.records < 0 || reader.info.counter < 0 ||
       reader.info.blocks < 0 || n < 0 { reader.status = bam::Status.Invalid; return reader.status; }
    if size > reader.max_container_bytes || n > reader.max_container_bytes / 4 {
        reader.status = bam::Status.TooLarge; return reader.status;
    }
    for i in n {
        let landmark = reader.field();
        if reader.status != bam::Status.Record { return reader.status; }
        if landmark < 0 || landmark >= size || (i > 0 && landmark <= reader.landmarks[i - 1]) {
            reader.status = bam::Status.Invalid; return reader.status;
        }
        reader.landmarks.push(landmark);
    }
    if !reader.required(word) { return reader.status; }
    if deflate::crc32(reader.container_header) != bgzf::le32(word, 0) {
        reader.status = bam::Status.Invalid; return reader.status;
    }
    reader.container.resize(size, 0);
    if !reader.required(reader.container) { return reader.status; }
    if reader.info.records == 0 && reader.info.tid == -1 && reader.info.start == 0x454f46 {
        reader.status = bam::Status.End;
    }
    reader.status
}

fn decode_block(reader: Reader&, c: cram_encoding::Cursor&) -> bam::Status {
    let start = c.at;
    let method = c.byte(); let kind = c.byte(); let id = c.itf8();
    let compressed = c.itf8() as i64; let size = c.itf8() as i64;
    if !c.ok { return bam::Status.Truncated; }
    if compressed < 0 || size < 0 || (kind != 0 && kind != 1 && kind != 2 && kind != 4 && kind != 5) {
        return bam::Status.Invalid;
    }
    let data = c.take(compressed); let end = c.at; let checksum = c.take(4);
    if !c.ok { return bam::Status.Truncated; }
    if deflate::crc32(c.bytes[start..end]) != bgzf::le32(checksum, 0) { return bam::Status.Invalid; }
    // Raw and empty blocks follow Seqair's size exceptions. Other codecs must
    // decode exactly the advertised length before any bytes are exposed.
    let count = if size == 0 { 0 } else if method == 0 { compressed } else { size };
    if size > reader.max_decoded_bytes || count > reader.max_decoded_bytes - reader.decoded.len {
        return bam::Status.TooLarge;
    }
    let begin = reader.decoded.len;
    reader.decoded.resize(begin + count, 0);
    let output = reader.decoded[begin..];
    if count != 0 {
        if method == 0 { output.move_from(data); }
        else if method == 1 {
            let result = reader.input.decoder.decompress(deflate::Format.Gzip, data, output);
            if result.status != deflate::Status.Success || result.written != count { return bam::Status.Invalid; }
        } else if method == 2 || method == 3 {
            let result = cram_arith::external_decode(method as i32, data, output, reader.max_decoded_bytes);
            if result < 0 { reader.unsupported_method = method as i32; }
            if result != 1 { return bam::Status.Invalid; }
        } else if method == 4 {
            if !cram_rans::decode(data, output, reader.rans_symbols, reader.rans_steps) { return bam::Status.Invalid; }
        } else if method == 5 {
            if !cram_nx16::decode(data, output, reader.nx_frequencies, reader.nx_cumulative, reader.nx_symbols,
                                 reader.max_decoded_bytes) { return bam::Status.Invalid; }
        } else if method == 6 {
            if !cram_arith::decode(data, output, reader.max_decoded_bytes) { return bam::Status.Invalid; }
        } else if method == 7 {
            if !cram_fqz::decode(data, output, reader.max_decoded_bytes) { return bam::Status.Invalid; }
        } else if method == 8 {
            if !cram_tok3::decode(data, output, reader.nx_frequencies, reader.nx_cumulative, reader.nx_symbols,
                                 reader.max_decoded_bytes) { return bam::Status.Invalid; }
        } else { reader.unsupported_method = method as i32; return bam::Status.Invalid; }
    }
    reader.blocks.push(cram_encoding::Block { id: id, kind: kind, begin: begin, size: count, .. });
    bam::Status.Record
}

// Open CRAM 3.0/3.1 and decode its SAM header. Record means ready, not a read.
// Validates container/block CRCs. Older CRAM versions return Invalid. Bzip2
// and XZ blocks require system libbz2/liblzma; missing codecs return Invalid
// with unsupported_method set. No reference is fetched during header parsing.
//
// ```goose
// let reader = cram::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Invalid);
// ```
fn open(reader: Reader&, path: const u8[:]) -> bam::Status {
    reader.close(); reader.header.clear(); reader.names.clear(); reader.references.clear();
    reader.groups.clear(); reader.group_names.clear(); reader.prefix.clear();
    reader.reference_bases.clear(); reader.reference_tid = -1;
    var basename = 0;
    for b, i in path { if b == '/' { basename = i + 1; } }
    reader.prefix.append(path[basename..]);
    reader.unsupported_method = -1;
    reader.status = bam::from_bgzf(reader.input.open(path));
    if reader.status != bam::Status.Record { return reader.status; }
    reader.status = reader.open_header();
    if reader.status != bam::Status.Record { reader.input.close(); }
    reader.status
}

fn open_header(reader: Reader&) -> bam::Status {
    var definition: u8[26] = [0; 26];
    if !reader.required(definition) { return reader.status; }
    if definition[..4] != "CRAM" || definition[4] != 3 || definition[5] > 1 { return bam::Status.Invalid; }
    if reader.next_container() != bam::Status.Record {
        return if reader.status == bam::Status.End { bam::Status.Truncated } else { reader.status };
    }
    let c = cram_encoding::Cursor { bytes: reader.container };
    let status = reader.decode_block(c);
    if status != bam::Status.Record { return status; }
    let b = reader.blocks[0];
    if b.kind != 0 || b.size < 4 { return bam::Status.Invalid; }
    let text_size = (bgzf::le32(reader.decoded, b.begin) as! i32) as i64;
    if text_size < 0 || text_size > b.size - 4 { return bam::Status.Invalid; }
    if text_size > reader.max_header_bytes { return bam::Status.TooLarge; }
    reader.header.append(reader.decoded[b.begin + 4..b.begin + 4 + text_size]);
    var at = 0;
    while at < reader.header.len {
        let newline = bstr::find_byte(reader.header[at..], '\n');
        var end = if newline < 0 { reader.header.len } else { at + newline };
        let next = end + 1;
        if end > at && reader.header[end - 1] == '\r' { end -= 1; }
        let parsed = sam::parse_header_line(reader.header[at..end], reader.references, reader.names, reader.max_references);
        if parsed != bam::Status.Record { return parsed; }
        let line = reader.header[at..end];
        if line.len >= 4 && line[..4] == "@RG\t" {
            bstr::split_str(line[4..], "\t") { field =>
                if field.len >= 3 && field[..3] == "ID:" {
                    reader.groups.push(bam::Reference { name: reader.group_names.len, name_len: field.len - 3, length: 0 });
                    reader.group_names.append(field[3..]);
                }
            };
        }
        at = next;
    }
    reader.landmarks.clear();
    bam::Status.Record
}

fn clear_slice(reader: Reader&) {
    reader.slice_records.clear(); reader.slice_mates.clear(); reader.slice_names.clear();
    reader.slice_bases.clear(); reader.slice_cigars.clear(); reader.slice_qualities.clear(); reader.slice_aux.clear();
    reader.next_record = 0;
}

fn output(reader: Reader&) -> cram_slice::Output {
    cram_slice::Output { records: reader.slice_records, mates: reader.slice_mates,
        names: reader.slice_names, bases: reader.slice_bases, cigars: reader.slice_cigars,
        qualities: reader.slice_qualities, auxiliary: reader.slice_aux, limit: reader.max_decoded_bytes }
}

fn prepare_container(reader: Reader&) -> bam::Status {
    if reader.next_container() != bam::Status.Record { return reader.status; }
    reader.decoded.clear(); reader.blocks.clear(); reader.encodings.clear(); reader.symbols.clear();
    reader.tags.clear(); reader.tag_uses.clear(); reader.tag_lines.clear(); reader.next_slice = 0;
    let c = cram_encoding::Cursor { bytes: reader.container };
    let status = reader.decode_block(c);
    if status != bam::Status.Record { return status; }
    let b = reader.blocks[0];
    if b.kind != 1 || (reader.landmarks.len > 0 && c.at > reader.landmarks[0]) { return bam::Status.Invalid; }
    let table = cram_encoding::Table { encodings: reader.encodings, symbols: reader.symbols };
    let header, ok = cram_encoding::compression(reader.decoded[b.begin..b.begin + b.size], table,
        reader.tags, reader.tag_uses, reader.tag_lines);
    reader.compression = header;
    if ok { bam::Status.Record } else { bam::Status.Invalid }
}

// Seqair's multi-reference pre-scan consumes all series before fetching one
// union window per referenced contig, then rewinds the shared cursors.
fn scan_windows(reader: Reader&, ds: cram_slice::Series&, slice: const cram_slice::Header&) -> bool {
    reader.windows.resize(reader.references.len, Window { .. });
    reader.windows.fill(Window { .. });
    var previous = slice.start;
    for i in slice.records {
        reader.scratch_name.clear(); reader.scratch_bases.clear(); reader.scratch_cigar.clear();
        reader.scratch_qualities.clear(); reader.scratch_aux.clear(); reader.scratch_payload.clear(); reader.scratch_md.clear();
        let scratch = cram_slice::Scratch { name: reader.scratch_name, bases: reader.scratch_bases,
            cigar: reader.scratch_cigar, qualities: reader.scratch_qualities, aux: reader.scratch_aux,
            payload: reader.scratch_payload, md: reader.scratch_md };
        let p = cram_slice::begin_record(ds, slice, previous, i, scratch, reader.tag_uses, reader.tag_lines, reader.references.len);
        if !ds.context.ok { return false; }
        let span = cram_slice::scan_span(ds, p, reader.scratch_payload);
        if !ds.context.ok { return false; }
        if p.tid >= 0 && p.flags & 4 == 0 {
            let length = reader.references[p.tid as i64].length;
            let w = reader.windows[p.tid as i64];
            let begin = min(p.pos, length); let end = min(p.pos + max(1, span), length);
            reader.windows[p.tid as i64] = Window { begin: if w.begin < 0 { begin } else { min(w.begin, begin) }, end: max(w.end, end), .. };
        }
    }
    for &b in ds.context.blocks { b.at = 0; }
    ds.context.bit = 0;
    true
}

// The reference source supplies fetch(name, begin, end, output) -> bool.
fn decode_slice<R>(reader: Reader&, source: R&) -> bam::Status {
    reader.clear_slice(); reader.decoded.clear(); reader.blocks.clear();
    reader.reference_bases.clear(); reader.reference_tid = -1;
    let at = reader.landmarks[reader.next_slice];
    let end = if reader.next_slice + 1 < reader.landmarks.len { reader.landmarks[reader.next_slice + 1] } else { reader.container.len };
    let c = cram_encoding::Cursor { bytes: reader.container[..end], at: at };
    var status = reader.decode_block(c);
    if status != bam::Status.Record { return status; }
    let b = reader.blocks[0];
    if b.kind != 2 { return bam::Status.Invalid; }
    let slice, ok = cram_slice::parse(reader.decoded[b.begin..b.begin + b.size]);
    if !ok || slice.records > reader.max_decoded_bytes / 32 || slice.blocks > reader.container.len / 9 {
        return bam::Status.Invalid;
    }
    for i in slice.blocks {
        status = reader.decode_block(c);
        if status != bam::Status.Record { return status; }
    }
    if c.at != end { return bam::Status.Invalid; }
    var core = -1; var embedded = -1;
    for &part, i in reader.blocks {
        if part.kind == 5 { core = i; }
        if part.kind == 4 && part.id == slice.embedded && slice.embedded >= 0 { embedded = i; part.kind = 0; }
    }
    if core < 0 { return bam::Status.Invalid; }
    // Resolve descriptors once per slice, not once per decoded value. Shared
    // content IDs still share the same cursor; unused missing blocks stay lazy.
    for &encoding in reader.encodings {
        encoding.external_slot = -1;
        if encoding.kind != 1 && encoding.kind != 5 { continue; }
        for part, i in reader.blocks {
            if part.kind == 4 && part.id == encoding.id { encoding.external_slot = i; break; }
        }
    }
    let no_reference = embedded < 0 && !reader.compression.reference;
    var previous = slice.start;
    let ctx = cram_encoding::Context { blocks: reader.blocks, data: reader.decoded, core: core, limit: reader.max_record_bytes, .. };
    let table = cram_encoding::Table { encodings: reader.encodings, symbols: reader.symbols };
    let ds = cram_slice::Series { header: reader.compression, table: table, context: ctx };
    let multi = slice.tid == -2 && !no_reference && embedded < 0;
    if multi {
        if !reader.scan_windows(ds, slice) { return bam::Status.Invalid; }
        for &w, tid in reader.windows {
            if w.begin < 0 { continue; }
            if w.end - w.begin > reader.max_decoded_bytes - reader.reference_bases.len { return bam::Status.TooLarge; }
            let ref = reader.references[tid];
            if !source.fetch(reader.names[ref.name..ref.name + ref.name_len], w.begin, w.end, reader.window_bases) ||
               reader.window_bases.len != w.end - w.begin { return bam::Status.Invalid; }
            w.offset = reader.reference_bases.len;
            for &base in reader.window_bases { if base >= 'a' && base <= 'z' { base -= 32; } }
            reader.reference_bases.append(reader.window_bases);
        }
    }
    for i in slice.records {
        reader.scratch_name.clear(); reader.scratch_bases.clear(); reader.scratch_cigar.clear();
        reader.scratch_qualities.clear(); reader.scratch_aux.clear(); reader.scratch_payload.clear(); reader.scratch_md.clear();
        let scratch = cram_slice::Scratch { name: reader.scratch_name, bases: reader.scratch_bases,
            cigar: reader.scratch_cigar, qualities: reader.scratch_qualities, aux: reader.scratch_aux,
            payload: reader.scratch_payload, md: reader.scratch_md };
        let pending = cram_slice::begin_record(ds, slice, previous, i, scratch, reader.tag_uses, reader.tag_lines, reader.references.len);
        if !ctx.ok { return bam::Status.Invalid; }
        var start: i64 = 0;
        if !multi && !no_reference && embedded < 0 && pending.tid >= 0 && pending.flags & 4 == 0 {
            let ref = reader.references[pending.tid as i64];
            start = if slice.tid >= 0 { max(0, slice.start - 1) } else { 0 };
            let stop = if slice.tid >= 0 && slice.span > 0 { min(ref.length, start + slice.span) } else { ref.length };
            if pending.tid != reader.reference_tid {
                if start > stop || stop - start > reader.max_decoded_bytes { return bam::Status.TooLarge; }
                if !source.fetch(reader.names[ref.name..ref.name + ref.name_len], start, stop, reader.reference_bases) ||
                   reader.reference_bases.len != stop - start { return bam::Status.Invalid; }
                for &base in reader.reference_bases { if base >= 'a' && base <= 'z' { base -= 32; } }
                if slice.tid >= 0 {
                    var nonzero = false;
                    for byte in slice.md5 { nonzero = nonzero || byte != 0; }
                    if nonzero {
                        var digest: u8[16] = [0; 16];
                        reference_md5(reader.reference_bases, digest);
                        if digest != slice.md5 { return bam::Status.Invalid; }
                    }
                }
                reader.reference_tid = pending.tid;
            }
        }
        let window = if multi && pending.tid >= 0 { reader.windows[pending.tid as i64] } else { Window { .. } };
        let reference = if embedded >= 0 { ctx.payload(embedded) } else if no_reference { "" }
            else if multi { if window.begin < 0 { "" } else { reader.reference_bases[window.offset..window.offset + window.end - window.begin] } }
            else { reader.reference_bases[..] };
        let out = reader.output();
        cram_slice::finish_record(ds, pending, scratch, out, reference,
            if embedded >= 0 { slice.start - 1 } else if multi { window.begin } else { start }, no_reference, reader.groups, reader.group_names);
        if !ctx.ok { return bam::Status.Invalid; }
    }
    let out = reader.output();
    if !out.resolve(reader.prefix, slice.counter, !reader.compression.names) { return bam::Status.TooLarge; }
    reader.next_slice += 1;
    bam::Status.Record
}

// Read one record, reconstructing its entire slice before resolving mates.
// source.fetch(name, begin, end, output) must replace output with exactly the
// requested zero-based half-open reference window and return true on success.
// An adapter can use goose-bio faidx::Reader; no FASTA dependency is imposed.
// Reference bases are uppercased and nonzero single-reference MD5s checked.
// Embedded/no-reference slices bypass source. Fetch failures and checked MD5
// mismatches fail; multi-reference/zero-MD5 slices cannot verify FASTA identity.
//
// ```goose
// let reader = cram::Reader { .. };
// let source = ExampleReference { .. };
// assert(reader.read(source) == bam::Status.Closed);
// ```
fn read<R>(reader: Reader&, source: R&) -> bam::Status {
    if reader.status != bam::Status.Record { return reader.status; }
    reader.bytes.clear();
    while reader.next_record >= reader.slice_records.len {
        if reader.next_slice >= reader.landmarks.len {
            reader.status = reader.prepare_container();
            if reader.status != bam::Status.Record { return reader.status; }
            if reader.landmarks.len == 0 { continue; }
        }
        reader.status = reader.decode_slice(source);
        if reader.status != bam::Status.Record { reader.clear_slice(); return reader.status; }
    }
    let f = reader.slice_records[reader.next_record];
    let a = bam_encode::Alignment {
        name: reader.slice_names[(f.name as i64)..f.name as i64 + f.name_len],
        tid: f.tid, pos: f.pos, flags: f.flags, mapq: f.mapq, mate_tid: f.mate_tid, mate_pos: f.mate_pos, template_len: f.template_len,
        cigar: reader.slice_cigars[(f.cigar as i64)..f.cigar as i64 + f.cigar_len],
        sequence: reader.slice_bases[(f.seq as i64)..f.seq as i64 + f.seq_len],
        qualities: reader.slice_qualities[(f.qual as i64)..f.qual as i64 + f.seq_len],
        aux: reader.slice_aux[(f.aux as i64)..f.aux as i64 + f.aux_len],
    };
    reader.status = bam_encode::encode_body(a, reader.bytes, reader.references.len, reader.max_record_bytes, bam_encode::SequenceEncoding.Codes);
    if reader.status != bam::Status.Record { reader.bytes.clear(); return reader.status; }
    let fields, ok = bam::decode(reader.bytes, reader.references.len);
    if !ok { reader.status = bam::Status.Invalid; reader.bytes.clear(); return reader.status; }
    reader.fields = fields; reader.next_record += 1;
    reader.status
}

// Borrow the current packed BAM record. Only valid after read returns Record;
// the next read, open or close invalidates it, including failed reads.
//
// ```goose
// let reader = cram::Reader { .. };
// assert(reader.status == bam::Status.Closed); // Read before calling record().
// ```
fn record(reader: const Reader&) -> bam::Record { bam::Record { bytes: reader.bytes, fields: reader.fields } }

// Resolve a SAM reference name; return -1 when absent.
//
// ```goose
// let reader = cram::Reader { .. };
// assert(reader.reference_id("missing") == -1);
// ```
fn reference_id(reader: const Reader&, name: const u8[:]) -> i32 { sam::find_reference(reader.references, reader.names, name) }

// Load a region through the shared pileup engine, with a reference provider.
// Requires sorted input and a freshly opened reader per independent query,
// like the BAM/SAM sequential path. Errors clear the store; success links mates.
//
// ```goose
// let reader = cram::Reader { .. };
// let source = ExampleReference { .. };
// let store = pileup::Store { .. };
// assert(cram::load(store, reader, source, 0, 0, 10) == bam::Status.Invalid);
// ```
fn load<R>(store: pileup::Store&, reader: Reader&, source: R&, tid: i32, begin: i64, end: i64,
           filter: pileup::Filter = pileup::Filter { .. }) -> bam::Status {
    store.clear(); store.tid = tid; store.begin = begin; store.end = end;
    if tid < 0 || tid >= reader.references.len || begin < 0 || end < begin || end > reader.references[tid as i64].length {
        return bam::Status.Invalid;
    }
    var status = bam::Status.Record;
    var last_tid: i32 = -1; var last_pos: i64 = -1;
    while begin != end {
        status = reader.read(source);
        if status == bam::Status.End { status = bam::Status.Record; break; }
        if status != bam::Status.Record { break; }
        let record = reader.record(); let f = record.fields;
        if f.flags & 4 != 0 || f.tid < 0 || f.pos < 0 { continue; }
        if f.tid < last_tid || (f.tid == last_tid && f.pos < last_pos) { status = bam::Status.Unsorted; break; }
        last_tid = f.tid; last_pos = f.pos;
        if f.tid > tid || (f.tid == tid && f.pos >= end) { break; }
        if f.tid != tid || f.end <= begin || f.flags & filter.exclude_flags != 0 ||
           f.mapq < filter.min_mapq || f.reference_len == 0 { continue; }
        status = pileup::retain(store, record);
        if status != bam::Status.Record { break; }
    }
    if status != bam::Status.Record { store.clear(); } else { store.link_mates(); }
    status
}
