// Reference-free CRAM 3.x encoding: inverse of the Seqair data-series model.
namespace cram_writer;
import std;
import bam;
import bam_encode;
import bgzf;
import sam;
import sam_writer;
import pileup;
import io from "github.com/sstadick/goose-io@v0.2.2";
import deflate from "github.com/sstadick/goose-deflate@8fcf5832a3fa680f4c6a5332e2436f8a6dd5c299";

enum Key { BF, CF, RI, RL, AP, RG, RN, MF, NS, NP, TS, TL, FN, FC, FP, DL, MQ, BA, QS, IN, SC, RS, PD, HC, BB, Tags }
inline fn index(key: Key) -> i64 {
    match key {
        BF => 0, CF => 1, RI => 2, RL => 3, AP => 4, RG => 5, RN => 6, MF => 7, NS => 8,
        NP => 9, TS => 10, TL => 11, FN => 12, FC => 13, FP => 14, DL => 15, MQ => 16,
        BA => 17, QS => 18, IN => 19, SC => 20, RS => 21, PD => 22, HC => 23, BB => 24, Tags => 25,
    }
}
struct TagLine { begin: i64, length: i64 }

// Owns a bounded, reference-free CRAM 3.0/3.1 output. Buffers one multi-reference
// slice per container, with independent reusable data-series buffers, gzip/raw
// blocks and detached mate fields. No reference fetching, sorting, lossy quality
// reduction or index creation. Set limits before open; check close for errors.
//
// ```goose
// let writer = cram_writer::Writer { .. };
// assert(writer.close() == bam::Status.Record);
// ```
root struct Writer {
    output: io::BufferedWriter = io::BufferedWriter { file: io::FileHandle { fd: -1 }, used: 0, buffer: [] },
    encoder: deflate::Compressor,
    status: bam::Status = bam::Status.Closed,
    references: i64 = 0, level: i32 = 6,
    records: i64 = 0, counter: i64 = 0, bases: i64 = 0, pending_bytes: i64 = 0,
    records_per_slice: i64 = 10000, max_slice_bytes: i64 = 67108864, max_record_bytes: i64 = 16777216,
    bf: u8[>..] = [], cf: u8[>..] = [], ri: u8[>..] = [], rl: u8[>..] = [], ap: u8[>..] = [], rg: u8[>..] = [],
    rn: u8[>..] = [], mf: u8[>..] = [], ns: u8[>..] = [], np: u8[>..] = [], ts: u8[>..] = [], tl: u8[>..] = [],
    fn_data: u8[>..] = [], fc: u8[>..] = [], fp: u8[>..] = [], dl: u8[>..] = [], mq: u8[>..] = [],
    ba: u8[>..] = [], qs: u8[>..] = [], in_data: u8[>..] = [], sc: u8[>..] = [], rs: u8[>..] = [],
    pd: u8[>..] = [], hc: u8[>..] = [], bb: u8[>..] = [], tag_data: u8[>..] = [],
    dictionary: u8[>..] = [], lines: TagLine[>..] = [], tags: u32[>..] = [],
    tag_line: u8[>..] = [], sequence: u8[>..] = [], payload: u8[>..] = [],
    header: u8[>..] = [], section: u8[>..] = [],
    container: u8[>..] = [], compressed: u8[>..] = [],
}

fn stream(writer: Writer&, i: i64) -> u8[>..]& {
    match i {
        0 => writer.bf, 1 => writer.cf, 2 => writer.ri, 3 => writer.rl, 4 => writer.ap, 5 => writer.rg,
        6 => writer.rn, 7 => writer.mf, 8 => writer.ns, 9 => writer.np, 10 => writer.ts, 11 => writer.tl,
        12 => writer.fn_data, 13 => writer.fc, 14 => writer.fp, 15 => writer.dl, 16 => writer.mq,
        17 => writer.ba, 18 => writer.qs, 19 => writer.in_data, 20 => writer.sc, 21 => writer.rs,
        22 => writer.pd, 23 => writer.hc, 24 => writer.bb, _ => writer.tag_data,
    }
}

fn itf8(out: u8[>..]&, value: i64) {
    let v = value as! u32;
    if v < 128 { out.push(v as u8); }
    else if v < 16384 { out.push((0x80 | (v >> 8)) as u8); out.push((v & 255) as u8); }
    else if v < 2097152 {
        out.push((0xc0 | (v >> 16)) as u8); out.push(((v >> 8) & 255) as u8); out.push((v & 255) as u8);
    } else if v < 268435456 {
        out.push((0xe0 | (v >> 24)) as u8);
        for shift in [16, 8, 0] { out.push(((v >> shift) & 255) as u8); }
    } else {
        out.push((0xf0 | (v >> 28)) as u8);
        for shift in [20, 12, 4] { out.push(((v >> shift) & 255) as u8); }
        out.push((v & 15) as u8);
    }
}

fn ltf8(out: u8[>..]&, value: i64) {
    let v = value as u64; var n = 0;
    while n < 8 && v >= ((1 as u64) << (7 * (n + 1))) { n += 1; }
    let first = if n == 8 { 255 } else { ((255 << (8 - n)) & 255) | (v >> (8 * n)) as i64 };
    out.push(first as u8);
    for i in n { out.push(((v >> (8 * (n - i - 1))) & 255) as u8); }
}

fn word(out: u8[>..]&, value: u32) { sam::append_word(out, value, 4); }
fn sized(out: u8[>..]&, data: const u8[:]) { itf8(out, data.len); out.append(data); }
fn put(writer: Writer&, key: Key, value: i64) { itf8(writer.stream(key.index()), value); }
fn array(writer: Writer&, key: Key, data: const u8[:]) { sized(writer.stream(key.index()), data); }

fn external(out: u8[>..]&, id: i64) {
    // Writer-assigned content IDs are 1..26, so their ITF8 parameter is one byte.
    out.append("\x01\x01"); out.push(id as u8);
}

fn array_encoding(out: u8[>..]&, id: i64) {
    // Length and byte values deliberately share the same external cursor.
    out.append("\x04\x06"); external(out, id); external(out, id);
}

fn encode_block(writer: Writer&, kind: u8, id: i64, data: const u8[:], compress: bool) {
    var n: i64 = 0;
    if compress && writer.level > 0 && data.len >= 64 {
        writer.compressed.resize(writer.encoder.bound(deflate::Format.Gzip, data.len), 0);
        n = writer.encoder.compress(deflate::Format.Gzip, data, writer.compressed);
        if n <= 0 { writer.status = bam::Status.Invalid; return; }
    }
    let gzip = n > 0 && n < data.len;
    let begin = writer.container.len;
    writer.container.push(if gzip { 1 } else { 0 }); writer.container.push(kind);
    itf8(writer.container, id); itf8(writer.container, if gzip { n } else { data.len }); itf8(writer.container, data.len);
    if gzip { writer.container.append(writer.compressed[..n]); } else { writer.container.append(data); }
    let crc = deflate::crc32(writer.container[begin..]); word(writer.container, crc);
}

fn container_header(writer: Writer&, tid: i64, start: i64, count: i64, counter: i64, bases: i64, blocks: i64, landmark: i64) {
    writer.header.clear(); word(writer.header, writer.container.len as u32);
    itf8(writer.header, tid); itf8(writer.header, start); itf8(writer.header, 0);
    itf8(writer.header, count); ltf8(writer.header, counter); ltf8(writer.header, bases); itf8(writer.header, blocks);
    itf8(writer.header, if landmark < 0 { 0 } else { 1 });
    if landmark >= 0 { itf8(writer.header, landmark); }
    let crc = deflate::crc32(writer.header); word(writer.header, crc);
}

fn send(writer: Writer&) {
    if writer.status != bam::Status.Record { return; }
    if io::write(writer.output, writer.header) != writer.header.len ||
       io::write(writer.output, writer.container) != writer.container.len { writer.status = bam::Status.IoError; }
}

fn reset_slice(writer: Writer&) {
    for i in 26 { writer.stream(i).clear(); }
    writer.dictionary.clear(); writer.lines.clear(); writer.tags.clear();
    writer.records = 0; writer.bases = 0; writer.pending_bytes = 0;
}

// Flush pending records into a complete container. Empty flush is a no-op.
// It also flushes the I/O buffer (not fsync); failures poison subsequent writes.
//
// ```goose
// let writer = cram_writer::Writer { .. };
// assert(writer.flush() == bam::Status.Closed);
// ```
fn flush(writer: Writer&) -> bam::Status {
    if writer.status != bam::Status.Record { return writer.status; }
    if writer.records > 0 {
        writer.container.clear(); writer.header.clear(); writer.section.clear();
        // Five preservation keys: names retained, absolute positions, no reference.
        itf8(writer.section, 5); writer.section.append("RN\x01AP\x00RR\x00SM\x1b\x1b\x1b\x1b\x1bTD");
        sized(writer.section, writer.dictionary); sized(writer.header, writer.section);
        writer.section.clear(); itf8(writer.section, 25);
        let keys = "BFCFRIRLAPRGRNMFNSNPTSTLFNFCFPDLMQBAQSINSCRSPDHCBB";
        for i in 25 {
            writer.section.append(keys[2 * i..2 * i + 2]);
            if i == Key.RN.index() || i == Key.IN.index() || i == Key.SC.index() {
                writer.section.append("\x05\x02\x00"); writer.section.push((i + 1) as u8);
            } else if i == Key.BB.index() { array_encoding(writer.section, i + 1); }
            else { external(writer.section, i + 1); }
        }
        sized(writer.header, writer.section);
        writer.section.clear(); itf8(writer.section, writer.tags.len);
        for tag in writer.tags { itf8(writer.section, tag as i64); array_encoding(writer.section, 26); }
        sized(writer.header, writer.section);
        writer.encode_block(1, 0, writer.header, false);
        let landmark = writer.container.len;
        writer.header.clear(); itf8(writer.header, -2); itf8(writer.header, 0); itf8(writer.header, 0);
        itf8(writer.header, writer.records); ltf8(writer.header, writer.counter);
        itf8(writer.header, 27); itf8(writer.header, 26);
        for i in 26 { itf8(writer.header, i + 1); }
        itf8(writer.header, -1);
        for i in 16 { writer.header.push(0); } // No external reference MD5.
        writer.encode_block(2, 0, writer.header, false);
        writer.encode_block(5, 0, "", false); // A core block exists even with no core bits.
        for i in 26 { writer.encode_block(4, i + 1, writer.stream(i), true); }
        writer.container_header(-2, 0, writer.records, writer.counter, writer.bases, 29, landmark);
        writer.send();
        writer.counter += writer.records; writer.reset_slice();
    }
    if writer.status == bam::Status.Record && !io::flush(writer.output) { writer.status = bam::Status.IoError; }
    writer.status
}

// Flush records and append the canonical 38-byte CRAM EOF container only on
// success. Always release the compressor/file. Repeated close succeeds.
//
// ```goose
// let writer = cram_writer::Writer { .. };
// assert(writer.close() == bam::Status.Record);
// assert(writer.close() == bam::Status.Record);
// ```
fn close(writer: Writer&) -> bam::Status {
    if writer.status == bam::Status.Record {
        writer.flush();
        if writer.status == bam::Status.Record {
            writer.container.clear(); writer.encode_block(1, 0, "\x01\x00\x01\x00\x01\x00", false);
            writer.container_header(-1, 0x454f46, 0, 0, 0, 1, -1); writer.send();
        }
    }
    var result = writer.status;
    if result != bam::Status.Record { writer.output.used = 0; }
    if !io::close(writer.output) && result == bam::Status.Record { result = bam::Status.IoError; }
    writer.encoder.close(); writer.status = bam::Status.Closed; writer.reset_slice();
    if result == bam::Status.Closed { bam::Status.Record } else { result }
}

// Create/truncate CRAM, writing the definition and SAM-header container eagerly.
// level 0 emits raw blocks; 1..12 use libdeflate gzip when smaller. minor_version
// is 0 or 1 (default 1). Header SQ names/lengths must match the target dictionary.
//
// ```goose
// let reader = bam::Reader { .. }; defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let writer = cram_writer::Writer { .. }; defer writer.close();
// let path = "tests/data/writer.cram"; defer delete_file(path);
// assert(writer.open(path, reader.header, reader.references, reader.names) == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(writer.write(record) == bam::Status.Record);
// assert(writer.close() == bam::Status.Record);
// let input = cram::Reader { .. }; defer input.close();
// let reference = ExampleReference {};
// assert(input.open(path) == bam::Status.Record && input.read(reference) == bam::Status.Record);
// let decoded = input.record();
// assert(decoded.name() == "match");
// ```
fn open(writer: Writer&, path: const u8[:], header: const u8[:], refs: const bam::Reference[:], names: const u8[:], level: i32 = 6, minor_version: i32 = 1) -> bam::Status {
    let previous = writer.close(); if previous != bam::Status.Record { return previous; }
    writer.status = bam::Status.Invalid;
    if level < 0 || level > 12 || minor_version < 0 || minor_version > 1 || writer.records_per_slice < 1 ||
       writer.records_per_slice > 2147483647 || writer.max_slice_bytes < 1 || writer.max_slice_bytes > 268435456 ||
       writer.max_record_bytes < 1 || writer.max_record_bytes > writer.max_slice_bytes { return writer.status; }
    writer.status = sam_writer::header_text(header, refs, names, writer.section);
    if writer.status != bam::Status.Record { return writer.status; }
    var parsed: bam::Reference[>..] = []; var parsed_names: u8[>..] = [];
    var at = 0;
    while at < writer.section.len {
        var end = at; while end < writer.section.len && writer.section[end] != '\n' { end += 1; }
        let line_end = if end > at && writer.section[end - 1] == '\r' { end - 1 } else { end };
        writer.status = sam::parse_header_line(writer.section[at..line_end], parsed, parsed_names, 1048576);
        if writer.status != bam::Status.Record { return writer.status; }
        at = end + 1;
    }
    if parsed.len != refs.len { writer.status = bam::Status.Invalid; return writer.status; }
    for r, i in refs {
        let p = parsed[i];
        if r.length != p.length || names[r.name..r.name + r.name_len] != parsed_names[p.name..p.name + p.name_len] {
            writer.status = bam::Status.Invalid; return writer.status;
        }
    }
    let encoder, allocated = deflate::compressor(level);
    if !allocated { writer.status = bam::Status.IoError; return writer.status; }
    writer.encoder = encoder;
    let file, ok = io::open_write(path);
    if !ok { writer.encoder.close(); writer.status = bam::Status.IoError; return writer.status; }
    writer.output = io::buffered_writer(file, 65536);
    writer.references = refs.len; writer.level = level; writer.counter = 0;
    var definition: u8[26] = [0; 26]; definition[..4].move_from("CRAM"); definition[4] = 3; definition[5] = minor_version as u8;
    if io::write(writer.output, definition) != definition.len { writer.status = bam::Status.IoError; return writer.status; }
    writer.header.clear(); word(writer.header, writer.section.len as u32); writer.header.append(writer.section);
    writer.container.clear(); writer.encode_block(0, 0, writer.header, true);
    writer.container_header(-1, 0, 0, 0, 0, 1, -1); writer.send();
    writer.status
}

fn encode_tags<T>(writer: Writer&, record: const T&) -> bam::Status {
    let aux = record.aux(); var at = 0; var skip = record.fields.cigar_len > 65535;
    writer.tag_line.clear();
    while at < aux.len {
        let end = bam::aux_end(aux, at); if end < 0 { return bam::Status.Invalid; }
        if skip && aux[at..at + 2] == "CG" { skip = false; at = end; continue; }
        // cF is reserved by CRAM for MD/NM suppression; double tags are not a
        // CRAM tag type. Reject rather than silently dropping/changing them.
        if aux[at..at + 2] == "cF" || aux[at + 2] == 'd' { return bam::Status.Invalid; }
        writer.tag_line.append(aux[at..at + 3]);
        let tag = (aux[at] as u32 << 16) | (aux[at + 1] as u32 << 8) | aux[at + 2];
        var found = false; for known in writer.tags { if tag == known { found = true; break; } }
        if !found { writer.tags.push(tag); }
        writer.array(Key.Tags, aux[at + 3..end]); at = end;
    }
    writer.tag_line.push(0);
    var line = writer.lines.len;
    for entry, i in writer.lines {
        if writer.dictionary[entry.begin..entry.begin + entry.length] == writer.tag_line { line = i; break; }
    }
    if line == writer.lines.len {
        writer.lines.push(TagLine { begin: writer.dictionary.len, length: writer.tag_line.len });
        writer.dictionary.append(writer.tag_line);
    }
    writer.put(Key.TL, line);
    bam::Status.Record
}

// Encode a validated BAM/slab record. CIGAR =/X normalize to M and adjacent
// equal operations coalesce, as CRAM decoders do. Names, bases, qualities, tags
// and detached mate fields are retained. Unsupported inputs (including mapped
// records without a usable CIGAR, unmapped CIGAR/nonzero MAPQ, cF/d tags and
// unpaired records carrying a mate target) fail instead of silently losing data.
//
// ```goose
// let writer = cram_writer::Writer { .. };
// let alignment = bam_encode::Alignment { flags: 4, .. };
// assert(writer.write(alignment) == bam::Status.Closed);
// ```
fn write<T>(writer: Writer&, record: const T&) -> bam::Status {
    if writer.status != bam::Status.Record { return writer.status; }
    let f = record.fields;
    let mapped = f.flags & 4 == 0;
    if f.tid < -1 || f.tid >= writer.references || f.mate_tid < -1 || f.mate_tid >= writer.references ||
       f.pos >= 2147483647 || f.mate_pos >= 2147483647 || record.name().len == 0 ||
       (mapped && (f.tid < 0 || f.pos < 0 || f.cigar_len == 0)) ||
       (!mapped && (f.cigar_len != 0 || f.mapq != 0)) || (f.flags & 1 == 0 && (f.mate_tid >= 0 || f.flags & 32 != 0)) {
        writer.status = bam::Status.Invalid; return writer.status;
    }
    if sam_writer::payload_size(record) > writer.max_record_bytes {
        writer.status = bam::Status.TooLarge; return writer.status;
    }
    var query = 0; var features = 0;
    for i in f.cigar_len {
        let op = record.cigar(i) & 15; let n = record.cigar(i) >> 4;
        if n == 0 { continue; }
        features += 1;
        if op == 0 || op == 1 || op == 4 || op == 7 || op == 8 { query += n; }
    }
    let length = if mapped && f.seq_len == 0 { query } else { f.seq_len };
    if mapped && f.seq_len > 0 && query != f.seq_len {
        writer.status = bam::Status.Invalid; return writer.status;
    }
    // Account for absent SEQ using its CIGAR-derived length before allocating
    // dummy bases/qualities. Reuse this upper bound instead of rescanning all
    // 26 stream lengths after every record.
    let estimate = 512 + length * 4 + f.cigar_len * 16 + record.aux().len * 3 + record.name().len;
    if length > writer.max_record_bytes || estimate > writer.max_slice_bytes {
        writer.status = bam::Status.TooLarge; return writer.status;
    }
    if writer.records > 0 && (writer.records >= writer.records_per_slice || estimate > writer.max_slice_bytes - writer.pending_bytes) {
        writer.flush(); if writer.status != bam::Status.Record { return writer.status; }
    }
    writer.sequence.clear(); sam_writer::sequence(record, writer.sequence);
    if f.seq_len == 0 { writer.sequence.resize(length, 'N'); }
    writer.put(Key.BF, f.flags); writer.put(Key.CF, 3 | (if f.seq_len == 0 { 8 } else { 0 }));
    writer.put(Key.RI, f.tid); writer.put(Key.RL, length); writer.put(Key.AP, f.pos + 1); writer.put(Key.RG, -1);
    writer.rn.append(record.name()); writer.rn.push(0);
    writer.put(Key.MF, (if f.flags & 32 != 0 { 1 } else { 0 }) | (if f.flags & 8 != 0 { 2 } else { 0 }));
    writer.put(Key.NS, f.mate_tid); writer.put(Key.NP, f.mate_pos + 1); writer.put(Key.TS, f.template_len);
    writer.status = writer.encode_tags(record); if writer.status != bam::Status.Record { return writer.status; }
    if mapped {
        writer.put(Key.FN, features);
        var pos = 0; var previous = 0;
        for i in f.cigar_len {
            let word = record.cigar(i); let n = (word >> 4) as i64; let op = word & 15;
            if n == 0 { continue; }
            let code: u8 = match op { 1 => 'I', 2 => 'D', 3 => 'N', 4 => 'S', 5 => 'H', 6 => 'P', _ => 'b' };
            writer.fc.push(code); writer.put(Key.FP, pos + 1 - previous); previous = pos + 1;
            if code == 'b' { writer.array(Key.BB, writer.sequence[pos..pos + n]); pos += n; }
            else if code == 'I' || code == 'S' {
                let key = if code == 'I' { Key.IN } else { Key.SC };
                writer.stream(key.index()).append(writer.sequence[pos..pos + n]); writer.stream(key.index()).push(0); pos += n;
            } else { writer.put(if code == 'D' { Key.DL } else if code == 'N' { Key.RS } else if code == 'H' { Key.HC } else { Key.PD }, n); }
        }
        writer.put(Key.MQ, f.mapq);
    } else { writer.ba.append(writer.sequence); }
    if f.seq_len == 0 { writer.qs.resize(writer.qs.len + length, 255); }
    else { writer.qs.append(sam_writer::qualities(record)); }
    writer.records += 1; writer.bases += length;
    writer.pending_bytes += estimate;
    writer.status
}

// Validate and encode caller-owned fields using reusable BAM scratch storage.
//
// ```goose
// let writer = cram_writer::Writer { .. };
// let alignment = bam_encode::Alignment { name: "r", flags: 4, sequence: "ACN", .. };
// assert(writer.write(alignment) == bam::Status.Closed);
// ```
fn write(writer: Writer&, alignment: const bam_encode::Alignment&) -> bam::Status {
    if writer.status != bam::Status.Record { return writer.status; }
    writer.status = bam_encode::encode(alignment, writer.payload, writer.references, writer.max_record_bytes);
    if writer.status != bam::Status.Record { return writer.status; }
    let fields, ok = bam::decode(writer.payload, writer.references);
    if !ok { writer.status = bam::Status.Invalid; return writer.status; }
    let record = bam::Record { fields: fields, bytes: writer.payload }; writer.write(record)
}
