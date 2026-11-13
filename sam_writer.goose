// Seqair sam/format.rs semantics, with reusable buffered output.
namespace sam_writer;
import std;
import binary;
import bam;
import bam_encode;
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";
import pileup;
import io from "github.com/sstadick/goose-io@v0.3.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

extern fn float_text(value: f64, array: bool, output: u8[:]) -> i64 = "goose_seqair_sam_float";
extern fn sequence_text(input: const u8[:], output: u8[:]) = "goose_seqair_sequence_text";

// Owns a plain SAM output and reusable buffers. Header/target names are copied
// at open; input record views are never retained. Check close for write errors.
// max_record_bytes bounds encoded BAM input; max_line_bytes bounds SAM output.
//
// ```goose
// let writer = sam_writer::Writer { .. };
// assert(writer.close() == bam::Status.Ok);
// ```
root struct Writer {
    output: io::BufferedWriter = io::BufferedWriter { file: io::FileHandle { fd: -1 }, used: 0, buffer: [] },
    status: bam::Status = bam::Status.Closed,
    references: bam::Reference[>..] = [], names: u8[>..] = [],
    line: u8[>..] = [], payload: u8[>..] = [],
    max_record_bytes: i64 = 16777216, max_line_bytes: i64 = 67108864,
}

// Flush and close. Errors poison writes until close/open; closing twice succeeds.
// A failed output may contain complete earlier records, never a success signal.
//
// ```goose
// let writer = sam_writer::Writer { .. };
// assert(writer.close() == bam::Status.Ok);
// assert(writer.close() == bam::Status.Ok);
// ```
fn close(writer: Writer&) -> bam::Status {
    var result = writer.status;
    if result != bam::Status.Ok { writer.output.used = 0; }
    if !io::close(writer.output) && result == bam::Status.Ok { result = bam::Status.IoError; }
    writer.status = bam::Status.Closed;
    if result == bam::Status.Closed { bam::Status.Ok } else { result }
}

// Seqair header formatting: stop at NUL, add a missing newline, synthesize SQ
// only when no SQ lines exist. Used by CRAM too, which embeds a SAM header.
fn header_text(header: const u8[:], refs: const bam::Reference[:], names: const u8[:], out: u8[>..]&) -> bam::Status {
    out.clear();
    if header.len > 16777216 || names.len > 16777216 || refs.len > 1048576 { return bam::Status.TooLarge; }
    for r in refs {
        if r.name < 0 || r.name_len < 1 || r.name > names.len || r.name_len > names.len - r.name ||
           r.length < 1 || r.length > 2147483647 { return bam::Status.Invalid; }
        for b in names[r.name..r.name + r.name_len] { if b <= 32 || b == 127 { return bam::Status.Invalid; } }
    }
    let nul = bstr::find_byte(header, 0);
    let text = header[..if nul < 0 { header.len } else { nul }];
    var sq = false; var at = 0;
    while at < text.len {
        let nl = bstr::find_byte(text[at..], '\n');
        let end = if nl < 0 { text.len } else { at + nl };
        if end == at || text[at] != '@' { return bam::Status.Invalid; }
        if end - at >= 4 && text[at..at + 4] == "@SQ\t" { sq = true; }
        at = end + 1;
    }
    out.append(text);
    if text.len > 0 && text[text.len - 1] != '\n' { out.push('\n'); }
    if !sq {
        for r in refs { format(out, "@SQ\tSN:", names[r.name..r.name + r.name_len], "\tLN:", r.length, "\n"); }
    }
    if out.len > 16777216 { return bam::Status.TooLarge; }
    bam::Status.Ok
}

// Create/truncate plain SAM and write its header eagerly. Target IDs in records
// must match this dictionary; they are not remapped. No output is opened until
// header validation succeeds. Header text follows Seqair's SAM formatter.
//
// ```goose
// let reader = bam::Reader { .. }; defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Ok);
// let writer = sam_writer::Writer { .. }; defer writer.close();
// let path = "tests/data/writer.sam"; defer delete_file(path);
// assert(writer.open(path, reader.header, reader.references, reader.names) == bam::Status.Ok);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(writer.write(record) == bam::Status.Ok);
// assert(writer.close() == bam::Status.Ok);
// let input = sam::Reader { .. }; defer input.close();
// assert(input.open(path) == bam::Status.Ok && input.read() == bam::Status.Record);
// let decoded = input.record();
// assert(decoded.name() == "match");
// ```
fn open(writer: Writer&, path: const u8[:], header: const u8[:], refs: const bam::Reference[:], names: const u8[:]) -> bam::Status {
    let previous = writer.close(); if previous != bam::Status.Ok { return previous; }
    if writer.max_record_bytes < 32 || writer.max_line_bytes < 1 { writer.status = bam::Status.Invalid; return writer.status; }
    writer.status = header_text(header, refs, names, writer.line);
    if writer.status != bam::Status.Ok { return writer.status; }
    writer.references.clear(); writer.references.append(refs);
    writer.names.clear(); writer.names.append(names);
    let file, ok = io::open_write(path);
    if !ok { writer.status = bam::Status.IoError; return writer.status; }
    writer.output = io::buffered_writer(file, 65536);
    if io::write(writer.output, writer.line) != writer.line.len { writer.status = bam::Status.IoError; }
    writer.status
}

fn target(refs: const bam::Reference[:], names: const u8[:], tid: i32) -> const u8[:] {
    if tid < 0 { return "*"; }
    let r = refs[tid as i64]; names[r.name..r.name + r.name_len]
}

fn number(out: u8[>..]&, bytes: const u8[:], at: i64, kind: u8, array: bool) -> bool {
    if kind == 'f' || kind == 'd' {
        let value = if kind == 'f' { binary::f32_of_bits(bgzf::le32(bytes, at)) as f64 }
            else { binary::f64_of_bits(binary::u64le(bytes, at)) };
        var text: u8[64] = [0; 64];
        let n = float_text(value, array, text);
        if n < 0 || n >= text.len { return false; }
        out.append(text[..n]);
    } else {
        let n: i64 = match kind {
            'c' => bytes[at] as! i8 as i64, 'C' => bytes[at] as i64,
            's' => bgzf::le16(bytes, at) as u16 as! i16 as i64,
            'S' => bgzf::le16(bytes, at) as i64,
            'i' => bgzf::le32(bytes, at) as! i32 as i64,
            _ => bgzf::le32(bytes, at) as i64,
        };
        format(out, n);
    }
    true
}

fn format_aux(aux: const u8[:], skip_cg: bool, out: u8[>..]&, limit: i64) -> bam::Status {
    var at = 0; var skip = skip_cg;
    while at < aux.len {
        let end = bam::aux_end(aux, at);
        if end < 0 { return bam::Status.Invalid; }
        if skip && aux[at..at + 2] == "CG" { skip = false; at = end; continue; }
        out.push('\t'); out.append(aux[at..at + 2]); out.push(':');
        let kind = aux[at + 2]; at += 3;
        if kind == 'Z' || kind == 'H' || kind == 'A' {
            out.push(kind); out.push(':'); out.append(aux[at..if kind == 'A' { end } else { end - 1 }]);
        } else if kind == 'B' {
            let sub = aux[at]; out.append("B:"); out.push(sub);
            let width = if sub == 'c' || sub == 'C' { 1 } else if sub == 's' || sub == 'S' { 2 } else { 4 };
            at += 5;
            while at < end {
                out.push(','); if !number(out, aux, at, sub, true) { return bam::Status.Invalid; }
                at += width;
                if out.len > limit { return bam::Status.TooLarge; }
            }
        } else {
            out.push(if kind == 'f' || kind == 'd' { kind } else { 'i' }); out.push(':');
            if !number(out, aux, at, kind, false) { return bam::Status.Invalid; }
        }
        if out.len > limit { return bam::Status.TooLarge; }
        at = end;
    }
    bam::Status.Ok
}

fn sequence(record: const bam::Record&, out: u8[>..]&) {
    let start = out.len; out.resize(start + record.fields.seq_len, 0);
    sequence_text(record.bytes[record.fields.seq..record.fields.qual], out[start..]);
}
fn sequence(record: const pileup::Record&, out: u8[>..]&) {
    for base in record.bases { out.push("=ACMGRSVTWYHKDBN"[base as i64]); }
}

fn payload_size(record: const bam::Record&) -> i64 { record.bytes.len }
fn payload_size(record: const pileup::Record&) -> i64 {
    32 + record.name().len + 1 + record.fields.cigar_len * 4 + (record.fields.seq_len + 1) / 2 +
        record.fields.seq_len + record.aux().len
}

fn qualities(record: const bam::Record&) -> const u8[:] { record.bytes[record.fields.qual..record.fields.aux] }
fn qualities(record: const pileup::Record&) -> const u8[:] { record.quality_bytes }

// Format and buffer one validated BAM or slab record in stored orientation.
// Long CIGARs are expanded and their first CG escape omitted. Aux order is
// retained; integer tag widths become SAM i, floats follow Seqair/htslib's
// six-significant-digit conventions. Invalid/oversized lines are not written.
//
// ```goose
// let writer = sam_writer::Writer { .. };
// let alignment = bam_encode::Alignment { flags: 4, .. };
// assert(writer.write(alignment) == bam::Status.Closed);
// ```
fn write<T>(writer: Writer&, record: const T&) -> bam::Status {
    if writer.status != bam::Status.Ok { return writer.status; }
    let f = record.fields;
    if f.tid < -1 || f.tid >= writer.references.len || f.mate_tid < -1 || f.mate_tid >= writer.references.len {
        writer.status = bam::Status.Invalid; return writer.status;
    }
    if payload_size(record) > writer.max_record_bytes {
        writer.status = bam::Status.TooLarge; return writer.status;
    }
    writer.line.clear();
    format(writer.line, record.name(), "\t", f.flags, "\t", target(writer.references, writer.names, f.tid),
        "\t", f.pos + 1, "\t", f.mapq, "\t");
    if f.cigar_len == 0 { writer.line.push('*'); }
    for i in f.cigar_len {
        let word = record.cigar(i); format(writer.line, word >> 4); writer.line.push("MIDNSHP=X"[(word & 15) as i64]);
        if writer.line.len >= writer.max_line_bytes { writer.status = bam::Status.TooLarge; return writer.status; }
    }
    format(writer.line, "\t", if f.mate_tid >= 0 && f.mate_tid == f.tid { "=" } else { target(writer.references, writer.names, f.mate_tid) },
        "\t", f.mate_pos + 1, "\t", f.template_len, "\t");
    if f.seq_len == 0 { writer.line.append("*\t*"); }
    else {
        sequence(record, writer.line); writer.line.push('\t');
        if record.quality(0) == 255 { writer.line.push('*'); }
        else {
            let start = writer.line.len; writer.line.resize(start + f.seq_len, 0);
            for q, i in qualities(record) { writer.line[start + i] = q + 33; }
        }
    }
    writer.status = format_aux(record.aux(), f.cigar_len > 65535, writer.line, writer.max_line_bytes);
    if writer.status != bam::Status.Ok { writer.line.clear(); return writer.status; }
    writer.line.push('\n');
    if writer.line.len > writer.max_line_bytes { writer.status = bam::Status.TooLarge; return writer.status; }
    if io::write(writer.output, writer.line) != writer.line.len { writer.status = bam::Status.IoError; }
    writer.status
}

// Serialize caller-owned alignment fields, validating through the BAM encoder.
// Uses a reusable payload buffer; no borrowed field survives this call.
//
// ```goose
// let writer = sam_writer::Writer { .. };
// let alignment = bam_encode::Alignment { name: "r", flags: 4, sequence: "ACN", .. };
// assert(writer.write(alignment) == bam::Status.Closed);
// ```
fn write(writer: Writer&, alignment: const bam_encode::Alignment&) -> bam::Status {
    if writer.status != bam::Status.Ok { return writer.status; }
    writer.status = bam_encode::encode(alignment, writer.payload, writer.references.len, writer.max_record_bytes);
    if writer.status != bam::Status.Ok { return writer.status; }
    let fields, ok = bam::decode(writer.payload, writer.references.len);
    if !ok { writer.status = bam::Status.Invalid; return writer.status; }
    let record = bam::Record { fields: fields, bytes: writer.payload };
    writer.write(record)
}
