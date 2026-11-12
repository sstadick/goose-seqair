// SAM text feeds the same packed records and decoded slabs as BAM.
namespace sam;
import std;
import binary;
import bam;
import bam_encode;
import bgzf;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

// Owns bounded line storage, the header dictionary and one reusable BAM
// payload. Construct with Reader { .. }; close explicitly. Plain and BGZF SAM
// are detected by content, not extension. Ordinary gzip is not supported yet.
//
// ```goose
// let reader = sam::Reader { .. };
// assert(reader.read() == bam::Status.Closed);
// reader.close();
// ```
root struct Reader {
    stream: bgzf::Stream,
    compressed: bool = false,
    status: bam::Status = bam::Status.Closed,
    header: u8[>..] = [], names: u8[>..] = [], references: bam::Reference[>..] = [],
    buffer: u8[>..] = [], start: i64 = 0, end: i64 = 0,
    pending: bool = false, line_start: i64 = 0, line_end: i64 = 0,
    line_number: i64 = 0, last_tid: i32 = -1,
    cigar: u32[>..] = [], aux: u8[>..] = [], qualities: u8[>..] = [],
    bytes: u8[>..] = [], fields: bam::Fields,
    max_line_bytes: i64 = 16777216, max_record_bytes: i64 = 16777216,
    max_header_bytes: i64 = 16777216, max_references: i64 = 1048576,
}

// Close input; record views expire, but the header remains until reopen.
//
// ```goose
// let reader = sam::Reader { .. };
// reader.close(); reader.close();
// assert(reader.status == bam::Status.Closed);
// ```
fn close(reader: Reader&) {
    reader.stream.close(); reader.status = bam::Status.Closed;
    reader.bytes.clear(); reader.start = 0; reader.end = 0; reader.pending = false;
}

// Open and read leading header lines, preserving the first alignment for
// read(). Empty/header-only files open successfully. Errors close the input.
//
// ```goose
// let path = "tests/data/doc.sam";
// assert(write_file(path, "@SQ\tSN:chrA\tLN:100\nr\t0\tchrA\t6\t40\t4M\t*\t0\t0\tACGT\tIIII\n"));
// defer delete_file(path);
// let reader = sam::Reader { .. };
// defer reader.close();
// assert(reader.open(path) == bam::Status.Record);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(record.name() == "r" && record.fields.pos == 5);
// assert(record.base(2) == 'G' && record.quality(2) == 40);
// assert(reader.open(path) == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 6, 8) == bam::Status.Record);
// var columns = 0;
// store.columns() { pos, entries => assert(entries.len == 1); columns += 1; };
// assert(columns == 2);
// ```
fn open(reader: Reader&, path: const u8[:]) -> bam::Status {
    reader.close(); reader.header.clear(); reader.names.clear(); reader.references.clear();
    reader.line_number = 0; reader.last_tid = -1;
    reader.status = bam::from_bgzf(reader.stream.open(path));
    if reader.status != bam::Status.Record { return reader.status; }
    var magic: u8[2] = [0; 2];
    let count = io::read(reader.stream.input, magic);
    if count < 0 || bgzf::seek_raw(reader.stream.input.fd, 0) != 0 {
        reader.stream.close(); reader.status = bam::Status.IoError; return reader.status;
    }
    reader.compressed = count == 2 && magic == "\x1f\x8b";
    reader.buffer.resize(min(65536, max(2, reader.max_line_bytes + 2)), 0);
    while true {
        let status = reader.next_line();
        if status == bam::Status.End { break; }
        if status != bam::Status.Record { reader.status = status; break; }
        let line = reader.buffer[reader.line_start..reader.line_end];
        if line.len == 0 || line[0] != '@' { reader.pending = true; break; }
        if line.len + 1 > reader.max_header_bytes - reader.header.len {
            reader.status = bam::Status.TooLarge; break;
        }
        reader.status = parse_header_line(line, reader.references, reader.names, reader.max_references);
        if reader.status != bam::Status.Record { break; }
        reader.header.append(line); reader.header.push('\n');
    }
    if reader.status != bam::Status.Record { reader.stream.close(); }
    reader.status
}

fn input(reader: Reader&, out: u8[:]) -> i64 {
    if !reader.compressed { return io::read(reader.stream.input, out); }
    if reader.stream.begin == reader.stream.end {
        let status = reader.stream.block_read();
        if status == bgzf::Status.End { return 0; }
        if status != bgzf::Status.Ready { reader.status = bam::from_bgzf(status); return -1; }
    }
    let count = min(out.len, reader.stream.end - reader.stream.begin);
    out[..count].move_from(reader.stream.decoded[reader.stream.begin..reader.stream.begin + count]);
    reader.stream.begin += count;
    count
}

fn next_line(reader: Reader&) -> bam::Status {
    var scan = reader.start;
    loop {
        let newline = bstr::find_byte(reader.buffer[scan..reader.end], '\n');
        if newline >= 0 {
            reader.line_start = reader.start; reader.line_end = scan + newline;
            reader.start = reader.line_end + 1;
            if reader.line_end > reader.line_start && reader.buffer[reader.line_end - 1] == '\r' { reader.line_end -= 1; }
            reader.line_number += 1;
            return if reader.line_end - reader.line_start <= reader.max_line_bytes { bam::Status.Record } else { bam::Status.TooLarge };
        }
        let remaining = reader.end - reader.start;
        // A CR at the limit may be followed by LF in the next BGZF block.
        // Neither byte belongs to the logical line length.
        if remaining > reader.max_line_bytes &&
           !(remaining == reader.max_line_bytes + 1 && reader.buffer[reader.end - 1] == '\r') {
            return bam::Status.TooLarge;
        }
        if reader.start > 0 { reader.buffer[..remaining].move_from(reader.buffer[reader.start..reader.end]); }
        reader.start = 0; reader.end = remaining; scan = remaining;
        if reader.end == reader.buffer.len { reader.buffer.resize(min(reader.buffer.len * 2, reader.max_line_bytes + 2), 0); }
        let count = reader.input(reader.buffer[remaining..]);
        if count < 0 { return if reader.status == bam::Status.Record { bam::Status.IoError } else { reader.status }; }
        if count == 0 {
            reader.line_start = 0; reader.line_end = remaining; reader.start = remaining;
            if remaining == 0 { return bam::Status.End; }
            if reader.buffer[remaining - 1] == '\r' { reader.line_end -= 1; }
            reader.line_number += 1; return bam::Status.Record;
        }
        reader.end += count;
    }
}

fn parse_header_line(line: const u8[:], refs: bam::Reference[>..]&, names: u8[>..]&, max_refs: i64) -> bam::Status {
    if bstr::find_byte(line, 0) >= 0 || !bstr::is_utf8(line) { return bam::Status.Invalid; }
    if line.len < 3 || line[..3] != "@SQ" { return bam::Status.Record; }
    var name = line[..0];
    var length: i64 = -1;
    var at: i64 = 4;
    while at <= line.len {
        let tab = bstr::find_byte(line[at..], '\t');
        let end = if tab < 0 { line.len } else { at + tab };
        let field = line[at..end];
        if field.len >= 3 && field[..3] == "SN:" { name = field[3..]; }
        if field.len >= 3 && field[..3] == "LN:" { length = digits(field[3..], 2147483647); }
        at = end + 1;
    }
    if name.len == 0 || length < 1 || find_reference(refs, names, name) >= 0 { return bam::Status.Invalid; }
    if refs.len >= max_refs { return bam::Status.TooLarge; }
    refs.push(bam::Reference { name: names.len, name_len: name.len, length: length });
    names.append(name); names.push(0);
    bam::Status.Record
}

fn find_reference(refs: const bam::Reference[:], names: const u8[:], name: const u8[:]) -> i32 {
    for ref, i in refs { if names[ref.name..ref.name + ref.name_len] == name { return i as i32; } }
    -1
}

// Find a reference by exact bytes, or -1 if absent.
//
// ```goose
// let reader = sam::Reader { .. };
// assert(reader.reference_id("missing") == -1);
// ```
fn reference_id(reader: const Reader&, name: const u8[:]) -> i32 { find_reference(reader.references, reader.names, name) }

// Decode the next line into reused BAM storage. Record views expire on the
// next read/open/close. Unlike Seqair's recoverable malformed-line error, all
// errors remain sticky here, matching bam::Reader's contract.
//
// ```goose
// let reader = sam::Reader { .. };
// assert(reader.read() == bam::Status.Closed);
// ```
fn read(reader: Reader&) -> bam::Status {
    if reader.status != bam::Status.Record { return reader.status; }
    reader.bytes.clear();
    let status = if reader.pending { bam::Status.Record } else { reader.next_line() };
    reader.pending = false;
    if status != bam::Status.Record { reader.status = status; return status; }
    reader.status = reader.encode(reader.buffer[reader.line_start..reader.line_end]);
    if reader.status != bam::Status.Record { reader.bytes.clear(); }
    reader.status
}

// Borrow the current validated packed alignment, using bam::Record accessors.
// Call only after read() returns Record, not merely after open().
//
// ```goose
// let reader = sam::Reader { .. };
// assert(reader.status == bam::Status.Closed); // No record is available yet.
// ```
fn record(reader: const Reader&) -> bam::Record {
    assert(reader.status == bam::Status.Record && reader.bytes.len >= 32);
    bam::Record { fields: reader.fields, bytes: reader.bytes }
}

fn digits(text: const u8[:], limit: i64) -> i64 {
    if text.len == 0 { return -1; }
    var n: i64 = 0;
    for byte in text {
        if byte < '0' || byte > '9' { return -1; }
        let d = byte - '0';
        if n > (limit - d) / 10 { return -1; }
        n = n * 10 + d;
    }
    return if n <= limit { n } else { -1 };
}

fn integer(text: const u8[:], low: i64, high: i64) -> i64, bool {
    let negative = text.len > 0 && text[0] == '-';
    let sign = text.len > 0 && (negative || text[0] == '+');
    let value = digits(text[if sign { 1 } else { 0 }..], if negative { max(0, -low) } else { high });
    let n = if negative { -value } else { value };
    return n, value >= 0 && n >= low && n <= high && !(negative && low == 0);
}

fn encode(reader: Reader&, line: const u8[:]) -> bam::Status {
    if bstr::find_byte(line, 0) >= 0 { return bam::Status.Invalid; }
    var fields: const u8[:][11] = [""; 11];
    var at: i64 = 0;
    for i in 11 {
        let tab = bstr::find_byte(line[at..], '\t');
        if tab < 0 && i < 10 { return bam::Status.Invalid; }
        let end = if tab < 0 { line.len } else { at + tab };
        fields[i] = line[at..end]; at = end + 1;
    }
    let flag = digits(fields[1], 65535);
    let pos = digits(fields[3], 2147483648);
    let mapq = digits(fields[4], 255);
    let mate_pos = digits(fields[7], 2147483648);
    let tlen, valid_tlen = integer(fields[8], -2147483648, 2147483647);
    if flag < 0 || (fields[1].len > 1 && fields[1][0] == '0') || pos < 0 || mapq < 0 || mate_pos < 0 || !valid_tlen {
        return bam::Status.Invalid;
    }
    var tid = reader.last_tid;
    if tid < 0 || reader.names[reader.references[tid as i64].name..reader.references[tid as i64].name + reader.references[tid as i64].name_len] != fields[2] {
        tid = reader.reference_id(fields[2]); reader.last_tid = tid;
    }
    if tid < 0 && fields[2] != "*" { return bam::Status.Invalid; }
    if pos == 0 { tid = -1; }
    var flags = flag as u16;
    if tid < 0 || fields[5] == "*" { flags |= 4; }
    let mate_tid = if mate_pos == 0 || fields[6] == "*" { -1 } else if fields[6] == "=" { tid } else { reader.reference_id(fields[6]) };
    reader.cigar.clear();
    let text = fields[5];
    if text != "*" {
        if text.len == 0 { return bam::Status.Invalid; }
        var start: i64 = 0;
        for ch, i in text {
            if ch >= '0' && ch <= '9' { continue; }
            let len = digits(text[start..i], 268435455);
            let op = bstr::find_byte("MIDNSHP=X", ch);
            if len < 0 || op < 0 { return bam::Status.Invalid; }
            reader.cigar.push((len as u32) << 4 | op as u32); start = i + 1;
        }
        if start != text.len { return bam::Status.Invalid; }
    }
    let seq = if fields[9] == "*" { "" } else { fields[9] };
    reader.qualities.clear();
    if fields[10] != "*" {
        if fields[10].len != seq.len { return bam::Status.Invalid; }
        for q in fields[10] { if q < 33 || q > 160 { return bam::Status.Invalid; } reader.qualities.push(q - 33); }
    }
    reader.aux.clear();
    var cg_start: i64 = -1;
    var cg_end: i64 = -1;
    while at < line.len {
        let tab = bstr::find_byte(line[at..], '\t');
        let end = if tab < 0 { line.len } else { at + tab };
        let tag_start = reader.aux.len;
        if !parse_aux(line[at..end], reader.aux) { return bam::Status.Invalid; }
        if cg_start < 0 && line[at..at + 2] == "CG" { cg_start = tag_start; cg_end = reader.aux.len; }
        at = end + 1;
    }
    // Seqair/htslib promote the first eligible CG array even when it is short.
    if tid >= 0 && pos > 0 && reader.cigar.len > 0 && reader.cigar[0] == ((seq.len as u32) << 4 | 4) && cg_start >= 0 {
        let tag = reader.aux[cg_start..cg_end];
        if tag.len >= 8 && tag[2] == 'B' && (tag[3] == 'I' || tag[3] == 'i') {
            let count = bgzf::le32(tag, 4) as i64;
            if count >= reader.cigar.len && count < 536870912 {
                reader.cigar.clear();
                for i in count { reader.cigar.push(bgzf::le32(tag, 8 + 4 * i)); }
                let remaining = reader.aux.len - cg_end;
                reader.aux[cg_start..cg_start + remaining].move_from(reader.aux[cg_end..]);
                reader.aux.resize(cg_start + remaining);
            }
        }
    }
    if seq.len > 0 && reader.cigar.len > 0 {
        var query: i64 = 0;
        for word in reader.cigar {
            let op = word & 15;
            if op == 0 || op == 1 || op == 4 || op == 7 || op == 8 { query += word >> 4; }
        }
        if query != seq.len { return bam::Status.Invalid; }
    }
    let a = bam_encode::Alignment { name: fields[0], tid: tid, pos: pos - 1, flags: flags, mapq: mapq as u8,
        mate_tid: mate_tid, mate_pos: mate_pos - 1, template_len: tlen as i32,
        cigar: reader.cigar, sequence: seq, qualities: reader.qualities, aux: reader.aux };
    let status = bam_encode::encode_body(a, reader.bytes, reader.references.len, reader.max_record_bytes, bam_encode::SequenceEncoding.Sam);
    if status != bam::Status.Record { return status; }
    let f, ok = bam::decode(reader.bytes, reader.references.len);
    if !ok { return bam::Status.Invalid; }
    reader.fields = f;
    bam::Status.Record
}

fn parse_float(text: const u8[:]) -> f64, bool {
    var lower: u8[>..] = [];
    let negative = text.len > 0 && text[0] == '-';
    let sign = text.len > 0 && (negative || text[0] == '+');
    if text.len >= 3 && text.len <= 9 {
        for ch in text[if sign { 1 } else { 0 }..] { lower.push(if ch >= 'A' && ch <= 'Z' { ch + 32 } else { ch }); }
        if lower == "inf" || lower == "infinity" || lower == "nan" {
            let bits: u64 = if lower == "nan" { 0x7ff8000000000000 } else { 0x7ff0000000000000 };
            return binary::f64_of_bits(bits | (if negative { (1 as u64) << 63 } else { 0 })), true;
        }
    }
    let number, ok = parse_flt(text);
    return number, ok;
}

fn element(kind: u8, text: const u8[:], out: u8[>..]&) -> bool {
    if kind == 'f' || kind == 'd' {
        let n, ok = parse_float(text); if !ok { return false; }
        // Goose's supported targets use little-endian IEEE floats.
        if kind == 'f' { let values: f32[1] = [n as f32]; out.append(bytes_of(values)); }
        else { let values: f64[1] = [n]; out.append(bytes_of(values)); }
        return true;
    }
    let low: i64 = match kind { 'c' => -128, 's' => -32768, 'i' => -2147483648, _ => 0 };
    let high: i64 = match kind { 'c' => 127, 'C' => 255, 's' => 32767, 'S' => 65535, 'i' => 2147483647, _ => 4294967295 };
    let n, ok = integer(text, low, high); if !ok { return false; }
    if kind == 'c' || kind == 'C' { out.push((n & 255) as u8); }
    else if kind == 's' || kind == 'S' { append_word(out, n as! u32, 2); }
    else { append_word(out, n as! u32, 4); }
    true
}

fn append_word(out: u8[>..]&, value: u32, size: i64) {
    for i in size { out.push(((value >> (8 * i)) & 255) as u8); }
}

fn parse_aux(field: const u8[:], out: u8[>..]&) -> bool {
    if field.len < 5 || field[2] != ':' || field[4] != ':' || field[0] < 33 || field[0] > 127 || field[1] < 33 || field[1] > 127 { return false; }
    out.append(field[..2]);
    let kind = field[3]; let value = field[5..];
    if kind == 'A' || kind == 'a' || kind == 'c' || kind == 'C' {
        if value.len != 1 || value[0] < 10 || value[0] > 127 { return false; }
        out.push('A'); out.append(value);
    } else if kind == 'Z' || kind == 'H' {
        if kind == 'H' {
            if value.len & 1 != 0 { return false; }
            for ch in value {
                if !((ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'f') || (ch >= 'A' && ch <= 'F')) { return false; }
            }
        }
        out.push(kind); out.append(value); out.push(0);
    } else if kind == 'i' || kind == 'I' {
        let n, ok = integer(value, -2147483648, 4294967295); if !ok { return false; }
        let negative = value.len > 0 && value[0] == '-';
        let code: u8 = if negative { if n >= -128 { 'c' } else if n >= -32768 { 's' } else { 'i' } }
            else { if n <= 255 { 'C' } else if n <= 65535 { 'S' } else { 'I' } };
        out.push(code); if !element(code, value, out) { return false; }
    } else if kind == 'f' || kind == 'd' {
        out.push(kind); if !element(kind, value, out) { return false; }
    } else if kind == 'B' {
        if value.len == 0 || bstr::find_byte("cCsSiIf", value[0]) < 0 || (value.len > 1 && value[1] != ',') { return false; }
        out.push('B'); out.push(value[0]);
        let count_at = out.len; append_word(out, 0, 4);
        var count: u32 = 0; var at: i64 = 2;
        while at <= value.len {
            let comma = bstr::find_byte(value[at..], ',');
            let end = if comma < 0 { value.len } else { at + comma };
            if !element(value[0], value[at..end], out) { return false; }
            count += 1; at = end + 1;
        }
        bgzf::put32(out, count_at, count);
    } else { return false; }
    true
}
