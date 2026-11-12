// Slice metadata and reconstruction helpers, following Seqair slice.rs.
namespace cram_slice;
import std;
import cram_encoding;
import bam;
import pileup;

struct Header {
    tid: i32, start: i64, span: i64, records: i64, counter: i64,
    blocks: i64, embedded: i32, md5: u8[16],
}

fn parse(data: const u8[:]) -> Header, bool {
    let c = cram_encoding::Cursor { bytes: data };
    let h = Header { tid: c.itf8(), start: c.itf8() as i64, span: c.itf8() as i64,
                     records: c.itf8() as i64, counter: c.ltf8(), blocks: c.itf8() as i64, .. };
    let n = c.itf8() as i64;
    if n < 0 || n > data.len { return h, false; }
    for i in n { c.itf8(); }
    h.embedded = c.itf8();
    let md5 = c.take(16);
    if c.ok { h.md5[..].move_from(md5); }
    return h, c.ok && h.tid >= -2 && h.records >= 0 && h.counter >= 0 && h.blocks >= 0;
}

// A compact handle to the compression header and the slice's shared streams.
struct Series {
    header: const cram_encoding::Compression&,
    table: const cram_encoding::Table&,
    context: cram_encoding::Context&,
}
type Key = cram_encoding::Series;

inline fn value(ds: Series&, key: Key) -> i32 {
    ds.table.integer(ds.header.series[key.index()], ds.context,
        key == cram_encoding::Series.FC || key == cram_encoding::Series.BA ||
        key == cram_encoding::Series.QS || key == cram_encoding::Series.BS)
}

fn array(ds: Series&, key: Key, out: u8[>..]&) {
    ds.table.array(ds.header.series[key.index()], ds.context, out);
}

fn bytes(ds: Series&, key: Key, n: i64, out: u8[>..]&) {
    ds.table.bytes(ds.header.series[key.index()], ds.context, n, out);
}

struct Reconstruction {
    bases: u8[>..]&, cigar: u32[>..]&, qualities: u8[>..]&, payload: u8[>..]&, md: u8[>..]&,
    reference: const u8[:], reference_start: i64, position: i64,
    no_reference: bool, generate: bool, quality_array: bool,
    read: i64 = 0, consumed: i64 = 0, distance: i64 = 0, nm: i64 = 0,
}

inline fn reference_at(r: const Reconstruction&, offset: i64) -> u8 {
    let at = r.position - r.reference_start + offset;
    if at < 0 || at >= r.reference.len { 'N' } else { r.reference[at] }
}

fn push_cigar(r: Reconstruction&, n: i64, op: u32, ctx: cram_encoding::Context&) {
    if n == 0 { return; }
    if n < 0 || n > 268435455 { ctx.ok = false; return; }
    if r.cigar.len > 0 && (r.cigar[r.cigar.len - 1] & 15) == op {
        let prior = (r.cigar[r.cigar.len - 1] >> 4) as i64;
        if n > 268435455 - prior { ctx.ok = false; return; }
        r.cigar[r.cigar.len - 1] = ((prior + n) as u32 << 4) | op;
    } else { r.cigar.push((n as u32 << 4) | op); }
}

fn flush_md(r: Reconstruction&) { if r.distance >= 0 { format(r.md, r.distance); } }

fn mismatch(r: Reconstruction&, b: u8) {
    format(r.md, r.distance); r.md.push(b); r.distance = 0; r.nm += 1;
}

fn matches(r: Reconstruction&, n: i64, ctx: cram_encoding::Context&) {
    if n <= 0 { return; }
    if n > ctx.limit - r.bases.len { ctx.ok = false; return; }
    let start = r.position - r.reference_start + r.consumed;
    if start >= 0 && n <= r.reference.len - start {
        let bases = r.reference[start..start + n];
        r.bases.append(bases);
        if r.generate {
            for b in bases { if b == 'N' { r.mismatch(b); } else { r.distance += 1; } }
        }
    } else {
        for i in n { r.bases.push(if r.no_reference { '=' } else { r.reference_at(r.consumed + i) }); }
        if r.generate && r.distance >= 0 { r.distance += n; }
    }
    r.push_cigar(n, 0, ctx); r.read += n; r.consumed += n;
}

fn quality(r: Reconstruction&, at: i64, data: const u8[:]) {
    if r.quality_array || r.qualities.len == 0 { return; }
    if r.qualities[0] == 255 { r.qualities.fill(30); }
    if at >= 0 && at < r.qualities.len {
        let n = min(data.len, r.qualities.len - at);
        r.qualities[at..at + n].move_from(data[..n]);
    }
}

fn features(r: Reconstruction&, ds: Series&, count: i64, length: i64) {
    if count < 0 || count > ds.context.limit || length < 0 || length > ds.context.limit {
        ds.context.ok = false; return;
    }
    // The caller clears scratch and initializes feature qualities before
    // constructing this borrowed view; no output buffer shrinks through it.
    if !r.quality_array && r.qualities.len != length { ds.context.ok = false; return; }
    var anchor: u32 = 0;
    for i in count {
        let code = ds.value(cram_encoding::Series.FC) as! u8;
        let delta = ds.value(cram_encoding::Series.FP) as! u32;
        anchor = (((anchor as u64) + delta as u64) & 4294967295) as u32;
        let target = max(0, anchor as i64 - 1);
        r.matches(max(0, min(target, length) - r.read), ds.context);
        let ref_index = r.position - r.reference_start + r.consumed;
        let have_reference = ref_index >= 0 && ref_index < r.reference.len;
        if code == 'X' {
            let substitution = ds.value(cram_encoding::Series.BS) & 3;
            let ref = r.reference_at(r.consumed);
            let row = if ref == 'A' { 0 } else if ref == 'C' { 1 } else if ref == 'G' { 2 } else if ref == 'T' { 3 } else { 4 };
            r.bases.push(ds.header.substitutions[row * 4 + substitution]);
            if r.generate {
                if have_reference { r.mismatch(ref); } else { r.flush_md(); r.distance = -1; }
            }
            r.push_cigar(1, 0, ds.context); r.read += 1; r.consumed += 1;
        } else if code == 'I' || code == 'S' || code == 'b' {
            let first = r.payload.len;
            ds.array(if code == 'I' { cram_encoding::Series.IN } else if code == 'S' { cram_encoding::Series.SC } else { cram_encoding::Series.BB }, r.payload);
            let n = r.payload.len - first;
            if n > ds.context.limit - r.bases.len { ds.context.ok = false; return; }
            r.bases.append(r.payload[first..]);
            if r.generate && code == 'I' { r.nm += n; }
            if code == 'b' {
                if r.generate {
                    r.flush_md();
                    for j in n {
                        if j > 0 { r.md.push('0'); }
                        if ref_index + j < 0 || ref_index + j >= r.reference.len { r.distance = -1; break; }
                        r.md.push(r.reference[ref_index + j]); r.nm += 1;
                    }
                    r.distance = 0;
                }
                r.consumed += n;
            }
            r.push_cigar(n, if code == 'I' { 1 } else if code == 'S' { 4 } else { 0 }, ds.context); r.read += n;
        } else if code == 'i' || code == 'B' {
            let base = ds.value(cram_encoding::Series.BA) as! u8; r.bases.push(base);
            if code == 'B' {
                let q = ds.value(cram_encoding::Series.QS) as! u8; var one: u8[1] = [q]; r.quality(target, one);
                if r.generate {
                    r.flush_md();
                    if have_reference { r.md.push(r.reference_at(r.consumed)); r.distance = 0; r.nm += 1; }
                    else { r.distance = -1; }
                }
                r.consumed += 1;
            } else if r.generate { r.nm += 1; }
            r.push_cigar(1, if code == 'i' { 1 } else { 0 }, ds.context); r.read += 1;
        } else if code == 'D' || code == 'N' || code == 'H' || code == 'P' {
            let key = if code == 'D' { cram_encoding::Series.DL } else if code == 'N' { cram_encoding::Series.RS } else if code == 'H' { cram_encoding::Series.HC } else { cram_encoding::Series.PD };
            let n = ds.value(key) as i64;
            let op: u32 = if code == 'D' { 2 } else if code == 'N' { 3 } else if code == 'H' { 5 } else { 6 };
            r.push_cigar(n, op, ds.context);
            if n < 0 { return; }
            if code == 'D' {
                if r.generate {
                    r.flush_md();
                    if ref_index >= 0 && ref_index <= r.reference.len {
                        let have = min(n, r.reference.len - ref_index);
                        r.md.push('^'); r.md.append(r.reference[ref_index..ref_index + have]); r.nm += have;
                        r.distance = 0;
                        if have < n { r.md.push('0'); r.distance = -1; }
                    } else { r.distance = -1; }
                }
                r.consumed += n;
            } else if code == 'N' { r.consumed += n; }
        } else if code == 'q' || code == 'Q' {
            let first = r.payload.len;
            if code == 'q' { ds.array(cram_encoding::Series.QQ, r.payload); } else { r.payload.push(ds.value(cram_encoding::Series.QS) as! u8); }
            r.quality(target, r.payload[first..]);
        } else { ds.context.ok = false; }
        if !ds.context.ok { return; }
    }
    r.matches(max(0, length - r.read), ds.context);
}

// Consume the same feature payloads as reconstruction, but retain only the
// reference span. Multi-reference slices use this before fetching windows.
fn scan_span(ds: Series&, p: const Pending&, payload: u8[>..]&) -> i64 {
    var read: i64 = 0; var span: i64 = 0; var anchor: u32 = 0;
    if p.flags & 4 != 0 { ds.bytes(cram_encoding::Series.BA, p.length, payload); }
    else {
        let count = ds.value(cram_encoding::Series.FN) as i64;
        if count < 0 || count > ds.context.limit { ds.context.ok = false; return 0; }
        for i in count {
            let code = ds.value(cram_encoding::Series.FC) as! u8;
            let delta = ds.value(cram_encoding::Series.FP) as! u32;
            anchor = ((anchor as u64 + delta as u64) & 4294967295) as u32;
            let fill = max(0, min(max(0, anchor as i64 - 1), p.length) - read);
            read += fill; span += fill;
            if code == 'X' { ds.value(cram_encoding::Series.BS); read += 1; span += 1; }
            else if code == 'i' || code == 'B' {
                ds.value(cram_encoding::Series.BA); read += 1;
                if code == 'B' { ds.value(cram_encoding::Series.QS); span += 1; }
            } else if code == 'I' || code == 'S' || code == 'b' || code == 'q' {
                let begin = payload.len;
                ds.array(if code == 'I' { cram_encoding::Series.IN } else if code == 'S' { cram_encoding::Series.SC }
                         else if code == 'b' { cram_encoding::Series.BB } else { cram_encoding::Series.QQ }, payload);
                let n = payload.len - begin;
                if code != 'q' { read += n; }
                if code == 'b' { span += n; }
            } else if code == 'Q' { ds.value(cram_encoding::Series.QS); }
            else if code == 'D' || code == 'N' || code == 'H' || code == 'P' {
                let n = ds.value(if code == 'D' { cram_encoding::Series.DL } else if code == 'N' { cram_encoding::Series.RS }
                                 else if code == 'H' { cram_encoding::Series.HC } else { cram_encoding::Series.PD });
                if n < 0 { ds.context.ok = false; }
                if code == 'D' || code == 'N' { span += n; }
            } else { ds.context.ok = false; }
            if !ds.context.ok { return 0; }
        }
        span += max(0, p.length - read);
        ds.value(cram_encoding::Series.MQ);
    }
    if p.cf & 1 != 0 { ds.bytes(cram_encoding::Series.QS, p.length, payload); }
    span
}

// Fixed metadata survives until the entire slice has been decoded. Mate links
// are slice-relative, not filtered-store indices; chains may have >2 records.
struct Mate {
    next: i64 = -1, end: i64 = 0, resolved: bool = false,
    explicit_tlen: bool = false, tlen: i32 = 0,
}

struct Scratch {
    name: u8[>..]&, bases: u8[>..]&, cigar: u32[>..]&, qualities: u8[>..]&,
    aux: u8[>..]&, payload: u8[>..]&, md: u8[>..]&,
}

struct Output {
    records: pileup::SlimRecord[>..]&, mates: Mate[>..]&,
    names: u8[>..]&, bases: u8[>..]&, cigars: u32[>..]&,
    qualities: u8[>..]&, auxiliary: u8[>..]&, limit: i64,
}

// Read the fields which precede reference-dependent reconstruction. Keeping
// this separate lets the reader fetch the appropriate reference in multi-ref
// slices without reading any series out of order.
struct Pending {
    tid: i32, pos: i64, flags: u16, cf: i32, length: i64, group: i32,
    mate_tid: i32 = -1, mate_pos: i64 = -1, tlen: i32 = 0,
    mate: Mate, stored_md: bool = false, stored_nm: bool = false,
}

fn begin_record(ds: Series&, slice: const Header&, previous: i64&, index: i64,
                scratch: Scratch&, uses: const cram_encoding::Tag[:],
                lines: const cram_encoding::TagLine[:], references: i64) -> Pending {
    let raw = ds.value(cram_encoding::Series.BF);
    let cf = ds.value(cram_encoding::Series.CF);
    let tid = if slice.tid == -2 { ds.value(cram_encoding::Series.RI) } else { slice.tid };
    let length = ds.value(cram_encoding::Series.RL) as i64;
    let ap = ds.value(cram_encoding::Series.AP) as i64;
    previous = if ds.header.delta { previous + ap } else { ap };
    let group = ds.value(cram_encoding::Series.RG);
    let p = Pending { tid: tid, pos: max(-1, previous - 1), flags: (raw & 65535) as u16,
                      cf: cf, length: length, group: group, .. };
    if raw < 0 || raw > 65535 || tid < -1 || tid >= references || length < 0 ||
       length > ds.context.limit || (raw & 4 == 0 && previous < 1) {
        ds.context.ok = false; return p;
    }
    if ds.header.names { ds.array(cram_encoding::Series.RN, scratch.name); }
    if cf & 2 != 0 {
        let flags = ds.value(cram_encoding::Series.MF);
        if !ds.header.names { ds.array(cram_encoding::Series.RN, scratch.name); }
        p.mate_tid = ds.value(cram_encoding::Series.NS);
        let np = ds.value(cram_encoding::Series.NP);
        p.mate_pos = np as i64 - 1;
        p.tlen = ds.value(cram_encoding::Series.TS);
        if flags & 1 != 0 { p.flags |= 33; }
        if flags & 2 != 0 { p.flags |= 8; }
        if p.flags & 1 == 0 { p.mate_tid = -1; }
        if np < 0 || p.mate_tid < -1 || p.mate_tid >= references { ds.context.ok = false; }
    } else {
        if cf & 4 != 0 { p.mate.next = index + 1 + ds.value(cram_encoding::Series.NF) as i64; }
        if cf & 16 != 0 { p.mate.explicit_tlen = true; p.mate.tlen = ds.value(cram_encoding::Series.TS); }
    }
    let line = ds.value(cram_encoding::Series.TL) as i64;
    if line < 0 { ds.context.ok = false; return p; }
    if line < lines.len {
        let range = lines[line];
        for tag in uses[range.begin..range.begin + range.count] {
            let first = scratch.payload.len;
            if tag.encoding < 0 { ds.context.ok = false; return p; }
            ds.table.array(tag.encoding, ds.context, scratch.payload);
            if tag.head >> 8 == 0x4d44 { p.stored_md = true; }
            if tag.head >> 8 == 0x4e4d { p.stored_nm = true; }
            if tag.head == 0x634643 && scratch.payload.len == first + 1 {
                let flags = scratch.payload[first];
                p.stored_md = p.stored_md || flags & 1 != 0;
                p.stored_nm = p.stored_nm || flags & 2 != 0;
            } else {
                scratch.aux.push(((tag.head >> 16) & 255) as u8);
                scratch.aux.push(((tag.head >> 8) & 255) as u8);
                scratch.aux.push((tag.head & 255) as u8);
                scratch.aux.append(scratch.payload[first..]);
            }
        }
    }
    p
}

fn finish_record(ds: Series&, p: Pending&, scratch: Scratch&, output: Output&,
                 reference: const u8[:], reference_start: i64, no_reference: bool,
                 groups: const bam::Reference[:], group_names: const u8[:]) {
    let quality_array = p.cf & 1 != 0;
    if !quality_array { for i in p.length { scratch.qualities.push(255); } }
    var span: i64 = 0; var mapq: u8 = 0;
    if p.flags & 4 == 0 {
        let generate = !no_reference && reference.len > 0 && p.cf & 8 == 0;
        let count = ds.value(cram_encoding::Series.FN) as i64;
        let r = Reconstruction { bases: scratch.bases, cigar: scratch.cigar, qualities: scratch.qualities,
            payload: scratch.payload, md: scratch.md, reference: reference, reference_start: reference_start,
            position: p.pos, no_reference: no_reference, generate: generate && (!p.stored_md || !p.stored_nm),
            quality_array: quality_array, .. };
        r.features(ds, count, p.length);
        span = r.consumed;
        if generate {
            r.flush_md();
            if !p.stored_md { scratch.aux.append("MDZ"); scratch.aux.append(scratch.md); scratch.aux.push(0); }
            if !p.stored_nm {
                if r.nm < 0 || r.nm > 4294967295 { ds.context.ok = false; return; }
                scratch.aux.append("NM");
                let n = if r.nm <= 255 { 1 } else if r.nm <= 65535 { 2 } else { 4 };
                scratch.aux.push(if n == 1 { 'C' } else if n == 2 { 'S' } else { 'I' });
                for i in n { scratch.aux.push(((r.nm >> (i * 8)) & 255) as u8); }
            }
        }
        mapq = ds.value(cram_encoding::Series.MQ) as! u8;
    } else { ds.bytes(cram_encoding::Series.BA, p.length, scratch.bases); }
    if quality_array { ds.bytes(cram_encoding::Series.QS, p.length, scratch.qualities); }
    let sequence_length = if p.cf & 8 != 0 { 0 } else { scratch.bases.len };
    if p.group >= 0 && p.group < groups.len {
        let group = groups[p.group as i64];
        let name = group_names[group.name..group.name + group.name_len];
        if p.group != groups.len - 1 || name != "UNKNOWN" {
            scratch.aux.append("RGZ"); scratch.aux.append(name); scratch.aux.push(0);
        }
    } else if p.group != -1 { ds.context.ok = false; }
    if !ds.context.ok { return; }
    if scratch.bases.len != scratch.qualities.len { ds.context.ok = false; return; }
    let used = output.names.len + output.bases.len + output.qualities.len + output.cigars.len * 4 + output.auxiliary.len;
    let added = scratch.name.len + sequence_length * 2 + scratch.cigar.len * 4 + scratch.aux.len;
    if added > min(output.limit, 4294967295) - used { ds.context.ok = false; return; }
    let item = pileup::SlimRecord {
        tid: p.tid, pos: p.pos, end: p.pos + max(span, 1), reference_len: span, flags: p.flags, mapq: mapq,
        mate_tid: p.mate_tid, mate_pos: p.mate_pos, template_len: p.tlen,
        name: output.names.len as u32, name_len: scratch.name.len as u32,
        cigar: output.cigars.len as u32, cigar_len: scratch.cigar.len as u32,
        seq: output.bases.len as u32, seq_len: sequence_length as u32,
        qual: output.qualities.len as u32, aux: output.auxiliary.len as u32, aux_len: scratch.aux.len as u32, .. };
    output.names.append(scratch.name); output.cigars.append(scratch.cigar);
    output.qualities.append(scratch.qualities[..sequence_length]); output.auxiliary.append(scratch.aux);
    let first_base = output.bases.len;
    output.bases.append(scratch.bases[..sequence_length]);
    for &base in output.bases[first_base..] {
        let upper = if base >= 'a' && base <= 'z' { base - 32 } else { base };
        let letters = "\x01\x0e\x02\x0d\x0f\x0f\x04\x0b\x0f\x0f\x0c\x0f\x03\x0f\x0f\x0f\x0f\x05\x06\x08\x0f\x07\x09\x0f\x0a\x0f";
        base = if upper >= 'A' && upper <= 'Z' { letters[(upper - 'A') as i64] } else if upper == '=' { 0 } else { 15 };
    }
    p.mate.end = item.end;
    output.records.push(item); output.mates.push(p.mate);
}

fn close_chain(records: pileup::SlimRecord[:], mates: Mate[:], first: i64) -> bool {
    var left = records[first].pos; var right = mates[first].end;
    var left_count = 0; var right_count = 0; var same_reference = true; var id = first;
    loop {
        let pos = records[id].pos; let end = mates[id].end;
        if left > pos { left = pos; left_count = 1; } else if left == pos { left_count += 1; }
        if right < end { right = end; right_count = 1; } else if right == end { right_count += 1; }
        let next = mates[id].next;
        if next < 0 { mates[id].next = first; break; }
        if next <= id || next >= mates.len { return false; }
        id = next; same_reference = same_reference && records[id].tid == records[first].tid;
    }
    let span = right - left;
    let f = records[first];
    let signed = if !same_reference { 0 }
        else if f.pos == left && (mates[first].end < right || left_count <= 1) { span }
        else if f.pos == left && mates[first].end == right && left_count > 1 && right_count > 1 {
            if f.flags & 64 != 0 { span } else { -span }
        } else { -span };
    id = first;
    loop {
        let tlen = if id == first { signed } else { -signed };
        records[id].template_len = if tlen < -2147483648 || tlen > 2147483647 { 0 } else { tlen as i32 };
        mates[id].resolved = true;
        id = mates[id].next;
        if id == first { break; }
    }
    true
}

fn resolve(output: Output&, prefix: const u8[:], counter: i64, names: bool) -> bool {
    for i in output.records.len {
        let next = output.mates[i].next;
        if next < 0 || next >= output.records.len { continue; }
        if !output.mates[i].resolved && !close_chain(output.records, output.mates, i) { continue; }
        let mate = output.records[output.mates[i].next];
        let own = output.records[i];
        own.mate_tid = mate.tid; own.mate_pos = max(-1, mate.pos); own.flags |= 1;
        if mate.flags & 4 != 0 { own.flags |= 8; }
        if mate.flags & 16 != 0 { own.flags |= 32; }
        if own.flags & 4 != 0 || mate.flags & 4 != 0 { own.template_len = 0; }
        output.records[i] = own;
    }
    for i in output.records.len {
        if output.mates[i].explicit_tlen { output.records[i].template_len = output.mates[i].tlen; }
        if !names || output.records[i].name_len != 0 { continue; }
        let next = output.mates[i].next;
        let have_mate = next >= 0 && next < output.records.len;
        let start = output.names.len;
        if have_mate && output.records[next].name_len > 0 {
            let mate = output.records[next];
            for j in mate.name_len as i64 { output.names.push(output.names[mate.name as i64 + j]); }
        } else {
            output.names.append(prefix[..min(prefix.len, 233)]);
            format(output.names, ":", counter + (if have_mate && next < i { next } else { i }) + 1);
        }
        if output.names.len > min(output.limit, 4294967295) { return false; }
        output.records[i].name = start as u32;
        output.records[i].name_len = (output.names.len - start) as u32;
    }
    true
}
