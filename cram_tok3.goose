// CRAM name tokenizer, following Seqair's stream descriptors and flat token/
// name arenas. Copied token streams share bytes but retain independent cursors.
namespace cram_tok3;
import std;
import cram_encoding;
import cram_nx16;
import cram_arith;

struct Stream { begin: i64 = 0, size: i64 = 0, at: i64 = 0, generated: i32 = -1 }
struct Token { value: u32, start: u32, length_kind: u32 }
struct Name { start: u32, length: u32, first: u32, count: u32 }

fn byte(streams: Stream[:], index: i64, data: const u8[:], ok: bool&) -> u8 {
    if !ok || index < 0 || index >= streams.len { ok = false; return 0; }
    let s = streams[index];
    if s.at >= s.size { ok = false; return 0; }
    streams[index].at += 1;
    if s.generated >= 0 { if s.at == 0 { s.generated as u8 } else { 10 } }
    else { data[s.begin + s.at] }
}

fn word(streams: Stream[:], index: i64, data: const u8[:], ok: bool&) -> u32 {
    var n: u32 = 0;
    for i in 4 { n |= (byte(streams, index, data, ok) as u32) << (8 * i); }
    n
}

fn copy(output: u8[:], at: i64&, begin: i64, length: i64) -> bool {
    if begin < 0 || length > at - begin || length > output.len - at { return false; }
    for i in length { output[at + i] = output[begin + i]; }
    at += length; true
}

fn digits(output: u8[:], at: i64&, value: u32, width: i64) -> bool {
    var bytes: u8[10] = [0; 10]; var n = value; var used = 0;
    loop {
        bytes[used] = '0' + (n % 10) as u8; used += 1; n /= 10;
        if n == 0 { break; }
    }
    let length = max(width, used);
    if length > output.len - at { return false; }
    output[at..at + length - used].fill('0');
    for i in used { output[at + length - used + i] = bytes[used - 1 - i]; }
    at += length; true
}

fn decode(input: const u8[:], output: u8[:], frequencies: u32[>..]&, cumulative: u32[>..]&, symbols: u8[>..]&, limit: i64) -> bool {
    let c = cram_encoding::Cursor { bytes: input };
    let length = cram_nx16::word(c) as i64; let count = cram_nx16::word(c) as i64;
    let arithmetic = c.byte() != 0;
    if !c.ok || length > limit || output.len > limit || output.len > length + 1024 || count > 10000000 || count > length + 1024 { return false; }
    var streams: Stream[1664] = [Stream { .. }; 1664];
    var positions = 0; var data: u8[>..] = [];
    while c.at < c.bytes.len {
        let descriptor = c.byte(); let kind = (descriptor & 63) as i64;
        if kind >= 13 { return false; }
        if descriptor & 128 != 0 {
            if positions == 128 { return false; }
            if kind != 0 { streams[positions * 13] = Stream { size: count, generated: kind as i32, .. }; }
            positions += 1;
        }
        if positions == 0 { return false; }
        let slot = (positions - 1) * 13 + kind;
        if descriptor & 64 != 0 {
            let position = c.byte() as i64; let ty = c.byte() as i64;
            if !c.ok || position >= positions || ty >= 13 { return false; }
            let source = streams[position * 13 + ty];
            // A decoded empty stream is set too; track it with a negative
            // generated marker distinct from the unset default.
            if source.generated == -1 { return false; }
            streams[slot] = source;
        } else {
            let compressed = cram_nx16::uint7(c) as i64;
            let bytes = c.take(compressed);
            let header = cram_encoding::Cursor { bytes: bytes };
            let flags = header.byte();
            let size = if flags & 16 == 0 || (arithmetic && flags & 8 != 0) { cram_nx16::uint7(header) as i64 } else { 0 };
            if !c.ok || !header.ok || size > limit - data.len { return false; }
            let begin = data.len; data.resize(begin + size, 0);
            let ok = if arithmetic { cram_arith::decode(bytes, data[begin..], limit) }
                else { cram_nx16::decode(bytes, data[begin..], frequencies, cumulative, symbols, limit) };
            if !ok { return false; }
            streams[slot] = Stream { begin: begin, size: size, generated: -2, .. };
        }
    }
    var names: Name[>..] = []; var tokens: Token[>..] = [];
    var at = 0; var ok = true;
    for i in count {
        if positions == 0 { return false; }
        let start = at;
        let kind = byte(streams, 0, data, ok) as i64;
        if kind >= 13 { return false; }
        let distance = word(streams, kind, data, ok) as i64;
        if !ok || distance > i { return false; }
        let prior = i - distance;
        let previous = if prior < names.len { names[prior] } else { Name { .. } };
        if kind == 5 {
            if prior >= names.len || !copy(output, at, previous.start as i64, previous.length as i64) || at >= output.len { return false; }
            output[at] = 0; at += 1;
            names.push(Name { start: start as u32, length: previous.length, first: previous.first, count: previous.count });
            continue;
        }
        let first = tokens.len;
        var position = 1;
        loop {
            if position >= positions { return false; }
            let base = position * 13;
            let ty = byte(streams, base, data, ok);
            let begin = at;
            let prev = if position <= previous.count { tokens[previous.first as i64 + position - 1] } else { Token { .. } };
            var value: u32 = 0; var width: i64 = 0; var token_kind: u32 = 0;
            if ty == 2 {
                let b = byte(streams, base + 2, data, ok);
                if at >= output.len { return false; }
                output[at] = b; at += 1; token_kind = 1;
            } else if ty == 1 {
                loop {
                    let b = byte(streams, base + 1, data, ok);
                    if !ok { return false; }
                    if b == 0 { break; }
                    if at >= output.len { return false; }
                    output[at] = b; at += 1;
                }
                token_kind = 1;
            } else if ty == 3 || ty == 7 {
                value = word(streams, base + ty, data, ok);
                if ty == 3 { width = byte(streams, base + 4, data, ok) as i64; }
                token_kind = if ty == 3 { 3 } else { 2 };
            } else if ty == 8 || ty == 9 {
                let delta = byte(streams, base + ty, data, ok);
                token_kind = if ty == 9 { 3 } else { 2 };
                if prev.length_kind >> 30 != token_kind { return false; }
                value = ((prev.value as u64 + delta as u64) & 4294967295) as u32;
                if ty == 9 { width = (prev.length_kind & 1073741823) as i64; }
            } else if ty == 10 {
                if prev.length_kind >> 30 == 0 || !copy(output, at, prev.start as i64, (prev.length_kind & 1073741823) as i64) { return false; }
            } else if ty != 11 && ty != 12 { return false; }
            if token_kind >= 2 && !digits(output, at, value, width) { return false; }
            if !ok || tokens.len >= limit / 12 { return false; }
            tokens.push(if ty == 10 { prev } else {
                Token { value: value, start: begin as u32, length_kind: ((at - begin) as u32) | (token_kind << 30) }
            });
            if ty == 12 { break; }
            position += 1;
        }
        if at >= output.len || names.len >= limit / 16 { return false; }
        names.push(Name { start: start as u32, length: (at - start) as u32, first: first as u32, count: (tokens.len - first) as u32 });
        output[at] = 0; at += 1;
    }
    ok && at == output.len
}
