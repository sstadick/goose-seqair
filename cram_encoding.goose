// CRAM v3 data-series decoding, following Seqair 8401961 encoding.rs.
// Each content ID has exactly one cursor, even when several series share it.
namespace cram_encoding;
import std;
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

struct Cursor { bytes: const u8[:], at: i64 = 0, ok: bool = true }

inline fn byte(c: Cursor&) -> u8 {
    if c.at >= c.bytes.len { c.ok = false; return 0; }
    let v = c.bytes[c.at]; c.at += 1; v
}

fn take(c: Cursor&, n: i64) -> const u8[:] {
    if n < 0 || n > c.bytes.len - c.at { c.ok = false; return c.bytes[..0]; }
    let begin = c.at; c.at += n; c.bytes[begin..c.at]
}

inline fn itf8(c: Cursor&) -> i32 {
    let b = c.byte();
    if b < 128 { return b as i32; }
    var v = b as u32;
    if b < 192 { v = ((v & 63) << 8) | c.byte() as u32; }
    else if b < 224 { v = ((v & 31) << 16) | ((c.byte() as u32) << 8) | c.byte() as u32; }
    else if b < 240 {
        v = ((v & 15) << 24) | ((c.byte() as u32) << 16) | ((c.byte() as u32) << 8) | c.byte() as u32;
    } else {
        v = ((v & 15) << 28) | ((c.byte() as u32) << 20) | ((c.byte() as u32) << 12) |
            ((c.byte() as u32) << 4) | ((c.byte() as u32) & 15);
    }
    v as! i32
}

fn ltf8(c: Cursor&) -> i64 {
    let b = c.byte();
    var n = 0;
    var mask: u8 = 128;
    while n < 8 && (b & mask) != 0 { n += 1; mask >>= 1; }
    var v: u64 = if n >= 7 { 0 } else { (b & (mask - 1)) as u64 };
    for i in n { v = (v << 8) | c.byte() as u64; }
    v as! i64
}

struct Block { id: i32, kind: u8, begin: i64, size: i64, at: i64 = 0 }
struct Context {
    blocks: Block[:], data: const u8[:], core: i64 = -1, bit: i64 = 0,
    ok: bool = true, limit: i64 = 16777216,
}

fn payload(ctx: const Context&, index: i64) -> const u8[:] {
    let b = ctx.blocks[index]; ctx.data[b.begin..b.begin + b.size]
}

fn bits(ctx: Context&, n: i64) -> u32 {
    if n < 0 || n > 32 || ctx.core < 0 || n > ctx.blocks[ctx.core].size * 8 - ctx.bit {
        ctx.ok = false; return 0;
    }
    var out: u32 = 0;
    var left = n;
    while left > 0 {
        let count = min(left, 8 - (ctx.bit & 7));
        let value = ctx.data[ctx.blocks[ctx.core].begin + (ctx.bit >> 3)] as u32;
        out = (out << count) | ((value >> (8 - (ctx.bit & 7) - count)) & (((1 as u32) << count) - 1));
        ctx.bit += count; left -= count;
    }
    out
}

struct Encoding {
    kind: i32 = 0, offset: i32 = 0, width: i64 = 0, id: i32 = 0,
    external_slot: i64 = -1,
    length: i64 = -1, value: i64 = -1, stop: u8 = 0,
    first: i64 = 0, count: i64 = 0, levels: i64[34] = [0; 34],
    max_bits: i64 = 0,
}

struct Symbol { value: i32, width: i64 = 0, code: u32 = 0 }
struct Table { encodings: Encoding[>..]&, symbols: Symbol[>..]& }

// Mode: 0 integer, 1 byte, 2 byte array. Descriptor trees are flattened into
// indices so container scratch can be reused without recursive allocations.
fn parse(table: Table&, c: Cursor&, mode: i32) -> i64 {
    if !c.ok { return -1; }
    let kind = c.itf8(); let size = c.itf8() as i64;
    let p = Cursor { bytes: c.take(size) };
    var enc: Encoding = Encoding { kind: kind, .. };
    if kind == 4 && mode == 2 {
        enc.length = table.parse_scalar(p, 0);
        enc.value = table.parse_scalar(p, 1);
    } else { enc = table.parameters(p, kind, mode); }
    c.ok = c.ok && p.ok;
    let index = table.encodings.len; table.encodings.push(enc); index
}

fn parse_scalar(table: Table&, c: Cursor&, mode: i32) -> i64 {
    let kind = c.itf8(); let size = c.itf8() as i64;
    let p = Cursor { bytes: c.take(size) };
    let enc = table.parameters(p, kind, mode);
    c.ok = c.ok && p.ok;
    let index = table.encodings.len; table.encodings.push(enc); index
}

fn parameters(table: Table&, p: Cursor&, kind: i32, mode: i32) -> Encoding {
    let enc = Encoding { kind: kind, .. };
    if kind == 0 { }
    else if kind == 1 { enc.id = p.itf8(); }
    else if kind == 3 && mode != 2 {
        let n = p.itf8() as i64;
        if n < 0 || n > p.bytes.len - p.at { p.ok = false; }
        else {
            enc.first = table.symbols.len; enc.count = n;
            for i in n { table.symbols.push(Symbol { value: p.itf8(), .. }); }
            let lengths = p.itf8() as i64;
            if lengths != n { p.ok = false; }
            let symbols = table.symbols[enc.first..enc.first + n];
            if p.ok {
                for i in n { let width = p.itf8() as i64; symbols[i].width = width; if width < 0 || width > 32 { p.ok = false; } }
            }
            if p.ok && n > 1 {
                // Alphabet sizes are small; sort only once per compression header.
                for i in n {
                    var j = i;
                    while j > 0 && (symbols[j].width < symbols[j - 1].width ||
                          (symbols[j].width == symbols[j - 1].width && symbols[j].value < symbols[j - 1].value)) {
                        var s: Symbol = symbols[j]; symbols[j] = symbols[j - 1]; symbols[j - 1] = s;
                        j -= 1;
                    }
                }
                var code: u64 = 0; var previous = symbols[0].width;
                for &symbol in symbols {
                    let width = symbol.width;
                    code <<= width - previous;
                    if width == 0 || code >= ((1 as u64) << width) { p.ok = false; break; }
                    symbol.code = code as u32; code += 1; previous = width;
                }
                enc.max_bits = symbols[n - 1].width;
                var index = 0;
                for level in 34 {
                    while index < n && symbols[index].width < level { index += 1; }
                    enc.levels[level] = index;
                }
            }
        }
    } else if kind == 5 && mode == 2 {
        enc.stop = p.byte(); enc.id = p.itf8();
    } else if kind == 6 && mode != 2 {
        enc.offset = p.itf8(); enc.width = p.itf8() as i64;
        if enc.width < 0 || enc.width > 32 { p.ok = false; }
    } else if kind == 7 && mode == 0 {
        enc.offset = p.itf8(); enc.width = (p.itf8() as! u32) as i64;
    } else if kind == 9 && mode == 0 { enc.offset = p.itf8(); }
    else { p.ok = false; }
    enc
}

fn huffman(table: const Table&, enc: const Encoding&, ctx: Context&) -> i32 {
    let symbols = table.symbols[enc.first..enc.first + enc.count];
    if enc.count == 1 { return symbols[0].value; }
    var code: u32 = 0;
    for i in enc.max_bits {
        let level = i + 1;
        code = (code << 1) | ctx.bits(1);
        let start = enc.levels[level]; let end = enc.levels[level + 1];
        if start < end && code >= symbols[start].code {
            let offset = (code - symbols[start].code) as i64;
            if offset < end - start { return symbols[start + offset].value; }
        }
    }
    ctx.ok = false; 0
}

fn integer(table: const Table&, index: i64, ctx: Context&, byte_value: bool = false) -> i32 {
    if !ctx.ok { return 0; }
    if index < 0 { return 0; } // Missing data series have Null encoding.
    let enc: const Encoding& = table.encodings[index];
    if enc.kind == 0 { return 0; }
    if enc.kind == 1 {
        let b = enc.external_slot;
        if b < 0 { ctx.ok = false; return 0; }
        let c = Cursor { bytes: ctx.payload(b), at: ctx.blocks[b].at };
        let v = if byte_value { c.byte() as i32 } else { c.itf8() };
        ctx.blocks[b].at = c.at; ctx.ok = c.ok; return v;
    }
    if enc.kind == 3 { return table.huffman(enc, ctx); }
    var raw: i64 = 0;
    if enc.kind == 6 { raw = (ctx.bits(enc.width) as! i32) as i64; }
    else if enc.kind == 7 {
        var unary = 0;
        while ctx.ok && ctx.bits(1) == 1 { unary += 1; if unary > 32 { ctx.ok = false; } }
        let width = if unary == 0 { enc.width } else { unary + enc.width - 1 };
        if !ctx.ok || width > 32 || (unary != 0 && width == 32) { ctx.ok = false; return 0; }
        raw = ctx.bits(width) as i64;
        if unary != 0 { raw += ((1 as i64) << width) - ((1 as i64) << enc.width); }
        raw = (raw as u32 as! i32) as i64;
    } else if enc.kind == 9 {
        var zeros = 0;
        while ctx.ok && ctx.bits(1) == 0 { zeros += 1; if zeros > 31 { ctx.ok = false; } }
        if !ctx.ok { return 0; }
        raw = ((((1 as u32) << zeros) | ctx.bits(zeros)) as! i32) as i64 - 1;
        if raw < -2147483648 { ctx.ok = false; return 0; }
    } else { ctx.ok = false; }
    let result = raw - enc.offset as i64;
    if !byte_value && (result < -2147483648 || result > 2147483647) { ctx.ok = false; return 0; }
    if byte_value { (result & 255) as i32 } else { result as i32 }
}

fn bytes(table: const Table&, index: i64, ctx: Context&, n: i64, out: u8[>..]&) {
    if !ctx.ok { return; }
    if n < 0 || n > ctx.limit - out.len { ctx.ok = false; return; }
    if n == 0 { return; }
    if index < 0 || table.encodings[index].kind == 0 { for i in n { out.push(0); } return; }
    let enc: const Encoding& = table.encodings[index];
    if enc.kind == 1 {
        let b = enc.external_slot; if b < 0 { ctx.ok = false; return; }
        let c = Cursor { bytes: ctx.payload(b), at: ctx.blocks[b].at };
        let data = c.take(n);
        if c.ok { out.append(data); }
        ctx.blocks[b].at = c.at; ctx.ok = c.ok;
    } else {
        for i in n { let value = table.integer(index, ctx, true); if !ctx.ok { break; } out.push(value as! u8); }
    }
}

fn array(table: const Table&, index: i64, ctx: Context&, out: u8[>..]&) {
    if !ctx.ok || index < 0 { return; }
    let enc: const Encoding& = table.encodings[index];
    if enc.kind == 0 { return; }
    if enc.kind == 4 {
        let n = table.integer(enc.length, ctx) as i64;
        table.bytes(enc.value, ctx, n, out);
    } else if enc.kind == 5 {
        let b = enc.external_slot; if b < 0 { ctx.ok = false; return; }
        let at = ctx.blocks[b].at;
        let data = ctx.payload(b)[at..];
        let n = bstr::find_byte(data, enc.stop);
        if n < 0 || n > ctx.limit - out.len { ctx.ok = false; return; }
        out.append(data[..n]); ctx.blocks[b].at += n + 1;
    } else { ctx.ok = false; }
}

enum Series { BF, CF, RI, RL, AP, RG, RN, MF, NS, NP, TS, NF, TL, FN, FC, FP, DL, MQ, BA, QS, BS, IN, SC, RS, PD, HC, BB, QQ }
inline fn index(key: Series) -> i64 {
    match key {
        BF => 0, CF => 1, RI => 2, RL => 3, AP => 4, RG => 5, RN => 6,
        MF => 7, NS => 8, NP => 9, TS => 10, NF => 11, TL => 12, FN => 13,
        FC => 14, FP => 15, DL => 16, MQ => 17, BA => 18, QS => 19, BS => 20,
        IN => 21, SC => 22, RS => 23, PD => 24, HC => 25, BB => 26, QQ => 27,
    }
}
const SERIES = "BFCFRIRLAPRGRNMFNSNPTSNFTLFNFCFPDLMQBAQSBSINSCRSPDHCBBQQ";
struct Tag { head: u32, encoding: i64 }
struct TagLine { begin: i64, count: i64 }
struct Compression {
    names: bool = true, delta: bool = true, reference: bool = true,
    substitutions: u8[20] = ['N'; 20], series: i64[28] = [-1; 28],
}

fn compression(data: const u8[:], table: Table&, tags: Tag[>..]&, uses: Tag[>..]&, lines: TagLine[>..]&) -> Compression, bool {
    let c = Cursor { bytes: data };
    let header = Compression { .. };
    let size = c.itf8() as i64;
    let p = Cursor { bytes: c.take(size) };
    let n = p.itf8() as i64;
    if n < 0 || n > p.bytes.len { return header, false; }
    for i in n {
        let key = p.take(2);
        if key == "RN" { header.names = p.byte() != 0; }
        else if key == "AP" { header.delta = p.byte() != 0; }
        else if key == "RR" { header.reference = p.byte() != 0; }
        else if key == "SM" {
            let alternatives = "CGTNAGTNACTNACGNACGT";
            for row in 5 {
                let permutation = p.byte();
                for column in 4 {
                    let slot = ((permutation >> (6 - column * 2)) & 3) as i64;
                    header.substitutions[row * 4 + slot] = alternatives[row * 4 + column];
                }
            }
        } else if key == "TD" {
            let len = p.itf8() as i64;
            let dictionary = Cursor { bytes: p.take(len) };
            var first = uses.len;
            while dictionary.at < dictionary.bytes.len && dictionary.ok {
                let a = dictionary.byte();
                if a == 0 { lines.push(TagLine { begin: first, count: uses.len - first }); first = uses.len; }
                else {
                    let b = dictionary.byte(); let kind = dictionary.byte();
                    uses.push(Tag { head: ((a as u32) << 16) | ((b as u32) << 8) | kind, encoding: -1 });
                }
            }
            if uses.len > first { lines.push(TagLine { begin: first, count: uses.len - first }); }
            p.ok = p.ok && dictionary.ok;
        } else { p.byte(); }
        if !p.ok { break; }
    }
    if !c.ok || !p.ok { return header, false; }
    let series_size = c.itf8() as i64;
    let ds = Cursor { bytes: c.take(series_size) };
    let count = ds.itf8() as i64;
    if count < 0 || count > ds.bytes.len { return header, false; }
    for i in count {
        let key = ds.take(2);
        let mode: i32 = if key == "RN" || key == "IN" || key == "SC" || key == "BB" || key == "QQ" { 2 }
                   else if key == "FC" || key == "BA" || key == "QS" || key == "BS" { 1 } else { 0 };
        let index = table.parse(ds, mode);
        for slot in 28 { if key == SERIES[slot * 2..slot * 2 + 2] { header.series[slot] = index; break; } }
        if !ds.ok { break; }
    }
    if !c.ok || !ds.ok { return header, false; }
    let tag_size = c.itf8() as i64;
    let ts = Cursor { bytes: c.take(tag_size) };
    let tag_count = ts.itf8() as i64;
    if tag_count < 0 || tag_count > ts.bytes.len { return header, false; }
    for i in tag_count {
        let head = ts.itf8() as! u32;
        let index = table.parse(ts, 2);
        tags.push(Tag { head: head, encoding: index });
        if !ts.ok { break; }
    }
    for &used in uses { for tag in tags { if tag.head == used.head { used.encoding = tag.encoding; } } }
    return header, c.ok && ts.ok;
}
