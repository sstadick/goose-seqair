// CRAM FQZComp version 5. Context selection, table coding, duplicate ordering
// and lazy quality models follow Seqair's fqzcomp.rs, including its first-
// parameter-block context/QMAP rule for multi-parameter streams.
namespace cram_fqz;
import std;
import cram_encoding;
import cram_nx16;
import cram_range;

struct Parameter {
    context: i64, flags: u8, max_symbol: i64, qshift: i64, qloc: i64, sloc: i64, qmask: u32,
    qmap: u8[256], qtab: u16[256], ptab: u16[1024], dtab: u16[256],
}

fn array(c: cram_encoding::Cursor&, output: u16[:]) {
    var runs: u8[1024] = [0; 1024];
    var count = 0; var total = 0; var last = -1;
    while total < output.len && c.at < c.bytes.len {
        let run = c.byte() as i64;
        runs[count] = run as u8; count += 1; total += run;
        if last == run {
            var copies = c.byte() as i64;
            total += run * copies;
            while copies > 0 && total <= output.len && count < 1024 {
                runs[count] = run as u8; count += 1; copies -= 1;
            }
        }
        if count >= 1024 { c.ok = false; return; }
        last = run;
    }
    var next = 0; var filled = 0; var value = 0;
    while filled < output.len {
        var length = 0;
        loop {
            if next >= count { c.ok = false; return; }
            let part = runs[next] as i64; next += 1; length += part;
            if part != 255 { break; }
        }
        let end = min(filled + length, output.len);
        output[filled..end].fill(value as u16); filled = end; value += 1;
    }
}

fn parameter(c: cram_encoding::Cursor&) -> Parameter {
    let p = Parameter { .. };
    if c.bytes.len - c.at < 7 { c.ok = false; return p; }
    let low = c.byte() as i64; p.context = low | (c.byte() as i64 << 8);
    p.flags = c.byte(); p.max_symbol = c.byte() as i64;
    let quality = c.byte(); let location = c.byte(); let position = c.byte();
    let qbits = (quality >> 4) as i64;
    p.qshift = (quality & 15) as i64; p.qloc = (location >> 4) as i64; p.sloc = (location & 15) as i64;
    p.qmask = ((1 as u32) << qbits) - 1;
    for i in 256 { p.qmap[i] = i as u8; p.qtab[i] = i as u16; }
    if p.flags & 16 != 0 {
        p.qmap.fill(255);
        for i in p.max_symbol { p.qmap[i] = c.byte(); }
    }
    if qbits > 0 && p.flags & 128 != 0 { array(c, p.qtab); }
    if p.flags & 32 != 0 {
        array(c, p.ptab);
        for &n in p.ptab { n = (((n as u32) << (position >> 4)) & 65535) as u16; }
    }
    if p.flags & 64 != 0 {
        array(c, p.dtab);
        for &n in p.dtab { n = (((n as u32) << (position & 15)) & 65535) as u16; }
    }
    p
}

struct Reverse { start: i64, length: i64 }

fn decode(input: const u8[:], output: u8[:], limit: i64) -> bool {
    let c = cram_encoding::Cursor { bytes: input };
    let length = cram_nx16::uint7(c) as i64;
    if !c.ok || length != output.len || length > limit || c.bytes.len - c.at < 10 { return false; }
    if c.byte() != 5 { return false; }
    let flags = c.byte();
    let count = if flags & 1 != 0 { c.byte() as i64 } else { 1 };
    if count == 0 { return false; }
    var max_selector = if count > 1 { count } else { 0 };
    var selectors: u16[256] = [0; 256];
    if flags & 2 != 0 { max_selector = c.byte() as i64; array(c, selectors); }
    var parameters: Parameter[>..] = []; var max_symbol = 0;
    for i in count {
        let p = parameter(c);
        if !c.ok || (p.flags & 8 != 0 && max_selector == 0) { return false; }
        max_symbol = max(max_symbol, p.max_symbol); parameters.push(p);
    }
    if length == 0 { return true; }
    let decoder = cram_range::start(c.bytes[c.at..]);
    var arena: u32[>..] = [];
    var lengths: i64[4] = [0; 4];
    for i in 4 { lengths[i] = cram_range::add_model(arena, 256); }
    let reverse_model = cram_range::add_model(arena, 2);
    let duplicate_model = cram_range::add_model(arena, 2);
    let selector_model = cram_range::add_model(arena, max_selector + 1);
    var models: u32[>..] = []; models.resize(65536, 0);
    var reversed: Reverse[>..] = [];
    let first = parameters[0];
    var last_length = -1; var at = 0;
    while at < output.len && decoder.input.ok {
        let selector = if first.flags & 8 != 0 { decoder.symbol(arena, selector_model, max_selector + 1) as i64 } else { 0 };
        let index = if flags & 2 != 0 { selectors[selector] as i64 } else { selector };
        if index >= parameters.len { return false; }
        let p = parameters[index];
        if last_length < 0 || p.flags & 4 == 0 {
            var n: u32 = 0;
            for i in 4 { n |= (decoder.symbol(arena, lengths[i], 256) as u32) << (i * 8); }
            last_length = n as i64;
        }
        if !decoder.input.ok || last_length == 0 || last_length > output.len - at { return false; }
        if flags & 4 != 0 && decoder.symbol(arena, reverse_model, 2) != 0 {
            reversed.push(Reverse { start: at, length: last_length });
        }
        if p.flags & 2 != 0 && decoder.symbol(arena, duplicate_model, 2) != 0 {
            if last_length > at { return false; }
            for i in last_length { output[at + i] = output[at - last_length + i]; }
            at += last_length; continue;
        }
        var context = p.context; var qcontext: u32 = 0; var previous: i64 = 0; var delta: i64 = 0;
        for i in last_length {
            if !decoder.input.ok { return false; }
            if models[context] == 0 { models[context] = cram_range::add_model(arena, max_symbol + 1) as u32; }
            let q = decoder.symbol(arena, models[context] as i64, max_symbol + 1) as i64;
            output[at + i] = first.qmap[q];
            qcontext = ((((qcontext as u64) << first.qshift) + first.qtab[q]) & 4294967295) as u32;
            context = ((((qcontext & first.qmask) as i64) << first.qloc) +
                first.ptab[min(last_length - i, 1023)] + first.dtab[min(delta, 255)] + (selector << first.sloc)) & 65535;
            if previous != q { delta += 1; }
            previous = q;
        }
        at += last_length;
    }
    if !decoder.input.ok { return false; }
    for r in reversed { reverse(output[r.start..r.start + r.length]); }
    true
}
