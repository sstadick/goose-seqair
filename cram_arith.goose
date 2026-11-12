// CRAM method 6: Seqair's adaptive arithmetic coder and transform ordering.
namespace cram_arith;
import std;
import cram_encoding;
import cram_nx16;
import cram_range;

extern fn external_decode(method: i32, input: const u8[:], output: u8[:], limit: i64) -> i32 = "goose_seqair_cram_external";

fn body(input: const u8[:], output: u8[:], flags: u8) -> bool {
    if input.len == 0 { return output.len == 0; }
    if flags & 32 != 0 {
        if input.len < output.len { return false; }
        output.move_from(input[..output.len]); return true;
    }
    if flags & 4 != 0 { return external_decode(2, input, output, 268435456) == 1; }
    if output.len == 0 { return true; }
    let n = if input[0] == 0 { 256 } else { input[0] as i64 };
    let order1 = flags & 3 == 1; let rle = flags & 64 != 0;
    let decoder = cram_range::start(input[1..]);
    var arena: u32[>..] = [];
    let contexts = if order1 { n } else { 1 };
    for i in contexts { cram_range::add_model(arena, n); }
    let runs = arena.len;
    if rle { for i in 258 { cram_range::add_model(arena, 4); } }
    var depth: u16[256] = [0; 256];
    var context = 0; var at = 0;
    while at < output.len && decoder.input.ok {
        let symbol = decoder.symbol(arena, context * (n + 1), n, depth[context..context + 1]) as i64;
        output[at] = symbol as u8; at += 1;
        if order1 { context = symbol; }
        if rle {
            var run = 0; var rctx = symbol;
            loop {
                let part = decoder.symbol(arena, runs + rctx * 5, 4) as i64;
                rctx = if rctx == symbol { 256 } else { 257 };
                run += part;
                if !decoder.input.ok || part != 3 || run >= output.len { break; }
            }
            let len = min(run, output.len - at);
            output[at..at + len].fill(symbol as u8); at += len;
        }
    }
    decoder.input.ok
}

fn leaf(input: const u8[:], output: u8[:]) -> bool {
    let c = cram_encoding::Cursor { bytes: input };
    let flags = c.byte();
    let len = if flags & 16 != 0 { output.len } else { cram_nx16::uint7(c) as i64 };
    if !c.ok || len != output.len { return false; }
    if flags & 128 != 0 {
        let n = c.byte() as i64; let map = c.take(n);
        let packed = cram_nx16::uint7(c) as i64;
        if !c.ok || n == 0 || n > 16 || packed > len { return false; }
        let needed = if n == 1 { 0 } else if n == 2 { (len + 7) / 8 } else if n <= 4 { (len + 3) / 4 } else { (len + 1) / 2 };
        if packed < needed { return false; }
        var bytes: u8[>..] = []; bytes.resize(packed, 0);
        if !body(c.bytes[c.at..], bytes, flags) { return false; }
        return cram_nx16::unpack(bytes, map, output);
    }
    body(c.bytes[c.at..], output, flags)
}

fn decode(input: const u8[:], output: u8[:], limit: i64) -> bool {
    var pending: cram_nx16::Stripe[>..] = [];
    pending.push(cram_nx16::Stripe { begin: 0, end: input.len, length: output.len, origin: 0, stride: 1, depth: 0 });
    var part: u8[>..] = [];
    while pending.len > 0 {
        let item = pending[pending.len - 1]; pending.resize(pending.len - 1, cram_nx16::Stripe { .. });
        let bytes = input[item.begin..item.end];
        let c = cram_encoding::Cursor { bytes: bytes };
        let flags = c.byte();
        // Arithmetic STRIPE stores its length even when NOSZ is set.
        let size = if flags & 8 != 0 || flags & 16 == 0 { cram_nx16::uint7(c) as i64 } else { item.length };
        if !c.ok || size != item.length || size > limit { return false; }
        if flags & 8 != 0 {
            if item.depth >= 4 { return false; }
            let n = c.byte() as i64; if n == 0 { return false; }
            var sizes: i64[256] = [0; 256];
            for i in n { sizes[i] = cram_nx16::uint7(c) as i64; }
            for i in n {
                let begin = c.at; c.take(sizes[i]);
                if !c.ok { return false; }
                pending.push(cram_nx16::Stripe { begin: item.begin + begin, end: item.begin + c.at,
                    length: size / n + (if i < size % n { 1 } else { 0 }),
                    origin: item.origin + i * item.stride, stride: item.stride * n, depth: item.depth + 1 });
            }
        } else if item.stride == 1 {
            if !leaf(bytes, output) { return false; }
        } else {
            part.resize(size, 0);
            if !leaf(bytes, part) { return false; }
            for byte, i in part { output[item.origin + i * item.stride] = byte; }
        }
    }
    true
}
