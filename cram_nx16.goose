// rANS Nx16 and PACK/RLE/STRIPE transforms, following Seqair 8401961.
// Shared order-1 scratch is reused across blocks and nested stripe streams.
namespace cram_nx16;
import std;
import cram_encoding;

extern fn packed_decode(input: const u8[:], output: u8[:], table: const u32[:],
                        bits: i32, lanes: i32, order: i32) -> i64 = "goose_seqair_nx16_packed";
extern fn unpack_nibbles(input: const u8[:], output: u8[:], mapping: const u8[:]) = "goose_seqair_unpack_nibbles";

// One load supplies symbol, frequency and bias. Non-exact tables retain the
// unpacked decoder's semantics rather than changing malformed-input behavior.
fn pack_row(frequencies: const u32[:], row: u32[:]) -> bool {
    var at = 0;
    for frequency, symbol in frequencies {
        if frequency == 0 { continue; }
        let end = at + frequency as i64;
        if end > row.len { return false; }
        let base = ((frequency - 1) << 20) | symbol as u32;
        for &slot, bias in row[at..end] { slot = base | ((bias as u32) << 8); }
        at = end;
    }
    at == row.len
}

fn packed(c: cram_encoding::Cursor&, output: u8[:], table: const u32[:], bits: i64, lanes: i64, order: i32) -> bool {
    if !c.ok { return false; }
    let consumed = packed_decode(c.bytes[c.at..], output, table, bits as i32, lanes as i32, order);
    if consumed < 0 { c.ok = false; return false; }
    c.at += consumed; true
}

fn uint7(c: cram_encoding::Cursor&) -> u32 {
    var n: u32 = 0;
    for i in 5 {
        let b = c.byte(); n = ((n as u64 << 7) & 4294967295) as u32 | (b & 127);
        if b & 128 == 0 { return n; }
    }
    c.ok = false; 0
}

fn word(c: cram_encoding::Cursor&) -> u32 {
    var n: u32 = 0;
    for i in 4 { n |= (c.byte() as u32) << (i * 8); }
    n
}

fn alphabet(c: cram_encoding::Cursor&, active: bool[:]) {
    active.fill(false);
    var symbol = c.byte() as i64;
    var previous = symbol;
    while c.ok {
        active[symbol] = true;
        symbol = c.byte() as i64;
        if symbol == 0 { break; }
        if symbol == ((previous + 1) & 255) {
            let n = c.byte() as i64;
            if symbol + n > 255 { c.ok = false; return; }
            for i in n { active[symbol] = true; symbol += 1; }
        }
        previous = symbol;
    }
}

fn normalize(freq: u32[:], bits: i64) -> bool {
    var sum: u64 = 0;
    for n in freq { sum += n; }
    if sum > 4294967295 { return false; }
    if sum == 0 { return true; }
    var shift = 0;
    while sum < ((1 as u64) << bits) { sum <<= 1; shift += 1; }
    for &n in freq { n <<= shift; }
    true
}

fn lookup(freq: const u32[:], cumulative: u32[:], symbols: u8[:]) {
    var sum: u64 = 0;
    symbols.fill(255);
    for n, i in freq {
        cumulative[i] = min(sum, 4294967295 as u64) as u32;
        let begin = min(sum, symbols.len as u64) as i64;
        sum += n;
        let end = min(sum, symbols.len as u64) as i64;
        symbols[begin..end].fill(i as u8);
    }
}

inline fn step(c: cram_encoding::Cursor&, state: u32, frequency: u32, cumulative: u32, bits: i64) -> u32 {
    let mask = ((1 as u32) << bits) - 1;
    var next = ((frequency as u64 * (state >> bits) as u64 + (state & mask) as u64 + 4294967296 - cumulative as u64) & 4294967295) as u32;
    while next < 32768 && c.ok {
        let low = c.byte() as u32; let high = c.byte() as u32;
        next = (next << 16) | low | (high << 8);
    }
    next
}

fn order0(c: cram_encoding::Cursor&, output: u8[:], lanes: i64) -> bool {
    var active: bool[256] = [false; 256]; alphabet(c, active);
    var frequencies: u32[256] = [0; 256]; var cumulative: u32[256] = [0; 256];
    var symbols: u8[4096] = [0; 4096];
    for on, i in active { if on { frequencies[i] = uint7(c); } }
    if !c.ok || !normalize(frequencies, 12) { return false; }
    var packed_row: u32[4096] = [0; 4096];
    if pack_row(frequencies, packed_row) { return packed(c, output, packed_row, 12, lanes, 0); }
    lookup(frequencies, cumulative, symbols);
    var states: u32[32] = [0; 32];
    for i in lanes { states[i] = word(c); }
    for &out, i in output {
        if !c.ok { return false; }
        let lane = i % lanes; let state = states[lane];
        let symbol = symbols[(state & 4095) as i64]; out = symbol;
        states[lane] = step(c, state, frequencies[symbol as i64], cumulative[symbol as i64], 12);
    }
    c.ok
}

fn frequencies1(c: cram_encoding::Cursor&, frequencies: u32[:], bits: i64, active: bool[:]) -> bool {
    alphabet(c, active);
    frequencies.fill(0);
    for on, row in active {
        if !on { continue; }
        var skipped = 0;
        for included, symbol in active {
            if !included { continue; }
            if skipped > 0 { skipped -= 1; continue; }
            let f = uint7(c); frequencies[row * 256 + symbol] = f;
            if f == 0 { skipped = c.byte() as i64; }
        }
        if !c.ok || !normalize(frequencies[row * 256..row * 256 + 256], bits) { return false; }
    }
    c.ok
}

fn order1(c: cram_encoding::Cursor&, output: u8[:], lanes: i64,
          frequencies: u32[>..]&, cumulative: u32[>..]&, symbols: u8[>..]&, limit: i64) -> bool {
    let header = c.byte(); let bits = (header >> 4) as i64;
    if bits > 12 { return false; }
    frequencies.resize(65536, 0);
    var active: bool[256] = [false; 256];
    if header & 1 != 0 {
        let size = uint7(c) as i64; let compressed = uint7(c) as i64;
        let input = cram_encoding::Cursor { bytes: c.take(compressed) };
        if !c.ok || size > limit { return false; }
        var decoded: u8[>..] = []; decoded.resize(size, 0);
        if !order0(input, decoded, 4) { return false; }
        let table = cram_encoding::Cursor { bytes: decoded };
        if !frequencies1(table, frequencies, bits, active) { return false; }
    } else if !frequencies1(c, frequencies, bits, active) { return false; }
    // Reuse cumulative scratch for packed rows; row stride follows precision,
    // avoiding the cache-set collisions of always spacing rows by 4096.
    let row_width = 1 << bits;
    cumulative.resize(max(65536, 256 * row_width), 0);
    var exact = active[0];
    for on, row in active {
        if on && !pack_row(frequencies[row * 256..row * 256 + 256], cumulative[row * row_width..(row + 1) * row_width]) {
            exact = false;
        }
    }
    if exact { return packed(c, output, cumulative, bits, lanes, 1); }
    symbols.resize(1048576, 0);
    for row in 256 {
        if !active[row] {
            cumulative[row * 256..row * 256 + 256].fill(0);
            symbols[row * 4096..row * 4096 + (1 << bits)].fill(255);
            continue;
        }
        lookup(frequencies[row * 256..row * 256 + 256], cumulative[row * 256..row * 256 + 256],
               symbols[row * 4096..row * 4096 + (1 << bits)]);
    }
    var states: u32[32] = [0; 32]; var previous: u8[32] = [0; 32];
    for i in lanes { states[i] = word(c); }
    let width = output.len / lanes; let mask = ((1 as u32) << bits) - 1;
    for i in width {
        for lane in lanes {
            if !c.ok { return false; }
            let row = previous[lane] as i64; let state = states[lane];
            let symbol = symbols[row * 4096 + (state & mask) as i64];
            output[lane * width + i] = symbol;
            states[lane] = step(c, state, frequencies[row * 256 + symbol], cumulative[row * 256 + symbol], bits);
            previous[lane] = symbol;
        }
    }
    let lane = lanes - 1;
    for &out in output[width * lanes..] {
        if !c.ok { return false; }
        let row = previous[lane] as i64; let state = states[lane];
        let symbol = symbols[row * 4096 + (state & mask) as i64]; out = symbol;
        states[lane] = step(c, state, frequencies[row * 256 + symbol], cumulative[row * 256 + symbol], bits);
        previous[lane] = symbol;
    }
    c.ok
}

// Seqair's 256-entry byte expansion table; fixed widths keep each decoded
// byte to one 8/4-byte store. Nibbles use the shared SIMD kernel instead.
inline fn unpack_table(input: const u8[:], mapping: const u8[:], output: u8[:], width: i64) {
    var table: u8[2048] = [0; 2048];
    let shift = 8 / width; let mask = (1 << shift) - 1;
    for byte in 256 {
        for i in width {
            let index = (byte >> (i * shift)) & mask;
            if index < mapping.len { table[byte * width + i] = mapping[index]; }
        }
    }
    let whole = min(input.len, output.len / width);
    for byte, i in input[..whole] {
        let begin = byte as i64 * width;
        output[i * width..(i + 1) * width].move_from(table[begin..begin + width]);
    }
    if whole < input.len {
        let begin = input[whole] as i64 * width; let remaining = output.len - whole * width;
        output[whole * width..].move_from(table[begin..begin + remaining]);
    }
}

fn unpack(input: const u8[:], mapping: const u8[:], output: u8[:]) -> bool {
    let n = mapping.len;
    if n < 1 || n > 16 { return false; }
    if n == 1 { output.fill(mapping[0]); return true; }
    output.fill(0);
    if output.len >= 1024 {
        if n == 2 { unpack_table(input, mapping, output, 8); }
        else if n <= 4 { unpack_table(input, mapping, output, 4); }
        else { unpack_nibbles(input, output, mapping); }
        return true;
    }
    if n == 2 {
        for &out, i in output[..min(output.len, input.len * 8)] {
            out = mapping[((input[i >> 3] >> (i & 7)) & 1) as i64];
        }
    } else if n <= 4 {
        for &out, i in output[..min(output.len, input.len * 4)] {
            let index = ((input[i >> 2] >> ((i & 3) * 2)) & 3) as i64;
            if index < n { out = mapping[index]; }
        }
    } else {
        for &out, i in output[..min(output.len, input.len * 2)] {
            let index = ((input[i >> 1] >> ((i & 1) * 4)) & 15) as i64;
            if index < n { out = mapping[index]; }
        }
    }
    true
}

// The caller supplies 16 bytes of writable slack, as in Seqair's apply_rle.
// Short runs/literal stretches use fixed-size stores; later writes overwrite
// their excess. Neither the slack nor unread input is exposed to the caller.
fn expand(input: const u8[:], metadata: const u8[:], output: u8[:], length: i64) -> bool {
    let c = cram_encoding::Cursor { bytes: metadata };
    let size = c.byte() as i64;
    var alphabet: bool[256] = [false; 256];
    for i in (if size == 0 { 256 } else { size }) { alphabet[c.byte() as i64] = true; }
    if !c.ok { return false; }
    var at = 0; var read = 0;
    while at < length {
        let room = length - at; let available = min(room, input.len - read);
        var literals = 0;
        while literals < available && !alphabet[input[read + literals] as i64] { literals += 1; }
        if literals <= 16 && input.len - read >= 16 {
            output[at..at + 16].move_from(input[read..read + 16]);
        } else { output[at..at + literals].move_from(input[read..read + literals]); }
        at += literals; read += literals;
        if literals == available { return at == length; }
        let byte = input[read]; read += 1;
        let extra = uint7(c) as i64;
        if !c.ok { return false; }
        let n = min(extra + 1, length - at);
        if n <= 16 { output[at..at + 16].fill(byte); }
        else { output[at..at + n].fill(byte); }
        at += n;
    }
    true
}

fn leaf(input: const u8[:], output: u8[:], frequencies: u32[>..]&, cumulative: u32[>..]&, symbols: u8[>..]&,
        limit: i64) -> bool {
    let c = cram_encoding::Cursor { bytes: input };
    let flags = c.byte(); let lanes = if flags & 4 != 0 { 32 } else { 4 };
    let size = if flags & 16 == 0 { uint7(c) as i64 } else { output.len };
    if !c.ok || size != output.len || size > limit { return false; }
    var length = size; var mapping: u8[>..] = [];
    if flags & 128 != 0 {
        let n = c.byte() as i64; mapping.append(c.take(n));
        length = uint7(c) as i64;
        if n == 0 || n > 16 || length > limit { return false; }
    }
    let expanded_size = length;
    var metadata: u8[>..] = [];
    if flags & 64 != 0 {
        let header = uint7(c); let n = (header >> 1) as i64;
        length = uint7(c) as i64;
        if n > limit || length > limit { return false; }
        if header & 1 == 0 {
            let compressed = uint7(c) as i64;
            let m = cram_encoding::Cursor { bytes: c.take(compressed) };
            metadata.resize(n, 0);
            if !c.ok || !order0(m, metadata, lanes) { return false; }
        } else { metadata.append(c.take(n)); }
    }
    if !c.ok { return false; }
    var body: u8[>..] = []; body.resize(length, 0);
    if flags & 32 != 0 {
        let bytes = c.take(length); if !c.ok { return false; } body[..].move_from(bytes);
    } else if flags & 1 == 0 {
        if !order0(c, body, lanes) { return false; }
    } else if !order1(c, body, lanes, frequencies, cumulative, symbols, limit) { return false; }
    var expanded: u8[>..] = [];
    if flags & 64 != 0 {
        expanded.resize(expanded_size + 16, 0);
        if !expand(body, metadata, expanded, expanded_size) { return false; }
        expanded.resize(expanded_size, 0);
    }
    let transformed = if flags & 64 != 0 { expanded[..] } else { body[..] };
    if flags & 128 != 0 { return unpack(transformed, mapping, output); }
    if transformed.len != output.len { return false; }
    output.move_from(transformed); true
}

// Explicit bounded traversal avoids recursive growable owners in Goose. A
// nested stripe composes its output origin/stride with its parent's mapping.
struct Stripe { begin: i64, end: i64, length: i64, origin: i64, stride: i64, depth: i64 }

fn decode(input: const u8[:], output: u8[:], frequencies: u32[>..]&, cumulative: u32[>..]&, symbols: u8[>..]&,
          limit: i64) -> bool {
    var pending: Stripe[>..] = [];
    pending.push(Stripe { begin: 0, end: input.len, length: output.len, origin: 0, stride: 1, depth: 0 });
    var part: u8[>..] = [];
    while pending.len > 0 {
        let item = pending[pending.len - 1]; pending.resize(pending.len - 1, Stripe { .. });
        let bytes = input[item.begin..item.end];
        let c = cram_encoding::Cursor { bytes: bytes };
        let flags = c.byte();
        let size = if flags & 16 == 0 { uint7(c) as i64 } else { item.length };
        if !c.ok || size != item.length || size > limit { return false; }
        if flags & 8 != 0 {
            if item.depth >= 4 { return false; }
            let n = c.byte() as i64;
            if n == 0 { return false; }
            var sizes: i64[256] = [0; 256];
            for i in n { sizes[i] = uint7(c) as i64; }
            for i in n {
                let begin = c.at; c.take(sizes[i]);
                if !c.ok { return false; }
                pending.push(Stripe { begin: item.begin + begin, end: item.begin + c.at,
                    length: size / n + (if i < size % n { 1 } else { 0 }),
                    origin: item.origin + i * item.stride, stride: item.stride * n, depth: item.depth + 1 });
            }
        } else if item.stride == 1 {
            if !leaf(bytes, output, frequencies, cumulative, symbols, limit) { return false; }
        } else {
            part.resize(size, 0);
            if !leaf(bytes, part, frequencies, cumulative, symbols, limit) { return false; }
            for byte, i in part { output[item.origin + i * item.stride] = byte; }
        }
    }
    true
}
