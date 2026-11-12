// Seqair's rANS 4x8 layout: reusable slot-to-symbol and symbol-to-step rows.
namespace cram_rans;
import std;
import cram_encoding;
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.0";

inline fn word(c: cram_encoding::Cursor&) -> u32 {
    let data = c.take(4); if !c.ok { return 0; } bgzf::le32(data, 0)
}

fn frequencies(c: cram_encoding::Cursor&, freq: u16[:]) {
    freq.fill(0);
    var sym = c.byte() as i64;
    var previous = sym;
    while c.ok {
        freq[sym] = (c.itf8() & 65535) as u16;
        sym = c.byte() as i64;
        if sym == 0 { break; }
        if sym == ((previous + 1) & 255) {
            let run = c.byte() as i64;
            if sym + run > 255 { c.ok = false; break; }
            for i in run { freq[sym] = (c.itf8() & 65535) as u16; sym += 1; }
        }
        previous = sym;
    }
}

fn row(freq: const u16[:], syms: u8[:], steps: u32[:]) {
    var total: u32 = 0;
    for f, i in freq { steps[i] = ((total & 65535) << 16) | f; total += f; }
    if total <= 4096 {
        var at = 0;
        for f, i in freq { syms[at..at + f].fill(i as u8); at += f; }
        syms[at..].fill(255);
    } else {
        var sym = 0;
        for i in 4096 {
            while sym < 255 && i >= (steps[sym + 1] >> 16) as i64 { sym += 1; }
            syms[i] = sym as u8;
        }
    }
}

inline fn advance(x: u32, entry: u32) -> u32 {
    // Explicit modulo arithmetic also works in assertion-enabled Goose builds.
    let value = ((entry & 65535) as u64) * ((x >> 12) as u64) + (x & 4095) as u64 +
                4294967296 - (entry >> 16) as u64;
    (value & 4294967295) as u32
}

inline fn renormalize(c: cram_encoding::Cursor&, state: u32&) {
    while state < 8388608 && c.ok { state = (state << 8) | c.byte(); }
}

fn decode(input: const u8[:], output: u8[:], symbols: u8[>..]&, steps: u32[>..]&) -> bool {
    if input.len < 9 || bgzf::le32(input, 1) as i64 != input.len - 9 ||
       bgzf::le32(input, 5) as i64 != output.len || input[0] > 1 { return false; }
    let c = cram_encoding::Cursor { bytes: input, at: 9 };
    let order = input[0];
    let rows = if order == 0 { 1 } else { 256 };
    symbols.resize(rows * 4096, 0); steps.resize(rows * 256, 0);
    var active: bool[256] = [false; 256];
    var freq: u16[256] = [0; 256];
    if order == 0 {
        frequencies(c, freq); row(freq, symbols[..4096], steps[..256]); active[0] = true;
    } else {
        var context = c.byte() as i64;
        var previous = context;
        while c.ok {
            frequencies(c, freq);
            row(freq, symbols[context * 4096..(context + 1) * 4096], steps[context * 256..(context + 1) * 256]);
            active[context] = true;
            context = c.byte() as i64;
            if context == 0 { break; }
            if context == ((previous + 1) & 255) {
                let run = c.byte() as i64;
                if context + run > 255 { c.ok = false; break; }
                for i in run {
                    frequencies(c, freq);
                    row(freq, symbols[context * 4096..(context + 1) * 4096], steps[context * 256..(context + 1) * 256]);
                    active[context] = true; context += 1;
                }
            }
            previous = context;
        }
    }
    if !c.ok { return false; }
    var states: u32[4] = [0; 4];
    for i in 4 { states[i] = word(c); }
    var prev: u8[4] = [0; 4];
    let quarter = output.len / 4;
    let rounds = if order == 0 { (output.len + 3) / 4 } else { quarter };
    for i in rounds {
        for lane in 4 {
            let at = if order == 0 { i * 4 + lane } else { lane * quarter + i };
            if at >= output.len { break; }
            let context = if order == 0 { 0 } else { prev[lane] as i64 };
            // Seqair's inactive rows behave like all-zero frequency tables.
            let sym = if active[context] { symbols[context * 4096 + (states[lane] & 4095) as i64] } else { 255 };
            let step = if active[context] { steps[context * 256 + sym] } else { 0 };
            output[at] = sym; states[lane] = advance(states[lane], step);
            renormalize(c, states[lane]); prev[lane] = sym;
            if !c.ok { return false; }
        }
    }
    if order == 1 {
        for tail in output.len - quarter * 4 {
            let at = quarter * 4 + tail;
            let context = prev[3] as i64;
            let sym = if active[context] { symbols[context * 4096 + (states[3] & 4095) as i64] } else { 255 };
            let step = if active[context] { steps[context * 256 + sym] } else { 0 };
            output[at] = sym; states[3] = advance(states[3], step); prev[3] = sym;
            if at + 1 < output.len { renormalize(c, states[3]); }
            if !c.ok { return false; }
        }
    }
    c.ok
}
