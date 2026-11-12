// Shared adaptive range coder; arithmetic uses Seqair's depth-selected search.
namespace cram_range;
import std;
import cram_encoding;
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import simd_bytes.simd_int;

struct Decoder { input: cram_encoding::Cursor, range: u32 = 4294967295, code: u32 = 0 }

struct Found { index: i64, cumulative: u32 }

// Seqair range_coder.rs: paired eight-lane inclusive sums, one exit test
// per pair, unsigned comparisons and reductions. Total remains below 65536.
simd fn find_vector(entries: const u32[:], cumulative: u32, target: u32) -> Found {
    let t = simd_int::splat(target as u16);
    var carry = simd_int::splat(cumulative as u16);
    var index = 4;
    while index + 16 <= entries.len {
        let p0 = simd_int::load_low_u16(entries[index..]).prefix_sum();
        let p1 = simd_int::load_low_u16(entries[index + 8..]).prefix_sum();
        let sums0 = p0.add(carry);
        let carry1 = carry.add(p0.broadcast_last());
        let sums1 = p1.add(carry1);
        let before1 = sums1.le(t); let bits1 = before1.to_bits();
        if bits1 != 255 {
            let before0 = sums0.le(t); let bits0 = before0.to_bits();
            let below = before0.select(sums0, carry).max(before1.select(sums1, carry));
            return Found { index: index + ctz((~(bits0 | (bits1 << 8))) & 65535), cumulative: below.reduce_max() };
        }
        carry = carry1.add(p1.broadcast_last()); index += 16;
    }
    if index + 8 <= entries.len {
        let sums = simd_int::load_low_u16(entries[index..]).prefix_sum().add(carry);
        let before = sums.le(t); let bits = before.to_bits();
        if bits != 255 {
            return Found { index: index + ctz((~bits) & 255), cumulative: before.select(sums, carry).reduce_max() };
        }
        carry = sums.broadcast_last(); index += 8;
    }
    var sum: u32 = carry.reduce_max();
    while index < entries.len && sum + (entries[index] & 65535) <= target {
        sum += entries[index] & 65535; index += 1;
    }
    Found { index: index, cumulative: sum }
}

fn start(input: const u8[:]) -> Decoder {
    let decoder = Decoder { input: cram_encoding::Cursor { bytes: input }, .. };
    for i in 5 { decoder.code = ((decoder.code as u64 << 8) & 4294967295) as u32 | decoder.input.byte(); }
    decoder
}

// One total followed by n packed (frequency:u16, symbol:u16) entries. Arena
// offsets let FQZ instantiate only the quality contexts actually encountered.
fn add_model(arena: u32[>..]&, n: i64) -> i64 {
    let first = arena.len; arena.push(n as u32);
    for i in n { arena.push(((i as u32) << 16) | 1); }
    first
}

inline fn symbol(decoder: Decoder&, arena: u32[:], first: i64, n: i64) -> i32 {
    var depth: u16[0] = [];
    return symbol(decoder, arena, first, n, depth);
}

fn symbol(decoder: Decoder&, arena: u32[:], first: i64, n: i64, depth: u16[:]) -> i32 {
    if !decoder.input.ok { return 0; }
    let total = arena[first];
    var target: u32 = 0;
    if total > 0 && decoder.range >= total {
        decoder.range /= total; target = decoder.code / decoder.range;
    }
    if target >= total { decoder.input.ok = false; return 0; }
    let entries = arena[first + 1..first + 1 + n];
    var cumulative: u32 = 0; var index = 0;
    if depth.len > 0 {
        while index < min(4, n) && cumulative + (entries[index] & 65535) <= target {
            cumulative += entries[index] & 65535; index += 1;
        }
        if index == 4 && n > 8 && depth[0] >= 512 {
            let found = find_vector(entries, cumulative, target);
            cumulative = found.cumulative; index = found.index;
        }
    }
    while index < n && cumulative + (entries[index] & 65535) <= target {
        cumulative += entries[index] & 65535; index += 1;
    }
    if index == n { decoder.input.ok = false; return 0; }
    if depth.len > 0 && index >= 4 {
        depth[0] = depth[0] - (depth[0] >> 4) + min(index, 255) as u16;
    }
    let frequency = entries[index] & 65535; let symbol = entries[index] >> 16;
    let product = ((cumulative as u64 * decoder.range as u64) & 4294967295);
    decoder.code = ((decoder.code as u64 + 4294967296 - product) & 4294967295) as u32;
    decoder.range = ((decoder.range as u64 * frequency as u64) & 4294967295) as u32;
    while decoder.range < 16777216 && decoder.input.ok {
        decoder.code = ((decoder.code as u64 << 8) & 4294967295) as u32 | decoder.input.byte();
        decoder.range = ((decoder.range as u64 << 8) & 4294967295) as u32;
    }
    if !decoder.input.ok { return 0; }
    entries[index] += 16; arena[first] += 16;
    if arena[first] > 65519 {
        var sum: u32 = 0;
        for &entry in entries {
            let freq = ((entry & 65535) + 1) >> 1;
            entry = (entry & 4294901760) | freq; sum += freq;
        }
        arena[first] = sum;
    }
    if index > 0 && (entries[index] & 65535) > (entries[index - 1] & 65535) {
        let previous = entries[index - 1]; entries[index - 1] = entries[index]; entries[index] = previous;
    }
    symbol as i32
}
