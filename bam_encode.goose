// BAM serialization shared by writers and text/container decoders.
namespace bam_encode;
import std;
import bam;
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

let EMPTY_CIGAR: u32[0] = [];
enum SequenceEncoding { Ascii, Codes, Sam }
extern fn pack_sam(input: const u8[:], output: u8[:]) -> i32 = "goose_seqair_pack_sam";
extern fn pack_codes(input: const u8[:], output: u8[:]) -> i32 = "goose_seqair_pack_codes";

// Caller-owned alignment data. Sequence is ASCII IUPAC in stored orientation;
// qualities are raw Phred bytes (empty means missing), CIGARs are native BAM
// words, and aux contains raw BAM tags. -1 means absent reference/position.
// No input is retained by encode or write. Empty sequence means absent SEQ.
//
// ```goose
// let cigar: u32[] = [2 << 4];
// let alignment = bam_encode::Alignment { name: "read", tid: 0, pos: 4,
//     cigar: cigar, sequence: "aN", .. };
// var bytes: u8[>..] = [];
// assert(bam_encode::encode(alignment, bytes, 1) == bam::Status.Record);
// assert(bytes.len > 32);
// ```
struct Alignment {
    name: const u8[:] = "*",
    tid: i32 = -1, pos: i64 = -1, flags: u16 = 0, mapq: u8 = 0,
    mate_tid: i32 = -1, mate_pos: i64 = -1, template_len: i32 = 0,
    cigar: const u32[:] = EMPTY_CIGAR, sequence: const u8[:] = "",
    qualities: const u8[:] = "", aux: const u8[:] = "",
}

// Encode a BAM payload (without block_size) into reusable caller storage.
// Check lengths, references, coordinates, bases, CIGAR and auxiliary framing.
// Supports CG:B:I long CIGARs; an existing CG tag must match the resolved
// long CIGAR and is retained. Output is cleared on failure.
//
// ```goose
// let alignment = bam_encode::Alignment { name: "unmapped", flags: 4,
//     sequence: "ACN", .. };
// var bytes: u8[>..] = [];
// assert(bam_encode::encode(alignment, bytes, 0) == bam::Status.Record);
// let fields, ok = bam::decode(bytes, 0);
// assert(ok && fields.seq_len == 3 && fields.flags == 4);
// ```
fn encode(a: const Alignment&, bytes: u8[>..]&, references: i64,
          max_bytes: i64 = 16777216) -> bam::Status {
    let status = encode_body(a, bytes, references, max_bytes, SequenceEncoding.Ascii);
    if status != bam::Status.Record { bytes.clear(); }
    status
}

fn bin(begin: i64, end: i64) -> u32 {
    if begin < 0 { return 4680; }
    let last = max(begin + 1, end) - 1;
    for i in 5 {
        let shift = 14 + i * 3;
        if begin >> shift == last >> shift { return ((((1 as u32) << (15 - i * 3)) - 1) / 7 + (begin >> shift) as u32) & 65535; }
    }
    0
}

fn encode_body(a: const Alignment&, bytes: u8[>..]&, references: i64,
               max_bytes: i64, encoding: SequenceEncoding) -> bam::Status {
    bytes.clear();
    if a.name.len > 254 || bstr::find_byte(a.name, 0) >= 0 ||
       a.tid < -1 || a.tid >= references || a.mate_tid < -1 || a.mate_tid >= references ||
       a.pos < -1 || a.pos > 2147483647 || a.mate_pos < -1 || a.mate_pos > 2147483647 ||
       (a.qualities.len != 0 && a.qualities.len != a.sequence.len) { return bam::Status.Invalid; }
    if a.sequence.len > max_bytes || a.cigar.len > max_bytes / 4 || a.aux.len > max_bytes {
        return bam::Status.TooLarge;
    }
    let cg, cg_len, valid_aux = bam::aux_cigar(a.aux, 0);
    if !valid_aux { return bam::Status.Invalid; }
    let long = a.cigar.len > 65535;
    var span: i64 = 0;
    var query: i64 = 0;
    for word in a.cigar {
        let op = word & 15;
        if op > 8 { return bam::Status.Invalid; }
        if op == 0 || op == 2 || op == 3 || op == 7 || op == 8 { span += word >> 4; }
        if op == 0 || op == 1 || op == 4 || op == 7 || op == 8 { query += word >> 4; }
    }
    if a.flags & 4 == 0 && a.cigar.len > 0 && a.sequence.len > 0 && query != a.sequence.len { return bam::Status.Invalid; }
    if long && (a.sequence.len > 268435455 || span > 268435455) { return bam::Status.TooLarge; }
    let ops = if long { 2 } else { a.cigar.len };
    let size = 32 + a.name.len + 1 + ops * 4 + (a.sequence.len + 1) / 2 + a.sequence.len +
               a.aux.len + (if long && cg == 0 { 8 + a.cigar.len * 4 } else { 0 });
    if size > max_bytes || size > 2147483647 { return bam::Status.TooLarge; }
    if long && cg > 0 {
        if cg_len != a.cigar.len { return bam::Status.Invalid; }
        for word, i in a.cigar { if bgzf::le32(a.aux, cg + 4 * i) != word { return bam::Status.Invalid; } }
    }
    bytes.resize(size, 0);
    bgzf::put32(bytes, 0, a.tid as! u32);
    bgzf::put32(bytes, 4, a.pos as! u32);
    bytes[8] = (a.name.len + 1) as u8;
    bytes[9] = a.mapq;
    bgzf::put16(bytes, 10, bin(a.pos, a.pos + (if a.flags & 4 != 0 { 1 } else { max(span, 1) })));
    bgzf::put16(bytes, 12, ops as u32);
    bgzf::put16(bytes, 14, a.flags);
    bgzf::put32(bytes, 16, a.sequence.len as u32);
    bgzf::put32(bytes, 20, a.mate_tid as! u32);
    bgzf::put32(bytes, 24, a.mate_pos as! u32);
    bgzf::put32(bytes, 28, a.template_len as! u32);
    bytes[32..32 + a.name.len].move_from(a.name);
    bytes[32 + a.name.len] = 0;
    var at = 33 + a.name.len;
    if long {
        bgzf::put32(bytes, at, (a.sequence.len as u32) << 4 | 4);
        bgzf::put32(bytes, at + 4, (span as u32) << 4 | 3);
    } else {
        for word, i in a.cigar { bgzf::put32(bytes, at + i * 4, word); }
    }
    at += ops * 4;
    let packed = bytes[at..at + (a.sequence.len + 1) / 2];
    if encoding == SequenceEncoding.Sam {
        if pack_sam(a.sequence, packed) == 0 { return bam::Status.Invalid; }
    } else if encoding == SequenceEncoding.Codes {
        if pack_codes(a.sequence, packed) == 0 { return bam::Status.Invalid; }
    } else {
        for base, i in a.sequence {
            let upper = if base >= 'a' && base <= 'z' { base - 32 } else { base };
            let code = bstr::find_byte("=ACMGRSVTWYHKDBN", upper);
            if code < 0 { return bam::Status.Invalid; }
            if i & 1 == 0 { packed[i / 2] = (code << 4) as u8; }
            else { packed[i / 2] |= code as u8; }
        }
    }
    at += (a.sequence.len + 1) / 2;
    if a.qualities.len == 0 {
        for i in a.sequence.len { bytes[at + i] = 255; }
    } else { bytes[at..at + a.sequence.len].move_from(a.qualities); }
    at += a.sequence.len;
    bytes[at..at + a.aux.len].move_from(a.aux);
    at += a.aux.len;
    if long && cg == 0 {
        bytes[at..at + 4].move_from("CGBI");
        bgzf::put32(bytes, at + 4, a.cigar.len as u32);
        for word, i in a.cigar { bgzf::put32(bytes, at + 8 + i * 4, word); }
    }
    // A caller may supply an already escaped CIGAR. Validate the resolved CG
    // too; ordinary and decoded long CIGARs were checked before serialization.
    if cg > 0 && !long {
        let fields, ok = bam::decode(bytes, references);
        if !ok { return bam::Status.Invalid; }
    }
    bam::Status.Record
}
