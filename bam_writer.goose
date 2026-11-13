// Seqair-style eager header, reusable record buffer and BGZF output.
namespace bam_writer;
import std;
import bam;
import bam_encode;
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";
import pileup;
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

// Owns one output and its compressor. Construct with Writer { .. }, open,
// write records in caller order, and check close (which flushes and writes EOF).
// Errors poison subsequent writes until close/open. No sorting/index building
// or parallel compression is implicit. max_record_bytes defaults to 16 MiB.
//
// ```goose
// let writer = bam_writer::Writer { .. };
// assert(writer.close() == bam::Status.Ok);
// ```
root struct Writer {
    stream: bgzf::Writer,
    status: bam::Status = bam::Status.Closed,
    buffer: u8[>..] = [],
    references: i64 = 0,
    max_record_bytes: i64 = 16777216,
}

// Flush/finish, close the OS handle and release the compressor. Check the
// result; deferred cleanup alone cannot report a failed final write. A poisoned
// output is closed without an EOF marker. Repeated close succeeds.
//
// ```goose
// let writer = bam_writer::Writer { .. };
// assert(writer.close() == bam::Status.Ok);
// assert(writer.close() == bam::Status.Ok);
// ```
fn close(writer: Writer&) -> bam::Status {
    let previous = writer.status;
    if previous != bam::Status.Ok && previous != bam::Status.Closed { writer.stream.status = bgzf::Status.Invalid; }
    let result = bam::from_bgzf(writer.stream.close());
    writer.status = bam::Status.Closed;
    if previous != bam::Status.Ok && previous != bam::Status.Closed { return previous; }
    result
}

fn bytes(writer: Writer&, input: const u8[:]) -> bam::Status {
    writer.status = bam::from_bgzf(writer.stream.write(input));
    writer.status
}

fn word(writer: Writer&, value: u32) -> bam::Status {
    var data: u8[4] = [0; 4];
    bgzf::put32(data, 0, value);
    writer.bytes(data)
}

// Create/truncate a BAM, writing its header eagerly. Reference name offsets
// address names exactly as in bam::Reader. Header/reference order must match
// the IDs in all subsequent records. Levels 0..12 use libdeflate (default 6).
// A header is bounded by the reader's default limits and checked before opening.
//
// ```goose
// let writer = bam_writer::Writer { .. };
// let refs: bam::Reference[] = [];
// assert(writer.open("", "", refs, "") == bam::Status.IoError);
// assert(writer.close() == bam::Status.IoError);
// ```
fn open(writer: Writer&, path: const u8[:], header: const u8[:],
        references: const bam::Reference[:], names: const u8[:], level: i32 = 6) -> bam::Status {
    let previous = writer.close();
    if previous != bam::Status.Ok { return previous; }
    writer.references = references.len;
    if header.len > 16777216 || names.len > 16777216 || references.len > 1048576 {
        writer.status = bam::Status.TooLarge; return writer.status;
    }
    for &ref in references {
        if ref.name < 0 || ref.name_len < 1 || ref.name > names.len || ref.name_len > names.len - ref.name ||
           ref.length < 0 || ref.length > 2147483647 {
            writer.status = bam::Status.Invalid; return writer.status;
        }
        if bstr::find_byte(names[ref.name..ref.name + ref.name_len], 0) >= 0 {
            writer.status = bam::Status.Invalid; return writer.status;
        }
    }
    writer.status = bam::from_bgzf(writer.stream.open(path, level));
    if writer.status != bam::Status.Ok { return writer.status; }
    writer.bytes("BAM\x01");
    writer.word(header.len as u32);
    writer.bytes(header);
    writer.word(references.len as u32);
    for &ref in references {
        writer.word((ref.name_len + 1) as u32);
        writer.bytes(names[ref.name..ref.name + ref.name_len]);
        writer.bytes("\0");
        writer.word(ref.length as u32);
    }
    writer.status
}

// Write an already validated reader record byte-for-byte (including aux/CG).
// Only target bounds and size are rechecked, as in Seqair write_raw. The two
// headers must name the same targets in the same order; IDs are not remapped.
//
// ```goose
// let writer = bam_writer::Writer { .. };
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Ok);
// assert(reader.read() == bam::Status.Record);
// let record = reader.record();
// assert(writer.write(record) == bam::Status.Closed);
// ```
fn write(writer: Writer&, record: const bam::Record&) -> bam::Status {
    if writer.status != bam::Status.Ok { return writer.status; }
    let f = record.fields;
    if f.tid < -1 || f.tid >= writer.references || f.mate_tid < -1 || f.mate_tid >= writer.references {
        writer.status = bam::Status.Invalid; return writer.status;
    }
    writer.payload(record.bytes)
}

fn payload(writer: Writer&, data: const u8[:]) -> bam::Status {
    if data.len > writer.max_record_bytes || data.len > 2147483647 {
        writer.status = bam::Status.TooLarge; return writer.status;
    }
    if data.len < 32 { writer.status = bam::Status.Invalid; return writer.status; }
    writer.status = bam::from_bgzf(writer.stream.reserve_record(data.len + 4));
    if writer.status != bam::Status.Ok { return writer.status; }
    writer.word(data.len as u32);
    writer.bytes(data)
}

// Serialize caller-owned fields using a reused payload buffer. No borrowed
// alignment data survives the call. Invalid input poisons the writer.
//
// ```goose
// let writer = bam_writer::Writer { .. };
// let alignment = bam_encode::Alignment { flags: 4, .. };
// assert(writer.write(alignment) == bam::Status.Closed);
// ```
fn write(writer: Writer&, alignment: const bam_encode::Alignment&) -> bam::Status {
    if writer.status != bam::Status.Ok { return writer.status; }
    writer.status = bam_encode::encode(alignment, writer.buffer, writer.references, writer.max_record_bytes);
    if writer.status != bam::Status.Ok { return writer.status; }
    writer.payload(writer.buffer)
}

// Repack a slab-backed record directly, without building an owned intermediate
// record or converting its 0..15 base codes into ASCII and back. Name, CIGAR,
// quality, aux, mate fields and stored orientation are preserved.
//
// ```goose
// let writer = bam_writer::Writer { .. };
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Ok);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Ok);
// let record = store.record(0);
// assert(writer.write(record) == bam::Status.Closed);
// ```
fn write(writer: Writer&, record: const pileup::Record&) -> bam::Status {
    if writer.status != bam::Status.Ok { return writer.status; }
    let f = record.fields;
    let alignment = bam_encode::Alignment { name: record.name_bytes, tid: f.tid, pos: f.pos,
        flags: f.flags, mapq: f.mapq, mate_tid: f.mate_tid, mate_pos: f.mate_pos,
        template_len: f.template_len, cigar: record.cigar_words, sequence: record.bases,
        qualities: record.quality_bytes, aux: record.aux_bytes };
    writer.status = bam_encode::encode_body(alignment, writer.buffer, writer.references, writer.max_record_bytes, bam_encode::SequenceEncoding.Codes);
    if writer.status != bam::Status.Ok { return writer.status; }
    writer.payload(writer.buffer)
}

// Virtual offset at the current output boundary. close/flush may move a
// boundary to the next block; this is not an index builder.
//
// ```goose
// let writer = bam_writer::Writer { .. };
// assert(writer.tell() == 0);
// ```
fn tell(writer: const Writer&) -> u64 { writer.stream.tell() }
