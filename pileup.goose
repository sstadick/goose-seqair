// Load alignment records for a genomic region, then visit each reference
// position together with the reads that overlap it. This is a pileup.
// Input records must be sorted by reference ID and alignment start.
namespace pileup;
import std;
import dictionary;
import bam;
import bam_index;
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";

// Choose which alignment records to keep when loading a region.
// exclude_flags rejects a record when any of these SAM flag bits are set.
// The default mask, 0x704, rejects unmapped, secondary, quality-control-failed,
// and duplicate records. Supplementary alignments are kept.
// min_mapq is the minimum mapping quality; it does not filter base qualities.
// Overlapping mates are not combined when producing pileup positions.
//
// Unmapped records are always skipped, even with exclude_flags set to zero.
// By default, records whose CIGAR consumes no reference positions are also
// skipped. Set retain_zero_span to keep those mapped records in input order,
// for example when reproducing another tool's read-count limit. They remain
// available through Store.record(), but do not contribute to pileup positions.
//
// ```goose
// let filter = pileup::Filter { min_mapq: 20, .. };
// assert(filter.exclude_flags == 0x704 && filter.min_mapq == 20);
// ```
struct Filter {
    exclude_flags:    u16  = 0x704,
    min_mapq:         u8   = 0,
    retain_zero_span: bool = false,
}

// Metadata for one stored alignment. tid is its reference's index in the BAM
// header. pos is the zero-based alignment start; end is the exclusive end.
//
// Variable-length data lives in separate arrays owned by Store. The name,
// cigar, seq, qual, and aux fields are offsets into those arrays, not file
// positions. Offsets and lengths count array elements: bytes for names, bases,
// qualities, and tags; u32 operations for CIGARs. Use Store.record() to access
// those values without calculating the array slices yourself.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// assert(store.records[0].pos == 5 && store.records[0].seq_len == 13);
// ```
struct SlimRecord {
    tid:           i32,
    pos:           i64,
    end:           i64,
    reference_len: i64,
    flags:         u16,
    mapq:          u8,
    mate_tid:      i32,
    mate_pos:      i64,
    template_len:  i32,
    name:          u32,
    name_len:      u32,
    cigar:         u32,
    cigar_len:     u32,
    seq:           u32,
    seq_len:       u32,
    qual:          u32,
    aux:           u32,
    aux_len:       u32,
    // Index in Store.records of the matched mate, or -1 if no mate was linked.
    mate_index: i64 = -1,
}

// Owns all records loaded for one genomic region. Loading another region
// replaces those records while reusing the allocated storage.
//
// Names, unpacked base codes, CIGAR operations, base qualities, and raw
// auxiliary tags each occupy a separate growable array. These arrays are also
// called slabs: each holds one kind of data for every record in the region.
// Each base code is a byte from 0 to 15; record.base() converts it to a letter.
// The original packed BAM records are not kept. Header data belongs to the
// reader, but the stored records remain usable after that reader closes.
//
// max_bytes limits the combined bytes in the five data arrays. It does not
// include record metadata, spare allocation capacity, or other readers and
// workers. Its default and hard upper bound are 2^32-1 bytes (4 GiB minus one),
// so the u32 offsets into these arrays remain representable. The arrays grow
// as needed; setting the limit does not allocate that much memory. Exceeding
// the limit fails the load instead of returning a partial set of records.
// Set max_bytes before loading; otherwise treat the fields as read-only.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// assert(store.records.len == 1);
// reader.close(); // Stored records own their bytes, independently of Reader.
// let record = store.record(0);
// assert(record.name() == "match");
// ```
root struct Store {
    records:    SlimRecord[>..] = [],
    names:      u8[>..]         = [],
    bases:      u8[>..]         = [],
    cigars:     u32[>..]        = [],
    qualities:  u8[>..]         = [],
    auxiliary:  u8[>..]         = [],
    mate_stats: MateLinkStats   = MateLinkStats { .. },
    tid:        i32             = -1,
    begin:      i64             = 0,
    end:        i64             = 0,
    max_bytes:  i64             = 4294967295,
}

// A view of one alignment owned by Store. Its slices refer to the store's
// arrays rather than copying their contents. Closing the BAM reader does not
// invalidate this view; loading, clearing, or modifying the store does.
// Copy any data you need to keep before reusing the store for another region.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// let record = store.record(0);
// assert(record.fields.seq_len == record.bases.len);
// ```
struct Record {
    fields:        SlimRecord,
    name_bytes:    const u8[:],
    bases:         const u8[:],
    cigar_words:   const u32[:],
    quality_bytes: const u8[:],
    aux_bytes:     const u8[:],
}

// Access a loaded alignment without copying it. index must be between zero
// and store.records.len - 1. You can iterate that range directly, or use the
// record_index from a pileup entry to find the read that contributed it.
// The returned view is valid until the store is cleared, loaded, or modified.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// let record = store.record(0);
// assert(record.base(1) == 'C');
// ```
fn record(store: const Store&, index: i64) -> Record {
    let f = store.records[index];
    Record {
        fields:        f,
        name_bytes:    store.names[(f.name as i64)..f.name as i64 + f.name_len],
        bases:         store.bases[(f.seq as i64)..f.seq as i64 + f.seq_len],
        cigar_words:   store.cigars[(f.cigar as i64)..f.cigar as i64 + f.cigar_len],
        quality_bytes: store.qualities[(f.qual as i64)..f.qual as i64 + f.seq_len],
        aux_bytes:     store.auxiliary[(f.aux as i64)..f.aux as i64 + f.aux_len],
    }
}

// Read name without a trailing NUL.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// let record = store.record(0);
// assert(record.name() == "match");
// ```
fn name(record: const Record&) -> const u8[:] { record.name_bytes }

// Return a resolved CIGAR word. index must be in [0, fields.cigar_len).
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// let record = store.record(0);
// assert(record.cigar(0) == 13 << 4);
// ```
fn cigar(record: const Record&, index: i64) -> u32 { record.cigar_words[index] }

// Convert an unpacked IUPAC code to its SAM byte, preserving all 16 codes.
// index must be in [0, fields.seq_len); sequence stays in stored orientation.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// let record = store.record(0);
// assert(record.base(1) == 'C');
// ```
fn base(record: const Record&, index: i64) -> u8 { "=ACMGRSVTWYHKDBN"[record.bases[index] as i64] }

// Raw Phred quality; 255 means absent. index must be in [0, fields.seq_len).
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// let record = store.record(0);
// assert(record.quality(1) == 20);
// ```
fn quality(record: const Record&, index: i64) -> u8 { record.quality_bytes[index] }

// Raw auxiliary bytes, including the original CG array for long CIGARs.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// let record = store.record(0);
// assert(record.aux() == "NMi\x03\x00\x00\x00RGZgroup\x00");
// ```
fn aux(record: const Record&) -> const u8[:] { record.aux_bytes }

fn clear(store: Store&) {
    store.records.clear();
    store.names.clear();
    store.bases.clear();
    store.cigars.clear();
    store.qualities.clear();
    store.auxiliary.clear();
    store.mate_stats = MateLinkStats { .. };
}

// Results of Seqair-style mate linking. pairs counts reciprocal links;
// ambiguous_qnames counts eligible records without a free reciprocal partner
// when another eligible record of that name was already seen. It is a record
// count, not a count of distinct names. Inspect nonzero ambiguity explicitly.
//
// ```goose
// let store = pileup::Store { .. };
// let stats = store.link_mates();
// assert(stats.pairs == 0 && stats.ambiguous_qnames == 0);
// ```
struct MateLinkStats { pairs: i64 = 0, ambiguous_qnames: i64 = 0 }

// Link primary mapped pairs already in the store, using equal name bytes,
// same-reference mate fields and reciprocal positions (not READ1/READ2 bits).
// Secondary/supplementary, unmapped/mate-unmapped and empty-name records are
// excluded. No outside-region mate fetch or overlap correction is performed.
// Both load paths call this once after loading. Calling again is idempotent;
// reordering/removing records requires relinking before reading mate_index.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// assert(store.link_mates().pairs == 0 && store.records[0].mate_index == -1);
// ```
fn link_mates(store: Store&) -> MateLinkStats {
    for &rec in store.records { rec.mate_index = -1; }
    var waiting = dictionary<const u8[:], i64> {};
    waiting.reserve(store.records.len);
    var next: i64[>..] = [];
    next.resize(store.records.len, -1);
    var stats = MateLinkStats { .. };
    for i in store.records.len {
        let rec = store.records[i];
        if rec.flags & 1 == 0 || rec.flags & 0x90C != 0 || rec.tid < 0 ||
           rec.mate_tid != rec.tid || rec.mate_pos < 0 || rec.name_len == 0 { continue; }
        let name = store.names[(rec.name as i64)..rec.name as i64 + rec.name_len];
        waiting.update(name, -1) { head =>
            var candidate = head;
            var partner: i64 = -1;
            while candidate >= 0 {
                let other = store.records[candidate];
                if other.mate_index < 0 && other.pos == rec.mate_pos && other.mate_pos == rec.pos {
                    partner = candidate;
                    break;
                }
                candidate = next[candidate];
            }
            if partner >= 0 {
                store.records[i].mate_index = partner;
                store.records[partner].mate_index = i;
                stats.pairs += 1;
            } else {
                if head >= 0 { stats.ambiguous_qnames += 1; }
                next[i] = head;
                head = i;
            }
        };
    }
    store.mate_stats = stats;
    stats
}

// Scan from the reader's current position for [begin, end) on tid. Requires
// coordinate-sorted input and a freshly opened reader for independent queries.
// Returns Record on success, including an empty region. Error clears the store.
// This is not an indexed query: it reads preceding records and consumes the
// first mapped record beyond the region. Reopen before a different query.
// Unsorted detects decreasing mapped coordinates only in the scanned prefix;
// it cannot validate the unvisited remainder of a file. Empty regions do no I/O.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, reader.reference_id("chrA"), 5, 7) == bam::Status.Record);
// assert(store.records.len == 1);
// assert(store.load(reader, 0, 9, 8) == bam::Status.Invalid);
// assert(store.records.len == 0 && store.bases.len == 0);
// ```
fn load<R>(store: Store&, reader: R&, tid: i32, begin: i64, end: i64,
        filter: Filter = Filter { .. }) -> bam::Status {
    store.clear();
    store.tid = tid;
    store.begin = begin;
    store.end = end;
    if tid < 0 || tid >= reader.references.len || begin < 0 || end < begin ||
       end > reader.references[tid as i64].length {
        return bam::Status.Invalid;
    }
    if begin == end { return bam::Status.Record; }
    let status = scan_stream(store, reader, filter);
    if status != bam::Status.Record {
        store.clear();
    } else { store.link_mates(); }
    status
}

// Load [begin, end) via a matching BAI/CSI index and link retained mates.
// Reuses the open reader, index and slab capacity across independent queries,
// including backward queries and queries after EOF. No prefix scan or reopen.
// Like Seqair, rejects non-overlapping/filtered raw records before full decode.
// BGZF and record framing are always checked; names/auxiliary/query lengths
// of discarded records are not validated. Retained records are fully checked.
// Invalidates the reader's record view; call read before borrowing it again.
// Failures clear the store. Sorting is checked only in visited chunks; indexes
// cannot prove that the BAM is sorted, complete, or unchanged since indexing.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let index = bam_index::Index { .. };
// assert(index.open("tests/data/example.bam.bai") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load_indexed(reader, index, 0, 95, 100) == bam::Status.Record);
// assert(store.load_indexed(reader, index, 0, 5, 7) == bam::Status.Record);
// {
//     let record = store.record(0);
//     assert(store.records.len == 1 && record.name() == "match");
// }
// assert(store.load_indexed(reader, index, 0, 0, 30,
//     pileup::Filter { retain_zero_span: true, .. }) == bam::Status.Record);
// var empty = 0;
// for record in store.records { if record.reference_len == 0 { empty += 1; } }
// assert(empty == 3); // Empty CIGAR, insertion-only, soft-clip-only.
// store.columns() { pos, entries =>
//     for entry in entries { assert(store.records[entry.record_index].reference_len > 0); }
// };
// ```
fn load_indexed(store: Store&, reader: bam::Reader&, index: const bam_index::Index&,
                tid: i32, begin: i64, end: i64, filter: Filter = Filter { .. }) -> bam::Status {
    store.clear();
    store.tid = tid; store.begin = begin; store.end = end;
    if reader.stream.input.fd < 0 { return bam::Status.Closed; }
    if index.status != bam::Status.Record { return index.status; }
    if index.references.len != reader.references.len || tid < 0 || tid >= reader.references.len ||
       begin < 0 || end < begin || end > reader.references[tid as i64].length { return bam::Status.Invalid; }
    var chunks: bam_index::Chunk[>..] = [];
    var status = index.query(tid, begin, end, chunks);
    if status == bam::Status.Record && chunks.len > 0 { status = scan(store, reader, filter, chunks); }
    if status != bam::Status.Record { store.clear(); }
    else { store.link_mates(); }
    status
}

// Load an indexed region and visit its covered positions. The callback gets
// a zero-based position and entries borrowed for that call; store.record()
// accesses contributing reads. Uses load_indexed's filters, mate linking and
// storage options. A failed load clears the store and invokes no callbacks.
// Returns Record on success, including empty regions. Do not reload or mutate
// the store inside the callback. Repeated queries reuse its storage.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// reader.open("tests/data/example.bam").check(bam::Status.Record);
// let index = bam_index::Index { .. };
// index.open_for("tests/data/example.bam").check(bam::Status.Record);
// let store = pileup::Store { .. };
// var count = 0;
// let status = store.pileup(reader, index, 0 as i32, 5, 7) { pos, entries =>
//     assert(pos == 5 + count && entries.len == 1);
//     let record = store.record(entries[0].record_index);
//     assert(record.name() == "match");
//     count += 1;
// };
// status.check(bam::Status.Record);
// assert(count == 2);
// ```
fn pileup<F: fn(i64, const Entry[:])>(store: Store&, reader: bam::Reader&,
        index: const bam_index::Index&, tid: i32, begin: i64, end: i64,
        filter: Filter = Filter { .. }) -> bam::Status {
    let status = store.load_indexed(reader, index, tid, begin, end, filter);
    if status == bam::Status.Record { store.columns(F); }
    status
}

// Indexed pileup by reference name. An unknown name returns Invalid and clears
// the store without invoking the callback; otherwise identical to the tid form.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// reader.open("tests/data/example.bam").check(bam::Status.Record);
// let index = bam_index::Index { .. };
// index.open_for("tests/data/example.bam").check(bam::Status.Record);
// let store = pileup::Store { .. };
// let status = store.pileup(reader, index, "chrA", 5, 7) { pos, entries =>
//     assert(entries.len == 1);
// };
// status.check(bam::Status.Record);
// ```
fn pileup<F: fn(i64, const Entry[:])>(store: Store&, reader: bam::Reader&,
        index: const bam_index::Index&, name: const u8[:], begin: i64, end: i64,
        filter: Filter = Filter { .. }) -> bam::Status {
    store.pileup(reader, index, reader.reference_id(name), begin, end, filter, F)
}

// Sequential pileup without an index. Uses load's sorted-input contract and
// scans forward from the current position; reopen before an independent query.
// Like the indexed form, callbacks run only after a successful complete load.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// reader.open("tests/data/example.bam").check(bam::Status.Record);
// let store = pileup::Store { .. };
// let status = store.pileup(reader, 0 as i32, 5, 7) { pos, entries =>
//     assert(entries.len == 1);
// };
// status.check(bam::Status.Record);
// ```
fn pileup<R, F: fn(i64, const Entry[:])>(store: Store&, reader: R&,
        tid: i32, begin: i64, end: i64, filter: Filter = Filter { .. }) -> bam::Status {
    let status = store.load(reader, tid, begin, end, filter);
    if status == bam::Status.Record { store.columns(F); }
    status
}

// Sequential pileup by reference name, with the same behavior as the tid form.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// reader.open("tests/data/example.bam").check(bam::Status.Record);
// let store = pileup::Store { .. };
// let status = store.pileup(reader, "chrA", 5, 7) { pos, entries =>
//     assert(entries.len == 1);
// };
// status.check(bam::Status.Record);
// ```
fn pileup<R, F: fn(i64, const Entry[:])>(store: Store&, reader: R&,
        name: const u8[:], begin: i64, end: i64, filter: Filter = Filter { .. }) -> bam::Status {
    store.pileup(reader, reader.reference_id(name), begin, end, filter, F)
}

fn scan_stream<R>(store: Store&, reader: R&, filter: Filter) -> bam::Status {
    var last_tid: i32 = -1;
    var last_pos: i64 = -1;
    loop {
        let status = reader.read();
        if status == bam::Status.End { return bam::Status.Record; }
        if status != bam::Status.Record { return status; }
        let record = reader.record();
        let f = record.fields;
        if f.flags & 4 != 0 || f.tid < 0 || f.pos < 0 { continue; }
        if f.tid < last_tid || (f.tid == last_tid && f.pos < last_pos) { return bam::Status.Unsorted; }
        last_tid = f.tid; last_pos = f.pos;
        if f.tid > store.tid || (f.tid == store.tid && f.pos >= store.end) { return bam::Status.Record; }
        if f.tid != store.tid || f.end <= store.begin || f.flags & filter.exclude_flags != 0 ||
           f.mapq < filter.min_mapq || (f.reference_len == 0 && !filter.retain_zero_span) { continue; }
        let kept = retain(store, record);
        if kept != bam::Status.Record { return kept; }
    }
}

fn scan(store: Store&, reader: bam::Reader&, filter: Filter,
        chunks: const bam_index::Chunk[:] = []) -> bam::Status {
    var last_tid: i32 = -1;
    var last_pos: i64 = -1;
    var chunk:    i64 = 0;
    if chunks.len > 0 {
        let status = reader.seek(chunks[0].begin);
        if status != bam::Status.Record { return if status == bam::Status.End { bam::Status.Truncated } else { status }; }
    }
    while true {
        if chunks.len > 0 {
            while reader.tell() >= chunks[chunk].end {
                chunk += 1;
                if chunk == chunks.len { return bam::Status.Record; }
                if reader.tell() < chunks[chunk].begin {
                    let status = reader.seek(chunks[chunk].begin);
                    if status != bam::Status.Record { return if status == bam::Status.End { bam::Status.Truncated } else { status }; }
                }
            }
        }
        let status = if chunks.len > 0 { reader.read_payload() } else { reader.read() };
        if status == bam::Status.End { return if chunks.len > 0 { bam::Status.Truncated } else { bam::Status.Record }; }
        if status != bam::Status.Record { return status; }
        if chunks.len > 0 {
            let raw = reader.record();
            let tid = bam::signed32(raw.bytes, 0);
            let pos = bam::signed32(raw.bytes, 4) as i64;
            let flags = bgzf::le16(raw.bytes, 14) as u16;
            if flags & 4 != 0 || tid < 0 || pos < 0 { continue; }
            if tid < last_tid || (tid == last_tid && pos < last_pos) { return bam::Status.Unsorted; }
            last_tid = tid;
            last_pos = pos;
            if tid > store.tid || (tid == store.tid && pos >= store.end) { break; }
            if tid != store.tid || flags & filter.exclude_flags != 0 ||
               raw.bytes[9] < filter.min_mapq || !bam::may_reach(raw.bytes, store.begin) { continue; }
            let decoded = reader.decode_record();
            if decoded != bam::Status.Record { return decoded; }
        }
        let record = reader.record();
        let f = record.fields;
        if f.flags & 4 != 0 || f.tid < 0 || f.pos < 0 { continue; }
        if f.tid < last_tid || (f.tid == last_tid && f.pos < last_pos) { return bam::Status.Unsorted; }
        last_tid = f.tid;
        last_pos = f.pos;
        if f.tid > store.tid || (f.tid == store.tid && f.pos >= store.end) { break; }
        if f.tid != store.tid || f.end <= store.begin || f.flags & filter.exclude_flags != 0 ||
           f.mapq < filter.min_mapq || (f.reference_len == 0 && !filter.retain_zero_span) { continue; }
        let kept = retain(store, record);
        if kept != bam::Status.Record { return kept; }
    }
    bam::Status.Record
}

fn retain(store: Store&, record: const bam::Record&) -> bam::Status {
    let f = record.fields;
    let name = record.name();
    let aux = record.aux();
    let used = store.names.len + store.bases.len + store.cigars.len * 4 +
               store.qualities.len + store.auxiliary.len;
    let added = name.len + 2 * f.seq_len + 4 * f.cigar_len + aux.len;
    if added > min(store.max_bytes, 4294967295) - used { return bam::Status.TooLarge; }
    let item = SlimRecord {
        tid:          f.tid, pos: f.pos, end: f.end, reference_len: f.reference_len,
        flags:        f.flags, mapq: f.mapq, mate_tid: f.mate_tid, mate_pos: f.mate_pos,
        template_len: f.template_len,
        name:         store.names.len as u32, name_len: name.len as u32,
        cigar:        store.cigars.len as u32, cigar_len: f.cigar_len as u32,
        seq:          store.bases.len as u32, seq_len: f.seq_len as u32,
        qual:         store.qualities.len as u32, aux: store.auxiliary.len as u32, aux_len: aux.len as u32,
        ..
    };
    store.names.append(name);
    store.qualities.append(record.bytes[f.qual..f.aux]);
    store.auxiliary.append(aux);
    store.cigars.resize(store.cigars.len + f.cigar_len, 0);
    let cigars = store.cigars[(item.cigar as i64)..];
    for i in f.cigar_len { cigars[i] = record.cigar(i); }
    // Expand BAM's high nibble first, directly into the bases slab.
    // The trailing low nibble of an odd-length read is padding, not a base.
    store.bases.resize(store.bases.len + f.seq_len, 0);
    unpack_bases(record.bytes[f.seq..f.qual], store.bases[(item.seq as i64)..]);
    store.records.push(item);
    bam::Status.Record
}

extern fn unpack_bases_raw(sequence: const u8[:], bases: u8[:]) = "goose_seqair_unpack_bases";

fn unpack_bases(sequence: const u8[:], bases: u8[:]) {
    assert(sequence.len == (bases.len + 1) / 2);
    unpack_bases_raw(sequence, bases);
}

// One read's contribution to a column. Deletions and reference skips are
// distinct entries with no query base; both count in raw depth.
//
// ```goose
// assert(pileup::Op.Deletion != pileup::Op.RefSkip);
// ```
enum Op { Base, Deletion, RefSkip }

// record_index indexes Store.records. qpos is -1 for deletion/skip or absent
// sequence; base is 'N' and quality 255 when unavailable. indel is anchored at
// this column: positive for inserted bases after it, negative for a deletion.
// Sequence remains in BAM's stored orientation, including reverse-strand reads.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 6) == bam::Status.Record);
// store.columns() { pos, entries =>
//     assert(pos == 5 && entries.len == 1);
//     let e = entries[0];
//     assert(e.qpos == 0 && e.base == 'A' && e.quality == 13);
//     assert(e.op == pileup::Op.Base && e.indel == 0);
// };
// ```
struct Entry {
    record_index: i64, qpos: i64, base: u8, quality: u8,
    op:           Op, indel: i64,
}

struct Cursor {
    record_index:  i64, index: i64, reference: i64, query: i64,
    reference_end: i64, code: u32,
    end:           i64, cigar_end: i64, seq: i64, qual: i64, seq_len: i64,
}

inline fn consumes_reference(code: u32) -> bool { code == 0 || code == 2 || code == 3 || code == 7 || code == 8 }
inline fn consumes_query(code: u32) -> bool { code == 0 || code == 1 || code == 4 || code == 7 || code == 8 }

// Decode only when crossing a CIGAR boundary, not at every covered base.
fn advance(cigars: const u32[:], cursor: Cursor&, pos: i64) {
    while cursor.index < cursor.cigar_end {
        let word = cigars[cursor.index];
        let code = word & 15;
        let len = (word >> 4) as i64;
        let ref_op = consumes_reference(code);
        if ref_op && pos < cursor.reference + len {
            cursor.reference_end = cursor.reference + len;
            cursor.code = code;
            return;
        }
        if ref_op { cursor.reference += len; }
        if consumes_query(code) { cursor.query += len; }
        cursor.index += 1;
    }
}

fn contribution(cigars: const u32[:], bases: const u8[:], qualities: const u8[:],
                cursor: Cursor&, pos: i64) -> Entry {
    if pos >= cursor.reference_end { advance(cigars, cursor, pos); }
    var entry = Entry { record_index: cursor.record_index, qpos: -1, base: 'N', quality: 255,
                        op: Op.Base, indel: 0 };
    if cursor.code == 2 { entry.op = Op.Deletion; }
    else if cursor.code == 3 { entry.op = Op.RefSkip; }
    else {
        let q = cursor.query + pos - cursor.reference;
        if q < cursor.seq_len {
            entry.qpos = q;
            entry.base = "=ACMGRSVTWYHKDBN"[bases[cursor.seq + q] as i64];
            entry.quality = qualities[cursor.qual + q];
        }
    }
    if pos == cursor.reference_end - 1 {
        var next = cursor.index + 1;
        var insertion = false;
        while next < cursor.cigar_end {
            let following = cigars[next];
            let n = (following >> 4) as i64;
            let op = following & 15;
            next += 1;
            if op == 6 { continue; }
            if op == 1 { entry.indel += n; insertion = true; continue; }
            if op == 2 && !insertion { entry.indel = -n; }
            break;
        }
    }
    entry
}

// Visit covered positions in [store.begin, store.end), in order. The callback
// receives a zero-based position and a borrowed slice valid for that call only.
// Column entries preserve input order. Scratch storage is reused; CIGAR cursors
// only move forward. No reference FASTA, BAQ, mate dedup or depth cap is applied.
// Only reference-consuming CIGAR operations produce entries; clips/insertions
// alone never create columns. Empty positions are omitted, not zero-filled.
//
// ```goose
// let reader = bam::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/data/example.bam") == bam::Status.Record);
// let store = pileup::Store { .. };
// assert(store.load(reader, 0, 5, 7) == bam::Status.Record);
// var count = 0;
// store.columns() { pos, entries =>
//     assert(pos == 5 + count && entries.len == 1);
//     count += 1;
// };
// assert(count == 2);
// ```
fn columns<F: fn(i64, const Entry[:])>(store: const Store&) {
    var active:  Cursor[>..] = [];
    var entries: Entry[>..]  = [];
    var next:    i64         = 0;
    var pos = store.begin;
    var expire = store.end;
    while pos < store.end {
        // Between read ends the active set cannot change; avoid a second
        // per-entry traversal just to copy every cursor back onto itself.
        if pos >= expire {
            var kept: i64 = 0;
            expire = store.end;
            for cursor in active {
                let end = cursor.end;
                if end > pos {
                    active[kept] = cursor;
                    kept += 1;
                    expire = min(expire, end);
                }
            }
            active.resize(kept);
        }
        if active.len == 0 && next < store.records.len {
            pos = max(pos, store.records[next].pos);
        }
        if pos >= store.end { break; }
        while next < store.records.len && store.records[next].pos <= pos {
            let fields = store.records[next];
            if fields.end > pos && fields.reference_len > 0 {
                active.push(Cursor {
                    record_index:  next, index: fields.cigar, reference: fields.pos, query: 0,
                    reference_end: fields.pos, code: 0, end: fields.end,
                    cigar_end:     fields.cigar as i64 + fields.cigar_len,
                    seq:           fields.seq, qual: fields.qual, seq_len: fields.seq_len,
                });
                expire = min(expire, fields.end);
            }
            next += 1;
        }
        if active.len == 0 && next == store.records.len { break; }
        entries.clear();
        for &cursor in active {
            entries.push(contribution(store.cigars, store.bases, store.qualities, cursor, pos));
        }
        if entries.len > 0 { F(pos, entries[..]); }
        pos += 1;
    }
}
