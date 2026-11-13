import std;
import os;
import seqair from "..";
import seqair.bam;
import seqair.bam_index;
import seqair.pileup;
import seqair.bam_writer;
import seqair.bam_encode;
import seqair.sam;
import seqair.sam_writer;
import seqair.cram;
import seqair.cram_index;
import seqair.cram_writer;
import io from "github.com/sstadick/goose-io@v0.3.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import simd_bytes.simd_int;
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";

fn status_name(status: bam::Status) -> const u8[:] {
    match status {
        Ok => "Ok", Record => "Record", End => "End", Invalid => "Invalid", Truncated => "Truncated",
        IoError => "IoError", Closed => "Closed", TooLarge => "TooLarge", Unsorted => "Unsorted",
    }
}

fn print_record<T>(record: const T&, line: u8[>..]&) {
    let f = record.fields;
    line.clear();
    format(line, "R\t", record.name(), "\t", f.tid, "\t", f.pos, "\t", f.end, "\t",
           f.flags, "\t", f.mapq, "\t", f.mate_tid, "\t", f.mate_pos, "\t", f.template_len, "\t");
    for i in f.cigar_len { format(line, record.cigar(i), ","); }
    line.push('\t');
    for i in f.seq_len { line.push(record.base(i)); }
    line.push('\t');
    for i in f.seq_len { format(line, record.quality(i) as i64, ","); }
    line.push('\t');
    for byte in record.aux() { format(line, byte as i64, ","); }
    line.push('\n');
    write_stdout(line);
}

// Deliberately simple in-memory FASTA provider. Tests can cap fetch windows
// independently of the CRAM reader's limits to catch whole-contig fetching.
root struct Reference { data: u8[>..] = [], max_window: i64 = 1000000 }
fn fetch(source: Reference&, name: const u8[:], begin: i64, end: i64, output: u8[>..]&) -> bool {
    output.clear();
    if begin < 0 || end < begin || end - begin > source.max_window { return false; }
    var selected = false; var position = 0;
    bstr::split_str(source.data, "\n") { line =>
        if line.len > 0 && line[0] == '>' {
            selected = line[1..] == name; position = 0;
        } else if selected {
            for base in line {
                if base == '\r' { continue; }
                if position >= begin && position < end { output.push(base); }
                position += 1;
            }
        }
    };
    output.len == end - begin
}

fn check_return(actual: bam::Status, expected: bam::Status, cleaned: bool&) -> i64, bool {
    defer cleaned = true;
    let status = actual.check(expected) { error =>
        assert(error == actual);
        return 17, false;
    };
    assert(status == actual);
    return 23, true;
}

fn check_named(actual: bam::Status) { assert(actual == bam::Status.Closed); }

fn first_record(reader: bam::Reader&) -> u8[], bool {
    reader.records() { record => return str(record.name()), true; };
    return "", false;
}

fn main() {
    let args = args();
    assert(args.len >= 3);
    if args[1] == "check-abort" {
        bam::Status.End.check(bam::Status.Record); return;
    }
    if args[1] == "check-abort-ok" {
        bam::Status.Record.check(bam::Status.Ok); return;
    }
    if args[1] == "check-abort-record" {
        bam::Status.Ok.check(bam::Status.Record); return;
    }
    if args[1] == "check-abort-message" {
        bam::Status.IoError.check(bam::Status.Ok, "Failed to open BAM file"); return;
    }
    if args[1] == "check-api" {
        let statuses = [bam::Status.Ok, bam::Status.Record, bam::Status.End, bam::Status.Invalid, bam::Status.Truncated,
            bam::Status.IoError, bam::Status.Closed, bam::Status.TooLarge, bam::Status.Unsorted];
        for actual in statuses {
            assert(actual.check(actual) == actual);
            assert(actual.check(actual, "matched") == actual);
            for expected in statuses {
                var called = 0;
                let status = actual.check(expected) { error =>
                    assert(error == actual); called += 1;
                };
                assert(status == actual && called == (if actual == expected { 0 } else { 1 }));
                var cleaned = false;
                let n, ok = check_return(actual, expected, cleaned);
                assert(cleaned && ok == (actual == expected) && n == (if ok { 23 } else { 17 }));
            }
        }
        assert(bam::Status.Closed.check(bam::Status.Record, check_named) == bam::Status.Closed);
        print("S\tOk"); return;
    }
    if args[1] == "records" || args[1] == "records-first" {
        let reader = bam::Reader { .. }; defer reader.close();
        var status = reader.open(args[2]);
        var line: u8[>..] = [];
        if status == bam::Status.Ok {
            assert(reader.status == bam::Status.Ok);
            if args[1] == "records-first" {
                let name, found = first_record(reader);
                assert(found && name == "match" && reader.status == bam::Status.Record);
            }
            status = reader.records() { record =>
                assert(reader.status == bam::Status.Record);
                print_record(record, line);
            };
            assert(reader.status == (if status == bam::Status.Ok { bam::Status.End } else { status }));
            let repeated = reader.records() { record => assert(false); };
            assert(repeated == status);
        }
        print("S\t", status_name(status)); return;
    }
    if args[1] == "open-for" {
        let index = bam_index::Index { .. };
        if args.len > 3 { index.open(args[3]).check(bam::Status.Ok); }
        let status = index.open_for(args[2]);
        if status != bam::Status.Ok { assert(index.references.len == 0); }
        print("S\t", status_name(status)); return;
    }
    if args[1] == "fetch" || args[1] == "fetch-tid" || args[1] == "fetch-tail" || args[1] == "columns" ||
       args[1] == "columns-tid" || args[1] == "columns-scan" || args[1] == "columns-scan-tid" {
        let reader = bam::Reader { .. }; defer reader.close();
        let index = bam_index::Index { .. };
        let store = pileup::Store { .. };
        let scan = args[1] == "columns-scan" || args[1] == "columns-scan-tid";
        let by_tid = args[1] == "fetch-tid" || args[1] == "columns-tid" || args[1] == "columns-scan-tid";
        var status = reader.open(args[2]);
        if status == bam::Status.Ok && !scan {
            status = if args[3] == "-" { index.open_for(args[2]) } else { index.open(args[3]) };
        }
        if status != bam::Status.Ok { print("S\t", status_name(status)); return; }
        // Exercise fetch/seek recovery after reaching EOF on this same reader.
        reader.records() { record => };
        var line: u8[>..] = [];
        fn column(pos: i64, entries: const pileup::Entry[:]) {
            for &entry in entries {
                let record = store.record(entry.record_index);
                let op = match entry.op { Base => 0, Deletion => 1, RefSkip => 2 };
                line.clear();
                format(line, pos, "\t", record.name(), "\t", entry.qpos, "\t", op, "\t",
                       entry.indel, "\t", entry.base as i64, "\t", entry.quality as i64, "\n");
                write_stdout(line);
            }
        }
        var at: i64 = 4;
        while at < args.len {
            let begin, b = parse_int(args[at + 1]); let end, e = parse_int(args[at + 2]); assert(b && e);
            let tid = reader.reference_id(args[at]);
            print("Q\t", args[at], "\t", begin, "\t", end);
            if scan { reader.open(args[2]).check(bam::Status.Ok); }
            if !scan && tid >= 0 && begin >= 0 && end >= begin && end <= reader.references[tid as i64].length {
                var chunks: bam_index::Chunk[>..] = [];
                assert(index.query(tid, begin, end, chunks) == bam::Status.Ok);
                if begin == end { assert(chunks.len == 0); }
            }
            if args[1] == "fetch" || args[1] == "fetch-tail" {
                status = index.fetch(reader, args[at], begin, end) { record => print_record(record, line); };
            } else if args[1] == "fetch-tid" {
                status = index.fetch(reader, tid, begin, end) { record => print_record(record, line); };
            } else if scan {
                status = if by_tid { store.pileup(reader, tid, begin, end, column) }
                    else { store.pileup(reader, args[at], begin, end, column) };
            } else {
                status = if by_tid { store.pileup(reader, index, tid, begin, end, column) }
                    else { store.pileup(reader, index, args[at], begin, end, column) };
            }
            if status != bam::Status.Ok { assert(store.records.len == 0); }
            if args[1] == "fetch-tail" {
                assert(status == bam::Status.Ok && reader.status == bam::Status.Record);
                assert(reader.read() == bam::Status.End && reader.status == bam::Status.End);
                assert(reader.read() == bam::Status.End);
            }
            print("S\t", status_name(status));
            at += 3;
        }
        reader.close();
        let closed = index.fetch(reader, 0, 0, 1) { record => assert(false); };
        assert(closed == bam::Status.Closed);
        return;
    }
    if args[1] == "cram-index" {
        let index = cram_index::Index { .. };
        if args.len > 3 { let limit, ok = parse_int(args[3]); assert(ok); index.max_bytes = limit; }
        print("S\t", status_name(index.open(args[2]))); return;
    }
    if args[1] == "cram-indexed" {
        let source = Reference { .. }; assert(read_file(args[3], source.data));
        let reader = cram::Reader { .. }; defer reader.close();
        let index = cram_index::Index { .. }; let cursor = cram_index::Cursor { .. };
        var status = reader.open(args[2]);
        if status == bam::Status.Ok {
            assert(reader.status == bam::Status.Ok);
            status = index.open(args[4]);
        }
        if status != bam::Status.Ok { print("S\t", status_name(status)); return; }
        var line: u8[>..] = []; var at: i64 = 5;
        while at < args.len {
            let begin, b = parse_int(args[at + 1]); let end, e = parse_int(args[at + 2]); assert(b && e);
            print("Q\t", args[at], "\t", begin, "\t", end);
            status = cursor.fetch(index, reader, reader.reference_id(args[at]), begin, end);
            assert(status != bam::Status.Record && status != bam::Status.End);
            if status == bam::Status.Ok {
                while true {
                    status = cursor.read(reader, source);
                    if status != bam::Status.Record { break; }
                    assert(reader.status == bam::Status.Record);
                    let record = reader.record(); print_record(record, line);
                }
                assert(cursor.read(reader, source) == status);
            }
            print("S\t", status_name(status)); at += 3;
        }
        return;
    }
    if args[1] == "to-cram" {
        let reader = bam::Reader { .. }; defer reader.close();
        assert(reader.open(args[2]) == bam::Status.Ok);
        let writer = cram_writer::Writer { .. }; defer writer.close();
        if args.len > 5 { let n, ok = parse_int(args[5]); assert(ok); writer.records_per_slice = n; }
        var level: i32 = 6;
        if args.len > 4 { let n, ok = parse_int(args[4]); assert(ok); level = n as i32; }
        var minor: i32 = 1;
        if args.len > 6 { let n, ok = parse_int(args[6]); assert(ok); minor = n as i32; }
        if args.len > 7 { let n, ok = parse_int(args[7]); assert(ok); writer.max_record_bytes = n; }
        if args.len > 8 { let n, ok = parse_int(args[8]); assert(ok); writer.max_slice_bytes = n; }
        var status = writer.open(args[3], reader.header, reader.references, reader.names, level, minor);
        while status == bam::Status.Ok && reader.read() == bam::Status.Record {
            let record = reader.record(); status = writer.write(record);
        }
        if status == bam::Status.Ok { status = writer.flush(); }
        let closed = writer.close();
        if status == bam::Status.Ok { status = if closed == bam::Status.Ok { reader.status } else { closed }; }
        print("S\t", status_name(status)); return;
    }
    if args[1] == "to-sam" {
        let reader = bam::Reader { .. }; defer reader.close();
        assert(reader.open(args[2]) == bam::Status.Ok);
        let writer = sam_writer::Writer { .. }; defer writer.close();
        if args.len > 4 { let n, ok = parse_int(args[4]); assert(ok); writer.max_line_bytes = n; }
        if args.len > 5 { let n, ok = parse_int(args[5]); assert(ok); writer.max_record_bytes = n; }
        var status = writer.open(args[3], reader.header, reader.references, reader.names);
        while status == bam::Status.Ok && reader.read() == bam::Status.Record {
            let record = reader.record(); status = writer.write(record);
        }
        let closed = writer.close();
        if status == bam::Status.Ok { status = if closed == bam::Status.Ok { reader.status } else { closed }; }
        print("S\t", status_name(status)); return;
    }
    if args[1] == "writer-api-sam" || args[1] == "writer-api-cram" {
        let reader = bam::Reader { .. }; defer reader.close();
        assert(reader.open(args[2]) == bam::Status.Ok);
        let store = pileup::Store { .. };
        assert(store.load(reader, 0, 0, 100, pileup::Filter { exclude_flags: 0, .. }) == bam::Status.Ok);
        if args[1] == "writer-api-sam" {
            let writer = sam_writer::Writer { .. };
            writer_api(writer, reader, store, args[3]);
        } else {
            let writer = cram_writer::Writer { records_per_slice: 2, .. };
            writer_api(writer, reader, store, args[3]);
        }
        print("S\tOk"); return;
    }
    if args[1] == "cram" || args[1] == "cram-write" || args[1] == "cram-pileup" || args[1] == "cram-limit" {
        let source = Reference { .. };
        assert(read_file(args[3], source.data));
        let reader = cram::Reader { .. }; defer reader.close();
        if args[1] == "cram-limit" { reader.max_container_bytes = 50; }
        var status = reader.open(args[2]);
        assert(status != bam::Status.Record && status != bam::Status.End);
        let writer = bam_writer::Writer { .. }; defer writer.close();
        if status == bam::Status.Ok {
            assert(reader.status == bam::Status.Ok);
            if args[1] == "cram-pileup" {
                let begin, b = parse_int(args[5]); let end, e = parse_int(args[6]); assert(b && e);
                let tid = reader.reference_id(args[4]);
                let store = pileup::Store { .. };
                status = cram::load(store, reader, source, tid, begin, end, pileup::Filter { exclude_flags: 0, .. });
                if status == bam::Status.Ok {
                    store.columns() { pos, entries =>
                        for &entry in entries {
                            let record = store.record(entry.record_index);
                            let op = match entry.op { Base => 0, Deletion => 1, RefSkip => 2 };
                            print(pos, "\t", record.name(), "\t", entry.qpos, "\t", op, "\t",
                                  entry.indel, "\t", entry.base as i64, "\t", entry.quality as i64);
                        }
                    };
                } else { assert(store.records.len == 0 && store.bases.len == 0); }
            } else {
                if args[1] == "cram-write" { assert(writer.open(args[4], reader.header, reader.references, reader.names) == bam::Status.Ok); }
                var line: u8[>..] = [];
                while true {
                    status = reader.read(source);
                    assert(reader.status == status);
                    if status != bam::Status.Record { break; }
                    let record = reader.record();
                    if args[1] == "cram-write" { assert(writer.write(record) == bam::Status.Ok); }
                    else { print_record(record, line); }
                }
                assert(reader.read(source) == status);
                if args[1] == "cram-write" { assert(writer.close() == bam::Status.Ok); }
            }
        }
        print("S\t", status_name(status)); return;
    }
    if args[1] == "sam-pileup" || args[1] == "sam-mates" {
        let reader = sam::Reader { .. };
        defer reader.close();
        let store = pileup::Store { .. };
        var status = reader.open(args[2]);
        let begin, b = parse_int(args[4]); let end, e = parse_int(args[5]);
        let mask, m = parse_int(args[6]); let mapq, q = parse_int(args[7]); assert(b && e && m && q);
        if status == bam::Status.Ok {
            status = store.load(reader, reader.reference_id(args[3]), begin, end,
                pileup::Filter { exclude_flags: mask as u16, min_mapq: mapq as u8 });
        }
        reader.close();
        if status == bam::Status.Ok {
            if args[1] == "sam-mates" {
                print("M\t", store.mate_stats.pairs, "\t", store.mate_stats.ambiguous_qnames);
                for i in store.records.len { print(i, "\t", store.records[i].mate_index); }
            } else {
                store.columns() { pos, entries =>
                    for &entry in entries {
                        let record = store.record(entry.record_index);
                        let op = match entry.op { Base => 0, Deletion => 1, RefSkip => 2 };
                        print(pos, "\t", record.name(), "\t", entry.qpos, "\t", op, "\t",
                              entry.indel, "\t", entry.base as i64, "\t", entry.quality as i64);
                    }
                };
            }
        } else { assert(store.records.len == 0 && store.bases.len == 0); }
        print("S\t", status_name(status)); return;
    }
    if args[1] == "sam" || args[1] == "sam-write" || args[1] == "sam-limit" {
        let reader = sam::Reader { .. };
        defer reader.close();
        if args[1] == "sam-limit" { reader.max_line_bytes = 100; }
        var status = reader.open(args[2]);
        assert(status != bam::Status.Record && status != bam::Status.End);
        let writer = bam_writer::Writer { .. };
        defer writer.close();
        var line: u8[>..] = [];
        if status == bam::Status.Ok {
            assert(reader.status == bam::Status.Ok);
            if args[1] == "sam-write" { assert(writer.open(args[3], reader.header, reader.references, reader.names) == bam::Status.Ok); }
            for ref, i in reader.references { print("H\t", i, "\t", reader.names[ref.name..ref.name + ref.name_len], "\t", ref.length); }
            while true {
                status = reader.read(); if status != bam::Status.Record { break; }
                assert(reader.status == bam::Status.Record);
                let record = reader.record();
                if args[1] == "sam-write" { assert(writer.write(record) == bam::Status.Ok); }
                else { print_record(record, line); }
            }
            if args[1] == "sam-write" { assert(writer.close() == bam::Status.Ok); }
            assert(reader.status == status);
            assert(reader.read() == status);
        }
        print("S\t", status_name(status)); return;
    }
    let reader = bam::Reader { .. };
    assert(reader.read() == bam::Status.Closed);
    defer reader.close();
    if args[1] == "limit" { reader.max_record_bytes = 64; reader.max_header_bytes = 64; }
    var status = reader.open(args[2]);
    assert(reader.status == status && status != bam::Status.Record && status != bam::Status.End);
    var line: u8[>..] = [];
    let lean = args[1] == "lean" || args[1] == "indexed-lean";
    let indexed = args[1] == "indexed" || args[1] == "indexed-store" || args[1] == "indexed-mates" || args[1] == "indexed-lean";
    let mates = args[1] == "mates" || args[1] == "indexed-mates";
    if args[1] == "write" || args[1] == "write-store" || args[1] == "write-new" || args[1] == "write-fail" {
        assert(status == bam::Status.Ok && args.len == 5);
        let level, ok = parse_int(args[4]); assert(ok);
        let writer = bam_writer::Writer { .. };
        defer writer.close();
        status = writer.open(args[3], reader.header, reader.references, reader.names, level as i32);
        if status == bam::Status.Ok {
            if args[1] == "write-new" {
                let cigar: u32[] = [3 << 4, 1 << 4 | 1, 2 << 4];
                let alignment = bam_encode::Alignment { name: "new", tid: 0, pos: 123,
                    flags: 16, mapq: 42, mate_tid: 1, mate_pos: 321, template_len: -99,
                    cigar: cigar, sequence: "=aCMnT", qualities: "\x01\x02\x03\x04\x05\x06",
                    aux: "NMi\x03\0\0\0", .. };
                var encoded: u8[>..] = [];
                assert(bam_encode::encode(alignment, encoded, reader.references.len) == bam::Status.Ok);
                status = writer.write(alignment);
                assert(status == bam::Status.Ok);
                let absent = bam_encode::Alignment { name: "absent", flags: 4, sequence: "NNN", .. };
                status = writer.write(absent);
            } else if args[1] == "write-fail" {
                let invalid = bam_encode::Alignment { name: "bad", flags: 4, sequence: "!", .. };
                assert(writer.write(invalid) == bam::Status.Invalid);
                let valid = bam_encode::Alignment { flags: 4, .. };
                assert(writer.write(valid) == bam::Status.Invalid);
                assert(writer.close() == bam::Status.Invalid);
                assert(writer.close() == bam::Status.Ok);
                status = bam::Status.Invalid;
            } else if args[1] == "write-store" {
                let store = pileup::Store { .. };
                status = store.load(reader, 0, 0, reader.references[0].length, pileup::Filter { exclude_flags: 0, .. });
                assert(status == bam::Status.Ok);
                for i in store.records.len {
                    let record = store.record(i);
                    status = writer.write(record);
                    if status != bam::Status.Ok { break; }
                }
            } else {
                while reader.read() == bam::Status.Record {
                    let record = reader.record();
                    status = writer.write(record);
                    if status != bam::Status.Ok { break; }
                }
                if status == bam::Status.Ok { assert(reader.status == bam::Status.End); }
            }
            if status == bam::Status.Ok { status = writer.close(); }
            else { writer.close(); }
        }
    } else if args[1] == "index" {
        let index = bam_index::Index { .. };
        if args.len == 4 { let limit, ok = parse_int(args[3]); assert(ok); index.max_bytes = limit; }
        status = index.open(args[2]);
        if status != bam::Status.Ok { assert(index.references.len == 0 && index.chunks.len == 0); }
    } else if args[1] == "seek" {
        assert(status == bam::Status.Ok && reader.status == bam::Status.Ok);
        var offsets: u64[>..] = [];
        while true {
            let offset = reader.tell();
            status = reader.read();
            if status == bam::Status.End { break; }
            assert(status == bam::Status.Record && reader.status == bam::Status.Record);
            offsets.push(offset);
        }
        assert(reader.read() == bam::Status.End);
        for n in offsets.len {
            assert(reader.seek(offsets[offsets.len - n - 1]) == bam::Status.Ok);
            assert(reader.status == bam::Status.Ok);
            assert(reader.read() == bam::Status.Record);
            assert(reader.status == bam::Status.Record);
            let record = reader.record();
            print_record(record, line);
            if n == 0 {
                // Seeking to the final alignment still delivers Record first.
                assert(reader.read() == bam::Status.End && reader.status == bam::Status.End);
                assert(reader.read() == bam::Status.End);
            }
        }
        assert(reader.seek(0) == bam::Status.Invalid);
        assert(reader.seek(offsets[0]) == bam::Status.Ok);
        assert(reader.status == bam::Status.Ok);
        status = reader.read();
    } else if args[1] == "pileup" || args[1] == "budget" || args[1] == "store" || indexed || mates || lean {
        assert(args.len == 7 || args.len == 8 || args.len == 9);
        let begin, a = parse_int(args[4]);
        let end, b = parse_int(args[5]);
        let mask, c = parse_int(args[6]);
        assert(a && b && c);
        let store = pileup::Store { .. };
        if lean {
            store.retain_auxiliary = false;
            store.auto_link_mates = false;
        }
        var min_mapq: i64 = 0;
        if args.len >= 8 {
            let value, ok = parse_int(args[7]);
            assert(ok);
            if args[1] == "store" { store.max_bytes = value; }
            else { min_mapq = value; }
        }
        if args[1] == "budget" { store.max_bytes = 150; }
        if indexed && status == bam::Status.Ok {
            let index = bam_index::Index { .. };
            status = index.open(args[8]);
            if status == bam::Status.Ok {
                status = store.load_indexed(reader, index, reader.reference_id(args[3]), begin, end,
                                    pileup::Filter { exclude_flags: mask as u16, min_mapq: min_mapq as u8 });
                if status == bam::Status.Ok {
                    // Query another reference and an empty window, then go back:
                    // all reuse the same file, index and store without reopening.
                    assert(store.load_indexed(reader, index, 1, 0, 1) == bam::Status.Ok);
                    assert(store.load_indexed(reader, index, 0, 0, 0) == bam::Status.Ok);
                    status = store.load_indexed(reader, index, reader.reference_id(args[3]), begin, end,
                                        pileup::Filter { exclude_flags: mask as u16, min_mapq: min_mapq as u8 });
                }
            }
        } else if status == bam::Status.Ok {
            status = store.load(reader, reader.reference_id(args[3]), begin, end,
                                pileup::Filter { exclude_flags: mask as u16, min_mapq: min_mapq as u8 });
        }
        if args[1] == "store" && status == bam::Status.Ok {
            // Reuse must replace, not append; access below must survive close.
            assert(reader.open(args[2]) == bam::Status.Ok);
            status = store.load(reader, reader.reference_id(args[3]), begin, end,
                                pileup::Filter { exclude_flags: mask as u16, .. });
        }
        reader.close();
        if status == bam::Status.Ok {
            if lean {
                assert(store.auxiliary.len == 0);
                assert(store.mate_stats.pairs == 0 && store.mate_stats.ambiguous_qnames == 0);
                for item in store.records {
                    assert(item.mate_index == -1);
                }
            }
            if mates {
                let stats = store.mate_stats;
                let again = store.link_mates();
                assert(stats.pairs == again.pairs && stats.ambiguous_qnames == again.ambiguous_qnames);
                print("M\t", stats.pairs, "\t", stats.ambiguous_qnames);
                for i in store.records.len { print(i, "\t", store.records[i].mate_index); }
            } else if args[1] == "store" || args[1] == "indexed-store" || lean {
                for i in store.records.len {
                    let record = store.record(i);
                    print_record(record, line);
                }
            } else {
                store.columns() { pos, entries =>
                    for &entry in entries {
                        let record = store.record(entry.record_index);
                        let op = match entry.op { Base => 0, Deletion => 1, RefSkip => 2 };
                        line.clear();
                        format(line, pos, "\t", record.name(), "\t", entry.qpos, "\t", op, "\t",
                               entry.indel, "\t", entry.base as i64, "\t", entry.quality as i64, "\n");
                        write_stdout(line);
                    }
                };
            }
        } else {
            assert(store.records.len == 0 && store.names.len == 0 && store.bases.len == 0 &&
                   store.cigars.len == 0 && store.qualities.len == 0 && store.auxiliary.len == 0);
        }
    } else {
        if status == bam::Status.Ok {
            for ref, i in reader.references {
                line.clear();
                format(line, "H\t", i, "\t", reader.names[ref.name..ref.name + ref.name_len], "\t", ref.length, "\n");
                write_stdout(line);
            }
            status = reader.read();
            while status == bam::Status.Record {
                assert(reader.status == bam::Status.Record);
                let record = reader.record();
                print_record(record, line);
                status = reader.read();
            }
        }
        assert(reader.status == status);
        assert(reader.read() == status);
    }
    line.clear();
    format(line, "S\t", status_name(status), "\n");
    write_stdout(line);
    reader.close();
    reader.close();
    assert(reader.read() == bam::Status.Closed);
}

fn writer_api<W>(writer: W&, reader: bam::Reader&, store: const pileup::Store&, path: const u8[:]) {
    let cigar: u32[] = [3 << 4];
    let alignment = bam_encode::Alignment { name: "owned", tid: 0, pos: 3, cigar: cigar,
        sequence: "aCn", qualities: "\x07\x16\x28", aux: "ZZZvalue\x00", .. };
    assert(writer.write(alignment) == bam::Status.Closed);
    assert(writer.open(path, reader.header, reader.references, reader.names) == bam::Status.Ok);
    assert(writer.write(alignment) == bam::Status.Ok);
    for i in store.records.len { let record = store.record(i); assert(writer.write(record) == bam::Status.Ok); }
    assert(writer.close() == bam::Status.Ok && writer.close() == bam::Status.Ok);
    assert(writer.write(alignment) == bam::Status.Closed);
    // The next output must not inherit records, tag dictionaries or counters.
    var other: u8[>..] = []; format(other, path, ".reopened");
    assert(writer.open(other, reader.header, reader.references, reader.names) == bam::Status.Ok);
    assert(writer.write(alignment) == bam::Status.Ok);
    assert(writer.close() == bam::Status.Ok);
    format(other, ".error");
    assert(writer.open(other, reader.header, reader.references, reader.names) == bam::Status.Ok);
    let invalid = bam_encode::Alignment { tid: 999, .. };
    assert(writer.write(invalid) == bam::Status.Invalid);
    assert(writer.write(alignment) == bam::Status.Invalid);
    assert(writer.close() == bam::Status.Invalid && writer.close() == bam::Status.Ok);
}
