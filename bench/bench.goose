import std;
import os;
import seqair from "..";
import seqair.bam;
import seqair.bam_index;
import seqair.pileup;
import seqair.bam_writer;
import seqair.sam;
import seqair.sam_writer;
import seqair.cram_writer;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.0";

fn main() {
    let args = args();
    if args.len == 4 && args[1] == "read-sam" {
        let iterations, ok = parse_int(args[3]); assert(ok && iterations > 0);
        let reader = sam::Reader { .. };
        defer reader.close();
        var elapsed: i64 = 0;
        var records: i64 = 0;
        var bases: i64 = 0;
        for iteration in iterations + 1 {
            let start = clock_ns();
            assert(reader.open(args[2]) == bam::Status.Record);
            records = 0; bases = 0;
            while reader.read() == bam::Status.Record {
                records += 1; bases += reader.record().fields.seq_len;
            }
            assert(reader.status == bam::Status.End);
            if iteration > 0 { elapsed += clock_ns() - start; }
        }
        print(elapsed, " ", records, " ", bases); return;
    }
    if args.len == 6 && (args[1] == "write" || args[1] == "write-sam" || args[1] == "write-cram") {
        write_benchmark(args[2], args[3], args[4], args[5], args[1]);
        return;
    }
    if args.len != 6 && args.len != 7 { abort("usage: bench input.bam contig begin end iterations [index]"); }
    let begin, b = parse_int(args[3]);
    let end, e = parse_int(args[4]);
    let iterations, n = parse_int(args[5]);
    assert(b && e && n && iterations > 0);
    let reader = bam::Reader { .. };
    defer reader.close();
    let store = pileup::Store { .. };
    let index = bam_index::Index { .. };
    var read_ns:          i64 = 0;
    var load_ns:          i64 = 0;
    var columns_ns:       i64 = 0;
    var expected_records: i64 = 0;
    var expected_entries: i64 = 0;
    var expected_sum:     i64 = 0;
    // Iteration zero warms the filesystem cache and reusable buffers.
    for iteration in iterations + 1 {
        let start = clock_ns();
        assert(reader.open(args[1]) == bam::Status.Record);
        var records: i64 = 0;
        var sum:     i64 = 0;
        var status = reader.read();
        while status == bam::Status.Record {
            let record = reader.record();
            records += 1;
            sum += record.fields.seq_len;
            status = reader.read();
        }
        assert(status == bam::Status.End);
        let read_end = clock_ns();
        assert(reader.open(args[1]) == bam::Status.Record);
        if args.len == 7 {
            // Charge index parsing to load, matching Seqair's indexed adapter.
            assert(index.open(args[6]) == bam::Status.Record);
            assert(store.load_indexed(reader, index, reader.reference_id(args[2]), begin, end) == bam::Status.Record);
        } else {
            assert(store.load(reader, reader.reference_id(args[2]), begin, end) == bam::Status.Record);
        }
        let load_end = clock_ns();
        var entries_count: i64 = 0;
        store.columns() { pos, entries =>
            entries_count += entries.len;
            for &entry in entries {
                let op = match entry.op { Base => 0, Deletion => 1, RefSkip => 2 };
                sum += pos + entry.record_index + entry.qpos + entry.base + entry.quality + entry.indel + op;
            }
        };
        let columns_end = clock_ns();
        if iteration == 0 {
            expected_records = records;
            expected_entries = entries_count;
            expected_sum = sum;
        } else {
            assert(records == expected_records && entries_count == expected_entries && sum == expected_sum);
            read_ns += read_end - start;
            load_ns += load_end - read_end;
            columns_ns += columns_end - load_end;
        }
    }
    print(read_ns, " ", load_ns, " ", columns_ns, " ", expected_records, " ",
          store.records.len, " ", expected_entries, " ", expected_sum);
}

fn emit<W>(writer: W&, fields: const bam::Fields[:], payloads: const u8[:], offsets: const i64[:]) {
    for f, i in fields {
        let record = bam::Record { fields: f, bytes: payloads[offsets[i]..offsets[i + 1]] };
        assert(writer.write(record) == bam::Status.Record);
    }
    assert(writer.close() == bam::Status.Record);
}

fn write_benchmark(input: const u8[:], output: const u8[:], level_arg: const u8[:], count_arg: const u8[:], mode: const u8[:]) {
    let level, l = parse_int(level_arg);
    let iterations, n = parse_int(count_arg);
    assert(l && n && iterations > 0);
    let reader = bam::Reader { .. };
    defer reader.close();
    assert(reader.open(input) == bam::Status.Record);
    var payloads: u8[>..] = [];
    var offsets: i64[>..] = [0];
    var fields: bam::Fields[>..] = [];
    while reader.read() == bam::Status.Record {
        let record = reader.record();
        payloads.append(record.bytes);
        offsets.push(payloads.len);
        fields.push(record.fields);
    }
    assert(reader.status == bam::Status.End);
    var elapsed: i64 = 0;
    for iteration in iterations + 1 {
        let start = clock_ns();
        if mode == "write-sam" {
            let writer = sam_writer::Writer { .. };
            assert(writer.open(output, reader.header, reader.references, reader.names) == bam::Status.Record);
            emit(writer, fields, payloads, offsets);
        } else if mode == "write-cram" {
            let writer = cram_writer::Writer { .. };
            assert(writer.open(output, reader.header, reader.references, reader.names, level as i32) == bam::Status.Record);
            emit(writer, fields, payloads, offsets);
        } else {
            let writer = bam_writer::Writer { .. };
            assert(writer.open(output, reader.header, reader.references, reader.names, level as i32) == bam::Status.Record);
            emit(writer, fields, payloads, offsets);
        }
        if iteration > 0 { elapsed += clock_ns() - start; }
    }
    print(elapsed, " ", fields.len, " ", payloads.len);
}
