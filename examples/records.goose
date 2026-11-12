import std;
import os;
import seqair from "..";
import seqair.bam;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import deflate from "github.com/sstadick/goose-deflate@8fcf5832a3fa680f4c6a5332e2436f8a6dd5c299";

fn main() {
    let args = args();
    if args.len != 2 { abort("usage: records input.bam"); }
    let reader = bam::Reader { .. };
    defer reader.close();
    if reader.open(args[1]) != bam::Status.Record { abort("cannot open BAM"); }
    var line: u8[>..] = [];
    var status = reader.read();
    while status == bam::Status.Record {
        let record = reader.record();
        line.clear();
        format(line, record.name(), "\t", record.fields.tid, "\t", record.fields.pos,
               "\t", record.fields.end, "\t", record.fields.seq_len, "\n");
        write_stdout(line);
        status = reader.read();
    }
    if status != bam::Status.End { abort("invalid or truncated BAM"); }
}
