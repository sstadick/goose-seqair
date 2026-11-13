import std;
import os;
import seqair from "..";
import seqair.bam;
import seqair.bam_index;
import io from "github.com/sstadick/goose-io@v0.3.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";

fn main() {
    let args = args();
    if args.len != 2 && args.len != 5 && args.len != 6 {
        abort("usage: records input.bam [contig begin end [index.bai|index.csi]] (0-based, half-open)");
    }
    let reader = bam::Reader { .. };
    defer reader.close();
    reader.open(args[1]).check(bam::Status.Ok, "opening BAM");
    var line: u8[>..] = [];
    fn emit(record: const bam::Record&) {
        line.clear();
        format(line, record.name(), "\t", record.fields.tid, "\t", record.fields.pos,
               "\t", record.fields.end, "\t", record.fields.seq_len, "\n");
        write_stdout(line);
    }
    if args.len == 2 {
        reader.records(emit).check(bam::Status.Ok, "reading BAM");
    } else {
        let begin, b = parse_int(args[3]); let end, e = parse_int(args[4]);
        if !b || !e { abort("invalid region coordinates"); }
        let index = bam_index::Index { .. };
        let opened = if args.len == 6 { index.open(args[5]) } else { index.open_for(args[1]) };
        opened.check(bam::Status.Ok, "opening BAM index");
        index.fetch(reader, args[2], begin, end, emit).check(bam::Status.Ok, "fetching BAM region");
    }
}
