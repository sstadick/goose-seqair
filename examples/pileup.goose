import std;
import os;
import seqair from "..";
import seqair.bam;
import seqair.bam_index;
import seqair.pileup;
import io from "github.com/sstadick/goose-io@v0.3.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";

fn main() {
    let args = args();
    if args.len != 5 && args.len != 6 { abort("usage: pileup input.bam contig begin end [index.bai|index.csi] (0-based, half-open)"); }
    let begin, start_ok = parse_int(args[3]);
    let end, end_ok = parse_int(args[4]);
    assert(start_ok && end_ok);
    let reader = bam::Reader { .. };
    defer reader.close();
    reader.open(args[1]).check(bam::Status.Ok, "opening BAM");
    let store = pileup::Store { .. };
    var line: u8[>..] = [];
    fn column(pos: i64, entries: const pileup::Entry[:]) {
        if store.mate_stats.ambiguous_qnames != 0 { abort("ambiguous template names in region"); }
        line.clear();
        format(line, pos, "\t", entries.len, "\n");
        write_stdout(line);
    }
    if args.len == 6 {
        let index = bam_index::Index { .. };
        index.open(args[5]).check(bam::Status.Ok, "opening BAM index");
        store.pileup(reader, index, args[2], begin, end, column).check(bam::Status.Ok, "loading pileup");
    } else {
        store.pileup(reader, args[2], begin, end, column).check(bam::Status.Ok, "loading pileup");
    }
    if store.mate_stats.ambiguous_qnames != 0 { abort("ambiguous template names in region"); }
}
