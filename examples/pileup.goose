import std;
import os;
import seqair from "..";
import seqair.bam;
import seqair.bam_index;
import seqair.pileup;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.0";

fn main() {
    let args = args();
    if args.len != 5 && args.len != 6 { abort("usage: pileup input.bam contig begin end [index.bai|index.csi] (0-based, half-open)"); }
    let begin, start_ok = parse_int(args[3]);
    let end, end_ok = parse_int(args[4]);
    assert(start_ok && end_ok);
    let reader = bam::Reader { .. };
    defer reader.close();
    assert(reader.open(args[1]) == bam::Status.Record);
    let store = pileup::Store { .. };
    if args.len == 6 {
        let index = bam_index::Index { .. };
        assert(index.open(args[5]) == bam::Status.Record);
        assert(store.load_indexed(reader, index, reader.reference_id(args[2]), begin, end) == bam::Status.Record);
    } else {
        assert(store.load(reader, reader.reference_id(args[2]), begin, end) == bam::Status.Record);
    }
    if store.mate_stats.ambiguous_qnames != 0 { abort("ambiguous template names in region"); }
    var line: u8[>..] = [];
    store.columns() { pos, entries =>
        line.clear();
        format(line, pos, "\t", entries.len, "\n");
        write_stdout(line);
    };
}
