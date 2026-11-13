import std;
import os;
import seqair from "..";
import seqair.cram;
import seqair.bam;
import bio from "github.com/sstadick/goose-bio@v0.4.1";
import bio.faidx;
import io from "github.com/sstadick/goose-io@v0.3.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import simd_bytes.simd_int;
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";

struct Reference { reader: faidx::Reader& }
fn fetch(source: Reference&, name: const u8[:], begin: i64, end: i64, output: u8[>..]&) -> bool {
    source.reader.fetch(name, begin, end, output) == faidx::Status.Ok
}

fn main() {
    let args = args();
    let iterations, ok = parse_int(args[3]); assert(ok && iterations > 0);
    var elapsed: i64 = 0; var records: i64 = 0; var bases: i64 = 0;
    for iteration in iterations + 1 {
        let start = clock_ns();
        let fasta = faidx::Reader { .. }; defer fasta.close();
        assert(fasta.open(args[2]) == faidx::Status.Ok);
        let source = Reference { reader: fasta };
        let reader = cram::Reader { .. }; defer reader.close();
        assert(reader.open(args[1]) == bam::Status.Ok);
        records = 0; bases = 0;
        while reader.read(source) == bam::Status.Record { records += 1; bases += reader.fields.seq_len; }
        assert(reader.status == bam::Status.End);
        if iteration > 0 { elapsed += clock_ns() - start; }
    }
    print(elapsed, " ", records, " ", bases);
}
