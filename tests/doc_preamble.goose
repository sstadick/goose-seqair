import std;
import os;
import seqair from "..";
import seqair.bam;
import seqair.bam_index;
import seqair.pileup;
import seqair.bam_encode;
import seqair.bam_writer;
import seqair.sam;
import seqair.cram;
import seqair.cram_index;
import seqair.sam_writer;
import seqair.cram_writer;
import io from "github.com/sstadick/goose-io@v0.3.0";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import simd_bytes.simd_int;
import deflate from "github.com/sstadick/goose-deflate@v0.1.0";
import bgzf from "github.com/sstadick/goose-bgzf@v0.1.2";

// A provider that refuses external references; sufficient for embedded/no-ref
// CRAM. Applications can adapt faidx::Reader as shown in README.md.
struct ExampleReference {}
fn fetch(source: ExampleReference&, name: const u8[:], begin: i64, end: i64, output: u8[>..]&) -> bool { false }
