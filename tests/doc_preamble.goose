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
import seqair.sam_writer;
import seqair.cram_writer;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.3.0";
import simd_bytes.simd_int;
import deflate from "github.com/sstadick/goose-deflate@8fcf5832a3fa680f4c6a5332e2436f8a6dd5c299";

// A provider that refuses external references; sufficient for embedded/no-ref
// CRAM. Applications can adapt faidx::Reader as shown in README.md.
struct ExampleReference {}
fn fetch(source: ExampleReference&, name: const u8[:], begin: i64, end: i64, output: u8[>..]&) -> bool { false }
