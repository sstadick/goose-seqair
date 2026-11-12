#!/usr/bin/env python3
"""Warm-cache baseline; builds timed C with Clang -O3 -march=native, including the codec."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam", type=Path)
    parser.add_argument("contig")
    parser.add_argument("begin", type=int)
    parser.add_argument("end", type=int)
    parser.add_argument("--cc", nargs="+", choices=["clang"], default=["clang"])
    parser.add_argument("--rounds", type=int, default=7)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--seqair", action="store_true", help="also build/run pinned Rust Seqair, sequential and indexed")
    parser.add_argument("--writer", action="store_true", help="also compare in-memory raw-record BAM writing at levels 0, 1, 6, 12")
    parser.add_argument("--format-writers", action="store_true", help="SAM vs Seqair and reference-free CRAM writing; requires pysam for output verification")
    parser.add_argument("--sam", type=Path, help="also compare sequential reading of this plain SAM")
    parser.add_argument("--index", type=Path, help="also measure Goose indexed loading; inferred beside BAM with --seqair")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.rounds < 1 or args.iterations < 1 or args.begin < 0 or args.end < args.begin:
        parser.error("invalid iteration count or region")
    if args.seqair and args.begin == args.end:
        parser.error("Seqair comparison requires a nonempty region")
    if args.seqair and args.index is None:
        args.index = next((p for p in (Path(str(args.bam) + ".bai"), args.bam.with_suffix(".bai"),
                                      Path(str(args.bam) + ".csi"), args.bam.with_suffix(".csi")) if p.exists()), None)
        if args.index is None:
            parser.error("--seqair needs a BAM index beside the input")
    goose = os.environ.get("GOOSE", "goose")
    with args.bam.open("rb") as file:
        digest = hashlib.file_digest(file, "sha256").hexdigest()
    report = {"platform": platform.platform(), "input_sha256": digest,
              "compressed_bytes": args.bam.stat().st_size,
              "affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
              "region": [args.contig, args.begin, args.end],
              "goose": subprocess.check_output([goose, "--version"], text=True).strip(),
              "goose_inlining": "-O2", "native_flags": ["-O3", "-march=native"],
              "codec_flags": ["-O3", "-march=native", "-DNDEBUG"],
              "rounds": args.rounds, "iterations": args.iterations, "results": []}
    programs = []
    if args.seqair:
        target = Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "bench/seqair/target")).resolve()
        env = dict(os.environ, CFLAGS="-O3 -march=native", CC="clang",
                   RUSTFLAGS="-C opt-level=3 -C target-cpu=native -C linker=clang",
                   CC_ENABLE_DEBUG_OUTPUT="1", CARGO_TARGET_DIR=str(target))
        command = ["cargo", "+1.98.1", "build", "--locked", "--release", "--manifest-path",
                   str(ROOT / "bench/seqair/Cargo.toml"), "-vv"]
        subprocess.run(command, env=env, check=True)
        for mode in ("sequential", "indexed"):
            row = {"implementation": f"seqair-{mode}",
                   "compiler": subprocess.check_output(["rustc", "+1.98.1", "--version"], text=True).strip(),
                   "build_command": command, "rustflags": env["RUSTFLAGS"], "cflags": env["CFLAGS"],
                   "cc": env["CC"],
                   "revision": "8401961de45883d4772226d91da95c9965104d59",
                   "samples_ns": []}
            extra = []
            if mode == "indexed":
                extra = ["indexed", str(args.index.resolve())]
                row["index_sha256"] = hashlib.sha256(args.index.read_bytes()).hexdigest()
            programs.append((target / "release/seqair-comparison", extra, row))
    with tempfile.TemporaryDirectory(prefix="goose-seqair-bench-") as temp:
        directory = Path(temp)
        source = directory / "bench.c"
        subprocess.run([goose, "-o", str(source), "--standalone", "-O2", "bench/bench.goose"],
                       cwd=ROOT, check=True)
        # Use the exact dependency checkouts selected by Goose, not sibling repos.
        headers = [(source.parent / p).resolve() for p in re.findall(r'^#include "([^"]+)"', source.read_text(), re.M)]
        codec, = [p.parent.parent for p in headers if p.name == "deflate.h"]
        io, = [p.parent.parent for p in headers if p.name == "goosio.h"]
        for index, cc in enumerate(args.cc):
            build = directory / f"codec-{index}"
            subprocess.run(["cmake", "-S", str(codec / "vendor/libdeflate"), "-B", str(build),
                            "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_C_COMPILER={shutil.which(cc) or cc}",
                            "-DCMAKE_C_FLAGS_RELEASE=-O3 -march=native -DNDEBUG", "-DLIBDEFLATE_BUILD_SHARED_LIB=OFF",
                            "-DLIBDEFLATE_BUILD_GZIP=OFF", "-DLIBDEFLATE_BUILD_TESTS=OFF"], check=True)
            subprocess.run(["cmake", "--build", str(build), "--target", "libdeflate_static", "--verbose"], check=True)
            archive, = build.rglob("libdeflate.a")
            exe = directory / f"bench-{index}"
            libraries = ["-lm", "-pthread"] + (["-ldl"] if platform.system() == "Linux" else [])
            command = [cc, "-O3", "-march=native", str(source), str(ROOT / "src/seek.c"), str(ROOT / "src/bases.c"), str(ROOT / "src/sam.c"),
                       str(io / "src/goosio.c"), str(archive), "-o", str(exe), *libraries]
            print("Native build:", command, flush=True)
            subprocess.run(command, check=True)
            row = {"implementation": f"goose-{cc}", "cc": cc,
                   "compiler": subprocess.check_output([cc, "--version"], text=True).splitlines()[0],
                   "build_command": command, "samples_ns": []}
            programs.append((exe, [], row))
            if args.index is not None:
                indexed = {**row, "implementation": f"goose-{cc}-indexed", "samples_ns": [],
                           "index_sha256": hashlib.sha256(args.index.read_bytes()).hexdigest()}
                programs.append((exe, [str(args.index.resolve())], indexed))
        reference_counts = None
        # Build everything first; alternate order to reduce thermal/drift bias.
        for round_ in range(args.rounds):
            for exe, extra, row in programs if round_ % 2 == 0 else reversed(programs):
                output = subprocess.check_output([str(exe), str(args.bam.resolve()), args.contig,
                                                  str(args.begin), str(args.end), str(args.iterations), *extra], text=True)
                values = list(map(int, output.split()))
                if reference_counts is None:
                    reference_counts = values[3:]
                assert values[3:] == reference_counts, (row["implementation"], "counts/checksum mismatch", values[3:], reference_counts)
                row["samples_ns"].append([v / args.iterations for v in values[:3]])
        for _, _, row in programs:
            samples = row["samples_ns"]
            medians = [statistics.median(s[i] for s in samples) for i in range(3)]
            row["median_ns"] = dict(zip(["read", "load", "columns"], medians))
            row["median_ns"]["region_total"] = statistics.median(s[1] + s[2] for s in samples)
            row["counts"] = dict(zip(["read_records", "retained_records", "entries", "checksum"], reference_counts))
            report["results"].append(row)
            print(f"{row['implementation']}: read {medians[0]/1e6:.3f} ms; load {medians[1]/1e6:.3f} ms; "
                  f"columns {medians[2]/1e6:.3f} ms; counts {reference_counts}", flush=True)
        if args.sam:
            report["sam_sha256"] = hashlib.sha256(args.sam.read_bytes()).hexdigest()
            sam_programs = [(exe, {"implementation": row["implementation"], "samples_ns": []})
                            for exe, extra, row in programs if not extra]
            counts = None
            for round_ in range(args.rounds):
                for exe, row in sam_programs if round_ % 2 == 0 else reversed(sam_programs):
                    result = subprocess.check_output([str(exe), "read-sam", str(args.sam.resolve()), str(args.iterations)], text=True)
                    elapsed, records, bases = map(int, result.split())
                    if counts is None: counts = [records, bases]
                    assert [records, bases] == counts
                    row["samples_ns"].append(elapsed / args.iterations)
            report["sam"] = []
            for _, row in sam_programs:
                row["median_ns"] = statistics.median(row["samples_ns"])
                row["records"], row["bases"] = counts
                report["sam"].append(row)
                print(f"{row['implementation']} SAM read: {row['median_ns']/1e6:.3f} ms; {counts}", flush=True)
        if args.writer:
            expected_payload = gzip.decompress(args.bam.read_bytes())
            report["writer"] = []
            writers = [(exe, row) for exe, extra, row in programs if not extra]
            for level in (0, 1, 6, 12):
                rows = [{"implementation": row["implementation"], "level": level, "samples_ns": []}
                        for _, row in writers]
                for round_ in range(args.rounds):
                    order = range(len(writers)) if round_ % 2 == 0 else reversed(range(len(writers)))
                    for i in order:
                        exe, _ = writers[i]
                        output = directory / f"writer-{i}.bam"
                        result = subprocess.check_output([str(exe), "write", str(args.bam.resolve()),
                                                          str(output), str(level), str(args.iterations)], text=True)
                        elapsed, records, payload_bytes = map(int, result.split())
                        assert records == reference_counts[0]
                        assert gzip.decompress(output.read_bytes()) == expected_payload
                        rows[i]["samples_ns"].append(elapsed / args.iterations)
                        rows[i]["records"] = records
                        rows[i]["payload_bytes"] = payload_bytes
                        rows[i]["output_bytes"] = output.stat().st_size
                for row in rows:
                    row["median_ns"] = statistics.median(row["samples_ns"])
                    print(f"{row['implementation']} write level {level}: {row['median_ns']/1e6:.3f} ms; "
                          f"{row['output_bytes']} bytes", flush=True)
                report["writer"].extend(rows)
        if args.format_writers:
            import pysam
            def normalized(path, normalize=False):
                with pysam.AlignmentFile(path) as reader:
                    result = []
                    for record in reader:
                        if normalize:
                            cigar = []
                            for op, count in record.cigartuples or []:
                                if not count: continue
                                op = 0 if op in (7, 8) else op
                                if cigar and cigar[-1][0] == op: cigar[-1] = (op, cigar[-1][1] + count)
                                else: cigar.append((op, count))
                            record.cigartuples = cigar
                        result.append(record.to_string())
                    return result
            expected = normalized(args.bam, True)
            report["format_writers"] = []
            sam_bytes = None
            for mode, level in [("sam", 0), ("cram", 0), ("cram", 1), ("cram", 6), ("cram", 12)]:
                writers = [(exe, {"implementation": row["implementation"], "format": mode, "level": level, "samples_ns": []})
                           for exe, extra, row in programs if not extra and
                           (mode != "cram" or not row["implementation"].startswith("seqair"))]
                output = directory / ("output." + mode)
                for round_ in range(args.rounds):
                    for exe, row in writers if round_ % 2 == 0 else reversed(writers):
                        result = subprocess.check_output([str(exe), "write-" + mode, str(args.bam.resolve()), str(output),
                                                          str(level), str(args.iterations)], text=True)
                        elapsed, records, payload_bytes = map(int, result.split())
                        assert records == len(expected)
                        row["samples_ns"].append(elapsed / args.iterations)
                        if mode == "sam":
                            data = output.read_bytes()
                            if sam_bytes is None: sam_bytes = data
                            assert data == sam_bytes, "SAM differs from matched formatter"
                        else:
                            assert normalized(output) == expected, "CRAM round-trip differs"
                        row.update(records=records, payload_bytes=payload_bytes, output_bytes=output.stat().st_size)
                for _, row in writers:
                    row["median_ns"] = statistics.median(row["samples_ns"])
                    print(f"{row['implementation']} {mode} level {level}: {row['median_ns']/1e6:.3f} ms; {row['output_bytes']} bytes", flush=True)
                    report["format_writers"].append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
