#!/usr/bin/env python3
"""Fresh-reader CRAM scans; Clang -O3 -march=native, alternating matched consumers."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cram", type=Path)
    p.add_argument("reference", type=Path)
    p.add_argument("--seqair-bin", type=Path, required=True)
    p.add_argument("--rounds", type=int, default=5)
    p.add_argument("--iterations", type=int, default=10)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.rounds < 1 or args.iterations < 1: p.error("counts must be positive")
    env = dict(os.environ, CC="clang", GOOSE_CC="clang")
    goose = env.get("GOOSE", "goose")
    flags = ["-O3", "-march=native"]
    report = {"native_flags": flags, "compiler": subprocess.check_output(["clang", "--version"], text=True),
              "goose": subprocess.check_output([goose, "--version"], text=True),
              "platform": platform.platform(), "rounds": args.rounds, "iterations": args.iterations,
              "affinity": sorted(os.sched_getaffinity(0)),
              "external_codec_library_path": env.get("LD_LIBRARY_PATH", "system libraries (build provenance unknown)"),
              "cram_sha256": hashlib.sha256(args.cram.read_bytes()).hexdigest(),
              "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest()}
    with tempfile.TemporaryDirectory(prefix="goose-cram-bench-") as temp:
        directory = Path(temp); source = directory / "bench.c"
        subprocess.run([goose, "-o", str(source), "--standalone", "-O2", "bench/cram.goose"], cwd=ROOT, env=env, check=True)
        headers = [(source.parent / name).resolve() for name in re.findall(r'^#include "([^"]+)"', source.read_text(), re.M)]
        codec, = [h.parent.parent for h in headers if h.name == "deflate.h"]
        io, = [h.parent.parent for h in headers if h.name == "goosio.h"]
        bio, = [h.parent.parent for h in headers if h.name == "faidx.h"]
        build = directory / "libdeflate"
        subprocess.run(["cmake", "-S", str(codec / "vendor/libdeflate"), "-B", str(build),
                        "-DCMAKE_C_COMPILER=clang", "-DCMAKE_BUILD_TYPE=Release",
                        "-DCMAKE_C_FLAGS_RELEASE=-O3 -march=native -DNDEBUG",
                        "-DLIBDEFLATE_BUILD_SHARED_LIB=OFF", "-DLIBDEFLATE_BUILD_GZIP=OFF",
                        "-DLIBDEFLATE_BUILD_TESTS=OFF"], env=env, check=True)
        subprocess.run(["cmake", "--build", str(build), "--target", "libdeflate_static"], env=env, check=True)
        archive, = build.rglob("libdeflate.a")
        binary = directory / "cram"
        command = ["clang", *flags, str(source), *map(str, [ROOT / "src/cram.c", ROOT / "src/bases.c", ROOT / "src/seek.c",
                   io / "src/goosio.c", bio / "src/faidx.c", archive]), "-lm", "-pthread", "-ldl", "-o", str(binary)]
        subprocess.run(command, env=env, check=True)
        report["build_command"] = command
        programs = {"goose": binary, "seqair": args.seqair_bin.resolve()}
        samples = {name: [] for name in programs}; counts = None
        for round in range(args.rounds):
            for name in list(programs) if round % 2 == 0 else reversed(programs):
                result = subprocess.check_output([str(programs[name]), str(args.cram.resolve()), str(args.reference.resolve()),
                                                  str(args.iterations)], env=env, text=True)
                elapsed, records, bases = map(int, result.split())
                if counts is None: counts = [records, bases]
                assert [records, bases] == counts, (name, records, bases, counts)
                samples[name].append(elapsed / args.iterations)
        report["counts"] = counts; report["samples_ns"] = samples
        report["median_ms"] = {name: statistics.median(values) / 1e6 for name, values in samples.items()}
        print(report["median_ms"])
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__": main()
