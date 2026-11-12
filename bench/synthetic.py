#!/usr/bin/env python3
"""Encode Seqair's DNA/spliced benchmark read shapes as BAM for matched I/O tests.

Reproduces synthetic_input's seed, RNG draws, starts, CIGARs, bases, qualities,
names and flags at Seqair 8401961d (crates/seqair/benches/bam.rs). Unlike that
in-memory benchmark, both consumers here read the same BAM and index. Requires
the pysam version in tests/requirements.txt. Generates 200,000 reads on chrA.

Adapted from Seqair, copyright (c) 2026 Pascal Hertleif, under the MIT license;
see ../LICENSE.
"""
import argparse
from pathlib import Path

import pysam


def generate(shape, output):
    state = 0x9E3779B97F4A7C15

    def below(n):
        nonlocal state
        state ^= (state << 13) & ((1 << 64) - 1)
        state ^= state >> 7
        state ^= (state << 17) & ((1 << 64) - 1)
        return state % n

    starts = sorted(below(1000000) for _ in range(200000))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrA", "LN": 1010001}, {"SN": "chrB", "LN": 400}]}
    output.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(output), "wb", header=header) as writer:
        for i, start in enumerate(starts):
            if shape == "dna":
                kind = below(50)
                if kind == 0:
                    cigar = [(0, 70), (1, 1 + below(3)), (0, 78)]
                elif kind == 1:
                    cigar = [(0, 70), (2, 1 + below(3)), (0, 80)]
                else:
                    cigar = [(0, 150)]
            elif below(3) == 0:
                cigar = [(0, 40), (3, 100 + below(4900)), (0, 70),
                         (3, 100 + below(4900)), (0, 40)]
            else:
                a = 20 + below(110)
                cigar = [(0, a), (3, 100 + below(4900)), (0, 150 - a)]
            length = sum(n for op, n in cigar if op in (0, 1))
            read = pysam.AlignedSegment(writer.header)
            read.query_name = f"r{i}"
            read.reference_id = 0
            read.reference_start = start
            read.flag = 0 if i % 2 == 0 else 16
            read.mapping_quality = 60
            read.cigartuples = cigar
            read.query_sequence = "".join("ACGT"[below(4)] for _ in range(length))
            read.query_qualities = [20 + below(20) for _ in range(length)]
            writer.write(read)
    pysam.index(str(output))
    pysam.index("-c", str(output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shape", choices=("dna", "spliced"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    generate(args.shape, args.output)
