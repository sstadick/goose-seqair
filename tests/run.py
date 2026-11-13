#!/usr/bin/env python3
"""BAM wire-format tests and per-read pileup comparisons with htslib via pysam."""
import argparse
import gzip
import os
from pathlib import Path
import struct
import subprocess
import tempfile

import pysam
from fixtures import REFS, bgzf, block, header, random_records, record, records

ROOT = Path(__file__).resolve().parents[1]


def oracle(path, begin, end, mask, min_mapq=0, contig="chrA"):
    result = []
    with pysam.AlignmentFile(path) as reader:
        for col in reader.pileup(contig, begin, end, truncate=True, stepper="samtools",
                                 flag_filter=mask, min_base_quality=0, min_mapping_quality=min_mapq,
                                 ignore_overlaps=False, ignore_orphans=False, compute_baq=False,
                                 max_depth=1000000):
            for p in col.pileups:
                r = p.alignment
                # Our qpos indexes stored SEQ, rather than hypothetical absent bases.
                q = p.query_position if r.query_sequence else None
                op = 2 if p.is_refskip else 1 if p.is_del else 0
                base = ord(r.query_sequence[q]) if q is not None and r.query_sequence else ord("N")
                qual = r.query_qualities[q] if q is not None and r.query_qualities is not None else 255
                result.append(f"{col.reference_pos}\t{r.query_name}\t{q if q is not None else -1}\t{op}\t{p.indel}\t{base}\t{qual}")
    return result


def suite(command, directory):
    def invoke(path, mode="read", *args):
        output = subprocess.check_output([*map(str, command), mode, str(path), *map(str, args)],
                                         cwd=ROOT, text=True)
        return output.splitlines()

    def write(name, data):
        path = directory / name
        Path(str(path) + ".bai").unlink(missing_ok=True)
        Path(str(path) + ".csi").unlink(missing_ok=True)
        path.write_bytes(data)
        return path

    expected_header = [f"H\t{i}\t{n}\t{l}" for i, (n, l) in enumerate(REFS)]
    fixtures = records()
    data = header() + b"".join(r for r, _ in fixtures)
    expected = expected_header + [e for _, e in fixtures] + ["S\tEnd"]
    assert invoke("unused", "check-api") == ["S\tRecord"]
    for mode, message in (("check-abort", "BAM status mismatch: expected Record, got End"),
                          ("check-abort-message", "Failed to open BAM file: expected Record, got IoError")):
        result = subprocess.run([*map(str, command), mode, "unused"], cwd=ROOT, text=True, capture_output=True)
        assert result.returncode != 0 and message in result.stdout + result.stderr, result
    for chunk in (1, 17, 65536):
        path = write(f"split-{chunk}.bam", bgzf(data, chunk, empty=True, extra=b"XY\x03\x00abc"))
        assert invoke(path) == expected, f"BGZF split {chunk}"
        assert invoke(path, "seek") == [e for _, e in reversed(fixtures)] + ["S\tRecord"]
        assert invoke(path, "records") == [e for _, e in fixtures] + ["S\tRecord"]
        assert invoke(path, "records-first") == [e for _, e in fixtures[1:]] + ["S\tRecord"]
    # Switch between borrowed block records and cross-block scratch, including
    # a split length word, a body ending exactly at the boundary, and a body
    # missing only its final byte. Seeking backward must make the same switch.
    for offset in (0, 1, 3, 4, 5, len(fixtures[0][0]) - 1, len(fixtures[0][0])):
        cut = len(header()) + offset
        path = write("boundary.bam", block(data[:cut]) + block(data[cut:]) + block(b""))
        assert invoke(path) == expected, ("record boundary", offset)
        assert invoke(path, "seek") == [e for _, e in reversed(fixtures)] + ["S\tRecord"]
    path = write("normal.bam", bgzf(data))
    assert invoke(path) == expected
    assert invoke(write("no-marker.bam", bgzf(data)[:-28])) == expected
    assert invoke(write("empty.bam", bgzf(header()))) == expected_header + ["S\tEnd"]
    assert invoke(directory / "missing.bam") == ["S\tIoError"]
    assert invoke(directory / "missing.bam", "records") == ["S\tIoError"]
    assert invoke(directory / "empty.bam", "records") == ["S\tRecord"]
    assert invoke(write("records-truncated.bam", bgzf(header() + fixtures[0][0] + fixtures[1][0][:-1])),
                  "records") == [fixtures[0][1], "S\tTruncated"]
    assert invoke(path, "limit")[-1] == "S\tTooLarge"
    assert invoke(path, "budget", "chrA", 0, 100, 0)[-1] == "S\tTooLarge"
    assert invoke(path, "pileup", "chrA", 20, 20, 0) == ["S\tRecord"]
    assert invoke(path, "pileup", "missing", 0, 10, 0) == ["S\tInvalid"]
    assert invoke(path, "pileup", "chrA", -1, 10, 0) == ["S\tInvalid"]
    assert invoke(path, "pileup", "chrA", 0, 200001, 0) == ["S\tInvalid"]

    long, long_expected = record("long", 1, [("M", 1)] * 65536)
    large = write("long.bam", bgzf(header() + long, 65536))
    assert invoke(large) == expected_header + [long_expected, "S\tEnd"]
    assert invoke(large, "seek") == [long_expected, "S\tRecord"]
    with pysam.AlignmentFile(large) as f:
        assert len(next(f).cigartuples) == 65536

    # Writer input/output are checked independently with zlib and htslib.
    eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    def check_bgzf(path):
        wire = path.read_bytes()
        assert wire.endswith(eof)
        at = 0
        while at < len(wire):
            size = struct.unpack_from("<H", wire, at + 16)[0] + 1
            assert 28 <= size <= 65536
            payload = gzip.decompress(wire[at:at + size])  # Checks CRC/ISIZE too.
            assert len(payload) <= 65280
            at += size
        assert at == len(wire)

    # Incompressible auxiliary bytes plus records larger than one block catch
    # a 64 KiB input buffer, missing tails, and mistaken record/block alignment.
    import random
    noise = random.Random(71).randbytes(200000)
    noisy = record("noise", 2, [("M", 17)], aux=b"ZZBC" + struct.pack("<i", len(noise)) + noise)
    sources = [write("writer-input.bam", bgzf(data)), large,
               write("noisy.bam", bgzf(header() + noisy[0])), write("header-only.bam", bgzf(header()))]
    output = directory / "written.bam"
    for source in sources:
        for level in (0, 1, 6, 12):
            assert invoke(source, "write", output, level) == ["S\tRecord"]
            check_bgzf(output)
            assert gzip.decompress(output.read_bytes()) == gzip.decompress(source.read_bytes())
            with pysam.AlignmentFile(source) as left, pysam.AlignmentFile(output) as right:
                assert [r.to_string() for r in left] == [r.to_string() for r in right]
    for source in (large, sources[2]):
        assert invoke(source, "write-store", output, 6) == ["S\tRecord"]
        check_bgzf(output)
        assert invoke(output) == invoke(source)
        pysam.index(str(output))
        with pysam.AlignmentFile(output) as f: assert len(list(f.fetch("chrA"))) == 1
    assert invoke(sources[0], "write-new", output, 0) == ["S\tRecord"]
    check_bgzf(output)
    with pysam.AlignmentFile(output) as f:
        a, b = list(f)
        assert (a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality,
                a.cigarstring, a.next_reference_id, a.next_reference_start, a.template_length,
                a.query_sequence, list(a.query_qualities), a.get_tag("NM")) == (
                    "new", 16, 0, 123, 42, "3M1I2M", 1, 321, -99, "=ACMNT", [1,2,3,4,5,6], 3)
        assert b.is_unmapped and b.reference_id == -1 and b.query_sequence == "NNN" and b.query_qualities is None
    assert invoke(sources[0], "write-fail", output, 6) == ["S\tInvalid"]
    assert not output.read_bytes().endswith(eof)
    output.write_bytes(b"untouched")
    assert invoke(sources[0], "write", output, 13) == ["S\tInvalid"]
    assert output.read_bytes() == b"untouched"
    if Path("/dev/full").exists():
        assert invoke(sources[0], "write", "/dev/full", 6) == ["S\tIoError"]

    # SAM has its own syntax/fixups but must produce the same BAM records as
    # htslib. Compare every field and the raw auxiliary encodings, not counts.
    sam_header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrA\tLN:2147483647\n@SQ\tSN:chrB\tLN:1000\n"
    def sam_line(**changes):
        fields = dict(name="r", flag="16", ref="chrA", pos="17", mapq="37", cigar="3M1I2M",
                      mate="=", mpos="31", tlen="-14", seq="aCMnTg", qual="!~ABCD")
        fields.update(changes)
        return "\t".join(fields.values()).encode()
    aux = (b"\tAA:A:Z\tNI:i:-129\tPI:i:65536\tUI:i:4294967295\tNZ:i:-0"
           b"\tFF:f:-1.25\tDD:d:3.125\tZZ:Z:spaces stay\tHH:H:0aFF"
           b"\tBC:B:c,-128,127\tBS:B:S,0,65535\tBI:B:I,0,4294967295\tBF:B:f,1.5,-2.75")
    lines = [sam_line() + aux,
             sam_line(flag="4", ref="*", pos="0", cigar="*", mate="*", mpos="0", seq="u0123?", qual="*"),
             sam_line(pos="0", mate="unknown"), sam_line(cigar="*"),
             sam_line(seq="*", qual="*"), sam_line(pos="2000000001"),
             sam_line(name="long", cigar="1M" * 65536, seq="A" * 65536, qual="*"),
             sam_line(name="escape", cigar="6S6N") + b"\tCG:B:I,48,17,32\tNM:i:1"]
    sam_text = sam_header + b"\n".join(lines) + b"\n"
    sam_path = write("input.sam", sam_text)
    sam_oracle = write("sam-oracle.bam", pysam.view("-b", "--no-PG", str(sam_path)))
    expected_sam = invoke(sam_oracle)
    for text in (sam_text, sam_text.replace(b"\n", b"\r\n"), sam_text.rstrip(b"\n")):
        actual_sam = invoke(write("input.sam", text), "sam")
        assert actual_sam == expected_sam, [(a[:250], b[:250]) for a, b in zip(actual_sam, expected_sam) if a != b][:3]
    for size in (17, 65536):
        assert invoke(write("input.sam.gz", bgzf(sam_text, size)), "sam") == expected_sam
    assert invoke(sam_path, "sam-write", output)[-1] == "S\tEnd"
    assert invoke(output) == expected_sam
    assert invoke(write("empty.sam", b""), "sam") == ["S\tEnd"]
    assert invoke(write("header.sam", sam_header), "sam") == expected_sam[:2] + ["S\tEnd"]
    assert invoke(write("unplaced.sam", b"r\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"), "sam")[-1] == "S\tEnd"
    assert invoke(write("limit.sam", sam_header + sam_line(name="x" * 101)), "sam-limit")[-1] == "S\tTooLarge"
    boundary = sam_line()
    boundary = sam_line(name="r" * (101 - len(boundary)))
    assert len(boundary) == 100
    for ending in (b"", b"\n", b"\r\n"):
        for padding in (b"", b"x"):
            text = sam_header + padding + boundary + ending
            for wire in (text, bgzf(text, 1)):
                result = invoke(write("line-boundary.sam", wire), "sam-limit")[-1]
                assert result == ("S\tTooLarge" if padding else "S\tEnd"), (ending, padding, result)
    all_bases = sam_line(cigar="94M", seq="".join(map(chr, range(33, 127))), qual="*")
    special = sam_header + all_bases + b"\tFF:f:nan\tFI:f:inf\tFN:f:-inf\n"
    special_path = write("all-ascii.sam", special)
    special_bam = write("all-ascii.bam", pysam.view("-b", "--no-PG", str(special_path)))
    assert invoke(special_path, "sam") == invoke(special_bam)
    invalid_sam = [sam_line(flag="016"), sam_line(flag="-1"), sam_line(ref="unknown"), sam_line(pos="2147483649"),
                   sam_line(mapq="256"), sam_line(tlen="2147483648"), sam_line(cigar="3M1"),
                   sam_line(cigar="268435456M"), sam_line(cigar="1Q"), sam_line(cigar="5M"),
                   sam_line(qual="abc"), sam_line(seq="AC GTN"), sam_line(name="bad\0name"), b"", b"r\t4"]
    invalid_sam += [sam_line() + b"\t" + tag for tag in
                    (b"NI:i:4294967296", b"NI:i:-2147483649", b"NI:i:12x", b"FF:f:1.5x", b"FF:f:0x1p3",
                     b"HH:H:GG", b"HH:H:ABC", b"BC:B:c,128", b"BC:B:C,-0", b"BC:B:i,1,", b"AA:A:ab", b"XX:q:3")]
    for bad in invalid_sam:
        assert invoke(write("bad.sam", sam_header + bad + b"\n"), "sam")[-1] == "S\tInvalid", bad[:100]
    for bad_header in (b"@SQ\tSN:chrA\n", b"@SQ\tLN:1\n", b"@SQ\tSN:chrA\tLN:bad\n"):
        assert invoke(write("bad.sam", bad_header), "sam")[-1] == "S\tInvalid"

    # Stored accessors must preserve decoded values at nonzero slab offsets,
    # across reuse and reader close, including odd padding and every IUPAC code.
    even = record("even", 2, [("M", 16)], seq="=ACMGRSVTWYHKDBN", aux=b"NMi\x03\0\0\0")
    odd = record("odd", 3, [("S", 1), ("M", 13), ("I", 1)], seq="NBDKHYWTGVSRCMA",
                 flags=16, missing_qual=True)
    stored = write("stored.bam", bgzf(header() + even[0] + odd[0], 17))
    assert invoke(stored, "store", "chrA", 0, 100, 0) == [even[1], odd[1], "S\tRecord"]
    assert invoke(large, "store", "chrA", 1, 2, 0) == [long_expected, "S\tRecord"]
    # Budget is name + decoded sequence + qualities + CIGAR words + raw aux.
    # The first record needs 4+16+16+4+7=47 bytes; both need 47+3+15+15+12=92.
    assert invoke(stored, "store", "chrA", 2, 3, 0, 46) == ["S\tTooLarge"]
    assert invoke(stored, "store", "chrA", 2, 3, 0, 47) == [even[1], "S\tRecord"]
    assert invoke(stored, "store", "chrA", 0, 100, 0, 91) == ["S\tTooLarge"]
    assert invoke(stored, "store", "chrA", 0, 100, 0, 92) == [even[1], odd[1], "S\tRecord"]

    # A region's decoded payload can exceed 256 MiB even though each record
    # is small enough. Keep the fixture compressed and emit only pileup rows.
    dense = directory / "large-payload.bam"
    with pysam.AlignmentFile(dense, "wb", header={"SQ": [{"SN": n, "LN": size} for n, size in REFS]}) as out:
        read = pysam.AlignedSegment()
        read.reference_id = 0; read.reference_start = 2; read.mapping_quality = 60
        read.cigartuples = [(0, 1)]; read.query_sequence = "A"; read.query_qualities = [30]
        read.set_tag("ZZ", "x" * (1024 * 1024), value_type="Z")
        for i in range(257):
            read.query_name = f"payload-{i}"
            out.write(read)
    expected = sorted(f"2\tpayload-{i}\t0\t0\t0\t65\t30" for i in range(257))
    actual = invoke(dense, "pileup", "chrA", 2, 3, 0)
    assert actual[-1] == "S\tRecord" and sorted(actual[:-1]) == expected

    # Every IUPAC code, asymmetric SIMD lanes, odd tails, and nonzero slab
    # offsets. A lane swap or treating the final low nibble as a base must fail.
    tails = []
    for i, n in enumerate((1, 2, 31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256, 257)):
        seq = "".join("=ACMGRSVTWYHKDBN"[(j * 7 + j // 16 + i) % 16] for j in range(n))
        tails.append(record(f"tail{i}", i, [("M", n)], seq=seq))
    path = write("tails.bam", bgzf(header() + b"".join(r for r, _ in tails), 83))
    assert invoke(path, "store", "chrA", 0, 1000, 0) == [e for _, e in tails] + ["S\tRecord"]

    # Every truncated prefix of a complete compressed header block is an error.
    encoded_header = block(header())
    for cut in range(len(encoded_header)):
        status = invoke(write("cut.bam", encoded_header[:cut]))[-1]
        assert status == "S\tTruncated", (cut, status)
    damaged = bytearray(block(header()))
    damaged[-8] ^= 1
    malformed = [
        (bytes(damaged), "Invalid"),
        (block(header(), extra=b"BC\x02\x00\x1b\x00"), "Invalid"),
        (block(header(), extra=b"XY\xff\xff"), "Invalid"),
        (block(b"BOOM"), "Invalid"),
        (block(b"BAM\x01" + struct.pack("<i", -1)), "Invalid"),
        (block(b"BAM\x01" + struct.pack("<i", 1 << 30)), "TooLarge"),
        (block(header() + struct.pack("<i", 31)), "Invalid"),
        (block(header() + struct.pack("<i", 1 << 30)), "TooLarge"),
        (block(header() + b"\x20"), "Truncated"),
        (block(header() + struct.pack("<i", 100) + bytes(10)), "Truncated"),
    ]
    bad = bytearray(fixtures[0][0]); bad[4 + 8] = 255
    malformed.append((bgzf(header() + bad), "Invalid"))
    bad = bytearray(fixtures[0][0]); struct.pack_into("<i", bad, 4 + 16, -1)
    malformed.append((bgzf(header() + bad), "Invalid"))
    bad, _ = record("badaux", 1, [("M", 1)], aux=b"ZZBf\xff\xff\xff\x7f")
    malformed.append((bgzf(header() + bad), "Invalid"))
    bad, _ = record("badcigar", 1, [("M", 3)], seq="AC")
    malformed.append((bgzf(header() + bad), "Invalid"))
    damaged = bytearray(block(header())); damaged[12:14] = b"XY"
    malformed.append((bytes(damaged), "Invalid"))
    damaged = bytearray(block(header())); struct.pack_into("<I", damaged, len(damaged) - 4, 65537)
    malformed.append((bytes(damaged), "Invalid"))
    for i, (data, status) in enumerate(malformed):
        assert invoke(write("bad.bam", data))[-1] == f"S\t{status}", i

    unsorted = write("unsorted.bam", bgzf(header() + fixtures[1][0] + fixtures[0][0]))
    assert invoke(unsorted, "pileup", "chrA", 0, 100, 0)[-1] == "S\tUnsorted"

    # Convenience queries preserve order, fields and filters, with name/tid
    # overloads, backward reuse, empty regions and recovery after invalid input.
    def helper_queries(path, alignments, queries, indexes=(".bai", ".csi")):
        by_name = {e.split("\t")[1]: e for _, e in alignments}
        fetched, columns = [], []
        with pysam.AlignmentFile(path) as reader:
            for contig, begin, end in queries:
                label = f"Q\t{contig}\t{begin}\t{end}"
                fetched.append(label); columns.append(label)
                valid = contig in reader.references and 0 <= begin <= end <= reader.get_reference_length(contig)
                if valid and begin < end:
                    fetched.extend(by_name[r.query_name] for r in reader.fetch(contig, begin, end))
                    columns.extend(oracle(path, begin, end, 0x704, contig=contig))
                status = "S\tRecord" if valid else "S\tInvalid"
                fetched.append(status); columns.append(status)
        args = [v for query in queries for v in query]
        for suffix in indexes:
            for mode in ("fetch", "fetch-tid"):
                assert invoke(path, mode, str(path) + suffix, *args) == fetched, (path, mode, suffix)
            for mode in ("columns", "columns-tid"):
                assert invoke(path, mode, str(path) + suffix, *args) == columns, (path, mode, suffix)
        for mode in ("columns-scan", "columns-scan-tid"):
            assert invoke(path, mode, "-", *args) == columns, (path, mode)

    # Compare individual records, qpos, base, quality, D/N and anchored indels,
    # not only depth. Window boundaries cut through reads and CIGAR operations.
    comparisons = 0
    for name, alignments in [("crafted", fixtures), ("random", random_records()), ("long", [(long, long_expected)])]:
        path = write(name + ".bam", bgzf(header() + b"".join(r for r, _ in alignments)))
        pysam.index(str(path))
        pysam.index("-c", str(path))
        helper_queries(path, alignments, [("chrA", 75, 100), ("chrA", 11, 18), ("chrA", 18, 19),
            ("chrA", 22, 25), ("chrA", 7, 7), ("chrA", 150000, 150010), ("missing", 0, 10),
            ("chrA", -1, 10), ("chrA", 10, 9), ("chrA", 0, 200001), ("chrA", 5, 7)])
        text_path = write(name + ".sam", pysam.view("-h", "--no-PG", str(path)).encode())
        for begin, end in [(0, 150), (11, 18), (18, 19), (75, 100)]:
            for mask, min_mapq in ((0, 0), (0x704, 0), (0x704, 20), (0, 42)):
                actual = invoke(path, "pileup", "chrA", begin, end, mask, min_mapq)
                assert actual[-1] == "S\tRecord"
                expected = oracle(path, begin, end, mask, min_mapq)
                assert actual[:-1] == expected, (name, begin, end, mask,
                    next(((a, b) for a, b in zip(actual, expected) if a != b), (len(actual)-1, len(expected))))
                assert invoke(text_path, "sam-pileup", "chrA", begin, end, mask, min_mapq) == expected + ["S\tRecord"]
                for suffix in (".bai", ".csi"):
                    assert invoke(path, "indexed", "chrA", begin, end, mask, min_mapq,
                                  str(path) + suffix) == expected + ["S\tRecord"], (name, suffix, begin, end)
                comparisons += 1
    second, second_expected = record("second_ref", 2, [("M", 4)], tid=1)
    unmapped, unmapped_expected = record("unmapped", -1, [], seq="ACG", tid=-1, flags=4)
    path = write("multiref.bam", bgzf(header() + fixtures[0][0] + second + unmapped))
    assert invoke(path) == expected_header + [fixtures[0][1], second_expected, unmapped_expected, "S\tEnd"]
    pysam.index(str(path))
    assert invoke(path, "pileup", "chrB", 0, 10, 0) == oracle(path, 0, 10, 0, contig="chrB") + ["S\tRecord"]
    assert invoke(path, "indexed", "chrB", 0, 10, 0, 0, str(path) + ".bai") == oracle(path, 0, 10, 0, contig="chrB") + ["S\tRecord"]
    helper_queries(path, [fixtures[0], (second, second_expected), (unmapped, unmapped_expected)],
                   [("chrB", 0, 10), ("chrA", 5, 18)], indexes=(".bai",))

    # Discovery: test every spelling independently, then precedence and a bad
    # preferred index alongside a valid fallback. Failures clear a prior index.
    discovered = write("discover.bam", bgzf(header() + b"".join(r for r, _ in fixtures)))
    pysam.index(str(discovered)); pysam.index("-c", str(discovered))
    bai, csi = (Path(str(discovered) + suffix).read_bytes() for suffix in (".bai", ".csi"))
    candidates = [Path(str(discovered) + ".csi"), discovered.with_suffix(".csi"),
                  Path(str(discovered) + ".bai"), discovered.with_suffix(".bai")]
    for candidate in candidates: candidate.unlink(missing_ok=True)
    saved = write("saved-index.bai", bai)
    assert invoke(discovered, "open-for", saved) == ["S\tIoError"]
    for candidate in candidates:
        candidate.write_bytes(csi if candidate.suffix == ".csi" else bai)
        assert invoke(discovered, "open-for") == ["S\tRecord"]
        assert invoke(discovered, "fetch", "-", "chrA", 5, 7) == ["Q\tchrA\t5\t7", fixtures[0][1], "S\tRecord"]
        candidate.unlink()
    for i, candidate in enumerate(candidates[:-1]):
        for fallback in candidates[i + 1:]: fallback.write_bytes(bai if fallback.suffix == ".bai" else csi)
        candidate.write_bytes(b"BAD!")
        assert invoke(discovered, "open-for", saved) == ["S\tInvalid"]
        for item in candidates: item.unlink(missing_ok=True)

    # Coarse bins, disjoint chunks, bin boundaries, long skips starting before
    # the window, and out-of-order reuse. htslib builds both kinds of index.
    spanning = [record("spanning", 1, [("M", 2), ("N", 180000), ("M", 2)])]
    spanning += [record(f"at{p}", p, [("M", 7)]) for p in (16380, 16384, 32768, 65534, 131072, 180000)]
    path = write("bins.bam", bgzf(header() + b"".join(r for r, _ in spanning), 83))
    pysam.index(str(path)); pysam.index("-c", "-m", "12", str(path))
    helper_queries(path, spanning, [("chrA", 179999, 180006), ("chrA", 16383, 16385),
                                    ("chrA", 65535, 65538), ("chrA", 2, 3)])
    for begin, end in ((16383, 16385), (65535, 65538), (131071, 131075), (179999, 180006), (2, 3), (190000, 190010), (8, 8)):
        expected = oracle(path, begin, end, 0) + ["S\tRecord"] if begin != end else ["S\tRecord"]
        for suffix in (".bai", ".csi"):
            assert invoke(path, "indexed", "chrA", begin, end, 0, 0, str(path) + suffix) == expected

    # A valid root-bin index, overlapping duplicate chunks and a pseudo-bin
    # whose count pair must not be interpreted as virtual offsets.
    raw = header() + fixtures[0][0] + fixtures[1][0]
    path = write("root-bin.bam", bgzf(raw))
    index_bytes = (b"BAI\x01" + struct.pack("<iiiiQQQQiiQQQQiQii", 2, 2, 0, 2,
        len(header()), len(raw), len(header()), len(raw), 37450, 2,
        len(header()), len(raw), 2, 0, 1, len(header()), 0, 0))
    index = write("root.bai", index_bytes)
    assert invoke(path, "indexed-store", "chrA", 5, 18, 0, 0, index) == [fixtures[0][1], fixtures[1][1], "S\tRecord"]
    assert invoke(path, "fetch", index, "chrA", 5, 18) == ["Q\tchrA\t5\t18", fixtures[0][1], fixtures[1][1], "S\tRecord"]
    assert invoke(unsorted, "fetch", index, "chrA", 0, 100)[-1] == "S\tUnsorted"
    assert invoke(index, "index", len(index_bytes)) == ["S\tRecord"]
    assert invoke(index, "index", len(index_bytes) - 1) == ["S\tTooLarge"]
    assert invoke(write("compressed.bai", bgzf(index_bytes, 17)), "index") == ["S\tRecord"]
    for cut in range(len(index_bytes)):
        assert invoke(write("cut.bai", index_bytes[:cut]), "index")[-1] != "S\tRecord", cut
    for bad in (b"BAI\x01" + struct.pack("<i", -1), b"CSI\x01" + struct.pack("<iii", 63, 5, 0),
                b"CSI\x01" + struct.pack("<iii", 14, -1, 0), b"CSI\x01" + struct.pack("<iii", 0, 11, 0),
                b"BAI\x01" + struct.pack("<iiiiQQ", 1, 1, 0, 1, 200, 100)):
        assert invoke(write("bad.idx", bad), "index")[-1] != "S\tRecord"
    for begin, end in ((-1, 5), (10, 9), (0, 200001)):
        assert invoke(path, "indexed", "chrA", begin, end, 0, 0, index) == ["S\tInvalid"]

    # Seqair's indexed path rejects non-overlaps before decoding auxiliary
    # tags. The sequential reader and a query retaining the bad read still fail.
    bad, _ = record("badaux", 1, [("M", 2)], aux=b"ZZQ")
    good, good_expected = record("after", 100, [("M", 5)])
    raw = header() + bad + good
    path = write("prefilter.bam", bgzf(raw))
    index = write("prefilter.bai", b"BAI\x01" + struct.pack("<iiiiQQiii", 2, 1, 0, 1,
                                                          len(header()), len(raw), 0, 0, 0))
    assert invoke(path)[-1] == "S\tInvalid"
    assert invoke(path, "indexed-store", "chrA", 2, 3, 0, 0, index) == ["S\tInvalid"]
    assert invoke(path, "indexed-store", "chrA", 3, 4, 0, 0, index) == ["S\tRecord"]
    assert invoke(path, "indexed-store", "chrA", 100, 105, 0, 0, index) == [good_expected, "S\tRecord"]

    # An index must really skip a prefix, not just wrap the sequential loader.
    # Keep the header valid, corrupt an irrelevant block, and index a later one.
    prefix = block(header())
    corrupt = bytearray(block(fixtures[0][0])); corrupt[-8] ^= 1
    late, late_expected = record("late", 100000, [("M", 5)])
    start = len(prefix) + len(corrupt)
    late_block = block(late)
    path = write("skip-prefix.bam", prefix + corrupt + late_block + block(b""))
    index = write("skip-prefix.bai", b"BAI\x01" + struct.pack("<iiiiQQiii", 2, 1, 0, 1,
        start << 16, (start + len(late_block)) << 16, 0, 0, 0))
    assert invoke(path)[-1] == "S\tInvalid"
    assert invoke(path, "indexed-store", "chrA", 100000, 100005, 0, 0, index) == [late_expected, "S\tRecord"]
    assert invoke(path, "fetch", index, "chrA", 100000, 100005) == ["Q\tchrA\t100000\t100005", late_expected, "S\tRecord"]
    # Nonexistent blocks and out-of-block uncompressed offsets must fail.
    for offset in ((start << 16) | 65535, (len(prefix + corrupt + late_block) + 1000) << 16):
        bad = write("bad-offset.bai", b"BAI\x01" + struct.pack("<iiiiQQiii", 2, 1, 0, 1, offset, offset + 1, 0, 0, 0))
        assert invoke(path, "indexed", "chrA", 100000, 100005, 0, 0, bad)[-1] in ("S\tInvalid", "S\tTruncated")
        assert invoke(path, "fetch", bad, "chrA", 100000, 100005)[-1] in ("S\tInvalid", "S\tTruncated")

    # CSI's configurable tree must handle coordinates BAI cannot represent.
    high, high_expected = record("high", 600000000, [("M", 8)])
    path = write("high.bam", bgzf(header(references=[("chrA", 700000000), REFS[1]]) + high))
    pysam.index("-c", str(path))
    csi = Path(str(path) + ".csi")
    raw_csi = write("raw.csi", gzip.decompress(csi.read_bytes()))
    for index in (csi, raw_csi):
        assert invoke(path, "indexed-store", "chrA", 600000003, 600000005, 0, 0, index) == [high_expected, "S\tRecord"]
        assert invoke(path, "fetch", index, "chrA", 600000003, 600000005) == ["Q\tchrA\t600000003\t600000005", high_expected, "S\tRecord"]
    for cut in range(len(raw_csi.read_bytes())):
        # CSI's optional final eight-byte count may be absent, but not partial.
        if cut == len(raw_csi.read_bytes()) - 8: continue
        assert invoke(write("cut.csi", raw_csi.read_bytes()[:cut]), "index")[-1] != "S\tRecord", cut

    # Mate links are store indexes, not qname matching in each column. READ1/
    # READ2 bits do not determine eligibility in Seqair. A third primary is
    # reported, while secondary/supplementary and unavailable mates are ignored.
    pairs = [record("pair", 5, [("M", 20)], flags=65, mate_tid=0, mate_pos=10),
             record("pair", 10, [("M", 20)], flags=129, mate_tid=0, mate_pos=5),
             record("pair", 11, [("M", 20)], flags=129, mate_tid=0, mate_pos=5),
             record("same_read_bit", 12, [("M", 20)], flags=65, mate_tid=0, mate_pos=14),
             record("same_read_bit", 14, [("M", 20)], flags=65, mate_tid=0, mate_pos=12)]
    for i, (flag, mate_tid, mate_pos) in enumerate(((0, 0, 50), (0x141, 0, 50), (0x841, 0, 50), (73, 0, 50), (65, 1, 50), (65, 0, -1))):
        pairs += [record(f"excluded{i}", 20 + i, [("M", 20)], flags=flag, mate_tid=mate_tid, mate_pos=mate_pos)]
    pairs += [record("", 30, [("M", 20)], flags=65, mate_tid=0, mate_pos=31),
              record("", 31, [("M", 20)], flags=129, mate_tid=0, mate_pos=30)]
    pairs += [record(f"excluded{i}", 50, [("M", 20)], flags=129, mate_tid=0, mate_pos=20 + i) for i in range(6)]
    path = write("mates.bam", bgzf(header() + b"".join(r for r, _ in pairs), 97))
    links = [1, 0, -1, 4, 3] + [-1] * (len(pairs) - 5)
    expected = ["M\t2\t1"] + [f"{i}\t{mate}" for i, mate in enumerate(links)] + ["S\tRecord"]
    assert invoke(path, "mates", "chrA", 0, 100, 0) == expected
    mate_sam = write("mates.sam", pysam.view("-h", "--no-PG", str(path)).encode())
    assert invoke(mate_sam, "sam-mates", "chrA", 0, 100, 0, 0) == expected
    pysam.index(str(path)); pysam.index("-c", str(path))
    for suffix in (".bai", ".csi"):
        assert invoke(path, "indexed-mates", "chrA", 0, 100, 0, 0, str(path) + suffix) == expected
        assert invoke(path, "indexed-mates", "chrA", 5, 6, 0, 0, str(path) + suffix) == ["M\t0\t0", "0\t-1", "S\tRecord"]
    # Opting out changes only auxiliary retention and automatic mate links.
    # The long-CIGAR record must still resolve its operations from the CG tag.
    for source in (path, large, sources[2]):
        if not Path(str(source) + ".bai").exists():
            pysam.index(str(source))
        normal = invoke(source, "store", "chrA", 0, 100000, 0)
        without_aux = [line.rsplit("\t", 1)[0] + "\t" if line.startswith("R\t") else line
                       for line in normal]
        assert invoke(source, "lean", "chrA", 0, 100000, 0) == without_aux
        assert invoke(source, "indexed-lean", "chrA", 0, 100000, 0, 0,
                      str(source) + ".bai") == without_aux
    print(f"PASS {command}: wire/index/seek/mate cases and {comparisons} htslib comparisons on each of scan/BAI/CSI", flush=True)
    from cram import suite as cram_suite
    cram_suite(command, directory)
    from writers import suite as writers_suite
    writers_suite(command, directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", nargs="*", default=["clang"])
    parser.add_argument("--native-only", action="store_true")
    args = parser.parse_args()
    goose = os.environ.get("GOOSE", "goose")
    source = ROOT / "tests/probe.goose"
    with tempfile.TemporaryDirectory(prefix="goose-seqair-tests-") as temp:
        directory = Path(temp)
        if not args.native_only:
            for opt in (0, 2):
                suite([goose, "run", source, "--debug", f"-O{opt}", "--"], directory)
        for cc in args.cc:
            for debug in (True, False):
                exe = directory / ("probe.exe" if os.name == "nt" else "probe")
                subprocess.run([goose, "build", str(source), "--cc", cc, "-O2",
                                *(["--debug"] if debug else []), "-o", str(exe)], cwd=ROOT, check=True)
                suite([exe], directory)


if __name__ == "__main__":
    main()
