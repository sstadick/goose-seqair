"""Independent htslib checks for SAM formatting and reference-free CRAM output."""
import gzip
import random
import struct
import subprocess
from pathlib import Path

import pysam
from fixtures import BASES, bgzf, header, record, records, random_records


def normalized(path, normalize=False):
    with pysam.AlignmentFile(path) as reader:
        result = []
        for r in reader:
            if normalize:
                cigar = []
                for op, n in r.cigartuples or []:
                    if not n: continue
                    op = 0 if op in (7, 8) else op
                    if cigar and cigar[-1][0] == op: cigar[-1] = (op, cigar[-1][1] + n)
                    else: cigar.append((op, n))
                r.cigartuples = cigar
            result.append(r.to_string())
        return result


def packed_aux(data, name):
    """Inspect independently emitted BAM, retaining NaN payloads and tag widths."""
    data = gzip.decompress(data)
    at = 8 + struct.unpack_from("<i", data, 4)[0]
    refs, = struct.unpack_from("<i", data, at); at += 4
    for _ in range(refs): at += 8 + struct.unpack_from("<i", data, at)[0]
    while at < len(data):
        size, = struct.unpack_from("<i", data, at); at += 4
        body = data[at:at + size]; at += size
        qname = body[8]; cigar, = struct.unpack_from("<H", body, 12)
        length, = struct.unpack_from("<i", body, 16)
        if body[32:32 + qname - 1] == name:
            return body[32 + qname + 4 * cigar + (length + 1) // 2 + length:]
    raise AssertionError("missing record")


def suite(command, directory):
    def invoke(mode, path, output, *args):
        return subprocess.check_output([*map(str, command), mode, str(path), str(output), *map(str, args)], text=True).splitlines()

    def make(name, rows, text=b"@HD\tVN:1.6\n@RG\tID:group\n"):
        path = directory / (name + ".bam")
        path.write_bytes(bgzf(header(text) + b"".join(r for r, _ in rows)))
        return path

    # Every scalar/array tag representation, including integer boundaries;
    # array floats around decimal bands, ties, signed zero, and IEEE extremes.
    aux = b"AAAxZZZvalue\0HHH01aB\0"
    for tag, fmt, values in [(b"a", "b", [-128, 127]), (b"b", "B", [0, 255]),
                             (b"c", "h", [-32768, 32767]), (b"d", "H", [0, 65535]),
                             (b"e", "i", [-2147483648, 2147483647]), (b"f", "I", [0, 4294967295])]:
        kind = {"b": b"c", "B": b"C", "h": b"s", "H": b"S", "i": b"i", "I": b"I"}[fmt]
        aux += tag + b"1" + kind + struct.pack("<" + fmt, values[0])
        aux += tag + b"2B" + kind + struct.pack("<I" + fmt * 2, 2, *values)
    floats = [-0.0, 0.0, 1.234565, 1.234575, -123.4565, 999999.0, 0.0001, float("inf"), -float("inf"), float("nan")]
    for exponent in range(-4, 7):
        bits, = struct.unpack("<I", struct.pack("<f", 10.0 ** exponent))
        floats.extend(struct.unpack("<f", struct.pack("<I", bits + d))[0] for d in (-1, 0, 1))
    rng = random.Random(773)
    # Keep raw NaN payloads too: casting float to Python and back can change them.
    float_bytes = b"".join(struct.pack("<f", f) for f in floats) + rng.randbytes(4000)
    aux += b"ffBf" + struct.pack("<I", len(float_bytes) // 4) + float_bytes
    aux += b"eeBi\0\0\0\0"
    tagged = record("tags", 1, [("M", 16)], seq=BASES, aux=aux)
    scalar_rows = [record(f"float{i}", i + 2, [("M", 1)], aux=b"fff" + float_bytes[i:i + 4])
                   for i in range(0, len(float_bytes), 4)]
    long = record("long", 1, [("M", 1)] * 65536)
    path = make("sam-writer", records() + [tagged, long] + scalar_rows + [record("double", 5, [("M", 1)], aux=b"ddd" + struct.pack("<d", -1.234567e199))])
    sam = directory / "written.sam"
    assert invoke("to-sam", path, sam) == ["S\tEnd"]
    assert sam.read_text() == pysam.view("-h", "--no-PG", str(path)), "SAM formatter differs from htslib"

    rows = [r for r in records() if not r[1].startswith("R\tempty_cigar\t")] + random_records() + [tagged]
    rows += [record("pair", 3, [("M", 8)], flags=99, mate_tid=1, mate_pos=1, template_len=12),
             record("pair", 1, [("M", 8)], flags=147, tid=1, mate_tid=0, mate_pos=3, template_len=-12),
             record("unmapped", -1, [], flags=4, tid=-1, seq="=ACN", mapq=0),
             record("missing", -1, [], flags=4, tid=-1, seq="", mapq=0),
             record("max-ap", 2147483646, [("M", 1)]), long]
    path = make("cram-writer", rows)
    expected = normalized(path, True)
    cram = directory / "written.cram"
    empty_ref = directory / "writer-absent.fa"; empty_ref.write_text("")
    bam = directory / "writer-roundtrip.bam"
    eof = bytes.fromhex("0f000000ffffffff0fe0454f4600000000010005bdd94f0001000606010001000100ee63014b")
    for minor, level, count in [(0, 0, 1), (1, 6, 7), (1, 12, 10000)]:
        assert invoke("to-cram", path, cram, level, count, minor) == ["S\tEnd"]
        assert cram.read_bytes()[:6] == b"CRAM\x03" + bytes([minor])
        assert cram.read_bytes().endswith(eof)
        assert normalized(cram) == expected, (minor, level, "htslib")
        assert packed_aux(pysam.view("-b", "--no-PG", str(cram)), b"tags") == aux
        assert invoke("cram-write", cram, empty_ref, bam) == ["S\tEnd"]
        assert normalized(bam) == expected, (minor, level, "Goose")
        assert packed_aux(bam.read_bytes(), b"tags") == aux

    # Empty output, owned alignment input, slab input, reopen, sticky errors,
    # and close failures are observable parts of both writers' public contract.
    empty = make("writer-empty", [])
    small = make("writer-small", [record("first", 5, [("M", 3)]), record("second", 8, [("M", 2)])])
    for kind in ("sam", "cram"):
        output = directory / ("api." + kind)
        assert invoke("to-" + kind, empty, output) == ["S\tEnd"]
        assert normalized(output) == []
        empty_read = invoke(kind, output, empty_ref)
        assert empty_read[-1] == "S\tEnd" and not any(r.startswith("R\t") for r in empty_read)
        assert invoke("writer-api-" + kind, small, output) == ["S\tOk"]
        result = normalized(output)
        assert result[0] == "owned\t0\tchrA\t4\t0\t3M\t*\t0\t0\tACN\t(7I\tZZ:Z:value"
        assert result[1:] == normalized(small)
        assert normalized(Path(str(output) + ".reopened")) == result[:1]
        # Reopened output contains exactly one alignment: the probe checks
        # Ok after open, Record on that read, then End and repeated End.
        single_read = invoke(kind, Path(str(output) + ".reopened"), empty_ref)
        assert single_read[-1] == "S\tEnd"
        assert [r.split("\t")[1] for r in single_read if r.startswith("R\t")] == ["owned"]
        assert invoke("to-" + kind, small, directory / "nonexistent-dir/out") == ["S\tIoError"]
        if Path("/dev/full").exists(): assert invoke("to-" + kind, small, "/dev/full") == ["S\tIoError"]

    one = make("writer-limit", [record("x", 5, [("M", 1)])])
    line_len = len(pysam.view("--no-PG", str(one)).encode())
    assert invoke("to-sam", one, sam, line_len) == ["S\tEnd"]
    assert invoke("to-sam", one, sam, line_len - 1) == ["S\tTooLarge"]
    assert invoke("to-sam", one, sam, 10000, 32) == ["S\tTooLarge"]
    assert invoke("to-cram", one, cram, 6, 100, 1, 32) == ["S\tTooLarge"]
    assert invoke("to-cram", small, cram, 6, 100, 1, 256, 1000) == ["S\tEnd"]
    assert normalized(cram) == normalized(small)
    absent = make("writer-absent-limit", [record("missing", 5, [("M", 1000)], seq="")])
    assert invoke("to-cram", absent, cram, 6, 100, 1, 1024, 1024) == ["S\tTooLarge"]
    for bad in [record("cigar", 5, [], seq="A"), record("mapq", -1, [], tid=-1, flags=4, mapq=30),
                record("reserved", 5, [("M", 1)], aux=b"cFC\x01"),
                record("double", 5, [("M", 1)], aux=b"ddd" + struct.pack("<d", 1.0)),
                record("mate", 5, [("M", 1)], mate_tid=0), record("coordinate", 2147483647, [("M", 1)])]:
        assert invoke("to-cram", make("writer-rejected", [bad]), cram) == ["S\tInvalid"]
    conflicting = make("writer-header", [], b"@SQ\tSN:wrong\tLN:10\n")
    assert invoke("to-cram", conflicting, cram) == ["S\tInvalid"]
    for text in (b"@HD\tVN:1.6", b"@HD\tVN:1.6\n\0ignored"):
        path = make("writer-text", [], text)
        assert invoke("to-sam", path, sam) == ["S\tEnd"]
        # Seqair explicitly truncates at NUL. htslib 1.24 warns but emits the
        # NUL and following garbage in this malformed-header case.
        assert sam.read_bytes() == b"@HD\tVN:1.6\n@SQ\tSN:chrA\tLN:200000\n@SQ\tSN:chrB\tLN:400\n"
    crlf = make("writer-crlf", [], b"@SQ\tSN:chrA\tLN:200000\r\n@SQ\tSN:chrB\tLN:400\r\n")
    assert invoke("to-cram", crlf, cram) == ["S\tEnd"]
    assert normalized(cram) == []
    print(f"PASS {command}: SAM/CRAM writers, htslib fidelity, floats, limits and lifecycle", flush=True)
