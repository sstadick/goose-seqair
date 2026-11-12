"""Small independently encoded BAM/BGZF fixtures; no Goose or htslib required."""
import random
import struct
import zlib

BASES = "=ACMGRSVTWYHKDBN"
OPS = "MIDNSHP=X"
REFS = [("chrA", 200000), ("chrB", 400)]


def block(data, extra=b""):
    encoder = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = encoder.compress(data) + encoder.flush()
    size = 18 + len(extra) + len(payload) + 8
    assert size <= 65536
    return (b"\x1f\x8b\x08\x04" + bytes(5) + b"\xff" + struct.pack("<H", 6 + len(extra))
            + extra + b"BC\x02\x00" + struct.pack("<H", size - 1) + payload
            + struct.pack("<II", zlib.crc32(data), len(data)))


def bgzf(data, chunk=60000, empty=False, extra=b""):
    chunks = [block(data[i:i + chunk], extra) for i in range(0, len(data), chunk)]
    return (block(b"") if empty else b"").join(chunks) + block(b"")


def header(text=b"@HD\tVN:1.6\tSO:coordinate\n", references=REFS):
    return (b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<i", len(references))
            + b"".join(struct.pack("<i", len(n) + 1) + n.encode() + b"\0" + struct.pack("<i", l)
                       for n, l in references))


def record(name, pos, cigar, seq=None, flags=0, tid=0, mapq=41, aux=b"", missing_qual=False,
           mate_tid=-1, mate_pos=-1, template_len=0):
    qlen = sum(n for op, n in cigar if op in "MIS=X")
    if seq is None:
        seq = ("ACGTNMR" * ((qlen + 6) // 7))[:qlen]
    qual = bytes([255] * len(seq) if missing_qual else [(i * 7 + 13) % 42 for i in range(len(seq))])
    words = [(n << 4) | OPS.index(op) for op, n in cigar]
    packed = bytes((BASES.index(seq[i]) << 4) | (BASES.index(seq[i + 1]) if i + 1 < len(seq) else 0)
                   for i in range(0, len(seq), 2))
    wire_words = words
    wire_aux = aux
    span = sum(n for op, n in cigar if op in "MDN=X")
    if len(words) > 65535:
        wire_words = [(len(seq) << 4) | 4, (span << 4) | 3]
        wire_aux += b"CGBI" + struct.pack("<i", len(words)) + struct.pack(f"<{len(words)}I", *words)
    core = struct.pack("<iiBBHHHiiii", tid, pos, len(name) + 1, mapq, 0, len(wire_words), flags,
                       len(seq), mate_tid, mate_pos, template_len)
    body = core + name.encode() + b"\0" + struct.pack(f"<{len(wire_words)}I", *wire_words) + packed + qual + wire_aux
    expected = (f"R\t{name}\t{tid}\t{pos}\t{pos + (1 if flags & 4 else max(span, 1))}\t{flags}\t{mapq}\t{mate_tid}\t{mate_pos}\t{template_len}\t"
                + "".join(f"{w}," for w in words) + "\t" + seq + "\t"
                + "".join(f"{v}," for v in qual) + "\t" + "".join(f"{v}," for v in wire_aux))
    return struct.pack("<i", len(body)) + body, expected


def records():
    specs = [
        ("match", 5, [("M", 13)]),
        ("reverse", 7, [("S", 2), ("M", 4), ("I", 3), ("=", 2), ("X", 2), ("H", 1)]),
        ("deletion", 8, [("M", 2), ("D", 3), ("M", 5)]),
        ("splice", 9, [("M", 2), ("N", 7), ("M", 3)]),
        ("complex", 10, [("M", 2), ("D", 2), ("I", 1), ("M", 3)]),
        ("pad", 11, [("M", 2), ("P", 1), ("I", 2), ("P", 2), ("M", 4)]),
        ("adjacent", 12, [("M", 2), ("I", 1), ("I", 2), ("M", 4)]),
        ("duplicate", 14, [("M", 3)]),
        ("missing_qual", 15, [("M", 4)]),
        ("secondary", 16, [("M", 6)]),
        ("supplementary", 17, [("M", 4)]),
        ("leading_deletion", 18, [("D", 2), ("M", 5)]),
        ("leading_skip", 19, [("N", 3), ("M", 2)]),
        ("leading_insertion", 20, [("I", 2), ("M", 5)]),
        ("missing_seq", 21, [("M", 5)]),
        ("empty_cigar", 22, []),
        ("insertion_only", 23, [("I", 3)]),
        ("soft_only", 24, [("S", 3)]),
        ("qc_failed", 25, [("M", 3)]),
        ("low_mapq", 26, [("M", 3)]),
        ("far", 95, [("M", 7)]),
    ]
    flags = {"reverse": 16, "duplicate": 1024, "secondary": 256, "supplementary": 2048, "qc_failed": 512}
    return [record(n, p, c, flags=flags.get(n, 0), missing_qual=n == "missing_qual",
                   seq="" if n == "missing_seq" else "ACG" if n == "empty_cigar" else None,
                   mapq=10 if n == "low_mapq" else 41,
                   aux=b"NMi" + struct.pack("<i", 3) + b"RGZgroup\0") for n, p, c in specs]


def random_records(seed=2026, count=70):
    rng = random.Random(seed)
    result = []
    for i in range(count):
        cigar = [("S", rng.randrange(4))]
        for j in range(rng.randrange(1, 9)):
            cigar.append((rng.choice("M=X"), rng.randrange(1, 8)))
            cigar.append((rng.choice("IDN"), rng.randrange(1, 6)))
        cigar.append(("M", rng.randrange(1, 8)))
        result.append(record(f"random{i}", i // 2 + 3, cigar, flags=16 if i % 3 == 0 else 0))
    return result


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    data = header() + b"".join(r for r, _ in records())
    args.output.write_bytes(bgzf(data))
    # A legal coarse root-bin index for the small, single-block documentation fixture.
    assert len(data) < 60000
    Path(str(args.output) + ".bai").write_bytes(
        b"BAI\x01" + struct.pack("<iiiiQQiQii", 2, 1, 0, 1, len(header()), len(data),
                                   1, len(header()), 0, 0))
