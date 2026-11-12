"""CRAM public-reader/region tests against independently encoded htslib files."""
from pathlib import Path
import subprocess

import pysam
from fixtures import bgzf, header, record


def suite(command, directory):
    def invoke(mode, path, reference, *args):
        return subprocess.check_output([*map(str, command), mode, str(path), str(reference),
                                        *map(str, args)], text=True).splitlines()

    reference = directory / "reference.fa"
    # Contigs exceed the probe provider's one-megabyte maximum fetch. Small
    # multi-reference slices must fetch span unions, never entire contigs.
    length = 1000100
    reference.write_text(">chrA\n" + "ACGTN" * (length // 5) + "\n>chrB\n" + "TGCAW" * (length // 5) + "\n")
    pysam.faidx(str(reference))
    empty = directory / "absent.fa"; empty.write_text("")
    wrong = directory / "wrong.fa"; wrong.write_text(reference.read_text().replace("ACGTN", "TTTTT"))
    text = (f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrA\tLN:{length}\n"
            f"@SQ\tSN:chrB\tLN:{length}\n@RG\tID:group\n").encode()
    reads = [record("pair", 5, [("M", 13)], flags=65, mate_tid=0, mate_pos=14, template_len=22),
             record("complex", 8, [("S", 2), ("M", 3), ("I", 2), ("M", 4), ("D", 3), ("M", 2), ("H", 1)]),
             record("splice", 9, [("M", 2), ("N", 7), ("M", 3)], aux=b"RGZgroup\0"),
             record("pair", 14, [("M", 13)], flags=129, mate_tid=0, mate_pos=5, template_len=-22),
             record("missing-quality", 15, [("M", 9)], missing_qual=True),
             record("missing-sequence", 20, [("M", 5)], seq=""),
             record("iupac", 24, [("M", 15)], seq="ACMGRSVTWYHKDBN"),
             record("other", 11, [("M", 15)], tid=1),
             record("unplaced", -1, [], flags=4, tid=-1, seq="NACGT", mapq=0)]
    bam = directory / "cram-input.bam"
    bam.write_bytes(bgzf(header(text, [("chrA", length), ("chrB", length)]) + b"".join(r for r, _ in reads)))
    output = directory / "cram-output.bam"

    def canonical(path, mode, ref=None):
        with pysam.AlignmentFile(path, mode, reference_filename=str(ref) if ref else None) as f:
            return [(r.query_name, r.flag, r.reference_id, r.reference_start, r.mapping_quality,
                     r.cigarstring, r.next_reference_id, r.next_reference_start, r.template_length,
                     r.query_sequence, list(r.query_qualities) if r.query_qualities is not None else None,
                     sorted(r.get_tags())) for r in f]

    cases = [("v30", ["version=3.0"]), ("v31", ["version=3.1"]),
             ("multi", ["version=3.1", "multi_seq_per_slice=1", "seqs_per_slice=100"]),
             ("slices", ["seqs_per_slice=3", "slices_per_container=2"]),
             ("embedded", ["embed_ref=2"]), ("no-ref", ["no_ref=1"])]
    for name, options in cases:
        path = directory / (name + ".cram")
        args = ["--no-PG", "-C", "-T", str(reference), "-o", str(path)]
        for option in options: args += ["--output-fmt-option", option]
        pysam.view(*args, str(bam), catch_stdout=False)
        ref = empty if name in ("embedded", "no-ref") else reference
        assert invoke("cram-write", path, ref, output) == ["S\tEnd"], name
        assert canonical(output, "rb") == canonical(path, "rc", reference), name
        pysam.index(str(output))
        # Compare every per-read column entry, including deletion/refskip,
        # query offsets, bases and missing qualities, not merely depth.
        from run import oracle
        for contig, begin, end in (("chrA", 9, 30), ("chrB", 12, 18), ("chrA", 80, 90)):
            expected = oracle(output, begin, end, 0, contig=contig)
            assert invoke("cram-pileup", path, ref, contig, begin, end) == expected + ["S\tRecord"], name
        assert invoke("cram-limit", path, ref)[-1] == "S\tTooLarge"
        if name in ("v30", "v31"):
            assert invoke("cram", path, empty) == ["S\tInvalid"]
            assert invoke("cram", path, wrong) == ["S\tInvalid"]
            assert invoke("cram-pileup", path, wrong, "chrA", 0, 100) == ["S\tInvalid"]

    path = directory / "v31.cram"
    raw = path.read_bytes()
    damaged = directory / "damaged.cram"
    for cut in (0, 3, 5, 25, 29, len(raw) // 2, len(raw) - 1):
        damaged.write_bytes(raw[:cut])
        assert invoke("cram", damaged, reference)[-1] in ("S\tTruncated", "S\tInvalid"), cut
    # Seqair recognizes EOF by its container header and does not decode the
    # marker's dummy blocks. Corrupt file headers and actual data blocks here.
    for offset in (0, 4, 5, 26, len(raw) // 2):
        corrupt = bytearray(raw); corrupt[offset] ^= 255; damaged.write_bytes(corrupt)
        assert invoke("cram", damaged, reference)[-1] != "S\tEnd", offset
    print(f"PASS {command}: CRAM 3.0/3.1 records, references, mates, regions and corruption", flush=True)
