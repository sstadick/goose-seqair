#!/usr/bin/env python3
"""Check unannotated seqair callees from checkfree code, and retained checks."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

POSITIVE = r'''
checkfree fn inspect(record: const pileup::Record&) {
    assert(record.name() == "match");
    assert(record.aux() == "NMi\x03\x00\x00\x00RGZgroup\x00");
    var called = false;
    bam::Status.End.check(bam::Status.Record) { actual =>
        assert(actual == bam::Status.End);
        called = true;
    };
    assert(called);
}
fn main() {
    let reader = bam::Reader { .. };
    defer reader.close();
    reader.open("tests/data/example.bam").check(bam::Status.Ok);
    let store = pileup::Store { .. };
    store.load(reader, 0, 5, 7).check(bam::Status.Ok);
    reader.close();
    let record = store.record(0);
    inspect(record);
}
'''

NEGATIVE = r'''
checkfree fn unchecked_quality(record: const pileup::Record&, index: i64) -> u8 {
    record.quality(index)
}
fn main() {
    let record = pileup::Record { .. };
    print(unchecked_quality(record, 0));
}
'''


def main():
    goose = os.environ.get("GOOSE", "goose")
    subprocess.run([sys.executable, "tests/fixtures.py", "tests/data/example.bam"],
                   cwd=ROOT, check=True)
    preamble = (ROOT / "tests/doc_preamble.goose").read_text()
    # Keep the preamble's relative root import valid, without changing sources.
    with tempfile.TemporaryDirectory(prefix=".checkfree-", dir=ROOT) as temp:
        directory = Path(temp)
        positive = directory / "positive.goose"
        negative = directory / "negative.goose"
        checked = directory / "checked.goose"
        positive.write_text(preamble + POSITIVE)
        negative.write_text(preamble + NEGATIVE)
        checked.write_text(preamble + NEGATIVE.replace("checkfree fn", "fn"))
        for opt in (0, 2):
            binary = directory / "probe"
            flags = ["--debug", f"-O{opt}", "--color", "never"]
            subprocess.run([goose, "build", str(positive), "--cc", "clang",
                            *flags, "-o", str(binary)], cwd=ROOT, check=True)
            subprocess.run([str(binary)], cwd=ROOT, check=True)
            rejected = subprocess.run([goose, "check", str(negative), *flags],
                                      cwd=ROOT, text=True, capture_output=True)
            diagnostic = rejected.stdout + rejected.stderr
            assert rejected.returncode != 0, diagnostic
            assert "checkfree fn unchecked_quality calls pileup::quality" in diagnostic, diagnostic
            assert "keeps a bounds check" in diagnostic, diagnostic
            subprocess.run([goose, "build", str(checked), "--cc", "clang",
                            *flags, "-o", str(binary)], cwd=ROOT, check=True)
            trapped = subprocess.run([str(binary)], cwd=ROOT, text=True, capture_output=True)
            assert trapped.returncode != 0, trapped
            assert "out of bounds" in trapped.stdout + trapped.stderr, trapped
        unverified = subprocess.run([goose, "check", str(positive), "--no-bce"],
                                    cwd=ROOT, text=True, capture_output=True)
        assert unverified.returncode == 0, unverified
        assert "checkfree is not verified under --no-bce" in unverified.stderr, unverified
    print("PASS checkfree callers, callee diagnostics, debug bounds checks (-O0/-O2), and --no-bce warning")


if __name__ == "__main__":
    main()
