#!/usr/bin/env python3
"""Run native rANS checks on the current CPU; require an ISA to prevent CI false coverage."""
import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-isa", choices=("auto", "avx2", "avx512", "neon"), default="auto")
    parser.add_argument("--sanitize", action="store_true")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="seqair-native-simd-") as temporary:
        for name, extra in (("scalar", ["-DGOOSIMD_SCALAR"]), ("native", [])):
            binary = Path(temporary) / name
            flags = ["-O3", "-march=native", "-Wall", "-Wextra"]
            if args.sanitize:
                flags += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
            libraries = ["-pthread"] + (["-ldl"] if sys.platform.startswith("linux") else [])
            subprocess.run(["clang", *flags, *extra, str(ROOT / "tests/native_simd.c"),
                            str(ROOT / "src/cram.c"), "-o", str(binary), *libraries], check=True)
            subprocess.run([str(binary), "scalar" if name == "scalar" else args.require_isa], check=True)


if __name__ == "__main__":
    main()
