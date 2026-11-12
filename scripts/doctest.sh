#!/usr/bin/env bash
set -euo pipefail
python3 tests/fixtures.py tests/data/example.bam
exec python3 "$(dirname "${GOOSE_DOC:-../goose-doc/main.goose}")/doctest.py" \
  --cc clang --symbols api.txt --preamble tests/doc_preamble.goose --require-examples \
  "$@" bam.goose bam_index.goose pileup.goose bam_encode.goose bam_writer.goose sam.goose cram.goose sam_writer.goose cram_writer.goose
