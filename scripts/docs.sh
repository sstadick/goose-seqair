#!/usr/bin/env bash
set -euo pipefail
exec "${GOOSE:-goose}" run "${GOOSE_DOC:-../goose-doc/main.goose}" -- \
  --title "goose-seqair API" --symbols api.txt --require-docs --require-examples \
  --output "${DOC_OUTPUT:-API.md}" "$@" bam.goose bam_index.goose pileup.goose bam_encode.goose bam_writer.goose sam.goose cram.goose sam_writer.goose cram_writer.goose
