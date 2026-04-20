#!/usr/bin/env bash

set -euox pipefail

INPUT="test/inputs/test_counts.csv"
SAMPLE_SHEET="test/inputs/test_sample_meta.csv"
STRAINS="test/inputs/test_strain_meta.csv"

bartab fit \
    "$INPUT" \
    --output "test/outputs/wls-results" \
    --sample-sheet "$SAMPLE_SHEET" \
    --barcode-sheet "$STRAINS" \
    --barcode-column "strain_id" \
    --culture-column "replicate" \
    --reference "wt" \
    --spike-name "spike" \
    --use-spike 

>&2 echo "[$(date)] Done!"
