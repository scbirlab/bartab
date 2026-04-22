#!/usr/bin/env bash

set -euox pipefail

INPUT="test/inputs/test_count.csv"
SAMPLE_SHEET="test/inputs/test_sample_meta.csv"
STRAINS="test/inputs/test_strain_meta.csv"

bartab fit \
    "$INPUT" \
    --output "test/outputs/wls-results.h5ad" \
    --sample-sheet "$SAMPLE_SHEET" \
    --barcode-sheet "$STRAINS" \
    --barcode-column "strain_id" \
    --culture-column "replicate" \
    --reference "wt" \
    --spike-name "spike" \
    --use-spike

bartab plot \
    "test/outputs/wls-results.h5ad" \
    --highlight mut_A mut_B mut_C mut_D \
    -o test/outputs/wls-plot

>&2 echo "[$(date)] Done!"
