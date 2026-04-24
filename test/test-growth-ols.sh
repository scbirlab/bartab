#!/usr/bin/env bash

set -euox pipefail

INPUT="test/inputs/test_count.csv"
SAMPLE_SHEET="test/inputs/test_sample_meta.csv"
STRAINS="test/inputs/test_strain_meta.csv"

bartab fit \
    "$INPUT" \
    --output "test/outputs/ols-growth-results.h5ad" \
    --sample-sheet "$SAMPLE_SHEET" \
    --barcode-sheet "$STRAINS" \
    --barcode-column "strain_id" \
    --culture-column "replicate" \
    --growth-column "growth" \
    --growth-type "density" \
    --reference "wt" \
    --spike-name "spike" \
    --model-type OLS

bartab plot \
    "test/outputs/ols-growth-results.h5ad" \
    --model-type OLS \
    --highlight mut_A mut_B mut_C mut_D \
    -o test/outputs/ols-growth-plot

>&2 echo "[$(date)] Done!"
