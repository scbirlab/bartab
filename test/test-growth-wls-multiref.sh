#!/usr/bin/env bash

set -euox pipefail

INPUT="test/inputs/test_count.csv"
SAMPLE_SHEET="test/inputs/test_sample_meta.csv"
STRAINS="test/inputs/test_strain_meta.csv"

bartab fit \
    "$INPUT" \
    --output "test/outputs/wls-growth-results-multiref.h5ad" \
    --sample-sheet "$SAMPLE_SHEET" \
    --barcode-sheet "$STRAINS" \
    --barcode-column "strain_id" \
    --culture-column "replicate" \
    --growth-column "growth" \
    --growth-type "density" \
    --reference "mut" \
    --spike-name "spike"

bartab plot \
    "test/outputs/wls-growth-results-multiref.h5ad" \
    --highlight mut_A mut_B mut_C mut_D \
    -o test/outputs/wls-growth-plot-multiref

>&2 echo "[$(date)] Done!"
