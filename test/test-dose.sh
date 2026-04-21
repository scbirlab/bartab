#!/usr/bin/env bash

set -euox pipefail

INPUT="test/inputs/dose-response_count.csv"
SAMPLE_SHEET="test/inputs/dose-response_sample_meta.csv"
STRAINS="test/inputs/dose-response_strain_meta.csv"

bartab fit \
    "$INPUT" \
    --output "test/outputs/dose-response-results.h5ad" \
    --sample-sheet "$SAMPLE_SHEET" \
    --barcode-sheet "$STRAINS" \
    --barcode-column "strain_id" \
    --culture-column "replicate" \
    --concentration-column "dose" \
    --reference "wt" \
    --spike-name "spike" \
    --use-spike

bartab plot \
    "test/outputs/dose-response-results.h5ad" \
    --model-type "HillFitnessModel" \
    -o test/outputs/dose-response-plot

>&2 echo "[$(date)] Done!"
