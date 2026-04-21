#/usr/bin/env bash

set -euox pipefail

N_CONTROLS=3 #0
N_MORE=0 #1000
TIME=5
TIMEPOINTS=10
INOCULUM=10000
K=10
N=3
SEQ_DEPTH=10000

bartab sim \
    "test/inputs/fitness.json" \
    --generate-controls "$N_CONTROLS" \
    --generate-more "$N_MORE" \
    --seed 42 \
    --max-time "$TIME" \
    --timepoints "$TIMEPOINTS" \
    --inoculum "$INOCULUM" \
    --carrying-capacity "$K" \
    --n-cultures "$N" \
    --reads-per-barcode "$SEQ_DEPTH" \
    --output "test/inputs/test"

bartab sim \
    "test/inputs/ic50s.json" \
    --generate-controls "$N_CONTROLS" \
    --generate-more "$N_MORE" \
    --seed 42 \
    --max-time "$TIME" \
    --timepoints "$TIMEPOINTS" \
    --inoculum "$INOCULUM" \
    --carrying-capacity "$K" \
    --n-cultures "$N" \
    --reads-per-barcode "$SEQ_DEPTH" \
    --n-dose 10 \
    --dose-max 50 \
    --dose-fold 2 \
    --output "test/inputs/dose-response"
