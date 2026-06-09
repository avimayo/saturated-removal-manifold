#!/bin/bash
WORK=/home/labs/alon/avimayo/sr_fits
PY=/home/labs/alon/navehr/.conda/envs/srtools/bin/python
LINE=$(sed -n "${LSB_JOBINDEX}p" $WORK/zims_groups.txt)
CLS=$(echo "$LINE"    | awk -F'\t' '{print $1}')
SPECIES=$(echo "$LINE" | awk -F'\t' '{print $2}')
SEX=$(echo "$LINE"    | awk -F'\t' '{print $3}')
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
$PY $WORK/scripts/zims/fit_zims_one.py "$CLS" "$SPECIES" "$SEX"
