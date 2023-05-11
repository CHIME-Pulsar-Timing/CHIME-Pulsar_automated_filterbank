#!/usr/bin/env bash
AFP="$(dirname $(readlink -f $0))"
for FIL in $@
do
    python $AFP/generate_rfi_mask.py $FIL &
done
