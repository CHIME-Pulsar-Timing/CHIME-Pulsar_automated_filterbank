#!/usr/bin/env bash
AFP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
for FIL in $@
do
    python $AFP/generate_rfi_mask.py $FIL &
done
