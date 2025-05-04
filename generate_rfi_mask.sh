#!/usr/bin/env bash
AFP="$(dirname $(readlink -f $0))"
for FIL in $@
do
    python $AFP/gwg_cand_search_pipeline.py --rfifind --sk_mask --kc_iqrm --fil $FIL &
done
