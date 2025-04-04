#!/bin/sh

for file in "$@"; do
    python ~/Documents/CHIME-Pulsar_automated_filterbank/gwg_cand_search_pipeline.py --sk_mask --kc_iqrm --rfifind --dedisp --dm 82.5 --one_dm --fil $file &
done
