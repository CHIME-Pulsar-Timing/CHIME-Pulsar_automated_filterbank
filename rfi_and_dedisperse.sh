#!/bin/sh

for file in "$@"; do
    SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

    echo $SCRIPT_DIR
    python $SCRIPT_DIR/gwg_cand_search_pipeline.py --sk_mask --kc_iqrm --rfifind --dedisp --dm 82.5 --one_dm --fil $file &
done
