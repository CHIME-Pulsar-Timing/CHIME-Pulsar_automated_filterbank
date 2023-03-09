#!/usr/bin/env python3
import subprocess
from gwg_cand_search_pipeline import run_rfifind
from gwg_cand_search_pipeline import edit_mask
import sys


if __name__=="__main__":
    observation_file = sys.argv[1]
    fname = observation_file.strip(".fil")
    ext = '.fil'
    run_rfifind(fname,ext)
    mask_name = '_rfifind.mask'
    edit_mask(fname,ext,mask_name)
