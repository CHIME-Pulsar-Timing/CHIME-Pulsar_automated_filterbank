#!/usr/bin/env python3

import os
import sys
import numpy as np
import argparse
from shutil import copy as cp
#this script will copy all pngs within a certain dm range.
#This is useful for cases where you know the DM but you're not sure that FETCH is getting everything
#(wide bursts etc)
def parse_arguments():
    parser = argparse.ArgumentParser(description='Grab all bursts within a certain dm')
    parser.add_argument('-dm', '--DM',
                        action='store',
                        default=0,
                        type=float)
    parser.add_argument('-r', '--dmrange',
                        action='store',
                        default=5,
                        type=float)
    parser.add_argument('-d', '--dir',
                        action='store',
                        default="./",)
    parser.add_argument('-od', '--outdir',
                        default="narrow_dm",
                        action='store',)

    return parser.parse_args()

parser = parse_arguments()
dm = parser.DM
dmr = parser.dmrange
bdir = parser.dir
od = parser.outdir
dm_max = dm+dmr
dm_min = dm-dmr
if not os.path.isdir(f"{od}_nsub0"):
    os.mkdir(f"{od}_nsub0")
if not os.path.isdir(f"{od}_nsub1"):
    os.mkdir(f"{od}_nsub1")
if not os.path.isdir(f"{od}_nsub_short"):
    os.mkdir(f"{od}_nsub_short")

#list all dirs in bdir
for my_dir in os.listdir(bdir):
    #check if the nsubs exist
    if os.path.isdir(f"{bdir}/{my_dir}/nsub_128_0"):
        #if it does, grab all the filesnames
        nsub0 = f"{bdir}/{my_dir}/nsub_128_0"
        nsub1 = f"{bdir}/{my_dir}/nsub_128_1"
        nsub_short = f"{bdir}/{my_dir}/nsub_128_0_short"
        nsubs = [nsub0,nsub1,nsub_short]
        for nsub in nsubs:
            fn = os.listdir(nsub)
            #only get the pngs
            for plot in fn:
                if 'png' in plot:
                    png_fp = f"{nsub}/{plot}"
                    #so this is the plot
                    splits = np.array(plot.split('_'))
                    dm_ind = np.argwhere(splits=='dm')
                    dm_val = float(splits[dm_ind+1])
                    #filter dm_val
                    if (dm_val>dm_min)&(dm_val<dm_max):
                        print(f"copying {plot} with dm:{dm_val}")
                        if ("nsub_128_0" in nsub) and (not ("short" in nsub)) :
                            cp(png_fp,f"{od}_nsub0")
                        elif "nsub_128_1" in nsub:
                            cp(png_fp,f"{od}_nsub1")
                        #elif "nsub_128_0_short" in nsub:
                        #    cp(png_fp,f"{od}_nsub_short")


