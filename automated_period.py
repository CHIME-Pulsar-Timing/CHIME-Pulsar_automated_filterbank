#!/usr/bin/env python3
import argparse
import numpy as np
from presto.psr_utils import rrat_period_multiday, rrat_period

def get_burst_dict(csvname):
    """From a csv of burst information named csv_name, extract burst info and
    compile a dictionary of burst data. The first column of the csv should look
    something like:

    [dir to filterbank file]/J0740+17_cand_59236_pow/0//nsub_128/J0740+17[cont.]
    _cand_59236_pow_153.1887616_sb_128_dm_44.85_snr_5.76

    Returns:
    Each key of the dictionary will be the string representation of the MJD the
    pulses were observed on. Each of these entries will be an array with the
    rows storing relevant information for each pulse.

    [TOA in s from observation start, subbanding, dm, S/N]

    'mean' and 'std' will also be keys representing the mean and std of each
    pulse characteristic over all days
    """
    # Open the positive bursts file
    positive_bursts_csv = open(csvname)
    burst_lines = positive_bursts_csv.readlines()

    # Get all of the ids and put them into a dictionary
    mjd = []
    burst_time = []
    subband = []
    dm = []
    snr = []
    for line in burst_lines:
        burst_str = line.split(',')[0]
        burst_numbers = [float(num) for num in burst_str.split('_')\
                         if num.replace('.', '1').isdigit()]
        burst_info = burst_numbers[-5:]
        mjd.append(burst_info[0])
        burst_time.append(burst_info[1])
        subband.append(burst_info[2])
        dm.append(burst_info[3])
        snr.append(burst_info[4])
    return mjd,burst_time,subband,burst_info

def format_burst(mjd,burst_time,min_bursts=3,min_seperation=0):
    umjd = set(mjd)
    mjd = np.array(mjd)
    burst_time = np.array(burst_time)
    p_array = []
    for u in umjd:
        #check if there's more than 2 bursts for a specific mjd
        if np.sum(u==mjd)>=min_bursts:
            #get indexes
            ind = u==mjd
            bt = burst_time[ind]
            diff_bt = np.diff(np.sort(bt))
            for d in diff_bt:
                if d<min_seperation:
                    print(u,':',np.sort(bt))
            p_array.append(bt)

    return p_array

import sys
mjd,burst_time,subband,burst_info = get_burst_dict(sys.argv[1])
p_array = format_burst(mjd,burst_time,min_seperation=2)

rrat_period_multiday(p_array)

for p in p_array:
    #get period individually
    print(rrat_period(p))
