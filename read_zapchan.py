#!/usr/bin/env python3
import sys
import numpy as np


def convert_to_presto(fn):
    #
    #fn is the file name to read
    #convert presto zapfile
    #
    with open(fn,'r') as zapfile:
        channels = zapfile.readlines()
        chans = (int(channel) for channel in channels)
        psrchive_chans = np.arange(0, 1024)[::-1]
        psrchive_chans = psrchive_chans[chans]
        chan_str = ",".join(psrchive_chans.astype(str))
    return chan_str


print(convert_to_presto(sys.argv[1]))
