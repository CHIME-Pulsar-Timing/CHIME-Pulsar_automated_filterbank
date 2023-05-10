#!/usr/bin/env python3

import numpy as np
from sigpyproc.readers import PFITSReader
import os
import argparse
from presto.filterbank import FilterbankFile, create_filterbank_file

parser = argparse.ArgumentParser()
parser.add_argument('-fil', type=str, help='Input filterbank file')
args = parser.parse_args()
fname = args.fil
fname_base = fname.replace(".fil","")
chunk_size = 1*1024*1024*1024 #Bytes each (first number is gigabytes)

filesize = os.path.getsize(fname)
total_files = np.ceil(filesize/chunk_size)
print("total files:",total_files)
#read the file
filfile = PFITSReader(fname)
nsubints = filfile.sub_hdr.nsubint
gulp_subints = nsubints//total_files
nsamps = filfile.header.nsamples
#gulp size
#this is a fits file, so it's split but subints.
subint_samples = filfile.sub_hdr.subint_samples
gulp = int(gulp_subints*subint_samples)

print("gulp size:",gulp)
current_samp = 0
i=0
filfile.header.foff=filfile.header.foff.value
filfile.header.fch1=filfile.header.fch1.value
#get the presto header
# filfile.header.azimuth=filfile.header.azimuth.value
# filfile.header.zenith=filfile.header.zenith.value
while current_samp<nsamps:
    if (nsamps-current_samp)<gulp:
        #set gulp to end of range
        gulp = nsamps-current_samp
    print("current sample",current_samp,"iteration",i)
    block = filfile.read_block(current_samp,gulp)
    # print(block.shape)
    out_fn = fname_base+f"_split_{i}.fil"
    block.to_file(filename=out_fn)
    #conver the split to presto format

    i += 1
    current_samp += gulp
