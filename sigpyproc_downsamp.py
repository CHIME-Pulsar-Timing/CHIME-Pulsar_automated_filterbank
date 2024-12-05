#!/usr/bin/env python3

import numpy as np
from sigpyproc.readers import PFITSReader
import os
import argparse
from presto.filterbank import FilterbankFile, create_filterbank_file
from presto import filterbank as fb
from presto import sigproc
from presto import psrfits
def translate_header(psrfits_file):
    fits_hdr = psrfits_file.header
    subint_hdr = psrfits_file.fits['SUBINT'].header
    subint_data = psrfits_file.fits['SUBINT'].data
    fil_header = {}

    if fits_hdr['TELESCOP'] in sigproc.telescope_ids:
        fil_header["telescope_id"] = \
                    sigproc.telescope_ids[fits_hdr['TELESCOP']]
    else:
        fil_header["telescope_id"] = -1
    if fits_hdr['BACKEND'] in sigproc.machine_ids:
        fil_header["machine_id"] = \
                    sigproc.machine_ids[fits_hdr['BACKEND']]
    else:
        fil_header["machine_id"] = -1

    fil_header["data_type"] = 1 # filterbank
    fn = psrfits_file.filename
    fil_header["rawdatafile"] = os.path.basename(fn)
    fil_header["source_name"] = fits_hdr['SRC_NAME']
    fil_header["barycentric"] = 0 # always not barycentered?
    fil_header["pulsarcentric"] = 0 # whats pulsarcentric?
    fil_header["az_start"] = subint_data[0]['TEL_AZ']
    fil_header["za_start"] = subint_data[0]['TEL_ZEN']
    fil_header["src_raj"] = float(fits_hdr['RA'].replace(':',''))
    fil_header["src_dej"] = float(fits_hdr['DEC'].replace(':',''))
    fil_header["tstart"] = fits_hdr['STT_IMJD'] + \
                           fits_hdr['STT_SMJD']/86400.0 + \
                           fits_hdr['STT_OFFS']/86400.0
    fil_header["tsamp"] = subint_hdr['TBIN']
    fil_header["nbits"] = None # set by user. Input should always be 4-bit.

    # first channel (fch1) in sigproc is the highest freq
    # foff is negative to signify this
    fil_header["fch1"] = fits_hdr['OBSFREQ'] + \
                         np.abs(fits_hdr['OBSBW'])/2.0 - \
                         np.abs(subint_hdr['CHAN_BW'])/2.0
    fil_header["foff"] = -1.0*np.abs(subint_hdr['CHAN_BW'])
    fil_header["nchans"] = subint_hdr['NCHAN']
    fil_header["nifs"] = subint_hdr['NPOL']

    return fil_header

parser = argparse.ArgumentParser()
parser.add_argument('-fil', type=str, help='Input filterbank file')
parser.add_argument('-downsamp', type=int, help='Downsampling factor')
args = parser.parse_args()
fname = args.fil
downsamp = args.downsamp
if fname.endswith(".fits"):
    fname_base = fname.replace(".fits","")
else:
    fname_base = fname.replace(".fil","")
chunk_size = 1024*1024*1024 #Bytes each (first number is gigabytes)

filesize = os.path.getsize(fname)
total_files = np.ceil(filesize/chunk_size)
print("total files:",total_files)
#read the file
fitsfile = PFITSReader(fname)
nsubints = fitsfile.sub_hdr.nsubint
gulp_subints = nsubints//total_files
nsamps = fitsfile.header.nsamples
#gulp size
#this is a fits file, so it's split but subints.
subint_samples = fitsfile.sub_hdr.subint_samples

gulp = int(gulp_subints*subint_samples)
gulp = gulp - (gulp % downsamp)
print("gulp size:",gulp)
current_samp = 0
i=0
#get the presto header
# filfile.header.azimuth=filfile.header.azimuth.value
# filfile.header.zenith=filfile.header.zenith.value
psrfits_file = psrfits.PsrfitsFile(fname)
fil_header = translate_header(psrfits_file)
fil_header['tsamp'] = fil_header['tsamp']*downsamp
fil_header['nbits'] = 32
my_filfile = create_filterbank_file(fname_base+".fil", fil_header, nbits=32)
while current_samp<nsamps:
    if (nsamps-current_samp)<gulp:
        #set gulp to end of range
        gulp = nsamps-current_samp
    print("current sample",current_samp,"iteration",i)
    block = fitsfile.read_block(current_samp,gulp)
    # print(block.shape)
    #conver the split to presto format
    #downsample the block
    block = block.downsample(tfactor=downsamp)
    my_filfile.append_spectra(block.T)




    i += 1
    current_samp += gulp
