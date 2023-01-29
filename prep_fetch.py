import numpy as np
import sys
import os
import csv
from presto.psrfits import PsrfitsFile as p
from presto import psrfits
from presto import filterbank
from presto import sigproc

DM_CONST = 4149.377593360996  # dispersion constant


def dispersion_delay(dm, f1, f2):
    """Return DM delay in seconds"""
    return DM_CONST * dm * (1.0 / f2 ** 2 - 1.0 / f1 ** 2)

def prep_fetch_csv(filfile,rank=5):
    #get spegid_python3 speg
    from SPEGID_Python3 import SinglePulseEventGroup
    spegs = np.load('spegs.npy',allow_pickle=1)
    #get only rank lower than the rank
    spegs = list([speg for speg in spegs if (speg.group_rank<=rank)&(speg.group_rank>0)])
    create_cands(spegs,filfile)

def create_cands(spegs,filfile):
    with open('cands.csv','w',newline='') as cands:
        writer_0=csv.writer(cands,delimiter=',')
        #write header
        writer_0.writerow(["file","snr","width","dm","label","stime","chan_mask_path","num_files"])
        for speg in spegs:
            if speg.peak_SNR>5.5:
                # define the width
                #the chunks are min size of 128 samples, this means that if we are less than 128, just round up to 128
                mint = speg.min_time
                maxt = speg.max_time
                SNR = float(speg.peak_SNR)
                dm = speg.peak_DM
                #fetch uses time calc is  timestamp-width- dispersion delay ---> timestamp+width+dispersion delay
                fn,tsamp,start = prep_fetch_scale_fil(filfile,mint,maxt,dm)
                #length has to be at least 0.5s
                #fetch takes log2 of the downfact
                print(f"Writing cands file!")
                writer_0.writerow([fn,SNR,1,dm,fn,start,"",1])

#copied from waterfaller.py
def maskfile(maskfn, data, start_bin, nbinsextra):
    from presto import rfifind
    print('loading mask')
    rfimask = rfifind.rfifind(maskfn)
    print('getting mask')
    mask = get_mask(rfimask, start_bin, nbinsextra)[::-1]
    print('get mask finished')
    masked_chans = mask.all(axis=1)
    #mask the data but set to the mean of the channel
    mask_vals = np.median(data,axis=1)
    for i in range(len(mask_vals)):
        _ = data[i,:]
        _m = mask[i,:]
        _[_m] = mask_vals[i]
        data[i,:] = _
    return data, masked_chans


def get_mask(rfimask, startsamp, N):
    """Return an array of boolean values to act as a mask
        for a Spectra object.
    
        Inputs:                                                                    
            rfimask: An rfifind.rfifind object                         
            startsamp: Starting sample
            N: number of samples to read

        Output:
            mask: 2D numpy array of boolean values. 
                True represents an element that should be masked.
    """
    sampnums = np.arange(startsamp, startsamp+N)
    blocknums = np.floor(sampnums/rfimask.ptsperint).astype('int')
    mask = np.zeros((N, rfimask.nchan), dtype='bool')
    for blocknum in np.unique(blocknums):
        blockmask = np.zeros_like(mask[blocknums==blocknum])  
        chans_to_mask = rfimask.mask_zap_chans_per_int[blocknum]
        if chans_to_mask.any():
            blockmask[:,chans_to_mask] = True
        mask[blocknums==blocknum] = blockmask
    return mask.T


def prep_fetch_scale_fil(filfile,min_burst_time,max_burst_time,dm):
    '''
    filfile: string input to filterbank filename
    filterbank_len: half the time length for filterbank file
    burst_time: time that the burst happens (optional)
    burst_sample: sample that the burst happens
    
    output
    filename: new string filename for the filterbank file created
    out_burst_time: time in the new filterbank file of the burst
    '''
    try:
        from presto.psrfits import PsrfitsFile
        fil = PsrfitsFile(filfile)
        tsamp = fil.tsamp
    except:
        from presto.filterbank import FilterbankFile
        fil = FilterbankFile(filfile)
        tsamp = fil.header['tsamp']
    # tsamp = float(fil.header.tsamp)
    # try:
    #     fch1 = float(fil.header.fch1._to_value("MHz"))
    #     foff = float(fil.header.foff._to_value("MHz"))
    # except:
    #     fch1 = float(fil.header.fch1)
    #     foff = float(fil.header.foff)
    # nchans = fil.header.nchans
    #get the number of samples at the bursts, i.e. how many bursts needed to get to sample
    bt = (max_burst_time+min_burst_time)/2
    return filfile,tsamp,bt
    
if __name__=='__main__':
    prep_fetch_csv(sys.argv[1],rank=5)
