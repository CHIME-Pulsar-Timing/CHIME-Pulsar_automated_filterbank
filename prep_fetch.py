import numpy as np
import sys
import os
import csv
def prep_fetch_csv(filfile,rank=1):
    #get spegid_python3 speg
    from SPEGID_Python3 import SinglePulseEventGroup
    spegs = np.load('spegs.npy',allow_pickle=1)
    #get only rank lower than the rank
    spegs = list([speg for speg in spegs if (speg.group_rank<=rank)&(speg.group_rank>0)])
    #create the subband 256 files
    #create_cands(spegs,256,filfile)
    #create the subband 128 files
    create_cands(spegs,128,filfile)

def create_cands(spegs,downsamp,filfile):
    with open('cands'+str(int(downsamp))+'.csv','w',newline='') as cands:
        writer=csv.writer(cands,delimiter=',')
        for speg in spegs:
            if speg.peak_SNR>5.5:
                #boxcar_w = np.around(np.log10(speg.peak_downfact)/np.log10(2))
                boxcar_w=0
                fn,peak_time=prep_fetch_scale_fil(filfile,speg.peak_time,float(speg.peak_DM),speg.peak_downfact,downsamp)
                #fetch takes log2 of the downfact
                writer.writerow([fn,speg.peak_SNR,peak_time,speg.peak_DM,boxcar_w,fn])

#copied from waterfaller.py
def maskfile(maskfn, data, start_bin, nbinsextra,extra_mask):    
    from presto import rfifind
    rfimask = rfifind.rfifind(maskfn)     
    mask = get_mask(rfimask, start_bin, nbinsextra)[::-1]    
    masked_chans = mask.all(axis=1)    
    # Mask data    
    if extra_mask:    
        masked_chans.append(extra_mask)    
    data = data.masked(mask, maskval='median-mid80')    
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

def prep_fetch_scale_fil(filfile,burst_time,dm,downsamp=32,subband=256):
    '''
    filfile: string input to filterbank filename
    filterbank_len: half the time length for filterbank file
    burst_time: time that the burst happens (optional)
    burst_sample: sample that the burst happens
    
    output
    filename: new string filename for the filterbank file created
    out_burst_time: time in the new filterbank file of the burst
    '''
    from presto.filterbank import FilterbankFile
    from presto import filterbank as fb
    from presto import rfifind
    #calculate the filterbank length required due to dispersion times 2 for plotting purposes
    filterbank_len=(4.15/1000)*dm*2
    fil = FilterbankFile(filfile,mode='read')
    tsamp = float(fil.header['tsamp'])
    burst_sample = burst_time/tsamp
    total_samples = fil.nspec
    nsamp = filterbank_len/tsamp
    if burst_sample<nsamp:
        #then there hasn't been enough time elapsed for this filterbank length
        nsamp=burst_sample
    if burst_sample+nsamp>total_samples:
        #we will over run
        nsamp=total_samples-burst_sample
    #reset new fb_len
    filterbank_len=nsamp*tsamp
    burst_sample=int(np.around(burst_sample))
    nsamp=int(np.around(nsamp))
    my_spec = fil.get_spectra(burst_sample-nsamp,nsamp*2)
    #mask the file
    maskfn = filfile.strip('.fil')+'_rfifind.mask'
    start_bin = burst_sample-nsamp
    nbinsextra = nsamp*2
    extra_mask=None
    data, masked_chans = maskfile(maskfn, my_spec, start_bin, nbinsextra,extra_mask)
    '''
    data_masked = np.ma.masked_array(data.data)
    data_masked[masked_chans] = np.ma.masked
    data.data = data_masked        
    '''
    #subband
    data.subband(subband,subdm=dm,padval='median')
    #downsample
    #find the highest value of power of 2
    data.downsample(int(downsamp))
    #may need to dedisperse

    my_spec = data
    
    my_spec = my_spec.scaled(False)
    #gotta resolve clipping issue
    #move all negative numbers to positive and scale if it's going to get clipped 
    my_spec.data = my_spec.data-np.min(my_spec.data)+1
    my_spec.data = my_spec.data*(255/np.max(my_spec.data))
    #modify the start time of the filterbank file
    fil.header['tstart'] = fil.header['tstart']+((burst_time-filterbank_len)/(60*60*24))
    fil.header['nchans'] = my_spec.numchans
    fil.header['tsamp'] = my_spec.dt
    fil.header['frequencies'] = my_spec.freqs
    fil.frequencies = my_spec.freqs

    #this is only because we are using int 8  bit, double check this!
    fil.bytes_per_spectrum = my_spec.numchans
    fil.nspec = my_spec.numspectra
    fil.dt = my_spec.dt
    fil.header['fch1'] = my_spec.freqs[0]
    fil.header['foff'] = np.diff(my_spec.freqs)[0]

    filename=filfile.rstrip('.fil')+'_'+str(float(burst_sample*tsamp))+'_sb_'+str(int(subband))+'.fil'
    fb.create_filterbank_file(filename,fil.header,spectra=my_spec.data.T,nbits=fil.header['nbits'])    
    return filename,filterbank_len
    
if __name__=='__main__':
    prep_fetch_csv(sys.argv[1])

 
