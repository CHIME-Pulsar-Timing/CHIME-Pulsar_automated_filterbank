import numpy as np
import sys
import os
import csv
def prep_fetch_csv(filfile,rank=2):
    #get spegid_python3 speg
    from SPEGID_Python3 import SinglePulseEventGroup
    spegs = np.load('spegs.npy',allow_pickle=1)
    #get only rank lower than the rank
    spegs = list([speg for speg in spegs if (speg.group_rank<=rank)&(speg.group_rank>0)])
    #create the subband 256 files
    #create_cands(spegs,256,filfile)
    #create the subband 128 files
    create_cands(spegs,128,filfile,1.5)
    create_cands(spegs,128,filfile,0)
    # create_cands(spegs,128,filfile,0.5)
    # create_cands(spegs,128,filfile,0.2)
def create_cands(spegs,subband,filfile,fb_len=0):
    with open('cands'+str(int(subband))+'_'+str(fb_len)+'.csv','w',newline='') as cands:
        writer=csv.writer(cands,delimiter=',')
        for speg in spegs:
            if speg.peak_SNR>5.5:
                #boxcar_w = np.around(np.log10(speg.peak_downfact)/np.log10(2))
                # boxcar_w=0
                fn,width,start=prep_fetch_scale_fil(filfile,speg.peak_time,float(speg.peak_DM),speg.peak_downfact,subband,downsamp=1,fb_len=fb_len)
                #fetch takes log2 of the downfact
                writer.writerow([fn,speg.peak_SNR,start,speg.peak_DM,np.log2(width),fn])

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

def prep_fetch_scale_fil(filfile,burst_time,dm,boxcar=32,subband=256,downsamp=1,fb_len=1):
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
    #start it a bit later to prevent errors
    # start_time = (4.15/1000)*(dm+50)*2

    fil = FilterbankFile(filfile,mode='read')
    tsamp = float(fil.header['tsamp'])
    if fb_len!=0:
        filterbank_len=((4.15/1000)*dm)+fb_len
    else:
        filterbank_len=((4.15/1000)*dm)+tsamp*20*boxcar
    burst_sample = burst_time/tsamp
    total_samples = fil.nspec
    #the downsamp is the amount of boxcar widths
    nsamp = (filterbank_len/tsamp)
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
    data, masked_chans = maskfile(maskfn, my_spec, start_bin, nbinsextra, extra_mask)
    #subband
    data.subband(subband,subdm=dm,padval='median')
    #add padding at start
    # my_spec_data = data.data
    # medians = np.median(my_spec_data,axis=1)
    # pad_samps = int(start_time/tsamp)+1
    # pad_data = np.tile(medians,(pad_samps,1)).T
    # new_data = np.concatenate((pad_data,my_spec_data),axis=1)
    # data.data = new_data
    # data.numspectra = new_data.shape[1]
    #downsample
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
    filename=filfile.rstrip('.fil')+'_'+str(float(burst_sample*tsamp))+'_sb_'+str(int(subband))+'_'+str(fb_len)+'.fil'
    fb.create_filterbank_file(filename,fil.header,spectra=my_spec.data.T,nbits=fil.header['nbits'])    
    # import pdb; pdb.set_trace()
    #the chunks are min size of 128 samples, this means that if we are less than 128, just round up to 128
    if (nsamp/downsamp) < 256:
        width = 2
    else:
        width = nsamp/downsamp/128
    #otherwise do nothing
    return filename,width,filterbank_len
    
if __name__=='__main__':
    prep_fetch_csv(sys.argv[1])

 
