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
    with open('cands.csv','w',newline='') as cands:
        writer=csv.writer(cands,delimiter=',')
        for speg in spegs:

            boxcar_w = np.around(np.log10(speg.peak_downfact)/np.log10(2))
            fn,peak_time=prep_fetch_scale_fil(filfile,speg.peak_time)
            #fetch takes log2 of the downfact
            writer.writerow([fn,speg.peak_SNR,peak_time,speg.peak_DM,boxcar_w])

def prep_fetch_scale_fil(filfile,burst_time,filterbank_len=5):
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
    my_spec = my_spec.scaled(False)
    #gotta resolve clipping issue
    #move all negative numbers to positive and scale if it's going to get clipped
    my_spec.data = my_spec.data-np.min(my_spec.data)
    if np.max(my_spec.data)>255:
        my_spec.data = my_spec.data*(255/np.max(my_spec.data))
    '''
    print(np.max(my_spec.data))
    print(np.min(my_spec.data))
    import matplotlib.pyplot as plt
    plt.imshow(my_spec.data)
    plt.show()
    plt.plot(np.sum(my_spec.data,1))
    plt.figure()
    plt.plot(np.sum(my_spec.data,1))
    '''

    #modify the start time of the filterbank file
    fil.header['tstart'] = fil.header['tstart']+(burst_time/(60*60*24))
    filename=filfile.rstrip('.fil')+'_'+str(int(burst_sample*tsamp))+'.fil' 
    fb.create_filterbank_file(filename,fil.header,spectra=my_spec.data.T,nbits=fil.header['nbits'])    
    return filename,filterbank_len
    

#prep_fetch_csv(sys.argv[1])
 
