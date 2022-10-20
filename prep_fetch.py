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
    create_cands(spegs,128,filfile)

def create_cands(spegs,subband,filfile):
    fetch_len_1 = 1
    fetch_len_0 = 0.5
    i=0
    with open('cands'+str(int(subband))+'_'+str(fetch_len_1)+'.csv','w',newline='') as cands_1:
        with open('cands'+str(int(subband))+'_'+str(fetch_len_0)+'.csv','w',newline='') as cands:
            writer_0=csv.writer(cands,delimiter=',')
            writer_1=csv.writer(cands_1,delimiter=',')
            writer_0.writerow(["file","snr","width","dm","label","stime","chan_mask_path","num_files"])
            for speg in spegs:
                if speg.peak_SNR>5.5:
                    # define the width
                    #the chunks are min size of 128 samples, this means that if we are less than 128, just round up to 128
                    mint = speg.min_time
                    maxt = speg.max_time
                    dm = float(speg.peak_DM)
                    snr = float(speg.peak_SNR)
                    fn,tsamp,start,delay = prep_fetch_scale_fil(filfile,mint,maxt,dm)

                    #try to create a width befitting of the width
                    width_box_0 = fetch_len_0/tsamp
                    #for 1 seconds of width, this gets the really long bursts
                    width_box_1 = fetch_len_1/tsamp
                    def get_width(width_box,delay):
                        if delay > width_box:
                            width = 10
                        else:
                            width = width_box-delay
                        return width
                    width_0 = get_width(width_box_0,delay)
                    width_1 = get_width(width_box_1,delay)
                    #fetch takes log2 of the downfact
                    print(f"Writing cands file!")
                    print([fn,snr,int(np.around(np.log2(width_box_0))),speg.peak_DM,fn,start,"",1])
                    writer_0.writerow([fn,snr,int(np.around(np.log2(width_box_0))),speg.peak_DM,fn,start,"",1])
                    # writer_1.writerow([fn,speg.peak_SNR,start,speg.peak_DM,int(np.around(np.log2(width_1))),fn])
                    i+=1
                    if i>1:
                        break



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


def prep_fetch_scale_fil(filfile,min_burst_time,max_burst_time,dm,ts=1e-3):
    '''
    filfile: string input to filterbank filename
    filterbank_len: half the time length for filterbank file
    burst_time: time that the burst happens (optional)
    burst_sample: sample that the burst happens

    output
    filename: new string filename for the filterbank file created
    out_burst_time: time in the new filterbank file of the burst
    '''
    from presto import rfifind
    from presto import filterbank as fb
    from sigpyproc import readers as r
    from presto.psrfits import PsrfitsFile as p
    #calculate the filterbank length required due to dispersion times 2 for plotting purposes
    #start it a bit later to prevent errors
    # start_time = (4.15/1000)*(dm+50)*2

    fil = p(filfile)
    tsamp = float(fil.tsamp)
    downsamp = int(ts/tsamp)
    #give it a generous filterbank length. but you only need to create one in this case

    #get the number of samples at the bursts, i.e. how many bursts needed to get to sample
    burst_sample = int(np.around((max_burst_time+min_burst_time)/(2*tsamp)))
    total_samples = fil.nspec
    #set the filterbank length to the dispersion delay+5s?
    # try:
    #     fch1 = float(fil.header.fch1._to_value("MHz"))
    #     foff = float(fil.header.foff._to_value("MHz"))
    # except:
    #     fch1 = float(fil.header.fch1)
    #     foff = float(fil.header.foff)
    # nchans = fil.header.nchans
    # freqs = [fch1,fch1+foff*nchans]
    freqs = fil.freqs
    #get the number of samples at the bursts, i.e. how many bursts needed to get to sample
    # bt = (max_burst_time+min_burst_time)/2
    delay = dispersion_delay(dm,max(freqs),min(freqs))

    filterbank_len = delay + 5
    #get the spectra
    nsamp = int(np.around(filterbank_len/tsamp))
    if burst_sample<(nsamp/2):
        #then there hasn't been enough time elapsed for this filterbank length
        start_samp = 0
        end_samp = nsamp

    elif (burst_sample+nsamp/2)>total_samples:
        #we will over run
        end_samp = total_samples
        start_samp = total_samples-nsamp

    else:
        start_samp = int(np.around(burst_sample-nsamp/2))
        end_samp = int(np.around(burst_sample+nsamp/2))

    bt = (burst_sample-start_samp)*tsamp
    filterbank_len=nsamp*tsamp
    try:
        print("reading block")
        nsamp = nsamp-nsamp%downsamp
        my_spec = fil.get_spectra(int(start_samp),nsamp)
        print("read block")
    except:
        import pdb; pdb.set_trace()

    #mask the file
    # maskfn = filfile.strip('.fits')+'_rfifind.mask'
    # my_spec, masked_chans = maskfile(maskfn, my_spec, start_samp, nsamp)
    #subband
    my_spec.downsample(downsamp)
    #my_spec = my_spec.normalise()
    #may need to dedisperse
    #my_spec = my_spec - np.min(my_spec)+1
    #my_spec = np.array(my_spec*(255/np.max(my_spec)))
    #modify the start time of the filterbank file
    #filheader = translate_header(fil)
    #import pdb; pdb.set_trace()
    filename=filfile.rstrip('.fits')+'_'+str(float(burst_sample*tsamp))+'.fil'
    filheader = translate_header(fil)

    # my_spec.header.foff = foff
    # my_spec.header.fch1 = fch1
    filheader['tsamp']=my_spec.dt
    filheader['nbits']=8
    fb.create_filterbank_file(filename,filheader,spectra=my_spec.data,nbits=8)
    # my_spec.to_file(filename=filename)
    # print(my_spec.header.tsamp)
    #otherwise do nothing

    return filename,my_spec.dt,bt,delay

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
if __name__=='__main__':
    prep_fetch_csv(sys.argv[1],rank=5)
