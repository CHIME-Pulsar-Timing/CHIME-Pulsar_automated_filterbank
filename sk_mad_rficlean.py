import numpy as np
import sys, os
import argparse
from matplotlib import pyplot as plt
import presto.filterbank as fb
from scipy import signal

from rfi_excision_tools import get_divisors, mad, mad2std, detect_rfi_mad, detect_rfi_sk
import pipeline_config

def plot_bandpass_fig(fname,spectra):

    if spectra.ndim == 2:

        spectra[spectra == 0] = np.nan 
        bandpass = np.nanmean(spectra,axis=1)

    else:

        print('Data dimension not correct, exiting...')
        return

    print('Plotting bandpass.')
    plt.plot(bandpass)
    plt.savefig(fname+'_bandpass.png')

def sk_mad_rfi_excision(fname,fil,apply_sk=True,apply_mad=True,apply_chanthresh=False,plot_bandpass=False):

    outopts = ''

    badchans = np.asarray(pipeline_config.zaplist.split(','),dtype=int)
    badchans = 1023-badchans

    print('Loading filterbank file '+str(fil))
    filfile = fb.FilterbankFile(fil)
    tsamp = filfile.dt
    nsamp = filfile.nspec
    nchan = filfile.nchan
    div = get_divisors(nsamp)
    samp_fac = int(round(0.00098304/tsamp)) # sampling time relative to 0.98304 ms

    #create a mask array from the loaded filterbank file
    print('Creating masked array')
    data = filfile.get_spectra(0,nsamp)
    data = np.ma.array(data.data.astype(float), mask=False)
    nsamp = np.size(data,axis=1)
 
    #create a badchan mask
    print('Creating bad channel mask')
    badchan_mask = np.zeros_like(data.mask, dtype=bool)
    badchan_mask[badchans, :] = True
    data.mask = np.logical_or(data.mask, badchan_mask)

    #sk filter
    if apply_sk:
        outopts += '_sk'
        print('Creating sk filter mask')
        min_div1 = np.argwhere(div >= 4 * samp_fac * 1024).min()
        min_div2 = np.argwhere(div >= samp_fac * 1024).min()
        min_div3 = np.argwhere(div >= samp_fac * 128).min()

        mset = [div[min_div1],div[min_div2],div[min_div3]]
        tset = [x * np.sqrt(samp_fac) for x in [7.0,4.0,3.0]]
        mindex = 1  

        for m, t in zip(mset, tset):

            for ichunk, chunk in enumerate(np.split(data, nsamp // m, axis=1)):
            
                chunk_mask = detect_rfi_sk(chunk, thresh=t, plot=False, ffact=16)
                chunk.mask = np.logical_or(chunk.mask, chunk_mask)
           
                #data.mask[:,ichunk*m:(ichunk+1)*m] = np.logical_or(chunk_mask, data.mask[:,ichunk*m:(ichunk+1)*m])
 
            print('sk mask no. %d done'%mindex)
            mindex+=1

    #mad filter
    if apply_mad:
        outopts += '_mad'
        print('Creating mad filter mask')
        mad_mask = detect_rfi_mad(data, stat_window=64, thresh=4, tfact=1)
        data.mask = np.logical_or(data.mask, mad_mask)

    #channel threshold
    if apply_chanthresh:
        outopts += '_chanthresh'
        print("flagging channels where >75% of time steps are masked")
        channel_masked_fraction = np.sum(data.mask, axis=1) / float(nsamp)
        threshold_mask = np.zeros_like(channel_masked_fraction, dtype=bool)
        threshold_mask[channel_masked_fraction > 0.75] = True
        data.mask = np.logical_or(data.mask, threshold_mask[:, np.newaxis])

    #apply the mask
    print('Applying the full mask')
    #data = np.ma.filled(data, fill_value=np.nan)  # data is no longer a masked array

    print(data.count(), data.size)
    masked_frac = (1.0 - (float(data.count()) / float(data.size)))
    print('Amount of data masked = %.3f' %masked_frac)

    nsamp_per_s_idx = np.argwhere(div >= samp_fac * 4096).min()
    nsamp_per_s = div[nsamp_per_s_idx]
    print(nsamp_per_s)
    for chan in range(data.shape[0]):
        for j in range(data.shape[1] // nsamp_per_s):
            _local = slice(j * nsamp_per_s, (j + 1) * nsamp_per_s)
            local_median = np.ma.median(data[chan, _local])
            data[chan, _local] = np.ma.filled(data[chan, _local], fill_value=local_median)
            data[chan, _local] = signal.detrend(data[chan, _local], type='linear')

    #for j in range(nchan):
        # for each channel, get the intensity and replace masked values with the local median of unflagged data
        #nsamp_per_s_idx = np.argwhere(div >= samp_fac * 1024).min()
        #nsamp_per_s = div[nsamp_per_s_idx]
        #for k in range(nsamp // nsamp_per_s):
            #local_time_slice = slice(k * nsamp_per_s, (k + 1) * nsamp_per_s)
            #intensity = data[j, local_time_slice]
            
            #med_intensity = np.nanmedian(intensity)
    
            #if np.isnan(med_intensity):
                # if the entire channel is flagged, set the median to 0
                #med_intensity = 0

            # do the replacement
            #data[j, local_time_slice][np.isnan(intensity)] = med_intensity

            # finally, remove the median intensity for the channel (even if that means subtracting 0)
            #data[j, local_time_slice] = data[j, local_time_slice] - med_intensity

    if plot_bandpass:
        plot_bandpass_fig(fname+outopts,data)
    # rescaling the data
    data = ((data - np.min(data)) / (np.max(data)- np.min(data))) * 256
    data = data.astype(int)

    print(np.max(data))

    print('Saving masked filterbank file')
    fb.create_filterbank_file(str(fname)+str(outopts)+'.fil',filfile.header,spectra=data.data.T,nbits=filfile.header['nbits'])    
   
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--fil', nargs=1, help='Input filterbank file')
    parser.add_argument("--no_sk", help='Do not run sk filter', action='store_false')
    parser.add_argument("--no_mad", help='Do not run mad filter', action='store_false')
    parser.add_argument("--chanthresh", help="Apply a mask to channel with >75%% of data masked", action='store_true')
    parser.add_argument("--bandpass", help="Plot the bandpass after cleaning", action='store_true')

    args = parser.parse_args()

    fil = args.fil[0]
    fname = fil.split('.')[0]
    apply_sk = args.no_sk
    apply_mad = args.no_mad
    apply_chanthresh = args.chanthresh
    plot_bandpass = args.bandpass

    print(apply_sk,apply_mad,apply_chanthresh,plot_bandpass)

    if os.path.islink(fil):
   
        fil = os.readlink(fil)
        if not os.path.isfile(fil):

            print('File does not exist')
            sys.exit()
    
    sk_mad_rfi_excision(fname,fil,apply_sk,apply_mad,apply_chanthresh,plot_bandpass) 
