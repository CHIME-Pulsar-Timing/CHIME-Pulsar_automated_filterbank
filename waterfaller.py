#!/usr/bin/env python

"""
waterfaller.py

Make waterfall plots to show frequency sweep of a single pulse.
Reads PSRFITS or SIGPROC filterbank format files.

Patrick Lazarus - Aug. 19, 2011
Paul Scholz - Nov 2015

"""

import optparse

import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np

from presto import psr_utils
from presto import rfifind
from presto import psrfits
from presto import filterbank


SWEEP_STYLES = ['r-', 'b-', 'g-', 'm-', 'c-']


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


def maskfile(maskfn, data, start_bin, nbinsextra,extra_mask):
    rfimask = rfifind.rfifind(maskfn) 
    mask = get_mask(rfimask, start_bin, nbinsextra)[::-1]
    masked_chans = mask.all(axis=1)
    # Mask data
    if extra_mask.any():
        masked_chans=np.append(masked_chans,extra_mask)
    data = data.masked(mask, maskval='median-mid80')

    #datacopy = copy.deepcopy(data)
    return data, masked_chans


def waterfall(rawdatafile, start, duration, dm=None, nbins=None, nsub=None,\
              subdm=None, zerodm=False, downsamp=1, scaleindep=False,\
              width_bins=1, mask=False, maskfn=None,extra_mask=None, bandpass_corr=False,
              ref_freq=None,freq_mask=None):
    """
    Create a waterfall plot (i.e. dynamic specrum) from a raw data file.
    Inputs:
       rawdatafile - a PsrfitsData instance.
       start - start time of the data to be read in for waterfalling.
       duration - duration of data to be waterfalled.
    Optional Inputs:
       dm - DM to use when dedispersing data.
             Default: Don't de-disperse
       nbins - Number of time bins to plot. This option overrides
                the duration argument. 
                Default: determine nbins from duration.
       nsub - Number of subbands to use. Must be a factor of number of channels.
               Default: Number of channels.
       subdm - DM to use when subbanding. Default: same as dm argument.
       zerodm - subtract mean of each time-sample from data before 
                 de-dispersing.
       downsamp - Factor to downsample in time by. Default: Don't downsample.
       scaleindep - Scale each channel independently.
                     Default: Scale using global maximum.
       width_bins - Smooth each channel/subband with a boxcar width_bins wide.
                     Default: Don't smooth.
       maskfn - Filename of RFIFIND mask to use for masking data.
                 Default: Don't mask data.
       bandpass_corr - Correct for the bandpass. Requires an rfifind
                        mask provided by maskfn keyword argument.
                        Default: Do not remove bandpass.
       ref_freq - Reference frequency to de-disperse to. 
                   If subbanding and de-dispersing the start time 
                   will be corrected to account for change in
                   reference frequency. 
                   Default: Frequency of top channel.
    Outputs:
       data - Spectra instance of waterfalled data cube.
       nbinsextra - number of time bins read in from raw data. 
       nbins - number of bins in duration.
       start - corrected start time. 
    """

    if subdm is None:
        subdm = dm

    # Read data
    if ref_freq is None:
        ref_freq = rawdatafile.freqs.max()

    if nsub and dm:
        df = rawdatafile.freqs[1] - rawdatafile.freqs[0]
        nchan_per_sub = rawdatafile.nchan/nsub
        top_ctrfreq = rawdatafile.freqs.max() - \
                      0.5*nchan_per_sub*df # center of top subband
        start += 4.15e3 * np.abs(1./ref_freq**2 - 1./top_ctrfreq**2) * dm

    start_bin = np.round(start/rawdatafile.tsamp).astype('int')
    dmfac = 4.15e3 * np.abs(1./rawdatafile.frequencies[0]**2 - 1./rawdatafile.frequencies[-1]**2)

    if nbins is None:
        nbins = np.round(duration/rawdatafile.tsamp).astype('int')

    if dm:
        nbinsextra = np.round((duration + dmfac * dm)/rawdatafile.tsamp).astype('int')
    else:
        nbinsextra = nbins

    # If at end of observation
    if (start_bin + nbinsextra) > rawdatafile.nspec-1:
        nbinsextra = rawdatafile.nspec-1-start_bin

    data = rawdatafile.get_spectra(start_bin, nbinsextra)

    if freq_mask:
        freq_lims = freq_mask.split(',')
        upper = float(freq_lims[1])
        lower = float(freq_lims[0])
        freqs = rawdatafile.frequencies
        extra_mask = np.squeeze(np.where((freqs>lower)&(freqs<upper)))
        print(np.squeeze(extra_mask))
    # Masking
    if mask and maskfn:
        data, masked_chans = maskfile(maskfn, data, start_bin, nbinsextra,extra_mask)
    else:
        masked_chans = np.zeros(rawdatafile.nchan,dtype=bool)

    # Bandpass correction
    if maskfn and bandpass_corr:
        bandpass = rfifind.rfifind(maskfn).bandpass_avg[::-1]
        #bandpass[bandpass == 0] = np.min(bandpass[np.nonzero(bandpass)])
        masked_chans[bandpass == 0] = True

        # ignore top and bottom 1% of band
        ignore_chans = np.ceil(0.01*rawdatafile.nchan)
        masked_chans[:ignore_chans] = True
        masked_chans[-ignore_chans:] = True


    data_masked = np.ma.masked_array(data.data)
    data_masked[masked_chans] = np.ma.masked
    data.data = data_masked

    if bandpass_corr:
       data.data /= bandpass[:, None]

    # Zerodm filtering
    if (zerodm == True):
        data.data -=  data.data.mean(axis=0)


    # Subband data
    if (nsub is not None) and (subdm is not None):
        data.subband(nsub, subdm, padval='mean')

    # Dedisperse
    if dm:
        data.dedisperse(dm, padval='mean')

    # Downsample
    data.downsample(downsamp)

    # scale data
    data = data.scaled(scaleindep)

    # Smooth
    if width_bins > 1:
        data.smooth(width_bins, padval='mean')

    return data, nbinsextra, nbins, start


def plot_waterfall(data, start, duration,
                   integrate_ts=False, integrate_spec=False, show_cb=False,
                   cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
                   ax_im=None, ax_ts=None, ax_spec=None, interactive=True,interactive_masking=False):
    """ I want a docstring too!
    """

    # Set up axes
    if interactive:
        fig = plt.figure()
        # fig.canvas.set_window_title("Frequency vs. Time")

    im_width = 0.6 if integrate_spec else 0.8
    im_height = 0.6 if integrate_ts else 0.8

    if not ax_im:
        ax_im = plt.axes((0.15, 0.15, im_width, im_height))
    if integrate_ts and not ax_ts:
        ax_ts = plt.axes((0.15, 0.75, im_width, 0.2),sharex=ax_im)

    if integrate_spec and not ax_spec:
        ax_spec = plt.axes((0.75, 0.15, 0.2, im_height),sharey=ax_im)

    # Ploting it up
    nbinlim = np.int(duration/data.dt)
    data_vals = data.data[..., :nbinlim]
    if interactive_masking:
        img = ax_im.imshow(data_vals,aspect = "auto",cmap='YlGnBu')
    else:
        img = ax_im.imshow(data_vals, aspect='auto',
                    cmap='YlGnBu',
                    interpolation='nearest', origin='upper',
                    extent=(data.starttime, data.starttime+ nbinlim*data.dt,
                            data.freqs.min(), data.freqs.max()))
    if show_cb:
        cb = ax_im.get_figure().colorbar(img)
        cb.set_label("Scaled signal intensity (arbitrary units)")

    #plt.axis('tight')
    # Sweeping it up
    for ii, sweep_dm in enumerate(sweep_dms):
        ddm = sweep_dm-data.dm
        delays = psr_utils.delay_from_DM(ddm, data.freqs)
        delays -= delays.min()

        if sweep_posns is None:
            sweep_posn = 0.0
        elif len(sweep_posns) == 1:
            sweep_posn = sweep_posns[0]
        else:
            sweep_posn = sweep_posns[ii]
        sweepstart = data.dt*data.numspectra*sweep_posn+data.starttime
        sty = SWEEP_STYLES[ii%len(SWEEP_STYLES)]
        ax_im.plot(delays+sweepstart, data.freqs, sty, lw=4, alpha=0.5)

    # Dressing it up
    ax_im.xaxis.get_major_formatter().set_useOffset(False)
    ax_im.set_xlabel("Time")
    ax_im.set_ylabel("Observing frequency (MHz)")

    # Plot Time series
    if integrate_ts:
        Data = np.array(data.data[..., :nbinlim])
        Dedisp_ts = Data.sum(axis=0)
        times = (np.arange(data.numspectra)*data.dt + start)[..., :nbinlim]
        ax_ts.plot(times, Dedisp_ts,"k")
        ax_ts.set_xlim([times.min(),times.max()])
        plt.setp(ax_ts.get_xticklabels(), visible = False)
        plt.setp(ax_ts.get_yticklabels(), visible = False)

    # Plot Spectrum                                                             
    if integrate_spec:                                                         
        spectrum_window = 0.05*duration
        window_width = int(spectrum_window/data.dt) # bins
        burst_bin = nbinlim//2
        on_spec = np.array(data.data[..., burst_bin-window_width:burst_bin+window_width])
        Dedisp_spec = on_spec.sum(axis=1)[::-1]                                 
                                                                                
        freqs = np.linspace(data.freqs.min(), data.freqs.max(), len(Dedisp_spec))           
        ax_spec.plot(Dedisp_spec,freqs,"k")                                       
        plt.setp(ax_spec.get_xticklabels(), visible = False)                   
        plt.setp(ax_spec.get_yticklabels(), visible = False)                    
        ax_spec.set_ylim([data.freqs.min(),data.freqs.max()])                   
        if integrate_ts:
            ax_ts.axvline(times[burst_bin]-spectrum_window,ls="--",c="grey")                  
            ax_ts.axvline(times[burst_bin]+spectrum_window,ls="--",c="grey")                  

    if interactive:
        global masked_channels
        masked_channels = []
        data_vals = data.data[..., :nbinlim]
        def key_press(event):
            global start_mask, end_mask, masked_channels
            # check if click is within plot boundaries
            if event.key == "n":
                if event.inaxes is not None:
                    # get y value at clicked location
                    y_val = event.ydata
                    #get evenly spaced y values between fch1 and fch1 + nchan*foff
                    y_vals = np.arange(data_vals.shape[0])
                    print(y_val,y_vals)
                    #figure out which y value is closest to the clicked y value
                    closest_y_val = np.argmin(np.abs(y_vals - y_val))
                    start_mask = closest_y_val
                    print(f"Start mask at {closest_y_val}")
            elif event.key == "m":
                if event.inaxes is not None:
                    # get y value at clicked location
                    y_val = event.ydata
                    #get evenly spaced y values between fch1 and fch1 + nchan*foff
                    y_vals = np.arange(data_vals.shape[0])
                    #figure out which y value is closest to the clicked y value
                    closest_y_val = np.argmin(np.abs(y_vals - y_val))
                    end_mask = closest_y_val
                    print(f"End mask at {closest_y_val}")
                if start_mask != -1 and end_mask != -1:
                    start_ = min([start_mask,end_mask])
                    end_ = max([start_mask,end_mask])
                    to_mask = np.arange(start_,end_)
                    print(f"Replacing rows {start_} to {end_} with median")
                    data_vals[to_mask,:] = np.median(data_vals)
                    ax_im.clear()
                    ax_im.imshow(data_vals,aspect = "auto",cmap='YlGnBu')
                    #plot the new data
                    plt.gcf().canvas.draw_idle()
                    start_mask = -1
                    end_mask = -1
                    masked_channels += list(to_mask)
                    print_mask_chans = 1023-np.array(masked_channels)
                    #print masked_channels as comma seperated values
                    csv_mask_chans = ",".join([str(x) for x in print_mask_chans])
                    print(f"Masked channels: {csv_mask_chans}")


        fig.suptitle("Frequency vs. Time")
        fig.canvas.mpl_connect('key_press_event',
                lambda ev: (ev.key in ('q','Q') and plt.close(fig)))
        cid = fig.canvas.mpl_connect('key_press_event', key_press)
        plt.show()

def main():
    fn = args[0]

    if fn.endswith(".fil"):
        # Filterbank file
        filetype = "filterbank"
        rawdatafile = filterbank.FilterbankFile(fn)
    elif fn.endswith(".fits"):
        # PSRFITS file
        filetype = "psrfits"
        rawdatafile = psrfits.PsrfitsFile(fn)
    else:
        raise ValueError("Cannot recognize data file type from "
                         "extension. (Only '.fits' and '.fil' "
                         "are supported.)")

    data, bins, nbins, start = waterfall(rawdatafile, options.start,
                            options.duration, dm=options.dm,
                            nbins=options.nbins, nsub=options.nsub,
                            subdm=options.subdm, zerodm=options.zerodm,
                            downsamp=options.downsamp,
                            scaleindep=options.scaleindep,
                            width_bins=options.width_bins, mask=options.mask,
                            maskfn=options.maskfile,
                            extra_mask=options.extra_mask,
                            bandpass_corr=options.bandpass_corr,
                            freq_mask=options.freq_mask)

    plot_waterfall(data, start, options.duration, integrate_ts=options.integrate_ts,
                   integrate_spec=options.integrate_spec, show_cb=options.show_cb, 
                   cmap_str=options.cmap, sweep_dms=options.sweep_dms, 
                   sweep_posns=options.sweep_posns,interactive_masking=options.interactive_masking)

if __name__=='__main__':
    parser = optparse.OptionParser(prog="waterfaller.py",
                        version="v0.9 Patrick Lazarus (Aug. 19, 2011)",
                        usage="%prog [OPTIONS] INFILE",
                        description="Create a waterfall plot to show the "
                                    "frequency sweep of a single pulse "
                                    "in psrFits data.")
    parser.add_option('--subdm', dest='subdm', type='float',
                        help="DM to use when subbanding. (Default: "
                                "same as --dm)", default=None)
    parser.add_option('--zerodm', dest='zerodm', action='store_true',
                        help="If this flag is set - Turn Zerodm filter - ON  (Default: "
                                "OFF)", default=False)
    parser.add_option('-s', '--nsub', dest='nsub', type='int',
                        help="Number of subbands to use. Must be a factor "
                                "of number of channels. (Default: "
                                "number of channels)", default=None)
    parser.add_option('-d', '--dm', dest='dm', type='float',
                        help="DM to use when dedispersing data for plot. "
                                "(Default: 0 pc/cm^3)", default=0.0)
    parser.add_option('--show-ts', dest='integrate_ts', action='store_true',
                        help="Plot the time series. "
                                "(Default: Do not show the time series)", default=False)
    parser.add_option('--show-spec', dest='integrate_spec', action='store_true',
                        help="Plot the spectrum. "
                                "(Default: Do not show the spectrum)", default=False)
    parser.add_option('--bandpass', dest='bandpass_corr', action='store_true',
                        help="Correct for the bandpass. Requires an rfifind "
                                "mask provided by --mask option."
                                "(Default: Do not remove bandpass)", default=False)
    parser.add_option('-T', '--start-time', dest='start', type='float',
                        help="Time into observation (in seconds) at which "
                                "to start plot.")
    parser.add_option('-t', '--duration', dest='duration', type='float',
                        help="Duration (in seconds) of plot.")
    parser.add_option('-n', '--nbins', dest='nbins', type='int',
                        help="Number of time bins to plot. This option takes "
                                "precedence over -t/--duration if both are "
                                "provided.")
    parser.add_option('--width-bins', dest='width_bins', type='int',
                        help="Smooth each channel/subband with a boxcar "
                                "this many bins wide. (Default: Don't smooth)",
                        default=1)
    parser.add_option('--sweep-dm', dest='sweep_dms', type='float',
                        action='append',
                        help="Show the frequency sweep using this DM. "
                                "(Default: Don't show sweep)", default=[])
    parser.add_option('--sweep-posn', dest='sweep_posns', type='float',
                        action='append',
                        help="Show the frequency sweep at this position. "
                                "The position refers to the high-frequency "
                                "edge of the plot. Also, the position should "
                                "be a number between 0 and 1, where 0 is the "
                                "left edge of the plot. "
                                "(Default: 0)", default=None)
    parser.add_option('--downsamp', dest='downsamp', type='int',
                        help="Factor to downsample data by. (Default: 1).",
                        default=1)
    parser.add_option('--maskfile', dest='maskfile', type='string',
                        help="Mask file produced by rfifind. Used for "
                             "masking and bandpass correction.",
                        default=None)
    parser.add_option('--extramask', dest='extra_mask', type='string',
                        help="The extra channels you want masked, on top of the rfifind mask files",
                        default=np.array([]))

    parser.add_option('--mask', dest='mask', action="store_true",
                        help="Mask data using rfifind mask (Default: Don't mask).",
                        default=False)
    parser.add_option('--scaleindep', dest='scaleindep', action='store_true',
                        help="If this flag is set scale each channel "
                                "independently. (Default: Scale using "
                                "global maximum.)",
                        default=False)
    parser.add_option('--show-colour-bar', dest='show_cb', action='store_true',
                        help="If this flag is set show a colour bar. "
                                "(Default: No colour bar.)",
                        default=False)
    parser.add_option('--colour-map', dest='cmap',
                        help="The name of a valid matplotlib colour map."
                                "(Default: gist_yarg.)",
                        default='gist_yarg')
    parser.add_option('--freq_mask', dest='freq_mask',
                        help="mask out certain frequencies",
                        default=None)
    parser.add_option('--interactive_masking', dest='interactive_masking',
                        help="mask out certain frequencies interactively",
                        action='store_true', default=False)
    options, args = parser.parse_args()
    
    if not hasattr(options, 'start'):
        raise ValueError("Start time (-T/--start-time) "
                            "must be given on command line!")
    if (not hasattr(options, 'duration')) and (not hasattr(options, 'nbins')):
        raise ValueError("One of duration (-t/--duration) "
                            "and num bins (-n/--nbins)"
                            "must be given on command line!")
    if options.subdm is None:
        options.subdm = options.dm
    main()
