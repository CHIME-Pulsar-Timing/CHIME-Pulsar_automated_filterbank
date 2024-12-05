from presto import FilterbankFile
import numpy as np
import matplotlib.pyplot as plt

def interactive_mask(freq_time,title = "Interactive Masking"):
    freq_time -= np.median(freq_time,axis=0)
    freq_time /= np.std(freq_time,axis=0)
    freq_time[np.isnan(freq_time)] = 0
    clip_percentiles = (2, 98)
    low_pc = np.percentile(freq_time.flatten(), clip_percentiles[0])
    high_pc = np.percentile(freq_time.flatten(), clip_percentiles[1])
    freq_time = freq_time.clip(low_pc,high_pc)

    fig,ax = plt.subplots()
    ax.imshow(freq_time.T,aspect = "auto",cmap='YlGnBu')
    ax.set_title(title)
    global masked_channels
    masked_channels = []
    def key_press(event):
        global start_mask, end_mask, masked_channels
        # check if click is within plot boundaries
        if event.key == "n":
            if event.inaxes is not None:
                # get y value at clicked location
                y_val = event.ydata
                #get evenly spaced y values between fch1 and fch1 + nchan*foff
                y_vals = np.arange(freq_time.shape[0])
                #figure out which y value is closest to the clicked y value
                closest_y_val = np.argmin(np.abs(y_vals - y_val))
                start_mask = closest_y_val
                print(f"Start mask at {closest_y_val}")
        elif event.key == "m":
            if event.inaxes is not None:
                # get y value at clicked location
                y_val = event.ydata
                #get evenly spaced y values between fch1 and fch1 + nchan*foff
                y_vals = np.arange(freq_time.shape[0])
                #figure out which y value is closest to the clicked y value
                closest_y_val = np.argmin(np.abs(y_vals - y_val))
                end_mask = closest_y_val
                print(f"End mask at {closest_y_val}")
        if start_mask != -1 and end_mask != -1:
            start_ = min([start_mask,end_mask])
            end_ = max([start_mask,end_mask])
            to_mask = np.arange(start_,end_)
            print(f"Replacing rows {start_} to {end_} with median")
            freq_time[:,to_mask] = np.median(freq_time)
            ax.clear()
            ax.imshow(freq_time.T,aspect = "auto",cmap='YlGnBu')
            #plot the new freq_time
            plt.gcf().canvas.draw_idle()
            start_mask = -1
            end_mask = -1
            masked_channels += list(to_mask)
            print(f"Masked channels: {masked_channels}")


    cid = fig.canvas.mpl_connect('key_press_event', key_press)
    plt.show()
    return freq_time,masked_channels

def maskfile(maskfn, start_bin, nbinsextra):
    from presto import rfifind
    print('loading mask')
    rfimask = rfifind.rfifind(maskfn)
    print('getting mask')
    mask = get_mask(rfimask, start_bin, nbinsextra)[::-1]
    print('get mask finished')
    masked_chans = mask.all(axis=1)
    return masked_chans

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

def grab_spectra_manual(
        gf, ts, te, mask_fn, dm, mask=True, downsamp=4, subband=256, manual=False, t_start = 4.1, t_dur = 1.8, fit_del = 100e-3,plot_name = ""
):
    # load the filterbank file
    g = r.FilReader(gf)
    if ts < 0:
        ts = 0
    tsamp = float(g.header.tsamp)
    total_time = float(g.header.nsamples) * tsamp
    if te>total_time:
        te = total_time
        #change t_dur to be the time from 4.1s to end
        t_right = te - ts - t_start
        if t_right<t_dur:
            t_dur = t_right
    if ts<0:
        #shift t_start back by however much time you need
        t_start = t_start + ts*tsamp
        ts = 0


    print("start and end times", ts, te)
    nsamps = int((te - ts) / tsamp)
    nsamps = nsamps - nsamps % downsamp
    ssamps = int(ts / tsamp)
    # sampels to burst
    nsamps_start_zoom = int(t_start / tsamp)
    nsamps_end_zoom = int((t_dur+t_start) / tsamp)
    try:
        spec = g.read_block(ssamps, nsamps)
    except Exception as e:
        print(e)
        import pdb; pdb.set_trace()

    # load mask
    if mask:
        print("masking data")
        data, masked_chans = maskfile(mask_fn, spec, ssamps, nsamps)
        print(sum(masked_chans))
    data = data.dedisperse(dm)
    waterfall_dat, dat_ts = extract_plot_data(data,masked_chans,dm,downsamp,nsamps_start_zoom,nsamps_end_zoom)

    if manual:
        # this gives you the location of the peak
        amp, std, loc, sigma_width = fit_SNR_manual(
            dat_ts,
            tsamp * downsamp,
            fit_del,
            nsamps=int(t_dur/2 / tsamp / downsamp),
            ds_data=waterfall_dat,
            downsamp=downsamp,
        )

        while amp==-1:
            #fit has failed, get a larger time window and try again
            # make a copy to plot the waterfall
            t_start = t_start - 1
            t_dur = t_dur + 2
            print(t_start,t_dur)
            if (t_start<0)|((t_start+t_dur)>(te-ts)):
                break
            nsamps_start_zoom = int(t_start / tsamp)
            nsamps_end_zoom = int((t_dur+t_start) / tsamp)
            print(nsamps_start_zoom,nsamps_end_zoom)
            waterfall_dat, dat_ts = extract_plot_data(data,masked_chans,dm,downsamp,nsamps_start_zoom,nsamps_end_zoom)
            amp, std, loc, sigma_width = fit_SNR_manual(
                dat_ts,
                tsamp * downsamp,
                fit_del,
                nsamps=int(t_dur/2 / tsamp / downsamp),
                ds_data=waterfall_dat,
                downsamp=downsamp,
            )
        if (amp!=-1)&((loc<(0.49*t_dur))|(loc>(t_dur*0.51))|(sigma_width>2e-2)):
            #repeat if initial loc guess is wrong
            amp, std, loc, sigma_width = fit_SNR_manual(
                dat_ts,
                tsamp * downsamp,
                fit_del,
                nsamps=int(loc / tsamp / downsamp),
                ds_data=waterfall_dat,
                downsamp=downsamp,
            )
        #scale the noise by the number of channels that are not masked
        std = std * np.sqrt(sum(~masked_chans)/len(masked_chans))
        SNR = amp / std
        if amp!=-1:
            loc = loc+t_start
            ts_no_ds_zoom_start = int(loc/tsamp - 0.9/tsamp)
            ts_no_ds_zoom_end = int(loc/tsamp + 0.9/tsamp)
            ts_no_ds = data[:, ts_no_ds_zoom_start : ts_no_ds_zoom_end]
            ts_no_ds = np.mean(ts_no_ds[~masked_chans, :], axis=0)
            #scale the ts with the std so everythin is in units of noise
            #FLUENCE = fit_FLUENCE(
            #    ts_no_ds/std,
            #    tsamp,
            #    3 * sigma_width,
            #    nsamp=int(loc / tsamp),
            #    ds_data=waterfall_dat,
            #    plot=False,
            #)
        else:
            FLUENCE = -1
    else:
        # fit using downsampled values
        # this is mostly used for the injections
        try:
            amp, std, loc, sigma_width = autofit_pulse(
                dat_ts,
                tsamp * downsamp,
                fit_del,
                nsamps=int(t_dur / 2 / tsamp / downsamp),
                ds_data=waterfall_dat,
                downsamp=downsamp,
                plot=False,
                plot_name=plot_name,
            )
            #refit with new initial params
            amp, std, loc, sigma_width = autofit_pulse(
                dat_ts,
                tsamp * downsamp,
                fit_del,
                nsamps=int(loc / tsamp / downsamp),
                ds_data=waterfall_dat,
                downsamp=downsamp,
                plot=True,
                plot_name=plot_name,
            )
        except Exception as e:
            print(e)
            amp, std, loc, sigma_width = -1, -1, -1, -1
        #scale std by the sqrt of non masked chans
        std = std * np.sqrt(sum(~masked_chans)/len(masked_chans))
        SNR = amp / std
        # because loc is predetermined set start and end a predifined spot
        ts_no_ds_zoom_start = int(4.1/tsamp)
        ts_no_ds_zoom_end = int(5.9/tsamp)
        ts_no_ds = data[:, ts_no_ds_zoom_start : ts_no_ds_zoom_end]
        ts_no_ds = np.mean(ts_no_ds[~masked_chans, :], axis=0)
        #FLUENCE = fit_FLUENCE(
        #    ts_no_ds/std,
        #    tsamp,
        #    fit_del,
        #    nsamp=int(loc / tsamp),
        #    ds_data=waterfall_dat,
        #    plot=False,
        #)
    FLUENCE = -1
    # recalculate the amplitude given a gaussian pulse shape
    gaussian_amp = FLUENCE / sigma_width / np.sqrt(2 * np.pi)
    print("filename:", gf, "downsample:", downsamp, "FLUENCE:", FLUENCE)
    approximate_toa = g.header.tstart + ((te+ts)/2)/86400
    return FLUENCE, std, amp, gaussian_amp, sigma_width, SNR, approximate_toa
