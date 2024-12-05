# ## Import stuff

import numpy as np

import matplotlib
matplotlib.use('pdf')

from matplotlib import pyplot as plt
import copy
from collections.abc import Iterable

from presto_without_presto import rfifind
from iqrm import iqrm_mask
from rfifind_numpy_tools import write_new_mask_from

import sys
import argparse
from scipy.ndimage import generic_filter, generic_filter1d
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap

from gen_utils import handle_exception
import logging

from scipy.signal import find_peaks
from itertools import groupby
from more_itertools import consecutive_groups
from math import ceil
from operator import itemgetter

import matplotlib.colors as mc
import colorsys

# catch uncaught exceptions and put them in log too
sys.excepthook = handle_exception

# ## Define functions

# ### General utils

def darken_colour(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def output_plot(fig, pdf=None):
    if pdf is None:
        plt.show()
    else:
        fig.savefig(pdf, format='pdf')
    plt.close(fig)

def masked_frac(mask):
    return mask.sum()/mask.size

def get_ignorechans_from_mask(mask):
    nint, nchan = mask.shape
    return np.where(mask.sum(axis=0) == nint)[0]

def get_ignoreints_from_mask(mask):
    nint, nchan = mask.shape
    return np.where(mask.sum(axis=1) == nchan)[0]

def np_ignorechans_to_presto_string(array1d):
    return ",".join(list(array1d.astype(str)))

def write_mask_and_ignorechans(mask, outname, rfifind_obj, infstats_too=True):
    ignorechans_fname = f"{outname[:-5]}.ignorechans"
    logging.info(f"Writing ignorechans to {ignorechans_fname}")
    with open(ignorechans_fname, "w") as fignore:
        fignore.write(np_ignorechans_to_presto_string(get_ignorechans_from_mask(mask)))
    logging.info(f"Writing mask to {outname}")
    write_new_mask_from(outname, mask, rfifind_obj, infstats_too=infstats_too)

def wrap_up(mask, mask_exstats, rfifind_obj, means, var, pdf, outfilename, infstats_too):
    logging.info(f"Fraction of data masked: {masked_frac(mask)}")
    write_mask_and_ignorechans(mask, outfilename, rfifind_obj, infstats_too=infstats_too)

    logging.info(f"Making summary plots")
    make_summary_plots(mask, mask_exstats, rfifind_obj, means, var, pdf, title_insert="final")

    if pdf is not None:
        logging.info("Writing pdf")
        pdf.close()
    logging.info("Done")


def make_summary_plots(mask, mask_exstats, rfifind_obj, means, var, pdf, title_insert=""):
    """Plot the mask, and the masked pow_stats, means, and var"""
    figtmp, axtmp = plt.subplots()
    plot_mask(mask, ax=axtmp)
    figtmp.suptitle(f"{title_insert} mask ({masked_frac(mask):.2f})")
    output_plot(figtmp, pdf=p)

    pow_stats_plot_mask = copy.deepcopy(mask)
    if (rfifind_obj.pow_stats[-1,:] == 1).all():
        logging.info("Weird final interval stats for rfifind, it will be masked in pow_stats plot, but not in the mask itself")
        pow_stats_plot_mask[-1,:] = True
    figtmp, axtmp, cbartmp = plot_map_plus_sums(rfifind_obj.pow_stats, mask=pow_stats_plot_mask, returnplt=True)
    figtmp.suptitle(f"{title_insert} pow_stats ({masked_frac(mask):.2f})")
    output_plot(figtmp, pdf=pdf)
    del pow_stats_plot_mask

    figtmp, axtmp, cbartmp = plot_map_plus_sums(means.data, mask=mask_exstats, returnplt=True)
    figtmp.suptitle(f"{title_insert} means ({masked_frac(mask_exstats):.2f})")
    output_plot(figtmp, pdf=pdf)

    figtmp, axtmp, cbartmp = plot_map_plus_sums(var.data, mask=mask_exstats, returnplt=True)
    figtmp.suptitle(f"{title_insert} var ({masked_frac(mask_exstats):.2f})")
    output_plot(figtmp, pdf=pdf)


# https://www.tutorialspoint.com/how-to-make-a-histogram-with-bins-of-equal-area-in-matplotlib
def equal_area(x, nbin):
   pow = 0.5
   dx = np.diff(np.sort(x))
   tmp = np.cumsum(dx ** pow)
   tmp = np.pad(tmp, (1, 0), 'constant')
   return np.interp(np.linspace(0, tmp.max(), nbin + 1), tmp, np.sort(x))


# ### Down/upsample masks

def tscrunch_mask(msk, fac):
    """
    mask of shape (int,chan)
    scrunch in time by <fac>, taking the logical or of the two rows
    returned mask will be shape (int/fac, chan)
    """
    remainder = msk.shape[0] % fac
    if remainder:
        tmp = msk[:-remainder,:].astype(int)
        excess = msk[-remainder,:].astype(int)
        return_nint = msk.shape[0] // fac + 1

        mout = np.zeros((return_nint, msk.shape[1]), dtype=bool)
        mout[:-1,:] = tmp.reshape(-1,fac,tmp.shape[-1]).sum(1) > 0
        mout[-1,:] = excess.sum(0) > 0
        return mout
    else:
        tmp = msk.astype(int)
        return tmp.reshape(-1,fac,tmp.shape[-1]).sum(1) > 0

def upsample_mask(msk, fac):
    """
    mask of shape (int, chan)
    upsample in time by <fac>
    returned mask will be shape (int*fac, chan)
    """
    return np.repeat(msk, fac, axis=0)

def reshape_extra_stats_mask(rfimask_shape, msk, fdp_gulp, ptsperint):
    """
    reshape mask derived from extra_stats to match <rfimask_shape>
    fdp_gulp = gulp used when running fdp
    ptsperint = rfimask's ptsperint
    msk = mask to be reshaped (shape is of form nint,nchan)

    will use tscrunch_mask if fdp_gulp < ptsperint
    will use upsample_mask if fdp_gulp > ptsperint
    """
    if msk.shape == rfimask_shape:
        return msk
    if fdp_gulp < ptsperint:
        logging.warning("WARNING ptsperint > fdp_gulp, unless ran into memory issues go redo fdp with higher gulp")
        if ptsperint % fdp_gulp:
            raise AttributeError(f"ptsperint ({ptsperint}) does not divide evenly into fdp_gulp ({fdp_gulp})")
        tscrunch_fac = int(ptsperint // fdp_gulp)
        return tscrunch_mask(msk, tscrunch_fac)
    else:
        if fdp_gulp % ptsperint:
            raise AttributeError(f"fdp_gulp ({fdp_gulp}) does not divide evenly into ptsperint ({ptsperint})")
        upsample_fac = int(fdp_gulp // ptsperint)
        tmp = upsample_mask(msk, upsample_fac)
        if tmp.shape != rfimask_shape:
            if tmp.shape[0] == rfimask_shape[0] + 1:
                return tmp[:-1,:]
            else:
                raise AttributeError(f"odd shape problem:\noriginal {msk.shape}, upsampled to {(upsample_fac)}, trying to match {rfimask_shape}")
        else:
            return tmp
                
def reshape_rfifind_mask(extra_stats_shape, msk, fdp_gulp, ptsperint):
    """
    reshape mask derived from rfifind stats to match <extra_stats_shape>
    fdp_gulp = gulp used when running fdp
    ptsperint = rfimask's ptsperint
    msk = mask to be reshaped (shape is of form nint,nchan)

    will use upsample_mask if fdp_gulp < ptsperint
    will use tscrunch_mask if fdp_gulp > ptsperint
    """
    if msk.shape == extra_stats_shape:
        return msk
    if fdp_gulp < ptsperint:
        logging.warning("WARNING ptsperint > fdp_gulp, unless ran into memory issues go redo fdp with higher gulp")
        if ptsperint % fdp_gulp:
            raise AttributeError(f"ptsperint ({ptsperint}) does not divide evenly into fdp_gulp ({fdp_gulp})")
        upsample_fac = int(ptsperint // fdp_gulp)
        tmp = upsample_mask(msk, upsample_fac)
        if tmp.shape != extra_stats_shape:
            if tmp.shape[0] == extra_stats_shape[0] + 1:
                return tmp[:-1,:]
            else:
                raise AttributeError(f"odd shape problem:\noriginal {msk.shape}, upsampled to {(upsample_fac)}, trying to match {extra_stats_shape}")
    else:
        if fdp_gulp % ptsperint:
            raise AttributeError(f"fdp_gulp ({fdp_gulp}) does not divide evenly into ptsperint ({ptsperint})")
        tscrunch_fac = int(fdp_gulp // ptsperint)
        return tscrunch_mask(msk, tscrunch_fac)


# ## Zapping functions
def get_zeros_mask_alt(var_stats, ignorechans=[], verbose=False, plot_diagnostics=True, ax=None, imshow=True):
    """
    Get a mask where the var_stats = 0
    (changed from std_stats as that seems to always be 0 for the final interval (always always or just with my data? tbd))
    Then shuffle it +- 1 interval on each side
    Often if something went wrong, like a node going down or a network error, intervals either side are affected also
    This, hopefully catches that but also doesn't throw out whole channel if only part of it has dropped out.

    verbose and plot_diagnostics both concern where std==0 in the data in places not covered by ignorechans
    """
    tmp = (var_stats == 0) | np.isnan(var_stats)
    working_mask = np.zeros_like(tmp, dtype=bool)
    if ignorechans:
        working_mask[:,np.array(ignorechans)] = True
    if plot_diagnostics:
        if ax is None:
            fig, ax = plt.subplots()
        if imshow:
            ax.imshow(np.ma.array(tmp, mask=working_mask).T, aspect='auto', origin='lower')
        else:
            ax.pcolormesh(np.ma.array(tmp, mask=working_mask).T)
        ax.set_xlabel("int")
        ax.set_ylabel("chan")
        ax.set_title("Plot of where var_stats==0 or NaN, masked by the ignorechans")

    # add ignorechans to mask
    for ii in ignorechans:
        tmp[:,ii] = 1

    if verbose:
        ignorechans = set(ignorechans)
        inv_tmp = ~tmp  # so inv_tmp is 0 anywhere data is zero

        whole_ints = set(np.where(inv_tmp.sum(axis=1)==0)[0])
        whole_chans = set(np.where(inv_tmp.sum(axis=0)==0)[0])
        additional_whole_chans = whole_chans.difference(ignorechans)

        if whole_ints:
            logging.info(f"Found whole interval/s where var_stats==0: {sorted(list(whole_ints))}")
            working_mask[np.array(list(whole_ints)),:] = 1

        if additional_whole_chans:
            logging.info(f"Found whole channel/s where var_stats==0 (not covered by ignorechans): {sorted(list(additional_whole_chans))}")
            working_mask[:,np.array(list(whole_chans))] = 1


        # assume we're dealing with partial channels only and NO partial intervals
        inv_tmp2 = np.ma.array(inv_tmp, mask=working_mask)
        partial_channels = list(np.where((inv_tmp2 == 0).any(axis=0))[0])
        if partial_channels:
            logging.info(f"Found partial channel/s where var_stats==0 (not covered by ignorechans): {partial_channels}")

    # shuffling mask +-1 int
    tmp[1:,:] = (tmp[1:,:] | tmp[:-1,:])
    tmp[:-1,:] = (tmp[:-1,:] | tmp[1:,:])

    return tmp.astype(bool)


def cut_off_high_fraction(existing_mask, hard_thresh_chans=0.2, hard_thresh_ints=0.3 ,cumul_threshold_chans=0.95, cumul_threshold_ints=0.95, plot_diagnostics=True, verbose=True, ax=None, axhard=None):  #, debug=False):
    """
    First cut off channels with a fraction masked above <hard_thresh_chans> (this is intfrac in rfifind)
    And cut off intervals with a fraction masked above <hard_thresh_ints> (this is chanfrac in rfifind)

    Then:
    For each channel, calculate the fraction of intervals which have been zapped (discounting entire zapped intervals)
    Make a cumulative distribution for these fractions
    Zap fractions above where the cumulative distribution is greater than <cumul_threshold_*>
    e.g. if have 10 channels with zapped fractions of 0.1,0.1,0.4,0.1,0.1,0.1,0.1,0.2,0.1,0.2
        and set cumul_threshold_chans=0.9
        the top 10% would be zapped, in this case just channel 2 with a fraction of 0.4

    Same thing for each interval

    <cumul_threshold_chans> is used to determine which channels to zap (per-channel fractions plot)
    <cumul_threshold_ints> is used to determine which ints to zap (per-interval fractions plot)

    ax, if passed in, is an axes array of shape (2,2) (with sharex='col')
    axhard, if passed in, is an axes array of shape (1,2) (this displays the plot for the hard cut)
    these two are different as want the sharex='col' for ax, and if included the hard cut in that the scale would be annoting
    """
    nint, nchan = existing_mask.shape
    zapped_chans = np.where(existing_mask.sum(axis=0) == nint)[0]
    nzapped_chans = zapped_chans.shape[0]
    zapped_ints = np.where(existing_mask.sum(axis=1) == nchan)[0]
    nzapped_ints = zapped_ints.shape[0]
    logging.debug(f"start nzapped_chans={nzapped_chans}, nzapped_ints={nzapped_ints}")

    frac_data_chans = (existing_mask.sum(axis=0) - nzapped_ints)/(nint-nzapped_ints)
    frac_data_ints = (existing_mask.sum(axis=1) - nzapped_chans)/(nchan-nzapped_chans)

    #if debug:
    #    for cc in range(nchan):
    #        logging.debug(f"{cc}:\t{frac_data_chans[cc]}")

    # cut off things above the hard threshold
    chans_to_zap_hard = np.where((frac_data_chans > hard_thresh_chans) & (frac_data_chans != 1))[0]
    if verbose:
        logging.info(f"Channels to zap from hard fraction threshold {hard_thresh_chans}: {chans_to_zap_hard}")
    ints_to_zap_hard = np.where((frac_data_ints > hard_thresh_ints) & (frac_data_ints != 1))[0]
    if verbose:
        logging.info(f"Intervals to zap from hard fraction threshold {hard_thresh_ints}: {ints_to_zap_hard}")

    if plot_diagnostics:
        show_plot = False
        if axhard is None:
            fighard, axhard = plt.subplots(1,2, figsize=(12,4))
            show_plot = True
        axhard[0].hist(frac_data_chans[frac_data_chans != 1], bins=40, density=True)
        axhard[0].axvline(hard_thresh_chans, c='red')
        axhard[1].hist(frac_data_ints[frac_data_ints != 1], bins=40, density=True)
        axhard[1].axvline(hard_thresh_ints, c='red')
        axhard[0].set_title("chans hard threshold cut")
        axhard[1].set_title("ints hard threshold cut")
        if show_plot:
            plt.show()
            plt.close()

    assert (zapped_chans == np.where(frac_data_chans == 1)[0]).all()
    assert (zapped_ints == np.where(frac_data_ints == 1)[0]).all()

    # remake frac_data_chans
    fracs_mask_hard = np.zeros((nint, nchan), dtype=bool)
    fracs_mask_hard[:,chans_to_zap_hard] = True
    fracs_mask_hard[ints_to_zap_hard, :] = True

    nzapped_ints = len(set(zapped_ints).union(set(ints_to_zap_hard)))
    nzapped_chans = len(set(zapped_chans).union(set(chans_to_zap_hard)))
    logging.debug(f"after hard nzapped_chans={nzapped_chans}, nzapped_ints={nzapped_ints}")

    # If don't need to do the cumulative cut, exit early
    if cumul_threshold_chans == 1 and cumul_threshold_ints == 1:
        logging.info("Not running cumulative cut as both thresholds are set to 1")
        return fracs_mask_hard
    
    frac_data_chans = ((existing_mask|fracs_mask_hard).sum(axis=0) - nzapped_ints)/(nint - nzapped_ints)
    frac_data_ints = ((existing_mask|fracs_mask_hard).sum(axis=1) - nzapped_chans)/(nchan - nzapped_chans)

    # make cumulative distributions and apply threshold
    sorted_fracs_chan = np.sort(frac_data_chans[frac_data_chans !=1])
    N_chans = sorted_fracs_chan.size
    cumul_chans = np.arange(N_chans)/N_chans
    try:
        chan_frac_threshold = sorted_fracs_chan[cumul_chans >= cumul_threshold_chans][0]
    except IndexError:
        chan_frac_threshold = 1
    chans_to_zap = np.where((frac_data_chans >= chan_frac_threshold) & (frac_data_chans != 1))[0]
    if verbose:
        logging.info(f"Channels to zap from cumulative threshold of {cumul_threshold_chans}: {sorted(list(set(chans_to_zap).difference(set(zapped_chans))))}")


    sorted_fracs_int = np.sort(frac_data_ints[frac_data_ints !=1])
    N_ints = sorted_fracs_int.size
    cumul_ints = np.arange(N_ints)/N_ints
    try:
        int_frac_threshold = sorted_fracs_int[cumul_ints >= cumul_threshold_ints][0]
    except IndexError:
        int_frac_threshold = 1
    ints_to_zap = np.where((frac_data_ints >= int_frac_threshold) & (frac_data_ints != 1))[0]
    if verbose:
        logging.info(f"Intervals to zap from cumulative threshold of {cumul_threshold_ints}: {sorted(list(set(ints_to_zap).difference(set(zapped_ints))))}")

    if plot_diagnostics:
        show_plot = False
        if ax is None:
            fig, ax = plt.subplots(2,2, sharex='col')
            show_plot = True
        ax[0,0].set_title(f"Per-channel fractions")
        ax[0,0].hist(frac_data_chans[frac_data_chans != 1], bins=40, density=True)
        if cumul_threshold_chans != 1:
            ax[0,0].axvline(chan_frac_threshold, c='orange')
        else:
            # the hard thresholds are not marked on the plot if cumul_threshold_chans != 1
            # because the fractions plotted are based on already having removed channels over the hard threshold
            # it's still a bit misleading to plot it when cumul_threshold_chans == 1, but a bit less so
            ax[0,0].axvline(hard_thresh_chans, c='red')

        ax[1,0].set_title(f"cumulative, threshold={cumul_threshold_chans}")
        ax[1,0].plot(sorted_fracs_chan, cumul_chans)
        ax[1,0].axhline(cumul_threshold_chans, c='orange')

        ax[0,1].set_title(f"Per-interval fractions")
        ax[0,1].hist(frac_data_ints[frac_data_ints != 1], bins=40, density=True)
        if cumul_threshold_ints != 1:
            ax[0,1].axvline(int_frac_threshold, c='orange')
        else:
            ax[0,1].axvline(hard_thresh_ints, c='red')
        ax[1,1].set_title(f"cumulative, threshold={cumul_threshold_ints}")
        ax[1,1].plot(sorted_fracs_int, cumul_ints)
        ax[1,1].axhline(cumul_threshold_ints, c='orange')

        if show_plot:
            plt.show()
            plt.close()

    fracs_mask = np.zeros((nint, nchan), dtype=bool)
    fracs_mask[:,zapped_chans] = True
    fracs_mask[zapped_ints,:] = True
    fracs_mask[:,chans_to_zap] = True
    fracs_mask[ints_to_zap,:] = True

    return fracs_mask|fracs_mask_hard

# helper function for get_step_chans
def rescale(a, b):
    """rescale a to b so that they have the same min and max"""
    mina = np.ma.min(a)
    dra = np.ma.max(a) - mina
    minb = np.ma.min(b)
    drb = np.ma.max(b) - minb

    return (a - mina)*drb/dra + minb

# ID channels with sharp steps in them
# based off https://stackoverflow.com/questions/48000663/step-detection-in-one-dimensional-data
def get_step_chans(stat, thresh=30, ignorechans=[], return_stats=False, return_plots=False, output_pdf=None, imshow=True):
    """stat of shape (nint, nchan) (numpy masked array)
    for each channel in stat, look for steps (via subtracting the mean and then taking the negative of np.cumsum)
    If the max of that is > thresh it gets zapped
    NB This means if a giant pulse or something happens will likely zap it

    Returns a list of channels to zap
    if return_stats returns
    <list of channels to zap> <channels> <max peak value>
    if return_plots will return a list of figs also (as last item returned)
    output_pdf is a bit weird. If return_plots is True then it does nothing
        otherwise it should be None (in which case plots are shown) or a PdfPages object
    """
    to_zap = []
    cs = []
    ms = []
    if return_plots:
        figs = []

    for c in np.arange(stat.shape[1]):
        if c not in ignorechans and not stat.mask[:,c].all():
            dary = copy.deepcopy(stat[:,c])
            dary -= np.average(stat[:,c])
            dary_step = -np.ma.cumsum(dary)
            m = dary_step.max()
            cs.append(c)
            ms.append(m)
            if m > thresh:
                figtmp, axtmp = plt.subplots(2,1)
                stat_tmp = stat[:,max(0,c-10):min(c+10,stat.shape[1])]
                if imshow:
                    axtmp[1].imshow(stat_tmp.T, vmin=stat_tmp.min(), vmax=stat_tmp.max(), aspect='auto', origin='lower')
                else:
                    axtmp[1].pcolormesh(stat_tmp.T, vmin=stat_tmp.min(), vmax=stat_tmp.max())
                axtmp[1].axhline(min(10,c), c='red')
                figtmp.suptitle(f"{c}: {dary_step.max()}")
                axtmp[0].plot(dary)
                axtmp[0].plot(rescale(dary_step, dary), c='orange')
                if return_plots:
                    figs.append(figtmp)
                else:
                    output_plot(figtmp, pdf=output_pdf)
                to_zap.append(c)
    figtmp, axtmp = plt.subplots()
    axtmp.hist(ms, equal_area(ms, 32), density=True)
    axtmp.axvline(thresh, c="orange")
    axtmp.set_title("distribution of stat used to cut on (max of a mean-subtracted cumulative sum)")
    if return_plots:
        figs.append(figtmp)
    else:
        output_plot(figtmp, pdf=output_pdf)

    if not return_stats and not return_plots:
        return to_zap
    else:
        to_return = [to_zap]
        if return_stats:
            to_return.extend([cs, ms])
        if return_plots:
            to_return.append(figs)
        return to_return


# ### iqrm

def run_iqrm_2D(data, mask, axis, r, size_fill=3, ignorechans=[], threshold=3):
    """Run iqrm on 1d slices of 2d data,
    axis means which direction to take the slice
    aka if data is shape (nint, nchan):
        axis=1 means for each interval, run iqrm looking for outlier channels
        axis=0 means for each channel, run iqrm looking for outlier intervals

    size_fill = factor to use for fill_value, med + <size_fill> * (max-med)

    ignorechans: channels to ignore, only used if axis=1
    
    tweaked so mask is used - votes from any masked channel/int intersection will be ignored
    """
    a0,a1 = data.shape
    out = np.zeros((a0,a1), dtype=bool)
    if mask is None:
        mask=np.zeros_like(data, dtype=bool)

    if mask.all() or np.isnan(data).all():
        return np.ones_like(data, dtype=bool)

    if axis == 1:
        for i in range(a0):
            ignore = list(set(np.where(mask[i,:])[0]).union(set(ignorechans)))
            out[i,:], v = iqrm_mask(data[i,:], r, ignorechans=ignore, threshold=threshold)

    if axis == 0:
        for j in range(a1):
            ignore = np.where(mask[:,j])[0]
            out[:,j], v = iqrm_mask(data[:,j].T, r, ignorechans=ignore, threshold=threshold)

    return out|mask

def get_fill_value(masked_statistic, factor, mode="max-med"):
    """Get fill value for filling in NaNs
    modes are
        "max" = Use factor*max
        "max-med" = Use median + factor*(max-median)
        "med" = Use median
    """
    if mode == "max":
        return factor * masked_statistic.max()
    if mode == "max-med":
        md = np.ma.median(masked_statistic)
        return md + factor * (masked_statistic.max() - md)
    if mode == "med":
        return np.ma.median(masked_statistic)
    else:
        logging.error("Invalid option for mode, must be one of ['max', 'max-med', 'med']")


def get_iqrm_chans(stat, mask, out='mask', rfac=8, flip=False, reduction_function=np.ma.median, size_fill=3, fill_mode="max-med", ignorechans=[]):
    """
    Mask channels based on (nint, nchan) stat passed in operated upon by <reduction_function> along the interval axis(=0)

    flip = do max - thing, to ID dips
    masked values filled in with 2*max

    The result can either be a set of channels (out='set')
    or the channel mask broadcast to shape (nint, nchan) (out='mask')

    the r used in iqrm is set to nchan/rfac
    IQRM documentation recommends 10 as a starting/generally-ok-for-most-setups option
    """
    if out not in ['mask', 'set']:
        raise AttributeError(f"out of {out} not recognised, must be one of 'mask', 'set'")
    if mask is None:
        masked_thing = stat
    else:
        masked_thing = np.ma.array(stat, mask=mask)
    reduced_thing = reduction_function(masked_thing, axis=0)
    non_nan = reduced_thing[np.where(np.invert(np.isnan(reduced_thing)))]
    med_non_nan = np.ma.median(non_nan)

    if flip:
        logging.debug("flipping")
        use = non_nan.max() - reduced_thing
        # setting the fill value to 3*(max-median) dist / same for flip, does improve things a bit
        # try that for other stuff too where it's still 2*max?
        fill_value1 = med_non_nan + size_fill*(med_non_nan - non_nan.min())
        fill_value = get_fill_value(use, size_fill, mode=fill_mode)
        logging.debug("fill value check", fill_value1, fill_value)
    else:
        logging.debug("not flipping")
        use = reduced_thing
        fill_value1 = med_non_nan + size_fill*(non_nan.max() - med_non_nan)
        fill_value = get_fill_value(use, size_fill, mode=fill_mode)
        logging.debug("fill value check", fill_value1, fill_value)

    r = stat.shape[1] / rfac

    #plt.plot(use.filled(fill_value))

    iqmask_stdavg, v = iqrm_mask(use.filled(fill_value), radius=r, ignorechans=ignorechans)

    chanset = set(np.where(iqmask_stdavg)[0])
    if out == 'set':
        return chanset

    iqmask_stdavg_broadcast = (iqmask_stdavg * np.ones(stat.shape)).astype(bool)
    if out == 'mask':
        return iqmask_stdavg_broadcast


def iqrm_of_median_of_means(means_data, mask, r, to_return="mask", plot=True, ax=None):
    """
    Find intervals to zap in means_data (shape (nint, nchan))
    by taking the median in each interval and running 1D iqrm to look for outliers

    r = radius to use with iqrm
    to_return governs the form of the output. Options are "array", "set", "mask", "chan_mask"
        array = array of interval indices to zap
        set = set of interval indices to zap
        mask = mask of same shape as means
        int_mask = mask of shape means.shape[0]

    plot=True, make a plot showing the median of the means with the masked version overlaid

    """
    if to_return not in ["array", "set", "mask", "int_mask"]:
        raise AttributeError(f"to_return ({to_return}) must be one of 'array', set', 'mask', 'int_mask'")
    thing = np.ma.median(np.ma.array(means_data, mask=mask), axis=1)
    q, v = iqrm_mask(thing, radius=r)
    if plot:
        if ax is None:
            figtmp, ax = plt.subplots()
        ax.plot(thing, "x")
        ax.plot(np.ma.array(thing.data, mask=(thing.mask|q)), "+")
        ax.set_xlabel("interval")
        ax.set_ylabel("median of the means")
    if to_return == "int_mask":
        return q
    if to_return == "array":
        return np.where(q)[0]
    if to_return == "set":
        return set(np.where(q)[0])
    if to_return == "mask":
        msk = np.zeros(means_data.shape, dtype=bool)
        msk[np.where(q)[0],:] = True
        return msk


# ### iterative +- sigma


# this works but iqrm is better
# and I was using a positive_only version of this which is literally what iqrm is designed for
# nope I take it back, iqrm zaps more than it needs to!

# this is less necessary now I have get_step_chans, but still good!
def reject_pm_sigma_iteration(arr1d, init_mask, thresh=5, plot=False, positive_only=False, iteration=0, prev_hits=np.array([9E9]), middle80=False):
    """
    Iteravely reject points more than <thresh>*std from the median

    positive_only=True: only reject points over median + <thresh>*std
    middle80=True: only do this process on the middle 80% of arr1d
        (i.e. ignoring 10% at each edge of the band)
    
    Don't pass in these arguments, they're so it iterates:
    iteration
    prev_hits
    """
    tmp0 = np.ma.array(arr1d, mask=init_mask)
    working_mask = copy.deepcopy(init_mask)

    nchans = tmp0.shape[0]
    if middle80:
        slc = slice(round(0.1*nchans), round(0.9*nchans)+1)
    else:
        slc = slice(None, None)

    md = np.ma.median(tmp0[slc])
    sig = np.ma.std(tmp0[slc])

    tmp = copy.deepcopy(tmp0)
    if middle80:
        tmp[:round(0.1*nchans)+1] = md
        tmp[round(0.9*nchans):] = md
    

    lo_lim = md - thresh*sig
    hi_lim = md + thresh*sig

    if positive_only:
        condit = (tmp > hi_lim)
    else:
        condit = ((tmp > hi_lim) | (tmp < lo_lim))

    all_hits = np.where(condit)[0]
    new_hits = np.ma.where(condit)[0]
    logging.debug(f"iteration {iteration}: std = {sig}, {condit.sum()} match/es : {np.ma.where(condit)[0]}" )

    if plot:
        plt.plot(np.ma.array(arr1d, mask=working_mask), "x")
        plt.axhline(hi_lim, c='orange')
        if not positive_only:
            plt.axhline(lo_lim, c='orange')
        plt.show()

    if len(new_hits) == 0:
        logging.info(f"channels zapped: {np.where(condit)[0]}")
        return working_mask

    for c in np.where(condit):
        working_mask[c] = True

    return reject_pm_sigma_iteration(arr1d, working_mask, thresh=thresh, plot=plot, positive_only=positive_only, iteration=iteration+1, prev_hits=all_hits, middle80=middle80)

    



# ### plotting

plt.rcParams['figure.figsize'] = [12, 8]

def plot_stat_map(stat, axis=None, cbar_axis=None, mask=None, imshow=True, **plot_kwargs):
    nint, nchan = stat.shape
    # grids = np.meshgrid(np.arange(nint + 1), np.arange(nchan + 1), indexing='ij')
    # for some WEIRD reason passing in the grids introduces a bunch of artifacts
    # sure not real as when zoom in they disappear/some are still there but much smaller than 1 int/chan
    # tried passing in centers and shading='nearest' and they're still there
    # tried a bunch of other things like snap=True, pcolor, saving the fig, etc but NADA! V confusing!!
    # doesn't show up is use pcolormesh without x and y though
    # BUT some bright single channel stuff you then can't see.
    # So I think he renderer is upping the resolution or something odd and it's producing weirdness
    # could not find a solution so using the thing that means I don't freak out that the data is terrible and lose half a day
    # UPDATE: also get with pcolormesh :( 
        # apparently this issue has a 10+year history and is really hard to fix 
        # https://stackoverflow.com/questions/27092991/white-lines-in-matplotlibs-pcolor
        # https://github.com/matplotlib/matplotlib/issues/1188
        # trying switching to imshow

    # cbar_axis only use if axis is not None

    if type(mask) != np.ndarray:
        to_plot = stat
    else:
        to_plot = np.ma.masked_array(data=stat, mask=mask)

    if axis == None:
        #im = plt.pcolormesh(grids[0], grids[1], to_plot, shading='flat', **plot_kwargs)
        if imshow:
            im = plt.imshow(to_plot.T, aspect='auto', origin='lower', **plot_kwargs)
        else:
            im = plt.pcolormesh(to_plot.T, **plot_kwargs)
        plt.xlabel="interval"
        plt.ylabel="channel"
        plt.colorbar(im)
        plt.show()
    else:
        #im = axis.pcolormesh(grids[0], grids[1], to_plot, shading='flat', **plot_kwargs)
        if imshow:
            im = axis.imshow(to_plot.T, aspect='auto', origin='lower', **plot_kwargs)
        else:
            im = axis.pcolormesh(to_plot.T, **plot_kwargs)
        if cbar_axis is None:
            cbar = plt.colorbar(im, ax=axis)
        else:
            cbar = plt.colorbar(im, cax=cbar_axis, orientation="horizontal")
        return axis, cbar

def plot_stat_v_nchan(stat, axis=None, mask=None, reduction_function=np.ma.median, **plot_kwargs):
    nint, nchan = stat.shape

    if not isinstance(reduction_function, Iterable):
        reduction_function = [reduction_function]

    if type(mask) != np.ndarray:
        to_plot = stat
    else:
        to_plot = np.ma.masked_array(data=stat, mask=mask)

    if axis == None:
        for red_func in reduction_function:
            plt.plot(red_func(to_plot, axis=0), np.arange(nchan), **plot_kwargs)
        plt.xlabel="channel"
        plt.show()
    else:
        z_min, z_max = None, None
        for red_func in reduction_function:
            z = red_func(to_plot, axis=0)
            axis.plot(z, np.arange(nchan), **plot_kwargs)
            if z_min is None:
                z_min = z.min()
            else:
                z_min = min(z_min, z.min())
            if z_max is None:
                z_max = z.max()
            else:
                z_max = max(z_max, z.max())
        return z_min, z_max
            

def plot_stat_v_nint(stat, axis=None, mask=None, reduction_function=np.ma.median, **plot_kwargs):
    nint, nchan = stat.shape

    if type(mask) != np.ndarray:
        to_plot = stat
    else:
        to_plot = np.ma.masked_array(data=stat, mask=mask)

    if not isinstance(reduction_function, Iterable):
        reduction_function = [reduction_function]

    x = np.arange(nint)

    if axis == None:
        for red_func in reduction_function:
            plt.plot(x, z, **plot_kwargs)

        plt.xlabel("interval")
        plt.show()
    else:
        z_min, z_max = None, None
        for red_func in reduction_function:
            z = red_func(to_plot, axis=1)
            axis.plot(x, red_func(to_plot, axis=1), **plot_kwargs)
            if z_min is None:
                z_min = z.min()
            else:
                z_min = min(z_min, z.min())
            if z_max is None:
                z_max = z.max()
            else:
                z_max = max(z_max, z.max())
        return z_min, z_max

def plot_map_plus_sums(stat, mask=None, reduction_function=np.ma.median, returnplt=False, fill=True, **plot_kwargs):
    """returnplt => return fig, ((ax0, ax1), (ax2, ax3))
    fill=True => if mask passed in or stat is a masked array, fill in masked values with the median"""
    nint, nchan = stat.shape
    #grids = np.meshgrid(np.arange(nint + 1), np.arange(nchan + 1), indexing='ij')
    widths = [3,1]
    heights = [1,3]

    figsize = None
    if "figsize" in plot_kwargs:
        figsize = plot_kwargs["figsize"]
    fig = plt.figure(figsize=figsize)#layout="constrained")
    # make outer grid to hold colourbar and data
    outer_grid = fig.add_gridspec(2, 1, wspace=0, hspace=0.3, height_ratios=[0.95,0.05])
    inner_grid = outer_grid[0].subgridspec(2, 2, wspace=0, hspace=0, width_ratios=widths, height_ratios=heights)

    ((ax0, ax1), (ax2, ax3)) = inner_grid.subplots()
    ax0.sharex(ax2)
    ax3.sharey(ax2)
    axcbar = fig.add_subplot(outer_grid[1])

    if mask is not None:
        statt = np.ma.masked_array(data=stat, mask=mask)
    else:
        statt = stat
    if isinstance(statt, np.ma.MaskedArray) and fill:
        ax2, cbar = plot_stat_map(statt.filled(np.ma.median(statt)), axis=ax2, cbar_axis=axcbar)
    else:
        plot_stat_map(statt, axis=ax2)
    x1, x2 = plot_stat_v_nchan(statt, axis=ax3, reduction_function=reduction_function)
    #ax3.set_xlim(0.9*x1, 1.1*x2)
    y1, y2 = plot_stat_v_nint(statt, axis=ax0, reduction_function=reduction_function)
    #ax0.set_ylim(0.9*y1, 1.1*y2)

    # colour bar throws things off, this rescales the v_nint plot
    # don't need this as cbar is in its own axis
    #pos_map = ax2.get_position()
    #pos_int = ax0.get_position()
    #ax0.set_position([pos_map.x0,pos_int.y0,pos_map.width,pos_int.height])

    for spine in ["top", "right"]:
        _ = ax1.spines[spine].set_visible(False)
    
    ax1.tick_params(labelbottom=False, labeltop=False)
    ax3.tick_params(labelleft=False)

    ax1.set_xticks([])
    ax1.set_yticks([])

    # labels
    ax2.set_xlabel("Interval")
    ax2.set_ylabel("Channel")
    if returnplt:
        return fig, ((ax0, ax1), (ax2, ax3)), cbar
    else:
        plt.show()

def plot_masked_channels_of_med(thing, channels, ax=None):  #, sig_lims=[3,3]):
    """
    plt the median of <thing> (shape (nint, nchan)) along axis 0
    and highlight <channels>
    sig_lims is for setting ylimits, limits come from med - sig_lims[0]*std and med + sig_lims[0]*std
        where med and std are calculated ignoring values in <channels>
        (if any non-<channels> values are greater or less than this that will be the limit instead)
    """
    nchan = thing.shape[1]
    med_thing = np.ma.median(thing, axis=0)
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(np.arange(nchan), med_thing, "+")
    ax.plot(np.array(list(channels)), med_thing[np.array(list(channels))], "x")
    #unmasked_channels = set(range(nchan)).difference(set(channels))
    #get_limits_from  = med_thing[np.array(list(unmasked_channels))]
    #md = np.ma.median(get_limits_from)
    #std = np.ma.std(get_limits_from)
    #lo = min([md-sig_lims[0]*std, get_limits_from.min()])
    #hi = max([md+sig_lims[1]*std, get_limits_from.max()])
    #ax.set_ylim(lo,hi)

def plot_mask(mask, ax=None, imshow=True):
    if ax is None:
        fig, ax = plt.subplots()
    if imshow:
        ax.imshow(mask.T, aspect='auto', origin='lower')
    else:
        ax.pcolormesh(mask.T)
    ax.set_ylabel("channel")
    ax.set_xlabel("interval")

def plot_mask_comparison(maska, maskb, title="", ax=None, returnplt=False, colorbar=False, ignorechans=None, imshow=True):
    """
    Plot maska - maskb
    If passed in, ignorechans are overplotted in black
    """
    imshow_kwargs =  dict(aspect='auto', origin='lower')
    if ax is None:
        fig, ax = plt.subplots()
    if imshow:
        im = ax.imshow((maska.astype(int) - maskb.astype(int)).T, **imshow_kwargs)
    else:
        im = ax.pcolormesh((maska.astype(int) - maskb.astype(int)).T)
    if colorbar:
        plt.colorbar(ax=ax)

    if ignorechans is not None and ignorechans != []:
        cmap2 = ListedColormap(['none', 'black'])
        ignorechans_mask = np.zeros_like(maska)
        for c in ignorechans:
            ignorechans_mask[:,c] = True
        if imshow:
            ax.imshow(ignorechans_mask.T, cmap=cmap2, **imshow_kwargs)
        else:
            ax.pcolormesh(ignorechans_mask.T, cmap=cmap2)


    
    if title:
        title += "\n"
    title += f"Newly masked data: {100*(masked_frac(maska) - masked_frac(maskb)):.2f}\%"
    ax.set_title(title)
    if returnplt:
        return fig, ax




def check_mask_and_continue(old_mask, old_mask_exstats, add_mask, add_mask_exstats, threshold, rfimask, means, var, pdf, stage=None, always_plot_summary=False, noupdate=False, time_lim_zapped_intervals_s=400):
    """
    Check if adding add_mask to old_mask means the masking fraction goes above the <threshold> 
    If it does:
        plot summary plots which will hopefully show what went wrong
        return the old masks
    If it doesn't:
        return (old_mask | add_mask), (old_mask_exstats | add_mask_exstats)

    if noupdate then don't update the maske even if the mask fraction is below the threshold

    if always_plot_summary then alway plot things even in masked fraction is over the threshold
    
    time_lim_zapped_intervals_s: Default 60. If a consecutive chunk of intervals is zapped > this many seconds in length, then there is probably a problem. Exit
    """
    zap_frac = masked_frac(old_mask | add_mask)
    zapped_ints = get_ignoreints_from_mask(old_mask | add_mask)
    zap_int_frac = len(zapped_ints) / old_mask.shape[0]

    # check for large chunks of zapped intervals
    interval_chunk_lengths = np.array([len(list(group)) for group in consecutive_groups(zapped_ints)])
    interval_chunk_length_limit = ceil(time_lim_zapped_intervals_s / rfimask.dtint)
    if (interval_chunk_lengths > interval_chunk_length_limit).any():
        logging.warning(f"{stage}: zaps a section of {max(interval_chunk_lengths)} intervals (={max(interval_chunk_lengths)*rfimask.dtint} s). This probably indicates a problem, plotting summary and exiting")
        make_summary_plots(add_mask, add_mask_exstats, rfimask, means, var, pdf, title_insert=f"ERROR stage {stage}")
        if pdf is not None:
            logging.info("Writing pdf")
            pdf.close()
        logging.info("Done")
        sys.exit(1)

    # check for large masking fractions
    if  zap_frac >= threshold:
        logging.warning(f"{stage}: zaps {zap_frac} of data, which is over the problem threshold, plotting summary and skipping")
        logging.info(f"{stage}: working mask unchanged")
        make_summary_plots(add_mask, add_mask_exstats, rfimask, means, var, pdf, title_insert=f"ERROR stage {stage}")
        return old_mask, old_mask_exstats
    elif zap_int_frac >= 0.3:  # Hopefully this is now obselete, leave in just in case
        logging.warning(f"{stage}: completely zaps {zap_int_frac} of the intervals, this probably indicates a problem, plotting summary and exiting")
        logging.info(f"{stage}: working msk unchanged")
        make_summary_plots(add_mask, add_mask_exstats, rfimask, means, var, pdf, title_insert=f"ERROR stage {stage}")
        if pdf is not None:
            logging.info("Writing pdf")
            pdf.close()
        logging.info("Done")
        sys.exit(1)
    else:
        logging.info(f"{stage}: zaps {zap_frac} of data (change = {zap_frac - masked_frac(old_mask)})")
        if always_plot_summary:
            make_summary_plots((add_mask|old_mask), (add_mask_exstats|old_mask_exstats), rfimask, means, var, pdf, title_insert=f"Summary stage {stage}")
        if noupdate:
            logging.info(f"{stage}: NOT updating working mask - ran with noupdate")
            return old_mask, old_mask_exstats
        else:
            logging.info(f"{stage}: updating mask")
            return (old_mask | add_mask), (old_mask_exstats | add_mask_exstats)

## New, mainly things to ID the stupid steps
def check_mean(arr, sign, slc_pre, slc_post, threshold, fn=np.ma.mean):
    """helper function for find_step"""
    stat_pre = fn(arr[slc_pre])
    stat_post = fn(arr[slc_post])
    stat_diff = stat_post - stat_pre
    stat_check = (sign*stat_diff >= threshold) and not (np.ma.is_masked(stat_diff))
    return stat_check, stat_diff 

def combine_posneg_negpos(posneg, negpos):
    """
    helper function for find_step
    
    return a sorted list combining found_pos_neg and found_neg_pos in find_step
    will return a list where each member has 3 elements
    the first 2 are the usual entries, the indices each side of the transition
    the third indicates whether the entry came from posneg or negpos
    posneg => -1
    negpos => +1
    (chosen from the sign of the expected difference in means if a step is present,
    this is also the expected sign of the peak in the gradient)
    """
    seq = [*posneg, *negpos]
    signs = [*[-1]*len(posneg), *[+1]*len(negpos)]
    to_sort = [[*x, signs[i]] for i,x in enumerate(seq)]
    return [y for x,y in sorted(enumerate(to_sort), key = lambda x: x[1][0][0])]

def get_peaks_iqrm(arr, r=10, threshold=5):
    if np.ma.is_masked(arr):
        ign = list(np.where(arr.mask)[0])
        stat = arr.filled(np.nan)
    else:
        ign = []
        stat=arr
    q, v = iqrm_mask(stat, radius=int(len(arr)/r), threshold=threshold, ignorechans=ign)
    pks = list(np.where(q)[0])
    return [pk for pk in pks if pk not in ign]

def split_into_consecutive(lst, startsends=True):
    out = []
    for k, g in groupby(enumerate(lst), lambda i_x: i_x[0] - i_x[1]):
        out.append(list(map(itemgetter(1), g)))
    if not startsends:
        return out
    starts = []
    ends = []
    for subsec in out:
        starts.append(subsec[0])
        ends.append(subsec[-1])
    return out, starts, ends

def any_indices_with_region(idx, other_idxes, region):
    """
    Returns True if there are indexes within +- region of idx within other_idxes
    (if idx is in other_idxes it is ignored)
    """
    tmp = other_idxes[other_idxes!=idx]
    diff = np.abs(tmp - idx)
    if (diff > region).all():
        return False
    return True

# try looking for places where iqrm switches from +ve to -ve with nothing or masked values in between?
def find_step(arr1d, debug=False, plots=False, return_plots=False, prom=1, mean_check_thresh=None, region_buffer=None):
    """
    Find a step in a 1d array
    run iqrm over +arr1d and -arr1d
    look for places where there's a section flagged in one directly followed by a section flagged in the other
    (or only separated by masked data is arr1d is a masked array)

    Then check using scipy.signal.find_peaks:
    Run find_peaks on +-arr1d using prominence=prom 
    if find +ve and -ve peaks in selection, indicates likely peak rather than step
    if find only -ve or only +ve peaks (for the correct kind of step) return True
    # region_buffer limits this check in the event that the length of the iqrm-flagged section before/after transition is > region_buffer

    mean_check_thresh
    if None then don't run
    otherwise look for a change in the mean of the correct size over the transition - this didn't work well, thus why default is None
    But if the baseline is flat might be worth checking out again

    debug = True
    print more info (if log level set to logging.DEBUG), and always make a plot
    plots = True
    Use if debug is False, but still want plots iff a step was found

    returns:
    False False if no step found, 
    True False if step found and a companion peak check is recommended
    True True if step found and companion peak check should be skipped.
        This is if there's only one peak. 
        If it only found one peak and found it of the right sign and in the correct region it's a pretty solid yes
        But it might be too short a peak to pass the companion check

    if return_plots, returns False/True, False/True, fig, ax

    """
    iqrm_pos = get_peaks_iqrm(arr1d)
    iqrm_neg = get_peaks_iqrm(-arr1d)

    superyes = False

    if np.ma.is_masked(arr1d):
        ign = list(np.where(arr1d.mask)[0])
    else:
        ign = []

    pos_split, pos_starts, pos_ends = split_into_consecutive(iqrm_pos)
    neg_split, neg_starts, neg_ends = split_into_consecutive(iqrm_neg)
    ign_split, ign_starts, ign_ends = split_into_consecutive(ign)


    logging.debug(f"positive sections flagged: {pos_split}")
    logging.debug(f"negative sections flagged: {neg_split}")

    # look for a positive section followed by a negative section or only separated by ignorechans
    found_pos_neg = []
    for i, pos in enumerate(pos_split):
        start = pos_starts[i]
        end = pos_ends[i]
        if end+1 in ign_starts:
            end = ign_ends[ign_starts.index(end+1)]

        # have to loop through in case of overlap
        for j, neg in enumerate(neg_split):
            if neg[0] <= end+1 and neg[-1] >= end+1:
                found_pos_neg.append([pos_split[i], neg])
                logging.debug("found pos-> neg", pos_split[i], neg)
        #if end+1 in neg_starts:
        #    found_pos_neg.append([pos_split[i], neg_split[neg_starts.index(end+1)]])

    found_neg_pos = []
    for i, pos in enumerate(neg_split):
        start = neg_starts[i]
        end = neg_ends[i]
        if end+1 in ign_starts:
            end = ign_ends[ign_starts.index(end+1)]
        # have to loop through in case of overlap
        for j, pos in enumerate(pos_split):
            if pos[0] <= end+1 and pos[-1] >= end+1:
                found_neg_pos.append([neg_split[i], pos])
                logging.debug("found neg-> pos", neg_split[i], pos)
        #if end+1 in pos_starts:
        #    found_neg_pos.append([neg_split[i], pos_split[pos_starts.index(end+1)]])

    
    if debug:
        make_plots = True
    elif plots and (found_neg_pos or found_pos_neg):
        make_plots = True
    else:
        make_plots = False

    if make_plots:
        fig, ax = plt.subplots(2,1)
        if debug:
            # plot the iqrm +ve and -ve outliers, in alpha=0.2 dodgerblue and red
            for x in iqrm_pos:
                ax[0].axvline(x, c='red', alpha=0.3, linewidth=0.5, zorder=1)
            for x in iqrm_neg:
                ax[0].axvline(x, c='dodgerblue', alpha=0.3, linewidth=0.5, zorder=1)

        ax[0].plot(arr1d, linewidth=0.5, zorder=2, c="k")
        ax[1].plot(np.gradient(arr1d), linewidth=0.5, zorder=2, c="k")


    if not found_pos_neg and not found_neg_pos:
        if make_plots:
            if return_plots:
                return False, False, fig, ax
            plt.show()
        return False, False, None, None  # []

    logging.debug("Running checks on potential steps")

    # setup to check means before and after potential steps:
    # was not helpful
    if mean_check_thresh is not None:
        if np.ma.is_masked(arr1d):
            masked_arr1d = copy.deepcopy(arr1d)
        else:
            masked_arr1d = np.ma.array(arr1d, mask=[False]*len(arr1d))
        masked_arr1d_wiqrm = copy.deepcopy(masked_arr1d)
        masked_arr1d_wiqrm.mask[pks_pos] = True
        masked_arr1d_wiqrm.mask[pks_neg] = True

    # check transition areas against scipy's find_peaks in the gradient
    grad = np.gradient(arr1d)
    pks_pos, _ = find_peaks(grad, prominence=prom)
    pks_neg, _ = find_peaks(-grad, prominence=prom)
    all_peaks = np.array([*pks_pos, *pks_neg])
    # plot all peaks found in magenta/green alpha=0.2
    if make_plots:
        for pk in pks_pos:
            ax[1].axvline(pk, c="darkred", alpha=1, linewidth=0.6, zorder=1)
        for pk in pks_neg:
            ax[1].axvline(pk, c="mediumblue", alpha=1, linewidth=0.6, zorder=1)

    # sorted combo has elements like:
    # [[iqrm, idx, before, transition], [iqrm, idx, after, transition], sign]
    # where sign is +1 for neg_pos and -1 for pos_neg
    # and is sorted based on the first index before transition
    # NB overlaps are possible because the +ve and -ve iqrm runs can overlap
    sorted_combo = combine_posneg_negpos(found_pos_neg, found_neg_pos)
    # Run some checks to see if it's likely this is actually a step

    sign_to_name = {+1: "neg->pos", -1: "pos->neg"}
    sign_to_col = {+1: "red", -1: "dodgerblue"}
    sign_to_want = {+1: "+ve", -1: "-ve"}

    overrule_region = int(len(arr1d)/10)
    step = False
    for i, [x,y,sign] in enumerate(sorted_combo):
        logging.debug(f"Checking {sign_to_name[sign]}: {x}, {y}")

        select_x = x
        select_y = y
        # if want to restrict or expand where look for peaks to a set distance around the step
        if region_buffer is not None:
            if len(x) > region_buffer:
                select_x = x[-region_buffer:]
            if len(y) > region_buffer:
                select_y = y[:region_buffer]


        # find_peaks check
        # look for a peak of the expected sign in the region around the step
        select = [*select_x, *select_y]
        pk_pos = [pk for pk in pks_pos if pk in select]
        pk_neg = [pk for pk in pks_neg if pk in select]

        if sign >0:
            want = pk_pos
            dontwant = pk_neg
        else:
            want = pk_neg
            dontwant = pk_pos
        
        if want and not dontwant:
            logging.debug(f"find_peaks: STEP: only {sign_to_want[sign]} peaks found at {want}")
            # check if there are no other close peaks and in which case overrule companion check
            for pk in pk_pos:
                if not any_indices_with_region(pk, all_peaks, overrule_region):
                    logging.info(f"Overruling companion test due to peak at {pk} with not others found within +-{overrule_region}")
                    superyes = True
            step = True
        elif pk_neg and pk_pos:
            logging.debug("find_peaks: NOT_STEP: +ve and -ve peaks found at", pk_pos, pk_neg)
        elif dontwant:
            logging.debug(f"find_peaks: WEIRD: only {sign_to_want[-sign]} peaks found at {dontwant}")
        else:
            logging.debug("find_peaks: NO PEAKS")


        # means check
        if mean_check_thresh is not None:
            slc_pre = slice(None, x[0])
            slc_post = slice(y[-1]+1, None)

            mn_check1, mn_diff1 = check_mean(masked_arr1d, sign, slc_pre, slc_post, mean_check_thresh, fn=np.ma.mean)
            mn_check2, mn_diff2 = check_mean(masked_arr1d_wiqrm, sign, slc_pre, slc_post, mean_check_thresh, fn=np.ma.mean)
            md_check1, md_diff1 = check_mean(masked_arr1d, sign, slc_pre, slc_post, mean_check_thresh, fn=np.ma.median)
            md_check2, md_diff2 = check_mean(masked_arr1d_wiqrm, sign, slc_pre, slc_post, mean_check_thresh, fn=np.ma.median)

            if len(sorted_combo) == 1:
                logging.debug(f"between mean/median checks invalid: i={i}, aka same as checks 1 and 2")
                mn_check3, mn_check4, md_check3, md_check4 = [None, None, None, None]
                mn_diff3, mn_diff4, md_diff3, md_diff4 = [None, None, None, None]
            else:
                if i == 0:
                    i_pre = 0
                else:
                    i_pre = sorted_combo[i-1][1][-1] + 1
                if i == len(sorted_combo) -1:
                    i_post = len(arr1d)
                else:
                    i_post = sorted_combo[i+1][0][0]

                if i_pre >= x[0]:
                    logging.debug("between mean/median checks invalid: no space before")
                    mn_check3, mn_check4, md_check3, md_check4 = [None, None, None, None]
                    mn_diff3, mn_diff4, md_diff3, md_diff4 = [None, None, None, None]
                elif i_post <=  y[-1] + 1:
                    logging.debug("between mean/median checks invalid: no space after")
                    mn_check3, mn_check4, md_check3, md_check4 = [None, None, None, None]
                    mn_diff3, mn_diff4, md_diff3, md_diff4 = [None, None, None, None]
                else:
                    slc_between_pre = slice(i_pre, x[0])
                    slc_between_post = slice(y[-1]+1, i_post)

                    mn_check3, mn_diff3 = check_mean(masked_arr1d, sign, slc_between_pre, slc_between_post, mean_check_thresh, fn=np.ma.mean)
                    mn_check4, mn_diff4 = check_mean(masked_arr1d_wiqrm, sign, slc_between_pre, slc_between_post, mean_check_thresh, fn=np.ma.mean)
                    md_check3, md_diff3 = check_mean(masked_arr1d, sign, slc_between_pre, slc_between_post, mean_check_thresh, fn=np.ma.median)
                    md_check4, md_diff4 = check_mean(masked_arr1d_wiqrm, sign, slc_between_pre, slc_between_post, mean_check_thresh, fn=np.ma.median)

            logging.info(f"MEAN/MED checks:")
            logging.info(mn_check1, md_check1, mn_check2, md_check2, mn_check3, md_check3, mn_check4, md_check4)
            diff_fstring = f"{mn_diff1:.2f}, {md_diff1:.2f}, {mn_diff2:.2f}, {md_diff2:.2f}"
            if [mn_check3, mn_check4, md_check3, md_check4] != [None, None, None, None]:
                diff_fstring += f", {mn_diff3:.2f}, {md_diff3:.2f}, {mn_diff4:.2f}, {md_diff4:.2f}"
            else:
                diff_fstring += f", {mn_diff3}, {md_diff3}, {mn_diff4}, {md_diff4}"
            logging.info(diff_fstring)

        if make_plots:
            # mark on gradient plot where looking for a peak and what sign
#            ax[1].axvline(x[-1], c=sign_to_col[sign])

            # mark on data plot the region where actually searched for a peak
            # go one alpha=0.2 more intense if region IDed as a positive followed by a negative or vice versa
            # go one alphs=0.2 more intense if it was the region where actually searched for peaks in gradient
            ax[0].axvspan(x[0], x[-1], color=sign_to_col[-sign], alpha=0.5, linewidth=0.5, zorder=1)
            ax[0].axvspan(select_x[0], select_x[-1], color=sign_to_col[-sign], alpha=0.5, linewidth=0.5, zorder=1)
            ax[0].axvspan(y[0], y[-1], color=sign_to_col[sign], alpha=0.5, linewidth=0.5, zorder=1)
            ax[0].axvspan(select_y[0], select_y[-1], color=sign_to_col[sign], alpha=0.5, linewidth=0.5, zorder=1)

            # mark region searched for peaks on gradient plot
            ax[1].axvspan(select_x[0], select_y[-1], color=sign_to_col[sign], alpha=0.5, linewidth=0.5, zorder=1)

#            for xx in select_x:
#                ax[0].axvline(xx, c=sign_to_col[-sign], alpha=0.2, linewidth=0.5, zorder=1)
#            for yy in select_y:
#                ax[0].axvline(yy, c=sign_to_col[sign], alpha=0.2, linewidth=0.5, zorder=1)
#            if region_buffer is not None:
#                if select_x != x:
#                    ax[0].axvline(select_x[0], c='black')
#                if select_y != y:
#                    ax[0].axvline(select_y[-1], c='black')
    
    # now defunct
    #if step and len(all_peaks) == 1:
    #    superyes = True

    if make_plots:
        if return_plots:
            return step, superyes, fig, ax
        plt.show()

    #if out_neg_pos or out_pos_neg:
    #    return out_pos_neg, out_neg_pos
    #else:
    #    return []
    return step, superyes

def is_there_a_peak(diffarr, peak, slc, rng, found_peaks, med, std, sign=1, debug_plot=False):
    """sign is the sign of peak, must be +-1
    so if sign is 1, peak is positive and looking for a negative companion peak in slc"""
    if sign not in [1,-1]:
        raise AttributeError(f"is_there_a_peak: sign ({sign}) must be +1 or -1")

    tmpx = np.arange(len(diffarr))
    if debug_plot:
        plt.plot(tmpx[slc], diffarr[slc])
        plt.plot(tmpx[peak], diffarr[peak])
        plt.axvline(tmpx[peak])

    # if a peak was found by found_peaks then this is easy
    for x in rng:
        if x in found_peaks:
            if debug_plot:
                plt.suptitle(f"{peak}: peak found through search")
                plt.axvline(tmpx[x], c='magenta')
            return True

    # if >half of the slice being searched is masked then you can't know
    # but do want to disregard peak so return true
    # (the half is a bit of an arbitrary threshold)
    if np.ma.is_masked(diffarr):
        if diffarr.mask[slc].sum() >= diffarr[slc].size/2:
            if debug_plot:
                plt.suptitle(f"{peak}: {diffarr.mask[slc].sum() / diffarr[slc].size} of slice masked")
            return True
    
    # look for a peak > 0.8*peak - std from mean
    peak_to_med = np.abs(diffarr[peak] - med)
    expect_height = 0.8*peak_to_med - std  # elsewhere ignore peaks under 3 std so this should always be +ve
    if sign > 0:
        if (diffarr[slc] < med - expect_height).any():
            if debug_plot:
                plt.axhline(med - expect_height)
                plt.suptitle(f"{peak}: found an expected value for a peak")
            return True
    else:
        if (diffarr[slc] > med + expect_height).any():
            if debug_plot:
                plt.axhline(med + expect_height)
                plt.suptitle(f"{peak}: found an expected value for a peak")
            return True
            
    return False
    
def get_search_indices(arr, pk, where=None):
    """
    Helper function for check_peaks_have_companions
    Gets range in which to search for companion peak

    returns range1, slice1
    where
    range1 = the range in which to search
    slice1 = the corresponding slice to range1
    or if it ran into an issue e.g. the peak is at the edge of the array, it returns None
    """
    if where == "after":
        zc = find_next_zero_crossing(arr, pk)
        if zc == len(arr) -1:
            return None
        w = zc - pk 
        #print(pk, zc, w)
        llim = pk  # zc
        ulim = min(len(arr), pk+3*w+1)
        return range(llim, ulim), slice(llim, ulim)

    elif where == "before":
        zc = find_prev_zero_crossing(arr, pk)
        if zc == 0:
            return None
        w = pk - zc
        #print(pk, zc, w)

        ulim = pk + 1  #zc + 1
        llim = max(0, pk-3*w)
        return range(ulim-1, llim-1, -1), slice(ulim-1, llim-1, -1)
    else:
        raise AttributeError(f"get_search_indices: <where> ({where}) must be either 'after' or 'before'")
    
def find_next_zero_crossing(arr, id):
    init = arr[id]

    if init > 0:
        for i in range(id+1, len(arr)):
            if arr[i] <= 0:
                return i
        return i

    if init < 0:
        for i in range(id+1, len(arr)):
            if arr[i] >= 0:
                return i
        return i

def find_prev_zero_crossing(arr, id):
    init = arr[id]

    if init > 0:
        for i in range(id-1, -1, -1):
            if arr[i] <= 0:
                return i
        return i

    if init < 0:
        for i in range(id-1, -1, -1):
            if arr[i] >= 0:
                return i
        return i
    
def check_peaks_have_companions2(diffarr1d, high_prom=1, low_prom=0.25, debug_plots=False, ax_debug=False, more_debug_plots=False, ignore_sig_thresh=3):
    """
    Helper function for finding steps in data.
    Meant to be run on a 1D array resulting from an np.diff
    diffarr1d = np.diff(arr1d)

    Idea is to check that all positive peaks are followed by a negative peak,
    and likewise that all negative peaks are rpeceded by a positive peak

    The above will be true for peaks in the original data but not true for a step.

    Does this by using scipy.signal.find_peaks to look for all +ve peaks with prominence > <high_prom>
    checks for nearby negative peaks in a search using prominence > <low_prom>
    Then vice versa for negative high-promimence peaks and positive low-prominence peaks
    
    "nearby" is determined by the distance to the  nearest zero-crossing

    doesn't look if a positive peak is located at len(diffarr1d) -1 or len(diffarr1d) -2
    or if a negative peak is located at 0 or 1

    ignore_sig_thresh=2, means don't look for companions for any peak <2 std from median
    (where median and std calculated using points within 3 std of the original data's median)

    debug_plots Plots 2 axes, showing a) the +ve peaks found and the -ve peaks being searched for companions and 
                                      b) the -ve peaks found and the +ve peaks being searched for companions
    ax_debug is the axes on which to plot these, (2,1) recommended, but at minimum ax_debug[0] and ax_debug[1] must be valid.

    more_debug_plots passes in to is_there_a_peak
    (. . . using debug_plots=True, more_debug_plots=True AND passing in someting to ax_debug will probably do weird things)


    If all peaks have companions this function returns:
        True, None, None
    If this function finds a peak without a companion it returns:
        False, peak, slc
        where peak is the index of the peak, and slc is the slice of diffarr1d containing the peak itself and the area searched for a companion
    """
    tmparr = copy.deepcopy(diffarr1d)
    tmparr -= tmparr.mean()

    # find large peaks
    peaks_pos, props = find_peaks(tmparr, prominence=high_prom)
    peaks_neg, props = find_peaks(-tmparr, prominence=high_prom)

    peaks_pos = list(peaks_pos)
    peaks_neg = list(peaks_neg)


    # find smaller peaks too
    peaks_pos_lowerprom, props = find_peaks(tmparr, prominence=low_prom)
    peaks_neg_lowerprom, props = find_peaks(-tmparr, prominence=low_prom)

    if debug_plots:
        if ax_debug is None:
            fig, ax = plt.subplots(2,1)
        else:
            ax = ax_debug
        ax[0].plot(tmparr, c="k", linewidth=0.5, zorder=3)
        for x in peaks_pos:
            ax[0].axvline(x, c='red', linewidth=0.5, zorder=2)
        for x in peaks_neg_lowerprom:
            ax[0].axvline(x, c='dodgerblue', linewidth=0.5, zorder=2)
        ax[0].set_title("Positive peaks, looking for companions in -ve peaks")

        ax[1].plot(tmparr, c="k", linewidth=0.5, zorder=3)
        for x in peaks_neg:
            ax[1].axvline(x, c='dodgerblue', linewidth=0.5, zorder=2)
        for x in peaks_pos_lowerprom:
            ax[1].axvline(x, c='red', linewidth=0.5, zorder=2)
        ax[1].set_title("Negative peaks, looking for companions in +ve peaks")
        if ax_debug is None:
            plt.show()

    # compute median and std only taking points within 3std of the median
    # use this to filter out some initial peaks likely don't care about
    medtmp = np.ma.median(tmparr)
    stdtmp = np.ma.std(tmparr)
    if np.ma.is_masked(tmparr):
        msk = tmparr.mask
    else:
        msk = np.zeros_like(tmparr, dtype=bool)
    msk2 = (tmparr < medtmp - 3*stdtmp) | (tmparr > medtmp + 3*stdtmp)
    med = np.ma.median(np.ma.array(tmparr, mask=(msk|msk2)))
    std = np.ma.std(np.ma.array(tmparr, mask=(msk|msk2)))

    for pk in peaks_pos:
        if pk == len(tmparr) -1 or pk ==len(tmparr) -2:
            continue
        if np.ma.is_masked(tmparr) and tmparr.mask[pk]:
            continue
        if tmparr[pk] < med + ignore_sig_thresh*std:  # probably don't care about peaks this small/baselines this noisy
            continue
        found_companion = False
        #print(f"--\n{pk}\t{tmparr[pk]}\n--")
        # use zero crossing point to find range in which to search for a positive peak
        rngslc = get_search_indices(tmparr, pk, 'after')
        if rngslc is None:
            continue
        else:
            rng, slc = rngslc

        found_companion = is_there_a_peak(tmparr, pk, slc, rng, peaks_neg_lowerprom, med, std, sign=1, debug_plot=more_debug_plots)
        if more_debug_plots:
            plt.show()

        if not found_companion:
            return False, pk, slc, (med, std)

    #print("---\nnegative peakes\n---")

    for pk in peaks_neg:
        if pk == 0 or pk ==1:
            continue
        if np.ma.is_masked(tmparr) and tmparr.mask[pk]:
            continue
        if tmparr[pk] > med - ignore_sig_thresh*std:  # probably don't care about peaks this small/baselines this noisy
            continue
        found_companion = False
        #print(f"--\n{pk}\t{tmparr[pk]}\n--")
        # use zero crossing point to find range in which to search for a positive peak
        rngslc = get_search_indices(tmparr, pk, 'before')
        if rngslc is None:
            continue
        else:
            rng, slc = rngslc
        
        found_companion = is_there_a_peak(tmparr, pk, slc, rng, peaks_pos_lowerprom, med, std, sign=-1, debug_plot=more_debug_plots)
        if more_debug_plots:
            plt.show()

        if not found_companion:
            return False, pk, slc, (med, std)
        

    return True, None, None, None




# ## Stuff that needs to be input to the code

#args_in = [
    #"test_data/SURVEYv0_point44_DM34_59831_pow_fdp_rfifind.mask",
    #"test_data/gulp76800/SURVEYv0_point44_DM34_59831_pow_fdp_stats.npz",
#    "test_data/SURVEYv0_point97_DM34_59634_pow_fdp_rfifind.mask",
#    "test_data/SURVEYv0_point97_DM34_59634_pow_fdp_stats.npz",
    #"--option", "0,1,2,3,4,5,6",
#    "--option", "0,1,2,3,4,6",
#    "--outfilename", "pipeline_rfifind.mask",
#    "--show",
#    "--ignorechans",
#    "1023,1000,970,965,957,941,928,917,914,910,909,908,906,905,904,903,902,901,900,899,898,897,896,895,894,893,892,891,890,889,888,887,886,885,882,881,880,879,878,877,876,875,874,873,872,871,870,869,868,866,865,862,861,860,859,858,857,856,855,854,853,852,851,850,849,848,847,846,845,843,840,835,832,826,825,821,813,757,755,750,731,728,712,710,706,704,702,699,696,691,689,688,685,682,672,664,662,661,658,656,645,617,616,611,610,606,600,594,592,586,585,584,579,578,577,576,575,574,573,572,571,570,569,568,567,566,564,563,562,561,560,559,558,556,544,536,533,517,512,509,501,493,485,475,471,469,468,467,466,465,464,462,461,460,459,458,457,456,455,451,448,440,439,438,436,434,433,432,430,429,428,425,424,421,397,385,384,380,379,378,368,363,345,339,333,332,325,323,317,309,302,301,285,282,280,272,269,267,259,258,257,256,253,251,246,241,235,234,232,225,224,222,221,215,206,189,181,175,141,125,109,104,101,85,83,77,69,45,24,18,13,5"
#]


# ## Parser
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Run some rfi mitigation, will write an rfifind-style mask and inf file. Will also write a .ignorechans file a list of (presto format) ignorechans to pass onto the next stage"
)

def _check_mask_name(string):
    print(string)
    if string[-12:] != "rfifind.mask":
        raise argparse.ArgumentTypeError(f'Invalid argument for maskfile/outfilename ({string}): must end in rfifind.mask')
    return string

def _check_m_frac_threshold(string):
    value = float(string)
    if value < 0 or value > 1:
        raise argparse.ArgumentTypeError(f'Invalid argument for m_frac_threshold ({value}): cannot be <0 or >1')
    return value

def _check_rfac(string):
    value = float(string)
    if value <= 0:
        raise argparse.ArgumentTypeError(f'Invalid argument for rfac ({value}): must be a positive number')
    return value

parser.add_argument(
    "maskfile",
    type=_check_mask_name,
    help=".mask file output from rfifind. Must also have a corresponding .stats and .inf file"
)
parser.add_argument(
    "--extra_stats_file",
    type=str,
    help="""npz file containing extra stats from the fdp process. Must contain n, num_unmasked_points, s1, s2, gulp.
    n = number of points in each interval, useful as the last one is often shorter
    num_unmasked_points = number of points that went into the s1 ans s2 calculation
    s1 = sum
    s1 = sum of the squares
    gulp = gulp used to make stats"""
)

parser.add_argument(
    "--option",
    type=str,
    default="1,2,3,4,5,6,8",
    help="""Which masking options to use. Comma separated string, will be run from 0 up regardless of order passed in.
    Default is '1,2,3,4,5,6,8'

    Alway included: basic processing: ignorechans, anywhere the std is 0, where the number of unmasked points is < set threshold,
    1: bad ints: do an interval high fraction cut on a 2D iqrm run on the means, aka flag intervals which are bad for many channels,
    2: find channels containing steps,
    3: find channels where the std of the means is a highly significant outlier (iqrm threshold=50),
    4: a 2D iqrm run on the means, both looking for bad intervals and channels,
    5: a 2D iqrm run on the variance, both looking for bad intervals and bad channels,
    6: first a 1D iqrm on the medians (along both axes) then a 2D iqrm run on the pow_stats,
    7: a 2D iqrm run on the generarlized spectral kurtosis statistic, looking for bad channels only,
    8: a high fraction cut. The threshold for the channels is varied so at most it masks an additional 0.1%% of the data,
    """
)

group = parser.add_mutually_exclusive_group(required=False)

group.add_argument(
    "--outfilename",
    type=_check_mask_name,
    help="filename to write output mask to"
)

group.add_argument(
    "-o",
    "--overwrite",
    action='store_true',
    help="Overwrite the mask with the new one"
)

parser.add_argument(
    "--m_frac_threshold",
    default=0.5,
    type=_check_m_frac_threshold,
    help="If num_unmasked_points in a block is < this fraction, mask this int,chan block"
)

parser.add_argument(
    "--rfac",
    default=10,
    type=_check_rfac,
    help="fraction of the band to use for iqrm window, e.g. if nchan=1024 and pass in 16, window will be 64 channels"
)

parser.add_argument(
    "--show",
    action='store_true',
    help="Show plots rather than saving them as a pdf"
)

parser.add_argument(
    "--dont_flip_band",
    action='store_true',
    help="""If this flag is passed in the channel order in <extra_stats_file> will NOT be flipped.
    If it came from an fdp script channels will be in sigproc order which is the opposite from presto order and thus most of the time they'll need to be flipped"""
)

parser.add_argument(
    "--ignorechans",
    type=str,
    default="",
    help="string of ignorechans in presto format, aka channel 0 is the lowest freq channel, and it's a string separated by commas"
)

parser.add_argument(
    "--problem_frac",
    type=float,
    default=0.6,
    help="If masking fraction goes above this threshold then there is a problem, skip whatever step did this"
)

parser.add_argument(
    "--include_rfifind",
    action='store_true',
    help="Include the rfifind mask in the final mask"
)

parser.add_argument(
    "--log", type=str, help="name of file to write log to", default=None
)

parser.add_argument(
    "-v",
    "--verbose",
    help="Increase logging level to debug",
    action="store_const",
    dest="loglevel",
    const=logging.DEBUG,
    default=logging.INFO,
)

# args = parser.parse_args(args_in)

if __name__ == "__main__":
    args = parser.parse_args()

    if args.log is not None:
        logging.basicConfig(
            filename=args.log,
            filemode="a",
            format="%(asctime)s %(levelname)s:%(message)s",
            datefmt="%d-%b-%y %H:%M:%S",
            level=args.loglevel,
        )
    else:
        logging.basicConfig(
            format="%(asctime)s %(levelname)s:%(message)s",
            datefmt="%d-%b-%y %H:%M:%S",
            level=args.loglevel,
            stream=sys.stdout,
        )

    logging.info("rfi_pipeline initialized with arguments:")
    logging.info(args)

    maskfile = args.maskfile
    extra_stats_fn = args.extra_stats_file
    rfac = args.rfac
    m_frac_threshold = args.m_frac_threshold

    if args.ignorechans == "":
        ignorechans = []
    else:
        ignorechans = [int(x) for x in args.ignorechans.split(",")]

    opts = [0]
    if args.option:
        opts.extend([int(x) for x in args.option.split(",") if int(x)>0])

    optstr="".join([f"{x}" for x in opts])
    basename_mask = maskfile[:maskfile.rfind("_rfifind.mask")]
    if args.overwrite:
        outfilename = maskfile
    elif args.outfilename is None:
        outfilename = basename_mask + "_" + optstr + "_rfifind.mask"
    else:
        outfilename = args.outfilename
    logging.info(f"New mask will be written to: {outfilename}")

    if args.show:
        p = None
    else:
        plotfname = f"rfipipeline_plots_{basename_mask}_{optstr}.pdf"
        p = PdfPages(plotfname)
        logging.info(f"Plots will be written to {plotfname}")

    opt_dict = {
        0: "basic processing: ignorechans, anywhere the std is 0, where the number of unmasked points is < set threshold",
        1: "bad ints: do an interval high fraction cut on a 2D iqrm run on the means, aka flag intervals which are bad for many channels",
        2: "find channels containing steps",
        3: "find channels where the std of the means is a highly significant outlier (iqrm threshold=50)",
        4: "a 2D iqrm run on the means, both looking for bad intervals and channels",
        5: "a 2D iqrm run on the variance, both looking for bad intervals and bad channels",
        6: "first a 1D iqrm on the medians (along both axes) then a 2D iqrm run on the pow_stats",
        7: "a 2D iqrm run on the generarlized spectral kurtosis statistic, looking for bad channels only",
        8: "a high fraction cut. The threshold for the channels is varied so at most it masks an additional 0.1% of the data",
    }
    logging.info(f"Options selected:")
    for x in opts:
        logging.info(f"\t{x}: {opt_dict[x]}")


    # ## Load files

    rfimask = rfifind.rfifind(maskfile)
    logging.info(f"loaded mask from {maskfile}")

    use_rfifind_meanstd = False
    if extra_stats_fn is None:
        if 7 in opts:
            raise RuntimeError(f"GSK requires extra_stats_fn as it contains s1 and s2")
        else:
            use_rfifind_meanstd = True

    if use_rfifind_meanstd:
        extra_stats_gulp = rfimask.ptsperint
        M = np.zeros_like(rfimask.mask)
        means = rfimask.avg_stats
        var = rfimask.std_stats**2
    else:
        extra_stats = np.load(extra_stats_fn, allow_pickle=True)
        if args.dont_flip_band:
            M = extra_stats["num_unmasked_points"][:,:]
            s1 = extra_stats["s1"][:,:]
            s2 = extra_stats["s2"][:,:]
        else:
            # flip everything to presto channel convention
            logging.info("Reversing channel order in extra stats to match presto convention")
            M = extra_stats["num_unmasked_points"][:,::-1]
            s1 = extra_stats["s1"][:,::-1]
            s2 = extra_stats["s2"][:,::-1]


        N = extra_stats["n"]
        extra_stats_gulp = extra_stats["gulp"]
        # something really weird happens with the means and var if you make them with masked arrays. Make and mask afterwards
        means = s1/M
        var = (s2 - s1**2/M)/M

    # ## make some other parameters
    r = rfimask.nchan/rfac
    r_int = rfimask.nint/rfac


    # ## Pipeline

    masks = {'rfifind': rfimask.mask}
    masks_exstats = {'rfifind': reshape_rfifind_mask(M.shape, rfimask.mask, extra_stats_gulp, rfimask.ptsperint)}

    # Stage 0 is non-optional
    # ### var==0 mask
    logging.info(f"Ignoring {len(ignorechans)}/{rfimask.nchan} channels")
    logging.info(f"Getting zeros mask")
    fig0, ax0 = plt.subplots()
    m0 = get_zeros_mask_alt(var, ignorechans=ignorechans, verbose=True, ax=ax0)
    output_plot(fig0, pdf=p)

    # ### M mask (where num_points_unmasked from fdp is under some threshold)
    # cut off anywhere where <0.5 of the gulp was nonzero
    if extra_stats_fn is None:
        logging.info("No extra stats file: skipping the mask based on a threshold number of unmasked points")
        mmask = np.zeros_like(m0)
    else:
        logging.info(f"Ignoring anywhere where the fraction of points used to calculate the stats was < {m_frac_threshold}")
        mmask = (M.T < m_frac_threshold*N).T
    # right the gulp is different for fdp and rfifind. damn. will have to record that in extra_stats


    base_mask_exstats = m0 | mmask
    # base_mask_exstats = mmask
    base_mask = reshape_extra_stats_mask(rfimask.mask.shape, base_mask_exstats, extra_stats_gulp, rfimask.ptsperint)
    logging.info(f"Reshaped base_mask_exstats from {base_mask_exstats.shape} to {base_mask.shape} for use with the rfifind mask")
    logging.info(f"ignorechans, std_stats=0 and lots of 0 data alone mask out {masked_frac(base_mask_exstats)} of data")

    masks[0] = copy.deepcopy(base_mask)
    masks_exstats[0] = copy.deepcopy(base_mask_exstats)
    # # ### Add initial rfifind mask and do a hard threshold cut
    # base_mask = base_mask | rfimask.mask
    # logging.info(f"+rfifind mask  masks out {masked_frac(base_mask)} of data")

    # fig01, ax01 = plt.subplots(1,2)
    # mcut = cut_off_high_fraction(base_mask, cumul_threshold_chans=1, cumul_threshold_ints=1, plot_diagnostics=True, ax=None, axhard=ax01) #, debug=True)
    # base_mask = base_mask | mcut
    # output_plot(fig01, pdf=p)

    # fig1, ax1 = plt.subplots()
    # ax1.imshow(base_mask.T, aspect='auto', origin='lower')
    # fig1.suptitle("base mask: std_stats==0, large fraction of masked data within an interval, rfifind mask")
    # output_plot(fig1, pdf=p)

    logging.info(f"0: base mask zaps {masked_frac(base_mask)} of data")


    # make gsk here for same reason if need it
    if 7 in opts:
        logging.info("Making the generalized spectral kurtosis statistic for the future, estimating d from np.ma.median(means**2/var, axis=channel)")
        all_deltas = means**2/var
        all_deltas[base_mask_exstats] = np.nan
        delta = np.nanmedian(all_deltas, axis=1)
        gsk_d_estimate = (((M.T * delta.T + 1) / (M.T - 1)) * (M * (s2 / s1**2) - 1).T).T
        gsk_d_estimate_masked = np.ma.array(gsk_d_estimate, mask=base_mask_exstats)
        #gsk_d_estimate_masked.mask[np.isnan(gsk_d_estimate)] = True
        # shouldn't need the above masking step, but let's check
        # (also it changes base_mask_exstats itself) (and base_mask if the shapes were the same too)
        assert not np.isnan(gsk_d_estimate_masked).any()

    # ### Now have base_mask should be using as minimum input for all other steps
    if not use_rfifind_meanstd:
        M = np.ma.array(M, mask=base_mask_exstats)
        s1 = np.ma.array(s1, mask=base_mask_exstats)
        s2 = np.ma.array(s2, mask=base_mask_exstats)
    means = np.ma.array(means, mask=base_mask_exstats)
    var = np.ma.array(var, mask=base_mask_exstats)

    if masked_frac(base_mask) >= args.problem_frac:
        logging.info("Something went wrong at stage 0, making summary plots and exiting")

        make_summary_plots(base_mask, base_mask_exstats, rfimask, means, var, p, title_insert="ERROR stage 0")
        if p is not None:
            logging.info("Writing pdf")
            p.close()
        logging.error("Something went horribly wrong at stage 0")
        sys.exit(1)

    base_ignorechans = get_ignorechans_from_mask(base_mask_exstats)


    # plot reference plots
    logging.info("Plotting reference plots, masked by the base mask must apply + original mask from rfifind")
    rfimask_mask_exstats = reshape_rfifind_mask(M.shape, rfimask.mask, extra_stats_gulp, rfimask.ptsperint)

    fig_ref_means, ax_ref_means, cbar_ref_means = plot_map_plus_sums(means.data, base_mask_exstats|rfimask_mask_exstats, returnplt=True)
    fig_ref_means.suptitle("Means rfimask + base mask")
    output_plot(fig_ref_means, pdf=p)

    fig_ref_var, ax_ref_var, cbar_ref_var = plot_map_plus_sums(var.data, base_mask_exstats|rfimask_mask_exstats, returnplt=True)
    fig_ref_var.suptitle("Var rfimask + base mask")
    output_plot(fig_ref_var, pdf=p)

    pow_stats_plot_mask = base_mask|rfimask.mask
    if (rfimask.pow_stats[-1,:] == 1).all():
        logging.info("Weird final interval stats for rfifind, it will be masked in pow_stats plot, but not in the mask itself")
        pow_stats_plot_mask[-1,:] = True
    fig_ref_pow, ax_ref_pow, cbar_ref_pow = plot_map_plus_sums(rfimask.pow_stats, pow_stats_plot_mask, returnplt=True)
    fig_ref_pow.suptitle("pow_stats rfimask + base mask")
    output_plot(fig_ref_pow, pdf=p)
    del pow_stats_plot_mask

    if 1 in opts:
        logging.info("1: Looking for bad intervals")

        mask_means_2diqrm_int1 = run_iqrm_2D(means.filled(np.nan), base_mask_exstats, 0,r, ignorechans=base_ignorechans, threshold=5)
        mask_means_2diqrm_int2 = run_iqrm_2D(-means.filled(np.nan), base_mask_exstats, 0,r, ignorechans=base_ignorechans, threshold=5)
        mask_means_2diqrm_int = mask_means_2diqrm_int1|mask_means_2diqrm_int2

        fig_iqrm_means_mask, ax_iqrm_means_mask = plt.subplots()
        plot_mask_comparison(mask_means_2diqrm_int, base_mask_exstats, title="iqrm_2D_means (bad ints) - base", ax=ax_iqrm_means_mask)

        fig_badint_hf, ax_badint_hf = plt.subplots(1,2)
        int_frac_msk = cut_off_high_fraction(mask_means_2diqrm_int, hard_thresh_chans=1, cumul_threshold_chans=1, cumul_threshold_ints=1, ax=None, axhard=ax_badint_hf)    
        

        bad_ints = get_ignoreints_from_mask(mask_means_2diqrm_int|base_mask_exstats|int_frac_msk)
        preexisting_bad_ints = get_ignoreints_from_mask(base_mask_exstats)
        bad_ints = sorted(list(set(bad_ints) - set(preexisting_bad_ints)))

        logging.info(f"Found new bad ints: {bad_ints}")
        for bad_int in bad_ints:
            ax_iqrm_means_mask.annotate('', 
                xy=(bad_int, 0), 
                xytext=(bad_int, -5), 
                arrowprops = dict(facecolor='red', edgecolor='red', arrowstyle='->')
            )

        output_plot(fig_iqrm_means_mask, pdf=p)
        output_plot(fig_badint_hf, pdf=p)

        mask_bad_ints = np.zeros_like(base_mask)
        mask_bad_ints_exstats = np.zeros_like(base_mask_exstats)
        mask_bad_ints[bad_ints,:] = True
        mask_bad_ints_exstats[bad_ints,:] = True

        masks[1] = mask_bad_ints
        masks_exstats[1] = mask_bad_ints_exstats
        
        logging.info("Testing and updating *base mask* this will be fed onto all following other stages")
        base_mask, base_mask_exstats = check_mask_and_continue(
            base_mask, base_mask_exstats, 
            mask_bad_ints, mask_bad_ints_exstats, 
            args.problem_frac, rfimask, means, var, p, 
            stage=1,
        )
        base_ignorechans = get_ignorechans_from_mask(base_mask_exstats)

        means.mask = base_mask_exstats
        var.mask = base_mask_exstats


    if 2 in opts:
        logging.info("2: Looking for steps")
        tmpx = np.arange(means.shape[0])

        chans = [cc for cc in range(means.shape[1]) if cc not in base_ignorechans]
        chans_w_step = []

        for c in chans:
            has_step, overrule_companion_check, fig_step_iqrm, ax_step_iqrm = find_step(means[:,c], plots=True, prom=0.25, region_buffer=5, return_plots=True)
            if fig_step_iqrm is not None:
                fig_title = f"{c}: iqrm:{has_step}"
                if has_step:
                    logging.info(f"{c}: iqrm found potential step")
                    grad = np.gradient(means[:,c])
                    fig_comp_debug, ax_comp_debug = plt.subplots(2,1)
                    condit, pk, slc, medstd = check_peaks_have_companions2(grad, high_prom=1, low_prom=0.25, debug_plots=True, ax_debug=ax_comp_debug, ignore_sig_thresh=3)

                    fig_title += f" companion_test:{not condit}"
                    if not condit:
                        logging.info("companion test says YES")
                        plt.close(fig_comp_debug)
                        # mark the actual peak the got IDed as corresponding to a step
                        ax_step_iqrm[0].plot(tmpx[slc], means[slc,c], c='black', linewidth=1, zorder=4)
                        ax_step_iqrm[0].axvline(pk, c='black', linestyle='--', linewidth=0.5, zorder=4)
                        ax_step_iqrm[1].plot(tmpx[slc], grad[slc], c='black', linewidth=1, zorder=4)
                        ax_step_iqrm[1].axvline(pk, c='black', linestyle='--', linewidth=0.5, zorder=4)
                        chans_w_step.append(c)
                        
                    else:
                        logging.info("companion test says NO")
                        if overrule_companion_check:
                            logging.info("Overruling companion check as only one close peak was found")
                            fig_title += ", OVERULE companion test"
                            chans_w_step.append(c)



                fig_step_iqrm.suptitle(fig_title)
                output_plot(fig_step_iqrm, pdf=p)
                if (has_step and condit) or overrule_companion_check:
                    fig_comp_debug.suptitle(f"{c}: companion test debug plot")
                    output_plot(fig_comp_debug, pdf=p)


                
        logging.info(f"Found channels with steps: {chans_w_step}")

        mask_step_chans = np.zeros_like(base_mask)
        mask_step_chans_exstats = np.zeros_like(base_mask_exstats)
        mask_step_chans[:,chans_w_step] = True
        mask_step_chans_exstats[:,chans_w_step] = True

        masks[2] = mask_step_chans
        masks_exstats[2] = mask_step_chans_exstats

        logging.info("Testing and updating *base mask* this will be fed onto all following other stages")
        base_mask, base_mask_exstats = check_mask_and_continue(
            base_mask, base_mask_exstats, 
            mask_step_chans, mask_step_chans_exstats, 
            args.problem_frac, rfimask, means, var, p, 
            stage=2,
        )
        base_ignorechans = get_ignorechans_from_mask(base_mask_exstats)

        means.mask = base_mask_exstats
        var.mask = base_mask_exstats

        # NB what the plots show is:
        # top:
        #   means
        #   -ve (dodgerblue) and +ve (red) ints flagged by iqrm
        #   low alpha if not a transition, and alph=1 if there is
        #   if there's a black line that shows region_buffer
        # bottom:
        #   np.gradient of means
        #   red and dodgerblue mark transitions IDed in the top plot
        #     dodgerblue = +ve to -ve and red = -ver to +ve
        #   green and magenta peaks are from find_peaks
        #   whole transition region or transition region limited by region_buffer is searched for find_peaks results
        #     if both signs of peak then iqrm gives False (there is no step)
        #     if onlt one sign of peak is present (and it's the correct sign) the iqrm gives True
        # if orange is present in either plot that shows the peak and region searched for a companion of the opposite sign in the companion check

    # krzy 12 colour
    # cset = '#9F0162', '#009F81', '#FF5AAF', '#00FCCF', '#8400CD', '#008DF9', '#00C2F9', '#FFB2FD', '#A40122', '#E20134', '#FF6E3A', '#FFC33B'[::-1]
    if 3 in opts:
        logging.info("3: Looking for channels where std of the means is a highly significant outlier")
        thresh=50
        logging.info(f"'highly significant' = iqrm threshold of {thresh}")

        thing = np.ma.std(means, axis=0)
        q,v = iqrm_mask(thing.filled(np.nan), radius=r, threshold=thresh, ignorechans=np.where(thing.mask)[0])
        chns = [c for c in np.where(q)[0] if c not in np.where(thing.mask)[0]]
        logging.info(f"found {len(chns)} channels with an unusually high std(means):")
        logging.info(f"{chns}")
        fig_outlier_std, ax_outlier_std = plt.subplots(2,1)
        ax_outlier_std[0].plot(thing, '-', c="k")
        ax_outlier_std[0].plot(np.ma.array(thing.data, mask=((~q)|thing.mask)), 'o', markersize=3, c="red", label="Outlier")
        ax_outlier_std[0].legend()
        ax_outlier_std[0].set_xlabel("Channel")
        ax_outlier_std[0].set_ylabel("Standard Deviation\nof Means")
        offset=5
        for i,c in enumerate(chns):
            ax_outlier_std[1].plot(means[:,c] - np.ma.median(means[:,c]) + i*offset)
        ax_outlier_std[1].set_xlabel("Interval")
        ax_outlier_std[1].set_ylabel("Means (offset)")
        ax_outlier_std[0].set_title(f"Channels with outlier std(means):\n{chns}")
        output_plot(fig_outlier_std, pdf=p)

        mask_outlier_std = np.zeros_like(base_mask)
        mask_outlier_std_exstats = np.zeros_like(base_mask_exstats)
        mask_outlier_std[:,chns] = True
        mask_outlier_std_exstats[:,chns] = True

        masks[3] = mask_outlier_std
        masks_exstats[3] = mask_outlier_std_exstats

        logging.info("Testing and updating *base mask* this will be fed onto all following other stages")
        base_mask, base_mask_exstats = check_mask_and_continue(
            base_mask, base_mask_exstats, 
            mask_outlier_std, mask_outlier_std_exstats, 
            args.problem_frac, rfimask, means, var, p, 
            stage=3,
        )
        base_ignorechans = get_ignorechans_from_mask(base_mask_exstats)

        means.mask = base_mask_exstats
        var.mask = base_mask_exstats

    working_mask_exstats = copy.deepcopy(base_mask_exstats)
    working_mask = copy.deepcopy(base_mask)


    if 4 in opts:
        # Run some iqrms
        logging.info("4: Running 2D iqrms on -means and +means, both int-wise and chan-wise")
        mask_means_2diqrm_chan1 = run_iqrm_2D(means.filled(np.nan), base_mask_exstats, 1,r, ignorechans=base_ignorechans, threshold=5)
        mask_means_2diqrm_chan2 = run_iqrm_2D(-means.filled(np.nan), base_mask_exstats, 1,r, ignorechans=base_ignorechans, threshold=5)
        mask_means_2diqrm_chan = mask_means_2diqrm_chan1|mask_means_2diqrm_chan2

        mask_means_2diqrm_int1 = run_iqrm_2D(means.filled(np.nan), base_mask_exstats, 0,r, ignorechans=base_ignorechans, threshold=5)
        mask_means_2diqrm_int2 = run_iqrm_2D(-means.filled(np.nan), base_mask_exstats, 0,r, ignorechans=base_ignorechans, threshold=5)
        mask_means_2diqrm_int = mask_means_2diqrm_int1|mask_means_2diqrm_int2

        mask_means_2diqrm = mask_means_2diqrm_chan|mask_means_2diqrm_int

        masks_exstats[4] = mask_means_2diqrm
        masks[4] =  reshape_extra_stats_mask(rfimask.mask.shape, mask_means_2diqrm, extra_stats_gulp, rfimask.ptsperint)

        fig_iqrm_means_mask, ax_iqrm_means_mask = plt.subplots(2,1)
        plot_mask_comparison(mask_means_2diqrm_chan, base_mask_exstats, title="iqrm_2D_means (bad chans) - base", ax=ax_iqrm_means_mask[0])
        plot_mask_comparison(mask_means_2diqrm_int, base_mask_exstats, title="iqrm_2D_means (bad ints) - base", ax=ax_iqrm_means_mask[1])
        output_plot(fig_iqrm_means_mask, pdf=p)


        fig_iqrm_means, ax_iqrm_means, cbar_iqrm_means = plot_map_plus_sums(means.data, mask_means_2diqrm|base_mask_exstats, returnplt=True)
        fig_iqrm_means.suptitle("Means post-2D-iqrm")
        output_plot(fig_iqrm_means, pdf=p)
        logging.info(f"masks {masked_frac(mask_means_2diqrm|base_mask_exstats)}")

        working_mask, working_mask_exstats = check_mask_and_continue(
            working_mask, working_mask_exstats, 
            masks[4], masks_exstats[4], 
            args.problem_frac, rfimask, means, var, p, 
            stage=4,
        )

        del mask_means_2diqrm_chan1
        del mask_means_2diqrm_int1
        del mask_means_2diqrm_chan2
        del mask_means_2diqrm_int2
        del mask_means_2diqrm_chan
        del mask_means_2diqrm_int


    if 5 in opts:
        logging.info("5: Running 2D iqrms")
        # On var (both + and -, both time- and chan-wise)
        logging.info("On -var and +var, both int-wise and chan-wise")
        mask_var_2diqrm_chan1 = run_iqrm_2D(var.filled(np.nan), base_mask_exstats, 1,r, ignorechans=base_ignorechans, threshold=5)
        mask_var_2diqrm_int1 = run_iqrm_2D(var.filled(np.nan), base_mask_exstats, 0,r, ignorechans=base_ignorechans, threshold=5)

        mask_var_2diqrm_chan2 = run_iqrm_2D(-var.filled(np.nan), base_mask_exstats, 1,r, ignorechans=base_ignorechans, threshold=5)
        mask_var_2diqrm_int2 = run_iqrm_2D(-var.filled(np.nan), base_mask_exstats, 0,r, ignorechans=base_ignorechans, threshold=5)

        mask_var_2diqrm_chan = mask_var_2diqrm_chan1|mask_var_2diqrm_chan2
        mask_var_2diqrm_int = mask_var_2diqrm_int1|mask_var_2diqrm_int2

        mask_var_2diqrm = mask_var_2diqrm_chan | mask_var_2diqrm_int

        masks_exstats[5] = mask_var_2diqrm
        masks[5] = reshape_extra_stats_mask(rfimask.mask.shape, mask_var_2diqrm, extra_stats_gulp, rfimask.ptsperint)

        fig_iqrm_var_mask, ax_iqrm_var_mask = plt.subplots(2,1)
        plot_mask_comparison(mask_var_2diqrm_chan, base_mask_exstats, title="iqrm_2D_var (bad chans) - base", ax=ax_iqrm_var_mask[0])
        plot_mask_comparison(mask_var_2diqrm_int, base_mask_exstats, title="iqrm_2D_var (bad ints) - base", ax=ax_iqrm_var_mask[1])
        output_plot(fig_iqrm_var_mask, pdf=p)


        fig_iqrm_var, ax_iqrm_var, cbar_iqrm_var = plot_map_plus_sums(var.data, mask_var_2diqrm|base_mask_exstats, returnplt=True)
        fig_iqrm_var.suptitle("Var post-2D-iqrm")
        output_plot(fig_iqrm_var, pdf=p)
        logging.info(f"masks {masked_frac(mask_var_2diqrm|base_mask_exstats)}")

        working_mask, working_mask_exstats = check_mask_and_continue(
            working_mask, working_mask_exstats, 
            masks[5], masks_exstats[5], 
            args.problem_frac, rfimask, means, var, p, 
            stage=5,
        )

        del mask_var_2diqrm_chan1
        del mask_var_2diqrm_int1
        del mask_var_2diqrm_chan2
        del mask_var_2diqrm_int2
        del mask_var_2diqrm_int
        del mask_var_2diqrm_chan


    if 6 in opts:
        # On pow
        logging.info("6: Running iqrm on median of pow_stats along each axis")
        pow_stats_mask = copy.deepcopy(base_mask)
        if (rfimask.pow_stats[-1,:] == 1).all():
            logging.info("Weird final interval stats for rfifind, will not include the final interval when computing pow_stats masks")
            pow_stats_mask[-1,:] = True
        pow_med_chans = np.ma.median(np.ma.array(rfimask.pow_stats, mask=pow_stats_mask), axis=0)
        pow_med_ints = np.ma.median(np.ma.array(rfimask.pow_stats, mask=pow_stats_mask), axis=1)
        pow_chans, v = iqrm_mask(pow_med_chans.filled(np.nan), threshold=5, radius=r, ignorechans=np.where(pow_med_chans.mask)[0])
        pow_ints, v = iqrm_mask(pow_med_ints.filled(np.nan), threshold=5, radius=r, ignorechans=np.where(pow_med_ints.mask)[0])

        iqrm_med_pow_mask = np.zeros_like(base_mask)
        for i, masked in enumerate(pow_chans):
            if masked:
                iqrm_med_pow_mask[:,i] = True
        for j, masked in enumerate(pow_ints):
            if masked:
                iqrm_med_pow_mask[j,:] = True
        pow_stats_plot_mask = (base_mask|iqrm_med_pow_mask)
        if (rfimask.pow_stats[-1,:] == 1).all():
            pow_stats_plot_mask[-1,:] = True
        fig_iqrm_med_pow, ax_iqrm_med_pow, cbar_iqrm_med_pow = plot_map_plus_sums(rfimask.pow_stats, mask=pow_stats_plot_mask, returnplt=True)
        fig_iqrm_med_pow.suptitle("pow_stats masked by base + iqrm 1D on medians along both axes")

        logging.info("5: Running 2D iqrm on pow_stats both int-wise and chan-wise")
        mask_pow_2diqrm_chan = run_iqrm_2D(rfimask.pow_stats, pow_stats_plot_mask, 1,r, ignorechans=base_ignorechans, threshold=5)
        mask_pow_2diqrm_int = run_iqrm_2D(rfimask.pow_stats, pow_stats_plot_mask, 0,r, ignorechans=base_ignorechans, threshold=5)
        if (rfimask.pow_stats[-1,:] == 1).all():
            tmp_mask = (base_mask[-1,:] | iqrm_med_pow_mask[-1,:])
            mask_pow_2diqrm_chan[-1,:] = tmp_mask
            mask_pow_2diqrm_int[-1,:] = tmp_mask
            del tmp_mask

        fig_iqrm_pow_mask, ax_iqrm_pow_mask = plt.subplots(2,1)
        plot_mask_comparison(mask_pow_2diqrm_chan, base_mask, title="iqrm_2D_pow (bad chans) - iqrm_1D_pow_on_medians", ax=ax_iqrm_pow_mask[0])
        plot_mask_comparison(mask_pow_2diqrm_int, base_mask, title="iqrm_2D_pow (bad ints) - iqrm_1D_pow_on_medians", ax=ax_iqrm_pow_mask[1])
        output_plot(fig_iqrm_pow_mask, pdf=p)

        mask_pow_iqrm_combo = mask_pow_2diqrm_chan|mask_pow_2diqrm_int|iqrm_med_pow_mask

        pow_stats_plot_mask = (mask_pow_iqrm_combo|base_mask)
        if (rfimask.pow_stats[-1,:] == 1).all():
            pow_stats_plot_mask[-1,:] = True
        fig_iqrm_pow, ax_iqrm_pow, cbar_iqrm_pow = plot_map_plus_sums(rfimask.pow_stats, pow_stats_plot_mask, returnplt=True)
        fig_iqrm_pow.suptitle("pow_stats post-1D-and-2D-iqrm")
        output_plot(fig_iqrm_pow, pdf=p)
        logging.info(f"masks {masked_frac(mask_pow_iqrm_combo|base_mask)}")

        mask_pow_iqrm_combo_exstats = reshape_rfifind_mask(M.shape, mask_pow_iqrm_combo, extra_stats_gulp, rfimask.ptsperint)
        masks[6] = mask_pow_iqrm_combo
        masks_exstats[6] = mask_pow_iqrm_combo_exstats

        working_mask, working_mask_exstats = check_mask_and_continue(
            working_mask, working_mask_exstats, 
            masks[6], masks_exstats[6], 
            args.problem_frac, rfimask, means, var, p, 
            stage=6,
        )
        
        del pow_med_chans
        del pow_med_ints
        del pow_chans
        del pow_ints
        del iqrm_med_pow_mask
        del mask_pow_2diqrm_chan
        del mask_pow_2diqrm_int


    if 7 in opts:
        logging.info("7: Precalculated gsk before, updating its mask")
        gsk_d_estimate_masked.mask[base_mask_exstats] = True
        # No idea why but iqrm_mask started failing on this when I switched to calculating gsk before the initial masking
        # If fill array and remask it works so that's the hacky thing I'm doing
        gsk_d_estimate_masked = np.ma.array(gsk_d_estimate_masked.filled(gsk_d_estimate_masked.mean()), mask=gsk_d_estimate_masked.mask)

        logging.info("Running 2D iqrm on -gsk and +gsk, chan-wise only")
        mask_gsk_2diqrm_chan1 = run_iqrm_2D(gsk_d_estimate_masked.filled(np.nan), gsk_d_estimate_masked.mask, 1,r, ignorechans=base_ignorechans, threshold=5)
        mask_gsk_2diqrm_chan2 = run_iqrm_2D(-gsk_d_estimate_masked.filled(np.nan), gsk_d_estimate_masked.mask, 1,r, ignorechans=base_ignorechans, threshold=5)
        mask_gsk_2diqrm_chan = mask_gsk_2diqrm_chan1|mask_gsk_2diqrm_chan2

        fig_iqrm_gsk_mask, ax_iqrm_gsk_mask = plt.subplots()
        plot_mask_comparison(mask_gsk_2diqrm_chan, base_mask_exstats, title="iqrm_2D_gsk (bad chans) - base", ax=ax_iqrm_gsk_mask)
        output_plot(fig_iqrm_gsk_mask, pdf=p)


        fig_iqrm_gsk, ax_iqrm_gsk, cbar_iqrm_gsk = plot_map_plus_sums(gsk_d_estimate_masked.data, mask_gsk_2diqrm_chan|base_mask_exstats, returnplt=True)
        fig_iqrm_gsk.suptitle("GSK post-2D-iqrm")
        output_plot(fig_iqrm_gsk, pdf=p)
        logging.info(f"masks {masked_frac(mask_gsk_2diqrm_chan|base_mask_exstats)}")

        masks_exstats[7] = mask_gsk_2diqrm_chan
        masks[7] = reshape_extra_stats_mask(rfimask.mask.shape, mask_gsk_2diqrm_chan, extra_stats_gulp, rfimask.ptsperint)

        working_mask, working_mask_exstats = check_mask_and_continue(
            working_mask, working_mask_exstats, 
            masks[7], masks_exstats[7], 
            args.problem_frac, rfimask, means, var, p, 
            stage=7,
        )

        del mask_gsk_2diqrm_chan1
        del mask_gsk_2diqrm_chan2

    if 8 in opts:
        logging.info("8: high fraction cut (with a limit of 0.1% on additional data masked)")

        hf = cut_off_high_fraction(working_mask_exstats, cumul_threshold_chans=1, cumul_threshold_ints=1, plot_diagnostics=False, verbose=False)
        diff = masked_frac(working_mask_exstats|hf) - masked_frac(working_mask_exstats)
        if diff > 0.001:
            threshes = np.linspace(0.2,1,101)[::-1]
            for i, thresh in enumerate(threshes):
                hf = cut_off_high_fraction(working_mask_exstats, hard_thresh_chans=thresh, cumul_threshold_chans=1, cumul_threshold_ints=1, plot_diagnostics=False, verbose=False)
                diff = masked_frac(working_mask_exstats|hf) - masked_frac(working_mask_exstats)
                if diff > 0.001:
                    logging.info(f"Over 0.1% threshold at thresh {thresh}")
                    if i == 0:
                        logging.warning(f"All extra data masked is coming from the int 0.3 fraction cut. This is a bit out of the ordinary, check the mask.")
                        thresh = 1
                    else:
                        thresh = threshes[i-1]
                    break
        else:
            logging.info(f"Under 0.1% limit at thresh 0.2")
            thresh = 0.2
        logging.info(f"Using thresh {thresh}")
        fighf, axhf = plt.subplots(1,2)
        hf = cut_off_high_fraction(working_mask_exstats, hard_thresh_chans=thresh, cumul_threshold_chans=1, cumul_threshold_ints=1, axhard=axhf)
        fighf.suptitle(f"Final high fraction cut using a hard theshold of {thresh} for the chans")
        output_plot(fighf, pdf=p)

        masks_exstats[8] = hf
        masks[8] = reshape_extra_stats_mask(rfimask.mask.shape, hf, extra_stats_gulp, rfimask.ptsperint)

        working_mask = working_mask | masks[8]
        working_mask_exstats = working_mask_exstats | masks_exstats[8]

    # Write out some info and make resulting mask

    bs_exstats = masks_exstats[0]
    bs_frac = masked_frac(bs_exstats)

    logging.info(f"Making base mask")
    logging.info(f"0: {bs_frac}")
    for opt in [1,2, 3]:
        if opt in opts:
            bs_exstats = bs_exstats|masks_exstats[opt]
            logging.info(f"+{opt}: {masked_frac(bs_exstats)} ({masked_frac(bs_exstats)-bs_frac})")
            bs_frac = masked_frac(bs_exstats)
    logging.info("")
    logging.info("Masks on top of base:")
    wm_exstats = copy.deepcopy(bs_exstats)
    for opt in [4,5,6]:
        if opt in opts:
            logging.info(f"base+{opt}: {masked_frac(bs_exstats|masks_exstats[opt])}")
            wm_exstats = wm_exstats | masks_exstats[opt]
    logging.info(f"all up to this point: {masked_frac(wm_exstats)}")
    if 7 in opts:
        logging.info(f"+7: {masked_frac(wm_exstats|masks_exstats[7])}")
        wm_exstats = wm_exstats | masks_exstats[7]
    if args.include_rfifind:
        logging.info(f"+rfifind: {masked_frac(wm_exstats|masks_exstats['rfifind'])}")

    logging.info("")
    logging.info(f"base+rfifind: {masked_frac(bs_exstats|masks_exstats['rfifind'])}")


    # ### Wrapping up

    logging.info("Wrapping up")
    if args.include_rfifind:
        logging.info("Adding in rfifind mask to the final mask")
        working_mask = working_mask | masks['rfifind']
        working_mask_exstats = working_mask_exstats | masks_exstats['rfifind']
    wrap_up(working_mask, working_mask_exstats, rfimask, means, var, p, outfilename, infstats_too=(not args.overwrite))

    sys.exit(0)
