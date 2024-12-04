from scipy.optimize import minimize
from scipy.stats import median_abs_deviation as mad
import numpy as np
from matplotlib import pyplot as plt
from presto.filterbank import FilterbankFile
from presto import filterbank as fb
from presto_without_presto import rfifind
from rfifind_numpy_tools import write_new_mask_from
import sys
import numpy as np
import logging
import matplotlib.pyplot as plt

from sigpyproc import readers as r
def get_filterbank_data_window(fn, s, samps):
    """
    Open the filterbank file, extract the spectra around
    the middle (+/- duration / 2) and return that data as a
    2D array.
    """

    print("getting filterbank data")

    filf = r.FilReader(fn)
    hdr = filf.header
    tsamp = hdr.tsamp
    fil_dur = hdr.nsamples * tsamp
    # start in the middle of the data
    start_bin = s
    nsamp = samps
    # get the data
    _ = filf.read_block(start_bin, nsamp)
    # read the block
    return _, hdr


def spectral_kurtosis(data, N=1, d=None):
    """
    Compute spectral kurtosis. See [Nita et al. (2016)](https://doi.org/10.1109/RFINT.2016.7833535) for details.

    Args:
        data (numpy.ndarray): 2D frequency time data
        N (int): Number of accumulations on the FPGA
        d (float): shape factor

    Returns:
         numpy.ndarray: Spectral Kurtosis along frequency axis

    """
    zero_mask = data == 0
    data = np.ma.array(data.astype(float), mask=zero_mask)
    S1 = data.sum(0)
    S2 = (data**2).sum(0)
    M = data.shape[0]
    if d is None:
        d = (np.nanmean(data.ravel()) / np.nanstd(data)) ** 2
    return ((M * d * N) + 1) * ((M * S2 / (S1**2)) - 1) / (M - 1)

def fit_model(x, a, A, freq, phase):
    # this will be a polynomial model with a sinusoidal with phase offset
    polynomial = np.polyval([a], x)
    sinusoidal = A * np.sin(2 * np.pi * freq * x + phase)
    return polynomial + sinusoidal

def loglikelihood(params, x, y):
    exponential_component = (y-fit_model(x, *params))**2/(2*0.04**2)
    return np.sum(exponential_component)-(0.5*len(y)*np.log(2*np.pi*0.04**2)) #std of 0.04 hardcoded via trial and error


def remove_rfi(gsk,threshold,fit = True, plot_orig = True,plot_name=1):
    mask = np.zeros(len(gsk),dtype=bool)
    #mask all gsk above 5
    mask[np.abs(gsk) > 10] = True
    gsk[mask] = np.median(gsk[~mask])
    gsk = gsk.data
    gsk -= np.median(gsk)
    gsk_mad = mad(gsk[~mask])
    x = np.arange(len(gsk))
    # plt.figure()
    # if plot_orig:
        # plt.plot(x,gsk,label="pre std cut")
    mask[np.abs(gsk) > threshold * gsk_mad] = True
    if sum(mask)>700:
        mask[:]=True
        print("removing all data")
        return x,gsk,mask
    gsk[mask] = np.median(gsk[~mask])
    subtracted_gsk = np.zeros_like(gsk)

    # plt.plot(x, gsk, label="post MAD cut")
    # plt.xlabel("Channel")
    # plt.ylabel("Spectral Kurtosis")
    # if not fit:
    #     plt.legend()
    #     plt.savefig(f"SK_secondcut_{plot_name}.png")
    #     print("saved plot")
    x_fit = np.zeros_like(x)
    y_fit = np.zeros_like(gsk)
    # lets see if we can fit the data with fit_model
    # do a chunked fit every 128 channels

    if fit:
        # gsk = gsk / np.max(gsk)
        for i in range(16):
            ind_start = i * 64
            ind_end = (i + 1) * 64
            # print(ind_start, ind_end)
            res = minimize(loglikelihood,[0,0.1,0.015,0],args=(x[ind_start:ind_end],gsk[ind_start:ind_end]),method="Nelder-Mead")
            x_fit[ind_start:ind_end] = x[ind_start:ind_end]
            y_fit[ind_start:ind_end] = fit_model(x[ind_start:ind_end], *res.x)
            subtracted_gsk[ind_start:ind_end] = gsk[ind_start:ind_end] - fit_model(x[ind_start:ind_end], *res.x)
        # plt.plot(x_fit, y_fit, label="fit",color="red")
        # plt.legend()
        # plt.savefig(f"SK_precut_{plot_name}.png")
        # plt.figure()
        # plt.plot(x, subtracted_gsk, label="subtracted",color="black")
        # plt.plot(x, gsk, label="original",color="blue",alpha=0.2)
        # plt.legend()
        # plt.savefig(f"SK_final_{plot_name}.png")
    else:
        subtracted_gsk = gsk

    # plt.close("all")
    return x,subtracted_gsk, mask


def perform_SK(fn,nints):
    fil = FilterbankFile(fn,mode='read')
    fn_clean = fn.strip('.fil')
    nspecs = fil.nspec
    loop_iters = nints
    chunk_size = nspecs//loop_iters
    # loop_iters = 5
    print(f"chunk size is {chunk_size}")
    rfimask = np.zeros((nints,int(fil.nchan)),dtype=bool)
    for i in range(loop_iters):
        s = i*chunk_size
        print(f"Processing chunk {i+1}/{loop_iters}")
        if i < loop_iters-1:
            e = (i+1)*chunk_size
        else:
            #if we're on the last chunk then take everything
            e = nspecs
        get_size = e-s
        _, hdr = get_filterbank_data_window(fn, s, get_size)
        gsk = spectral_kurtosis(_.T,N=512,d=0.5)
        if chunk_size < 45:
            x,gsk,mask_first_pass = remove_rfi(gsk,threshold=2,fit=False,plot_orig=True,plot_name=f"{i}_{chunk_size}")
            # x,gsk,mask_second_pass = remove_rfi(gsk,threshold=2,fit=False,plot_orig=False,plot_name=f"{i}_{chunk_size}")
            mask = mask_first_pass
        else:
            x,gsk,mask_first_pass = remove_rfi(gsk,threshold=3,fit=True,plot_orig=False,plot_name=f"{i}_{chunk_size}")
            x,gsk,mask_second_pass = remove_rfi(gsk,threshold=3,fit=False,plot_orig=False,plot_name=f"{i}_{chunk_size}")
            mask = mask_first_pass | mask_second_pass
        #for some reason mask is reversed
        rfimask[i,:] = mask[::-1]
    fil.close()
    return np.array(rfimask)

def load_rfi_mask(mask_fn,plot_name="initial_mask.png"):
    rfimask = rfifind.rfifind(mask_fn)
    nints = rfimask.nint
    plt.figure()
    plt.imshow(rfimask.mask,aspect="auto")
    plt.colorbar()
    plt.savefig(plot_name)

    return rfimask,nints

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compute the spectral kurtosis of a data file")
    parser.add_argument("fn", type=str, help="File name of the data")
    args = parser.parse_args()
    #parse the basename of fn
    basename = args.fn.split(".")[0]
    #construct the mask filename
    rfimask_fn = f"{basename}_rfifind.mask"
    rfimask,nints = load_rfi_mask(rfimask_fn)
    import os
    #copy the original mask to a new file
    os.system(f"cp {rfimask_fn} {basename}_original_rfifind.mask")
    #plot the rfimask
    plt.figure()
    plt.imshow(rfimask.mask,aspect="auto")
    plt.colorbar()
    plt.savefig("initial_mask.png")
    mask = perform_SK(args.fn,nints)
    write_new_mask_from(basename+"_gsk_rfifind.mask", mask, rfimask, include_old=True, infstats_too=True)
