#!/usr/bin/env python3
import numpy as np
from scipy import stats
from your.candidate import Candidate
import your
from presto import rfifind
from multiprocessing import Pool
def maskfile(maskfn, data, start_bin, nbinsextra):
    print('loading mask')
    rfimask = rfifind.rfifind(maskfn)
    print('getting mask')
    mask = get_mask(rfimask, start_bin, nbinsextra)[::-1]
    print('get mask finished')
    masked_chans = mask.all(axis=1)
    #mask the data but set to the mean of the channel
    mask_vals = np.median(data,axis=1)
    data[masked_chans,:] = mask_vals[masked_chans,np.newaxis]
    print(f"Total_masked chan length {sum(masked_chans)}")
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

def calculate_merge_sk(X):
    chunk_sz = X['chunk_sz']
    samp = X['samp']
    fil = X['fil']
    rfifind_mask = X['rfifind_mask']
    your_object = your.Your(fil)
    i = X['i']
    your_data = your_object.get_data(samp, chunk_sz)
    print(f"starting sk filter on data of shape {your_data.shape}")
    sk_mask = your.utils.rfi.sk_filter(
        your_data,
        your_object.your_header.foff,
        your_object.your_header.nchans,
        your_object.your_header.tsamp,
        sigma=5,
    )
    print("finished sk filter")
    rfimask = rfifind.rfifind(rfifind_mask)
    mask_arr = rfimask.mask_zap_chans_per_int
    int_mask = mask_arr[i]
    sk_chan = np.where(sk_mask)
    for skc in sk_chan:
        if not (skc in int_mask):
            int_mask = np.append(int_mask,skc)
    int_mask = np.sort(int_mask)
    return int_mask

def merge_mask(fil,rfifind_mask,presto_block = 8):
    #chunk size in seconds
    your_object = your.Your(fil)
    #presto uses 2400 as a standard block size
    chunk_sz = presto_block*2400
    total_ints = np.ceil(your_object.your_header.nspectra/chunk_sz)
    print(f"total ints {total_ints}")
    samp = 0
    i = 0
    new_mask_arr = []
    pool_arr = []
    while samp<your_object.your_header.nspectra:
        pool_arr.append({"chunk_sz":chunk_sz,"samp":samp,"fil":fil,"rfifind_mask":rfifind_mask,"i":i})
        samp += chunk_sz
        i+=1
    p = Pool(20)
    new_mask_arr = p.map(calculate_merge_sk,pool_arr)
    rfimask = rfifind.rfifind(rfifind_mask)
    maskarr = np.full((rfimask.nint,rfimask.nchan),False)
    for i,m in enumerate(new_mask_arr):
        maskarr[i,m] = True
    write_mask_file(rfifind_mask.strip(".mask")+"_SK", np.array(maskarr), rfimask.__dict__)


def mask_params_from_array(maskarr):
    """
    Input: 2D mask array of shape (nint, nchan), where the channels follow rfifind convention (index 0 corresponds to the lowest frequency channel)
    (Masking convention: 1/True means mask, 0/False means don't mask)
    Output: mask_zap_chans, mask_zap_ints, mask_zap_chans_per_int which can be used to write a presto .mask file
    """
    nint, nchan = maskarr.shape

    chans_to_mask = np.arange(nchan)[maskarr.sum(axis=0) == nint].astype(np.int32)
    ints_to_mask = np.arange(nint)[maskarr.sum(axis=1) == nchan].astype(np.int32)

    chans_per_int = []
    for i in range(nint):
        chans_per_int.append(np.where(maskarr[i, :] == 1)[0].astype(np.int32))

    return chans_to_mask, ints_to_mask, chans_per_int


def write_mask_file(filename, maskarr, header):
    """Write a mask numpy array as a rfifind .mask file
    filename: filename to write mask to (will add a .mask extension if not present)
    maskarr: 2D mask array of shape (nint, nchan)
             where the channels follow rfifind convention (index 0 corresponds to the lowest frequency channel)
             and 1/True in the array means mask, 0/False means don't mask
    header: dictionary which must contain keys:
        'time_sig' (float) - from rfifind options, PRESTO default is 10
        'freq_sig' (float) - from rfifind options, PRESTO default is 4
        'MJD' (float) - starting MJD of the observation
        'dtint' (float) - length of one time interval in seconds
        'lofreq' (float) - center frequency of lowest channel
        'df' (float) - width of one frequency channel
        'nchan' (int) - number of frequency channels in the mask
        'nint' (int) - number of time intervals in the mask
        'ptsperint' (int) - number of time samples per interval
    """
    # Massage inputs
    if filename[-5:] != ".mask":
        filename += ".mask"

    header_params = [
        np.array(header["time_sig"], dtype=np.float64),
        np.array(header["freq_sig"], dtype=np.float64),
        np.array(header["MJD"], dtype=np.float64),
        np.array(header["dtint"], dtype=np.float64),
        np.array(header["lofreq"], dtype=np.float64),
        np.array(header["df"], dtype=np.float64),
        np.array(header["nchan"], dtype=np.int32),
        np.array(header["nint"], dtype=np.int32),
        np.array(header["ptsperint"], dtype=np.int32),
    ]
    # Check maskarr shape matches nint and nchan in header
    if header["nint"] != maskarr.shape[0]:
        raise ValueError(
            f"nint in header ({header['nint']}) does not match maskarr shape {maskarr.shape}"
        )
    if header["nchan"] != maskarr.shape[1]:
        raise ValueError(
            f"nchan in header ({header['nchan']}) does not match maskarr shape {maskarr.shape}"
        )

    # Write to file
    with open(filename, "wb") as fout:
        for variable in header_params:
            variable.tofile(fout)

        zap_chans, zap_ints, zap_chans_per_int = mask_params_from_array(maskarr)

        nzap_chan = np.asarray(zap_chans.size, dtype=np.int32)
        nzap_chan.tofile(fout)
        if nzap_chan:
            zap_chans.tofile(fout)

        nzap_int = np.asarray(zap_ints.size, dtype=np.int32)
        nzap_int.tofile(fout)
        if nzap_int:
            zap_ints.tofile(fout)

        nzap_per_int = []
        for an_arr in zap_chans_per_int:
            nzap_per_int.append(an_arr.size)
        if len(nzap_per_int) != header["nint"]:
            raise AttributeError("BUG: nzap_per_int should be of length nint!")
        nzpi = np.asarray(nzap_per_int, dtype=np.int32)
        nzpi.tofile(fout)

        # rfifind.py only calls fromfile if nzap != 0 and nzap != nchan
        for i in range(header["nint"]):
            if nzap_per_int[i]:
                if nzap_per_int[i] != header["nchan"]:
                    zap_chans_per_int[i].tofile(fout)

if __name__ == '__main__':
    import sys

    merge_mask(sys.argv[1],sys.argv[2])
