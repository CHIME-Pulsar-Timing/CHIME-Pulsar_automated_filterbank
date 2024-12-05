#!/usr/bin/env python
import numpy as np
from presto_without_presto import infodata
import copy
import os


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

def mask_params_from_array_as_sets(maskarr):
    """
    Input: 2D mask array of shape (nint, nchan), where the channels follow rfifind convention (index 0 corresponds to the lowest frequency channel)
    (Masking convention: 1/True means mask, 0/False means don't mask)
    Output: mask_zap_chans, mask_zap_ints, mask_zap_chans_per_int which can be used to write a presto .mask file
    """
    nint, nchan = maskarr.shape

    chans_to_mask = set(np.arange(nchan)[maskarr.sum(axis=0) == nint].astype(np.int32))
    ints_to_mask = set(np.arange(nint)[maskarr.sum(axis=1) == nchan].astype(np.int32))

    chans_per_int = []
    for i in range(nint):
        chans_per_int.append(set(np.where(maskarr[i, :] == 1)[0].astype(np.int32)))

    return chans_to_mask, ints_to_mask, chans_per_int


def write_new_mask_from(filename, maskarr, rfifind_obj, include_old=False, infstats_too=False):
    """
    Write a mask numpy array as a rfifind .mask file, using attributes from rfifind_obj
    filename: filename to write mask to (will add a .mask extension if not present)
    maskarr: 2D mask array of shape (nint, nchan)
             where the channels follow rfifind convention (index 0 corresponds to the lowest frequency channel)
             and 1/True in the array means mask, 0/False means don't mask
    rfifind_obj: a rfifind object, e.g. a mask file read in using presto's rfifind.rfifind
    if include_old=True then take the logical or of maskarr and rfifind_obj.mask
    infstats_too: bool, False by default. Whether to also write out an inf file and a stats file
        NB this will not pass along the Notes of the inf file since those never get read in
    """
    if include_old:
        maskarr = (maskarr | rfifind_obj.mask)
    header = dict(
        timesig=rfifind_obj.time_sig,
        freqsig=rfifind_obj.freq_sig,
        MJD=rfifind_obj.MJD,
        dtint=rfifind_obj.dtint,
        lofreq=rfifind_obj.lofreq,
        df=rfifind_obj.df,
        nchan=rfifind_obj.nchan,
        nint=rfifind_obj.nint,
        ptsperint=rfifind_obj.ptsperint,
    )
    if filename[-5:] != ".mask":
        filename += ".mask"
    write_mask_file(filename, maskarr, header)
    if infstats_too:
        inf_filename = filename[:-5] + ".inf"
        stats_filename = filename[:-5] + ".stats"
        infdata = copy.deepcopy(rfifind_obj.idata)
        infdata.basenm = os.path.basename(filename)[:-5]
        infdata.to_file(inf_filename)
        stats_header = dict(
            nchan=rfifind_obj.nchan,
            nint=rfifind_obj.nint,
            ptsperint=rfifind_obj.ptsperint,
            lobin=rfifind_obj.lobin,
            numbetween=rfifind_obj.numbetween,
        )
        write_stats_file(
            stats_filename,
            rfifind_obj.pow_stats,
            rfifind_obj.avg_stats,
            rfifind_obj.std_stats,
            stats_header
        )


def write_stats_file(filename, pow_stats, avg_stats, std_stats, header):
    """
    Write pow_atats, avg_stats and std_stats arrays as a rfifind .stats file
    filename: filename to write stats to (will add a .stats extension if not present)
    The arrays should all be of shape (nint, nchan)
    header: dictionary which must contain keys:
        'nchan' (int) - number of frequency channels in the mask
        'nint' (int) - number of time intervals in the mask
        'ptsperint' (int) - number of time samples per interval
        'lobin' (int) - no idea what this is but it gets read in from the stats file in rfifind.read_stats
        'numbetween' (int) - no idea what this is but it gets read in from the stats file in rfifind.read_stats
    """
    if filename[-6:] != ".stats":
        filename += ".stats"

    header_params = [
        np.array(header["nchan"], dtype=np.int32),
        np.array(header["nint"], dtype=np.int32),
        np.array(header["ptsperint"], dtype=np.int32),
        np.array(header["lobin"], dtype=np.int32),
        np.array(header["numbetween"], dtype=np.int32)
    ]

    # check shapes match what's in header
    shape_should_be = (header["nint"], header["nchan"])
    if pow_stats.shape != shape_should_be:
        raise AttributeError(f"pow_stats of shape {pow_stats.shape}, from header expect {shape_should_be}")
    if avg_stats.shape != shape_should_be:
        raise AttributeError(f"avg_stats of shape {avg_stats.shape}, from header expect {shape_should_be}")
    if std_stats.shape != shape_should_be:
        raise AttributeError(f"std_stats of shape {std_stats.shape}, from header expect {shape_should_be}")

    with open(filename, "wb") as fout:
        for variable in header_params:
            variable.tofile(fout)

        pow_stats.astype(np.float32).tofile(fout)
        avg_stats.astype(np.float32).tofile(fout)
        std_stats.astype(np.float32).tofile(fout)

    fout.close()

def write_mask_file(filename, maskarr, header):
    """Write a mask numpy array as a rfifind .mask file
    filename: filename to write mask to (will add a .mask extension if not present)
    maskarr: 2D mask array of shape (nint, nchan)
             where the channels follow rfifind convention (index 0 corresponds to the lowest frequency channel)
             and 1/True in the array means mask, 0/False means don't mask
    header: dictionary which must contain keys:
        'timesig' (float) - from rfifind options, PRESTO default is 10
        'freqsig' (float) - from rfifind options, PRESTO default is 4
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
        np.array(header["timesig"], dtype=np.float64),
        np.array(header["freqsig"], dtype=np.float64),
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
