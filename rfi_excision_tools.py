import numpy as np

def get_divisors(n):
    d = []
    i = 1
    while i <= n:
        if n % i == 0:
            d.append(i)

        i += 1

    return np.array(d)

def mad(data, axis=None):
    median = np.ma.median(data, axis=axis)
    return np.ma.median(np.abs(data - median), axis=axis)

def mad2std(m):
    return 1.4826 * m

def detect_rfi_mad(spectra, stat_window=64, tfact=1, thresh=5.0):
    """ 
    This function detects RFI in frequency by computing the median and median absolute deviation
    of channels from adjacent channels. This helps to recognize dropouts due to packet loss or 
    GPU node failures, while also picking up on some RFI that spectral kurtosis may miss.
    Adapted from Bradley Meyers code developed for slow pulsar search

    This function is based on functions used to estimate calibration information for CHIME/FRB, 
    written by Bridget Anderson.

    Parameters
    ----------
    spectra : np.array of floats
            The dynamic spectrum with a shape (nchan, M) where M is the number 
            of time steps to be combined.
    stat_window : int
            The size (in chunks) of the sliding window in frequency over which 
            statistics are calculated.
    thresh : float
            The number of standard deviations the median of a channel
            must differ from the median to be flagged.
    
    Returns
    -------
    final_mask : np.array of bool
            The output RFI mask. It is the same shape as the input dynamic spectrum.
            Indices that are flagged are denoted by a value of True.
    """
    old_np_settings = np.seterr(divide='ignore', invalid='ignore')

    final_mask = np.zeros_like(spectra, dtype=bool)

    nchan, ntimes = spectra.shape
    step = 1024 // tfact
    # do this for every time step (i.e. every ~1 ms)
    for i in range(ntimes // step):
        # predetermine the time-index slice
        tsl = slice(i * step, (i + 1) * step)
        intensity = np.ma.mean(spectra[:, tsl], axis=1)

        for j in range(stat_window, nchan + stat_window, stat_window):
            # predetermine the freqeuncy-index slice based on stat_window
            upper_index = j if j < nchan else nchan
            lower_index = j - stat_window if (j - stat_window) > 0 else 0
            fsl = slice(lower_index, upper_index)

            spec_chunk = intensity[fsl]

            # calculate the median and the median absolute deviation
            med = np.ma.median(spec_chunk)
            mad = np.ma.median(np.abs(spec_chunk - med))
            # convert he MAD to an equivalent standard deviation
            std = mad2std(mad)

            # filter out samples outside the nominal threshold defined by thresh
            pthresh = med + thresh * std
            nthresh = med - thresh * std
            filt = np.logical_or(spec_chunk > pthresh, spec_chunk < nthresh)

            # update the entries of the final data mask
            final_mask[fsl, tsl] = np.tile(filt, (step, 1)).T

    # also do the median filtering based on the average spectrum
    mspec = np.ma.mean(spectra, axis=-1)
    mspec_med = np.ma.median(mspec)
    mspec_mad = np.ma.median(np.abs(mspec - mspec_med))
    mspec_std = mad2std(mspec_mad)
    pthresh = mspec_med + 3 * mspec_std
    nthresh = mspec_med - 1.5 * mspec_std
    mspec_mask = np.logical_or(mspec > pthresh, mspec < nthresh)

    # combine the mean spectrum mask with the individual time step mask
    final_mask = np.logical_or(final_mask, mspec_mask[:, np.newaxis])

    np.seterr(**old_np_settings)
    return final_mask

def detect_rfi_sk(spectra, thresh=3.0, ffact=1, plot=False):
    """
    Calculate the generalized spectral kurtosis from a dynamic spectrum in order to flag RFI.
    This function operates under the assumption that the `spectra` (which can be a masked array)
    passed has shape (nchan, M), which means that the organization of data needs to be done prior
    to calling this function.
    Adapted from Bradley Meyers code developed for slow pulsar search

    Parameters
    ----------
    spectra : np.array of floats
            The dynamic spectrum with a shape (nchan, M) where M is the number 
            of time steps to be combined.
    thresh : float
            The number of standard deviations from the mean (E[SK] = 1) above 
            and below which channels will be masked. The lower threshold is 
            multiplied by 0.75.
    
    Returns
    -------
    combined_mask : np.array of bool
            The output RFI mask. It is the same shape as the input dynamic spectrum.
            Indices that are flagged are denoted by a value of True.

    """
    old_np_settings = np.seterr(divide='ignore', invalid='ignore')
    nchan, m = spectra.shape

    # perform a summation over the m power measurements and their square
    s1_sq = (np.ma.sum(spectra, axis=1).astype(float))**2
    s2 = np.ma.sum(spectra**2, axis=1)

    # generalized spectral kurtosis estimator (Nita & Gary 2010, eq. 7)
    spec_sk = ((m + 1) / float(m - 1)) * (m * (s2 / s1_sq) - 1)

    # calcualte the mean of the SK estimator in clean parts of the spectrum
    spec_sk_norm = np.ma.mean([
        np.ma.median(spec_sk[3100//ffact:3800//ffact]),
        np.ma.median(spec_sk[6300//ffact:7000//ffact])])

    # calcualte the appropriate d so that the mean value of SK is 1
    #d = (np.ma.mean(s1_sq) * ((m - 1) * spec_sk_norm + 1) - m * np.ma.mean(s2)) / (m * (m * np.ma.mean(s2) - np.ma.mean(s1_sq))))
    d = 1.0 / spec_sk_norm # approximately
    #print("recalculated d = {0}".format(d))
    # recalculate the SK with the d correction. Solve for X in:
    #
    #   X * ((m + 1) / (m - 1)) * (m * (s2 / s1_sq) - 1) = ((d * m + 1) / (m - 1)) * (m * (s2 / s1_sq) - 1)
    #   X = (d * m + 1) / (m + 1)
    spec_sk = ((d * m + 1) / float(m + 1)) * spec_sk

    # the expectation values for a correctly scaled SK measurement (Nita & Gary 2010, eq. 9 & 10)
    expectation_mean = 1.0
    expectation_std = np.sqrt((2.0 / m) * (1.0 + 1.0 / d))

    # based on the desired threshold (equivalent to Gaussian sigma) for RFI removal
    nthresh = expectation_mean - thresh * expectation_std
    pthresh = expectation_mean + thresh * expectation_std

    # create RFI mask from SK estimator
    sk_mask = np.ma.logical_or(spec_sk > pthresh, spec_sk < nthresh)
    #sk_mask = np.tile(sk_mask.astype(bool), (m, 1))        

    # combine the SK mask with the original spectra mask (from the message packets themselves)
    combined_mask = np.logical_or(spectra.mask, sk_mask[:, np.newaxis])

    np.seterr(**old_np_settings)
    return combined_mask

