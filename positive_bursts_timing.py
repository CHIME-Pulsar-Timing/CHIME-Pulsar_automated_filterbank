import argparse
import numpy as np
from presto.psr_utils import rrat_period_multiday

"""
Script for extracting rrat TOA's from positive_bursts.csv. The id string should 
look something like:

[dir to filterbank files]/J2325-0530_59266_pow/0//nsub_128/cand_tstart_59266.889074767183_tcand_0.5816484_dm_14.95000_snr_6.67000

And it should be in the first column of the csv
"""

# Default filename, can be changed
FNAME = 'positive_bursts.csv'


def get_burst_ids(csvname):
    """From a csv file named csvname extract a list of burst ids. Generally this
    function expects the file positive_bursts.csv and first column values that
    look like: 

    [dir to filterbank files]/J2325-0530_59266_pow/0//nsub_128/cand_tstart_
    59266.889074767183_tcand_0.5816484_dm_14.95000_snr_6.67000
    
    Returns: 
    burst_ids: A list of strings identifying each burst. For the example above
    this function would extract:
    cand_tstart_59266.889074767183_tcand_0.5816484_m_14.95000_snr_6.67000

    burst_ids will be returned in sorted order.
    """
    # Open the positive bursts file
    positive_bursts_csv = open(csvname)
    burst_lines = positive_bursts_csv.readlines()

    # Get all of the ids and sort them for TOA extraction
    burst_ids = []
    for line in burst_lines:
        burst_str = line.split(',')[0]
        id_index = burst_str.find('cand_tstart_')
        burst_ids.append(burst_str[id_index:])
    burst_ids.sort()
    
    return burst_ids


def extract_pulse_info(burst_ids, sigma=1):
    """From a list of burst ids generated from a positive_bursts.csv with
    get_burst_ids, extract the day and DM or each burst. Pulses whose DM's are 
    more than sigma standard deviations away from the mean DM are not included

    Returns:
    day_array: An array of the MJD TOA of each pulse
    dm_array: An array of the DM of each pulse
    """
    day_list = []
    dm_list = []
    
    for i in range(len(burst_ids)):
        day_list.append(burst_ids[i][12:burst_ids[i].find('_tcand_')])
        dm_list.append(burst_ids[i][burst_ids[i].find('_dm_')+4:\
                       burst_ids[i].find('_snr_')])
    
    ay_array = np.array(day_list, dtype=np.float128)
    dm_array = np.array(dm_list, dtype=np.float128)

    # Eliminate everything that has a DM to far away
    dm_mean = np.mean(dm_array)
    dm_std = np.std(dm_array)

    bad_indeces = np.where(np.logical_or(dm_array < dm_mean - sigma*dm_std,
                                         dm_array > dm_mean + sigma*dm_std))
    day_array = np.delete(day_array, bad_indeces)
    dm_array = np.delete(dm_array, bad_indeces)
    return day_array, dm_array


def build_multiday(day_array, min_day=1):
    """From an iterable day_array of floating point MJD TOA's for each pulse, 
    construct a nested list multiday_times to be passed into presto's
    rrat_period_multiday

    Returns:
    multiday_times: A nested list with each entry containing the TOA of each
    pulse relative to the start of each observation day
    """
    # Build the list that rrat_period multiday wants
    multiday_times = []
    multiday_times.append([86400*(day_array[0] % 1)])

    # iterate through the timestamps
    for i in range(1, len(day_array)):

        # add a new entry if we need to
        if int(day_array[i-1]) != int(day_array[i]):
            multiday_times.append([])
    
        multiday_times[-1].append(86400*(day_array[i] % 1))
    
    # Set Each day to 0 seconds for the first pulse
    for day in multiday_times:
        first_time = day[0]
        for i in range(len(day)):
            day[i] -= first_time

    # remove days with too few pulses
    multiday_times = [dy for dy in multiday_times if len(dy) >= min(2, min_day)]

    return multiday_times


def take_good_days(multiday_times, n=2):
    """From a nested TOA list multiday_times, take only the n best days and
    return a nested list of the same form

    Returns:
    multiday_times_cut: a nested list with the same structure as multiday_times
    containing only the n best observing days
    """
    multiday_times_cut = []
    multiday_count = np.array([len(day) for day in multiday_times])
    sort_array = np.argsort(multiday_count)[::-1]
    for i in range(0, min(len(sort_array), n)):
        multiday_times_cut.append(multiday_times[sort_array[i]])
    
    return multiday_times_cut


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--keep_best', type=int, 
                        help='Optional: Keep the best n days of observation',
                        default=None)
    parser.add_argument('--filename', type=str,
                        help="Optional: Process a file that is not "\
                        +"positive_bursts.csv", default=None)
    parser.add_argument('--sigma', type=float,
                        help='Optional: Exclude pulses with DM sigma standard '\
                        +'deviations away from the mean DM, default=2', 
                        default=2)
    parser.add_argument('--min_count', type=int,
                        help='Optional: Only include days with at least this '\
                        +'many pulses', default=2)
    parser.add_argument('-v', action='store_true',
                        help='Set Flag to output timing info')

    args = parser.parse_args()
    if args.filename is not None:
        FNAME = args.filename
    
    try:
        burst_ids = get_burst_ids(FNAME)
    except FileNotFoundError:
        print(f'No {FNAME} found')
        exit()

    if len(burst_ids) == 0:
        print('No Bursts!')
        exit()

    day_array, dm_array = extract_pulse_info(burst_ids, sigma=args.sigma)

    print(f'DM = {round(np.mean(dm_array), 2)}({round(np.std(dm_array), 2)})',
            '\n')
    if args.v:
        print('Pulse MJD Times:\n', [day for day in day_array], '\n')

    multiday_times = build_multiday(day_array, min_day=args.min_count)
    
    if args.keep_best is not None:
        multiday_times = take_good_days(multiday_times, n=args.keep_best)

    if args.v:
        print('TOA on each day in seconds relative to the first pulse:')
        [print(day) for day in multiday_times]
        print()

    if len(multiday_times) < 1:
        print("Insufficient Data")
    else:
        rrat_period_multiday(multiday_times)
    print()

