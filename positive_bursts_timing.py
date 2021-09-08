#!/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/python/3.8.2/bin/python

import argparse
import numpy as np
from presto.psr_utils import rrat_period_multiday, rrat_period

"""
Script for extracting rrat TOA's from positive_bursts.csv. The id string should 
look something like:

[dir to filterbank file]/J0740+17_cand_59236_pow/0//nsub_128/J0740+17_cand_59236_pow_153.1887616_sb_128_dm_44.85_snr_5.76

And it should be in the first column of the csv
"""

# Default filename, can be changed
FNAME = 'positive_bursts.csv'

def get_burst_dict(csvname):
    """From a csv of burst information named csv_name, extract burst info and 
    compile a dictionary of burst data. The first column of the csv should look
    something like:

    [dir to filterbank file]/J0740+17_cand_59236_pow/0//nsub_128/J0740+17[cont.]
    _cand_59236_pow_153.1887616_sb_128_dm_44.85_snr_5.76
    
    Returns: 
    Each key of the dictionary will be the string representation of the MJD the
    pulses were observed on. Each of these entries will be an array with the 
    rows storing relevant information for each pulse.

    [TOA in s from observation start, subbanding, dm, S/N]

    'mean' and 'std' will also be keys representing the mean and std of each
    pulse characteristic over all days
    """
    # Open the positive bursts file
    positive_bursts_csv = open(csvname)
    burst_lines = positive_bursts_csv.readlines()

    # Get all of the ids and put them into a dictionary
    burst_dict = {}
    for line in burst_lines:
        burst_str = line.split(',')[0]
        burst_numbers = [float(num) for num in burst_str.split('_')\
                         if num.replace('.', '1').isdigit()]
        burst_info = burst_numbers[-5:]

        if str(int(burst_info[0])) not in burst_dict:
            burst_dict[str(int(burst_info[0]))] = [burst_info[1:]]
        else: burst_dict[str(int(burst_info[0]))].append(burst_info[1:])

    # Parse the data to sort the entries and get the mean and std for each value
    full_array = np.zeros((len(burst_lines), 4))
    i = 0

    for key in burst_dict.keys():
        day_array = np.array(burst_dict[key])
        burst_dict[key] = day_array[np.argsort(day_array[:,0])]
        full_array[i:i+len(burst_dict[key])] += burst_dict[key]
        i += len(burst_dict[key])

    burst_dict['mean'] = np.mean(full_array, axis=0)
    burst_dict['std'] = np.std(full_array, axis=0)

    return burst_dict

def print_burst_dict(burst_dict):
    """Print info about pulses extracted from the dictionary
    """
    # Take the mean and STD to flag outlier DMs
    dm_mean = burst_dict['mean'][2]
    dm_std = burst_dict['std'][2]
    print(26*'='+'\n'+'=  Detected Pulse Info:  =\n'+26*'=')
    
    # Loop through every day, which are keys in burst_dict
    for key in sorted(burst_dict):
        if key.isnumeric():
            # Extract DM and TOAs
            print(f'Pulse info for MJD {key}:')
            day_array, dm_array = burst_dict[key][:,0], burst_dict[key][:,2]
            last_day = -10

            for day, dm in zip(day_array, dm_array):
                # Start with basic info
                print_str = 'Time: {:10.6f}  Pulse DM: {:5.2f}  '\
                            .format(day, dm)
                # Append warnings if the dm is irregular of the pulse is too
                # close to the previous one
                if dm_mean-dm_std > dm or dm_mean+dm_std < dm:
                    print_str += 'Irregular DM! '
                if day - last_day < 0.5:
                    print_str += 'Pulses very close!'
                last_day = day
                print(print_str)

            print(f'{len(day_array)} pulses detected'+'\n')


def trim_burst_dict(burst_dict, min_count=2, min_time=0, sigma=1):
    """
    From a disctionary returned by get_burst_dict, remove data points according
    to the given constraints.

    Arguments:
    burst_dict: Dictionary of burst_information returned by get_burst_dict
    min_day: minimum number of pulses in a day, default 2
    min_time: minimum time in seconds between pulses, useful for elmiminating
    double hits
    sigma: Pulses with dm more than sigma standard deviations away from the mean
    will not be counted

    Returns:
    burst_dict_trimmed: An edited version of burst_dict with data that conforms
    to the input parameters
    """
    # Start a new dictionary to be populated
    burst_dict_trimmed = {}
    dm_mean, dm_std = burst_dict['mean'][2], burst_dict['std'][2]
    
    for key in sorted(burst_dict.keys()):
        if key.isnumeric():
            data_array = burst_dict[key].copy()
            # Find rows with irregular DM and delete them
            bad_dm =  np.logical_or(data_array[:,2] < dm_mean - sigma * dm_std,
                                    data_array[:,2] > dm_mean + sigma * dm_std)
            data_array = np.delete(data_array, bad_dm, axis=0)
            # Find rows where the pulse is too close to the preceding pulse
            if len(data_array) >= 2:
                double_hits = data_array[1:,0] - data_array[:-1,0] < min_time
                double_hits = np.append([False], double_hits).flatten()
                data_array = np.delete(data_array, double_hits, axis=0)
            
            # Only add to the dict if the remaining data is enough
            if len(data_array) >= max(2, min_count):
                burst_dict_trimmed[key] = data_array
    
    return burst_dict_trimmed


def build_multiday_from_dict(burst_dict, min_count=2, min_time=0, sigma=1):
    """From a dictionary returned by get_burst_dict, construct a nested list
    multiday_times to be passed into presto's rrat_period_multiday.

    Arguments:
    burst_dict: Dictionary of burst_information returned by get_burst_dict
    min_day: minimum number of pulses in a day, default 2
    min_time: minimum time in seconds between pulses, useful for elmiminating
    double hits
    sigma: Pulses with dm more than sigma standard deviations away from the mean
    will not be counted

    Returns:
    multiday_times: A nested list with each entry containing the TOA of each
    pulse relative to the start of each observation day
    """
    # Call trim_burst_dict to extract the bad data points
    burst_dict_trimmed = trim_burst_dict(burst_dict, min_count, min_time, sigma)
    # Extact the TOAs to multiday_times
    multiday_times = []
    for key in sorted(burst_dict_trimmed.keys()):
        if key.isnumeric():
            multiday_times.append(list(burst_dict_trimmed[key][:,0]))
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
    parser.add_argument('--keep_best', type=int, default=None,
                        help='Optional: Only process the best n days of'\
                        +'observation. Incompatible with --individual')
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
    parser.add_argument('--min_time', type=float, default=0,
                        help='Minimum time between pulses, useful for removing'\
                        +' double hits')
    parser.add_argument('-v', action='store_true',
                        help='Set Flag to output timing info')
    parser.add_argument('--individual', action='store_true',
                        help='Process each valid day individually with '\
                        +'rrat_period instead of rrat_period_multiday.'\
                        +'Incompatible with --keep_best')

    args = parser.parse_args()
    if args.filename is not None:
        FNAME = args.filename

    # Check and make sure the file is present
    try:
        # Extract burst info
        burst_dict = get_burst_dict(FNAME)
        if len(burst_dict) == 0:
            print('No Bursts!')
            exit()
        dm_mean, dm_std = burst_dict['mean'][2], burst_dict['std'][2]
    except FileNotFoundError:
        print(f'No {FNAME} found')
        exit()
    
    # Print DM info
    print(f'DM = {round(dm_mean, 2)}({round(dm_std, 2)})'+'\n')
    
    # Verbose output of pulse info. Very useful for census work
    if args.v:
        print_burst_dict(burst_dict)

    # Either run everything through multiday or process each day individually
    if not args.individual:
        multiday_times = build_multiday_from_dict(burst_dict,
                                                  min_count=args.min_count,
                                                  min_time=args.min_time,
                                                  sigma=args.sigma)
        if args.keep_best is not None:
            multiday_times = take_good_days(multiday_times, n=args.keep_best)

        if args.v:
            print('TOA on each day in seconds relative to the first pulse:')
            [print(day) for day in multiday_times]
            print()

        if len(multiday_times) < 1:
            print("Insufficient Data")

        rrat_period_multiday(multiday_times)

    else:
        # I could make them both work, but it would be tricky
        if args.individual and args.keep_best:
            print('--individual and --keep_best options are incompatible!')
            exit()
        
        # Trim burst_dict
        burst_dict_trimmed = trim_burst_dict(burst_dict, args.min_count, 
                                             args.min_time, args.sigma)
        # Different way of getting the days
        keys = [key for key in sorted(burst_dict_trimmed.keys())\
                if key.isnumeric]
        
        # Track periods for processing
        period_array = np.zeros(len(keys))
        # Go through each day
        for i in range(len(keys)):
            times = burst_dict_trimmed[keys[i]][:,0]
            print(f'Period Computation for MJD {keys[i]}')
            period_array[i] = rrat_period(times, numperiods=100)
            print()
        # Print final statistics
        print('Computed periods: \n'+str(period_array))
        print(f'Mean:   {np.mean(period_array)}')
        print(f'Median: {np.median(period_array)}')
        print(f'STD:    {np.std(period_array)}')

# Newline for readability            
print()

