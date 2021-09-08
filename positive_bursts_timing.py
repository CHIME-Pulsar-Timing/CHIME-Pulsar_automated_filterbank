import argparse
import numpy as np
from presto.psr_utils import rrat_period_multiday, rrat_period

"""
Script for extracting rrat TOA's from positive_bursts.csv. The id string should 
look something like:

[dir to filterbank file]/J0740+17_cand_59236_pow/0//nsub_128/J0740+17_cand_59236_pow_153.1887616_sb_128_dm_44.85_snr_5.76

And it should be in the first column of the csv

Using the --mjd flag, the script can also process old id strind which look
something like this:

[dir to filterbank files]/J2325-0530_59266_pow/0//nsub_128/cand_tstart_59266.889074767183_tcand_0.5816484_dm_14.95000_snr_6.67000
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
    
    # Make sure everything is an array
    master_array = np.zeros((len(burst_lines), 4))
    i = 0
    
    for key in burst_dict.keys():
        day_array = np.array(burst_dict[key])
        burst_dict[key] = day_array[np.argsort(day_array[:,0])]
        master_array[i:i+len(burst_dict[key])] += burst_dict[key]
        i += len(burst_dict[key])

    burst_dict['mean'] = np.mean(master_array, axis=0)
    burst_dict['std'] = np.std(master_array, axis=0)
   
    return burst_dict
    

def extract_pulse_info(burst_ids):
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
    
    day_array = np.array(day_list, dtype=np.float128)
    dm_array = np.array(dm_list, dtype=np.float128)

    return day_array, dm_array


def print_burst_dict(burst_dict):
    """Print info about pulses extracted from the dictionary
    """
    dm_mean = burst_dict['mean'][2]
    dm_std = burst_dict['std'][2]
    print(26*'='+'\n'+'=  Detected Pulse Info:  =\n'+26*'=')

    for key in sorted(burst_dict):
        if key.isnumeric():
            print(f'Pulse info for MJD {key}:')
            day_array, dm_array = burst_dict[key][:,0], burst_dict[key][:,2]
            last_day = -10

            for day, dm in zip(day_array, dm_array):
                print_str = 'Time: {:10.6f}  Pulse DM: {:5.2f}  '\
                            .format(day, dm)
                if dm_mean-dm_std > dm or dm_mean+dm_std < dm:
                    print_str += 'Irregular DM! '
                if day - last_day < 0.5:
                    print_str += 'Pulses very close!'
                last_day = day
                print(print_str)
            
            print(f'{len(day_array)} pulses detected'+'\n')


def pulse_print(day_array, dm_array):
    """Print info about pulses
    """

    print("Detected Pulse Info:")
    MJD = int(day_array[0])
    print(f'Pulse info for MJD {MJD}:')
    
    last_day = -100
    dm_mean, dm_std = np.mean(dm_array), np.std(dm_array)
    count = 0
    for day, dm in zip(day_array, dm_array):
        if int(day) != MJD:
            MJD = int(day)
            print(f'{count} pulses detected'+2*'\n'+\
                  f'Pulse info for MJD {MJD}:')
            count = 0
        print_str = 'Time: {:10.6f}  Pulse DM: {:5.2f}  '.format(day, dm)
        if dm_mean-dm_std > dm or dm_mean+dm_std < dm:
            print_str += 'Irregular DM! '
        if (day - last_day)*86400 < 0.5:
            print_str += 'Pulses very close!'
        last_day = day
        print(print_str)
        count += 1
    print(f'{count} pulses detected')
    print('\n')


def build_multiday_array(day_array, dm_array, min_day=2, min_time=0, sigma=1):
    """From an iterable day_array of floating point MJD TOA's for each pulse, 
    construct a nested list multiday_times to be passed into presto's
    rrat_period_multiday

    Arguments:
    day_array: iterable day_array of floating point MJD TOA's for each pulse, 
    construct a nested list multiday_times
    dm_array: iterable of dispersion measures corresponding to the days in
    day_array
    min_day: minimum number of pulses in a day, default 2
    min_time: minimum time in seconds between pulses, useful for elmiminating
    double hits
    sigma: Pulses with dm more than sigma standard deviations away from the mean
    will not be counted

    Returns:
    multiday_times: A nested list with each entry containing the TOA of each
    pulse relative to the start of each observation day
    """
    # Eliminate everything that has a DM to far away
    dm_mean = np.mean(dm_array)
    dm_std = np.std(dm_array)
    bad_indeces = np.where(np.logical_or(dm_array < dm_mean - sigma*dm_std,
                                         dm_array > dm_mean + sigma*dm_std))
    day_array = np.delete(day_array, bad_indeces)
    dm_array = np.delete(dm_array, bad_indeces)

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
    # Remove pulses that are too close together
    for day in multiday_times:
        first_time = day[0]
        pop_indeces = []
        for i in range(len(day)):
            day[i] -= first_time
            if i >= 1 and day[i] - day[i-1] < min_time:
                pop_indeces.append(i)
        [day.pop(i) for i in pop_indeces[::-1]]

    # remove days with too few pulses
    multiday_times = [dy for dy in multiday_times if len(dy) >= max(2, min_day)]

    return multiday_times


def build_multiday_from_dict(burst_dict, min_day=2, min_time=0, sigma=1):
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
    multiday_times = []
    for key in sorted(burst_dict):
        if key.isnumeric():
            pulse_times = burst_dict[key][:,0]
            dm_array = burst_dict[key][:,2]
            good_times = np.where(np.logical_and(dm_array > dm_mean \
                                                 -sigma*dm_std,
                                                 dm_array < dm_mean \
                                                 +sigma*dm_std))
            multiday_times.append(list(pulse_times[good_times]))
    
    # Remove pulses that are too close together
    for day in multiday_times:
        pop_indeces = []
        for i in range(len(day)):
            if i >= 1 and day[i] - day[i-1] < min_time:
                pop_indeces.append(i)
        [day.pop(i) for i in pop_indeces[::-1]]

    # remove days with too few pulses
    multiday_times = [dy for dy in multiday_times if len(dy) >= max(2, min_day)]

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
                        help='Optional: Keep the best n days of observation')
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
                        +'rrat_period instead of rrat_period_multiday')
    parser.add_argument('--mjd', action='store_true',
                        help='Run script with the old MJD id strings instead')

    args = parser.parse_args()
    if args.filename is not None:
        FNAME = args.filename
    
    try:
        if args.mjd:
            burst_ids = get_burst_ids(FNAME)
            if len(burst_ids) == 0:
                print('No Bursts!')
                exit()
            day_array, dm_array = extract_pulse_info(burst_ids)
            dm_mean, dm_std = np.mean(dm_array), np.std(dm_array)
        else:
            burst_dict = get_burst_dict(FNAME)
            if len(burst_dict) == 0:
                print('No Bursts!')
                exit()
            dm_mean, dm_std = burst_dict['mean'][2], burst_dict['std'][2]
    except FileNotFoundError:
        print(f'No {FNAME} found')
        exit()

    print(f'DM = {round(dm_mean, 2)}({round(dm_std, 2)})'+'\n')
    if args.v:
        if args.mjd:
            pulse_print(day_array, dm_array)
        else:
            print_burst_dict(burst_dict)
    
    if args.mjd:
        multiday_times = build_multiday_array(day_array, dm_array,
                                              min_day=args.min_count,
                                              min_time=args.min_time, 
                                              sigma=args.sigma)
    else:
        multiday_times = build_multiday_from_dict(burst_dict,
                                                  min_day=args.min_count,
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
    else:
        if args.individual:
            periods = np.zeros(len(multiday_times))
            for i in range(len(multiday_times)):
                periods[i] = rrat_period(multiday_times[i], numperiods=1000)
                print()
            print('Computed periods: \n'+str(periods))
            print(f'Mean:   {np.mean(periods)}')
            print(f'Median: {np.median(periods)}')
            print(f'STD:    {np.std(periods)}')
        else:                
            rrat_period_multiday(multiday_times)
    print()
