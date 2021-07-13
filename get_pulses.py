#!/usr/bin/env python3
import numpy as np
import csv
import positive_bursts_timing as pt
import sys
from presto.filterbank import FilterbankFile



fn = sys.argv[1]
basename = sys.argv[2]
burst_dict = pt.get_burst_dict(fn)
if len(burst_dict) == 0:
    print('No Bursts!')
    sys.exit()
#breakpoint()
with open('extracted_bursts.csv','w') as csv_file:
    writer = csv.writer(csv_file,delimiter = ',')
    for key in burst_dict.keys():
        if (key=='mean')| (key=='std'):
            continue
        if isinstance(burst_dict[key][0],np.ndarray):
            for burst in burst_dict[key]:
                # writer.writerow([key,burst[0]])
                # the key is the day, the burst[0] is the timestamp,
                # first gotta find the file and load up filterbank
                fb_file = '%s_%s_pow.fil'%(basename,key)
                fb_file = FilterbankFile(fb_file)
                mjd_pulse = fb_file.tstart+float(burst[0])/60/60/24
                obs_length = fb_file.nspec*fb_file.dt
                writer.writerow([mjd_pulse,obs_length])
        else:
            # writer.writerow([key,burst_dict[key][0]])
            fb_file = '%s_%s_pow.fil'%(basename,key)
            fb_file = FilterbankFile(fb_file)
            mjd_pulse = fb_file.tstart+float(burst[key][0])/60/60/24
            obs_length = fb_file.nspec*fb_file.dt
            writer.writerow([mjd_pulse,obs_length])



#multiday_times = pt.build_multiday_from_dict(burst_dict,min_day=1,min_time=0,sigma=3)
