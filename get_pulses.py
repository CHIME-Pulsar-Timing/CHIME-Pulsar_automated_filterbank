#!/usr/bin/env python3
import numpy as np
import csv
import positive_bursts_timing as pt
import sys
fn = sys.argv[1]
burst_dict = pt.get_burst_dict(fn)
if len(burst_dict) == 0:
    print('No Bursts!')
    exit()
#breakpoint()
with open('extracted_bursts.csv','w') as csv_file:
    writer = csv.writer(csv_file,delimiter = ',')
    for key in burst_dict.keys():
        if (key=='mean')| (key=='std'):
            continue
        if isinstance(burst_dict[key][0],np.ndarray):
            for burst in burst_dict[key]:
                writer.writerow([key,burst[0]])
                print(burst[0])
        else:
            writer.writerow([key,burst_dict[key][0]])

#multiday_times = pt.build_multiday_from_dict(burst_dict,min_day=1,min_time=0,sigma=3)
