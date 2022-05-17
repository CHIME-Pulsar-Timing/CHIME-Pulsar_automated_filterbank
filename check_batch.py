#!/usr/bin/env python3
import sys
import csv
import numpy as np
#first argument is speg file
SPEG_file = sys.argv[1]
#argv are all the fil files
filfiles = sys.argv[3:]
rank = float(sys.argv[2])
success = True
#create filfiles list
times = []
for fil in filfiles:
    try:
        splits = np.array(fil.split('_'))
        sb_ind = np.where(splits=='sb')
        time_ind = int(sb_ind[0][0]-1)
        t = np.around(float(splits[time_ind]),2)
        times.append(t)
    except:
        #for some reason can't read fil file
        print(fil)

with open(SPEG_file,'r') as speg:
    reader = csv.reader(speg,delimiter=',')
    for i,row in enumerate(reader):
        if i>0:
            #first line is a header
            #4 is group rank
            #16 is the peak_downfact
            #12 is the peak_DM
            #13 is peak_time
            #14 is peak_SNR
            #19 is min_time
            #20 is max_time

            group_rank = float(row[4])
            #we shouldn't be searching for peak time
            min_time = float(row[19])
            max_time = float(row[20])
            av_time = np.mean([min_time,max_time])
            av_time = np.around(av_time,2)
            if (group_rank<=rank) & (group_rank>0):
                if av_time in times:
                    pass
                else:
                    success = False
                    print(av_time)
                    print(group_rank)
if not success:
    sys.exit(1)
