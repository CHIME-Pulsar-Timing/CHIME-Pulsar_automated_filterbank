#!/usr/bin/env python3
import numpy as np
import csv
import positive_bursts_timing as pt
import sys



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
                fb_folder = '%s_%s_pow/0'%(basename,key)
                SPEG_file = '%s/0_SPEG.csv'%(fb_folder)
                success=False
                with open(SPEG_file,'r') as speg:
                    reader = csv.reader(speg,delimiter=',')
                    for i,row in enumerate(reader):
                        if i>0:
                            #first line is a header
                            #11 is the peak_downfact
                            #7 is the DM
                            #8 is peak_time
                            peak_time = row[8]
                            if float(peak_time)==float(burst[0]):
                                success=True
                                DM = row[7]
                                peak_downfact = row[11]
                                break
                            success=False
                if success:
                    writer.writerow([key,burst[0],DM,peak_downfact])
                else:
                    print("failed on burst"+str(burst))


#multiday_times = pt.build_multiday_from_dict(burst_dict,min_day=1,min_time=0,sigma=3)
