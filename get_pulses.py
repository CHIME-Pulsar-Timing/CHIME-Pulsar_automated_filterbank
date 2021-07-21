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
                fb_folder = '%s_%s_pow/0'%(basename,key)
                SPEG_file = '%s/0_SPEG_all.csv'%(fb_folder)
                mask_file = '%s/%s_%s_pow_rfifind.mask'%(fb_folder,basename,key)
                success=False

                #grab information needed for Bradley
                fb = FilterbankFile(fb_file)
                mjd_pulse = fb.tstart+float(burst[0])/60/60/24
                obs_length = fb.nspec*fb_file.dt

                with open(SPEG_file,'r') as speg:
                    reader = csv.reader(speg,delimiter=',')
                    for i,row in enumerate(reader):
                        if i>0:
                            #first line is a header
                            #16 is the peak_downfact
                            #12 is the peak_DM
                            #13 is peak_time
                            #14 is peak_SNR
                            peak_time = row[13]
                            SNR = row[14]
                            #import pdb; pdb.set_trace()
                            if (int(float(peak_time))==int(burst[0])) & (burst[3]==float(SNR)):
                                success=True
                                DM = row[12]
                                peak_downfact = int(float(row[16]))*2
                                break
                            success=False
                if success:
                    writer.writerow([fb_file,burst[0],60,DM,mask_file,mjd_pulse,obs_length])
                else:
                    print("day " + str(key))
                    print("failed on burst "+str(burst))


#multiday_times = pt.build_multiday_from_dict(burst_dict,min_day=1,min_time=0,sigma=3)
