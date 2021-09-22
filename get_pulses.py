#!/usr/bin/env python3
import numpy as np
import csv
import positive_bursts_timing as pt
import sys
from presto.filterbank import FilterbankFile
import glob

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
                fb_file = glob.glob('%s_*%s_pow.fil'%(basename,key))
                if len(fb_file)>1:
                    print('error globbing')
                    sys.exit(1)
                else:
                    fb_file=fb_file[0]
                fb_basename = fb_file.rstrip('.fil')



                fb_folders = '%s/0'%(fb_basename)
                print(fb_folders)
                SPEG_file = '%s/0_SPEG_all.csv'%(fb_folders)
                mask_file = glob.glob('%s/*.mask'%(fb_folders))

                if len(mask_file)>1:
                    print('error globbing mask')
                    sys.exit(1)
                else:
                    print(mask_file)
                    mask_file = mask_file[0]
                success=False


                #grab information needed for Bradley
                fb = FilterbankFile(fb_file)
                mjd_fb_start = fb.tstart#+float(burst[0])/60/60/24
                obs_length = fb.nspec*fb.dt
                pointing_ra = fb.src_raj
                pointing_dec = fb.src_dej
                rah = int(pointing_ra/10000)
                ram = int(np.mod(pointing_ra,10000)/100)
                ras = np.mod(np.mod(pointing_ra,10000),100)
                pointing_ra_h = '%d:%d:%f'%(rah,ram,ras)

                dech = int(pointing_dec/10000)
                decm = int(np.mod(pointing_dec,10000)/100)
                decs = np.mod(np.mod(pointing_dec,10000),100)
                pointing_dec_h = '%d:%d:%f'%(dech,decm,decs)

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
                    writer.writerow([fb_file,burst[0],60,DM,mask_file,mjd_pulse,obs_length,pointing_ra_h,pointing_dec_h])
                else:
                    print("day " + str(key))
                    print("failed on burst "+str(burst))


#multiday_times = pt.build_multiday_from_dict(burst_dict,min_day=1,min_time=0,sigma=3)
