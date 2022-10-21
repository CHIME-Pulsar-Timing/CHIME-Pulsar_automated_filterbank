#!/usr/bin/env python3
import numpy as np
import csv
import positive_bursts_timing as pt
import sys
from presto.filterbank import FilterbankFile
import glob

fn = sys.argv[1]
basename = sys.argv[2]
if len(sys.argv)>3:
    nom_DM = float(sys.argv[3])
else:
    nom_DM = -1
burst_dict = pt.get_burst_dict(fn)
if len(burst_dict) == 0:
    print('No Bursts!')
    sys.exit()

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
                fb_file = glob.glob('%s_*%s_pow_fdp.fil'%(basename,key))
                # print(basename,key)
                if len(fb_file)>1:
                    print('error globbing')
                    sys.exit(1)
                else:
                    fb_file=fb_file[0]
                fb_basename = fb_file.rstrip('.fil')



                fb_folders = '%s/'%(fb_basename)
                SPEG_file = glob.glob('%s/*_SPEG_all.csv'%(fb_folders))[0]
                mask_file = glob.glob('%s/*.mask'%(fb_folders))
                if len(mask_file)>1:
                    print('error globbing mask')
                    sys.exit(1)
                else:
                    # print(mask_file)
                    mask_file = mask_file[0]


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

                success=False
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
                            min_time = float(row[19])
                            max_time = float(row[20])

                            av_time = np.mean([min_time,max_time])
                            av_time = np.around(av_time,decimals=1)
                            burst[0] = np.around(burst[0],decimals=1)

                            SNR = row[14]
                            #import pdb; pdb.set_trace()
                            if (av_time==burst[0]) & (burst[3]==float(SNR)):
                                success=True
                                if nom_DM >0:
                                    DM = nom_DM
                                else:
                                    DM = row[12]
                                peak_downfact = int(float(row[16]))*2
                                break
                            success=False
                if success:
                    print('wrote successfully')
                    writer.writerow([fb_file,burst[0],60,DM,mask_file,mjd_fb_start,obs_length,pointing_ra_h,pointing_dec_h,SNR,''])
                else:
                    writer.writerow([fb_file,burst[0],60,DM,mask_file,mjd_fb_start,obs_length,pointing_ra_h,pointing_dec_h,SNR,''])

                    print("Never found burst in SPEG file")
                    print("day " + str(key))
                    print("failed on burst "+str(burst))
                    print("Wrote in extracted anything with previous downfact")

#multiday_times = pt.build_multiday_from_dict(burst_dict,min_day=1,min_time=0,sigma=3)
