import csv
import sys
import os
import numpy as np
from presto.filterbank import FilterbankFile
from presto import infodata
'''
gotta write the .csv and also the inf file
The inf file takes the following format
object_obs,RA,Dec,central_freq_low_chan,total_bandwidth
a,b,c,d,e
'''
def prep_speg(inffile):
    folder=os.getcwd()
    all_files= np.array(os.listdir(folder))
    sp_files =all_files[list(['.singlepulse' in file for file in all_files])]
    #find the postion of the 'DM' string
    dtype =[('fn','U200'),('dm',float)]
    unsorted=[]
    for file in sp_files:
        #parse the fn
        index=file.find('DM')
        DM = file[index+2:]
        DM = DM.split('.')
        DM = float(DM[0]+'.'+DM[1])
        unsorted.append((file,DM))
    
    #writes the csv file
    sorted = np.sort(np.array(unsorted,dtype=dtype),order='dm')
    last_folder=folder.split('/')
    last_folder=last_folder[-1]
    with open(last_folder+'singlepulses.csv','w',newline='') as spegscsv:
        writer=csv.writer(spegscsv,delimiter=',')
        writer.writerow(['DM','SNR','time','sample','downfact'])
        for file,DM in sorted:
            with open(file) as csvfile:
                reader=csv.reader(csvfile,delimiter=' ')
                for i,row in enumerate(reader):
                    row = list(filter(None,row))
                    if i>0:
                        writer.writerow(row)

    #write the inf file
    with open(last_folder+'_inf.txt','w',newline='') as speginf:
        writer=csv.writer(speginf,delimiter=',')
        writer.writerow(['object_obs','RA','Dec','central_freq_low_chan','total_bandwidth'])
        info = infodata.infodata(inffile)
        obj_n = info.basenm.split('_')[0] 
        ra=info.RA
        dec=info.DEC
        fch = info.lofreq
        bw=info.BW
        writer.writerow([obj_n,ra,dec,fch,bw])
        
#prep_speg(sys.argv[1])
