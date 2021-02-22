import sys
from presto import waterfaller
from presto.filterbank import FilterbankFile
import csv

#lets first parse the single pulse file
fil_file = sys.argv[1]
sp_file = sys.argv[3:]
fil = fil_file.strip('.fil')
snr_default = float(sys.argv[2])

#this is code to select time ranges to make waterfalls
search_file = fil+'_search.csv'
import os.path as path
search =path.exists(search_file)
if search:
    with open(search_file) as file:
        line = file.read()
        l=line.split(',')
else:
    l=[0,999999999]

for file in sp_file:
    times=[]
    with open(file) as csvfile:
        r=csv.reader(csvfile,delimiter=' ')
        for i,row in enumerate(r):
            if i>0:
                new_row = []
                for item in row:
                    if item:
                        new_row.append(item)
                if float(new_row[1])>snr_default:
                    times.append(float(new_row[2])-0.25)
                    dm=float(new_row[0])

    for time in times:
        for i in range(0,len(l),2):
            if (time>float(l[i])) & (time<float(l[i+1])):
                maskfn=fil+'_rfifind.mask'
                raw = FilterbankFile(fil_file)
                data,nbinsextra,nbins,start=waterfaller.waterfall(raw,time,0.5,dm=dm,nsub=256,downsamp=16,mask=True,maskfn=maskfn)
                waterfaller.plot_waterfall(data,start,0.5,integrate_ts=True,integrate_spec=True,cmap_str='binary',save_path=fil)
        
