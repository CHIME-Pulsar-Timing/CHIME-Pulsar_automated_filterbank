#!/usr/bin/env python3
from presto.filterbank import FilterbankFile
from presto import filterbank as fb
import sys
import numpy as np
import logging
import matplotlib.pyplot as plt
from presto import filterbank as fb
log = logging.getLogger(__name__)

fn = sys.argv[1]
fil = FilterbankFile(fn,mode='read')
nspecs = fil.nspec
#number of samples for 1 sec just like sps
#channel wide cleaning
chunk_size = 30674
thresh_sig = 4.5
#loop through chunks
nspecs = fil.nspec
loop_iters = int(nspecs/chunk_size)
fn_clean = fn.strip('.fil')
fdp_fn = f"{fn_clean}_fdp.fil"
fb.create_filterbank_file(fdp_fn,fil.header,nbits=fil.header['nbits'])
new_fil = FilterbankFile(fdp_fn,mode='append')

for i in range(loop_iters):
    print(i,'/',loop_iters)
    # gsv_arr = []
    s = i*chunk_size
    if i < loop_iters-1:
        e = (i+1)*chunk_size
    else:
        #if we're on the last chunk then take everything
        e = nspecs
    get_size = e-s
    spectra = fil.get_spectra(s,get_size)
    spec = spectra.data
    med = np.median(spec,axis=1)
    mean = np.mean(spec,axis=1)
    std = np.std(spec,axis=1)
    #set the thresold
    threshold = med - thresh_sig*std
    tile_count = e-s
    # for m,st in zip(mean,std):
    #     gsv = np.random.normal(loc=m,scale=st,size=tile_count)
    #     gsv_arr.append(gsv)
    # gsv_arr = np.array(gsv_arr)

    thresh_mat = np.tile(threshold,(tile_count,1)).T
    med_mat = np.tile(med,(tile_count,1)).T
    #generate an array of gaussian random values
    #find values below threshold and replace with the median
    mask = spec < thresh_mat
    spec[mask] = med_mat[mask]
    # spec[mask] = gsv_arr[mask]
    new_fil.append_spectra(spec.T)

fil.close()
new_fil.close()
# debugging plots
# nspecs = new_fil.nspec
# spectra = new_fil.get_spectra(0,nspecs)
# spectra = spectra.data
# plt.figure()
# plt.hist(spectra[701,:],bins=100)
# plt.figure()
# plt.plot(spectra[701,:])
# plt.xlabel('sample number')
# plt.ylabel('sample value')
# plt.show()
