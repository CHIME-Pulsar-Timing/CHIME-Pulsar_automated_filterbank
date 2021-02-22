from multiprocessing import Pool
import sys
from gwg_cand_search_pipeline import run_rfifind,run_prepsubband,run_sp
from gwg_cand_search_pipeline import run_sk_mad
from presto.filterbank import FilterbankFile
import pipeline_config

def run_prep_subband(fname,prepsub,dm):
    #fname = fil.rstrip('.fil')
    #get some header details from the filterbank file
    print(fname)
    filfile = FilterbankFile(fname+'.fil')
    tsamp = filfile.dt
    nsamp = filfile.nspec
    nchan = filfile.nchan

    sk_mad=False
    coherent=False
    dmlist = [i for i in pipeline_config.dm_set if (i < dm+20)]
    if prepsub:
        #start mp
        pool = Pool(len(dmlist))
        #don't need coherent dm value.
        pool.starmap(run_prepsubband,((fname,tsamp,nsamp,dm,123,sk_mad,coherent) for dm in dmlist))
        pool.close()
    return dmlist 


rfi=True
prepsub=True
sp=True

fil_files = sys.argv[2:]
fil_files = list(fil.rstrip('.fil') for fil in fil_files)
if rfi:
    pool=Pool(len(fil_files))
    pool.map(run_rfifind,fil_files)
    pool.close()
else:
    #if not rfifind run skmad
    for i,fil in enumerate(fil_files):
        fil_files[i]=run_sk_mad(fil,fil+'.fil')

dm_lists=[]
dm=float(sys.argv[1])
for fil in fil_files:
    dm_lists.append(run_prep_subband(fil,prepsub,dm))
step=1
if sp:
    for fname in fil_files:
        run_sp(fname)


