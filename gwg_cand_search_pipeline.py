import numpy as np
import subprocess
import os, glob
import argparse
import re

from presto.filterbank import FilterbankFile
import pipeline_config
from sk_mad_rficlean import sk_mad_rfi_excision
import sys
#original GWG pipeline written by Chiamin
def run_rfifind(fname):

    rfifind_command = 'rfifind -blocks %d -zapchan %s -o %s %s.fil' %(pipeline_config.rfiblocks,pipeline_config.zaplist,fname,fname)
    run_rfifind_cmd = subprocess.Popen([rfifind_command], shell=True)
    run_rfifind_cmd.wait()

def run_sk_mad(fname,fil):

    sk_mad_rfi_excision(fname,fil)
    fnamenew = str(fname)+'_sk_mad'

    return fnamenew

def run_prepsubband(fname,tsamp,nsamp,dm,coherent_dm,coherent=True):

    if coherent:
        dms, ds, sb = pipeline_config.coherent_ddplan(tsamp, dm, coherent_dm)
    else:
        dms, ds, sb = pipeline_config.ddplan(tsamp, dm)

    numout = (((nsamp / (2400*256))+1)*(2400*256))/ds
    if numout % 7 == 0:
        numout += (numout/7)

    prepsubband_command = 'prepsubband -lodm %.2f -dmstep %.2f -numdms 100 -downsamp %d -nsub %d -mask %s_rfifind.mask -o %s %s.fil' %(dm,dms,numout,ds,sb,fname,fname,fname)
    run_prepsubband_cmd = subprocess.Popen([prepsubband_command],shell=True)
    run_prepsubband_cmd.wait()

def run_realfft(fname,fil,rednoise=True,zaplist=None):

    datfiles = sorted(glob.glob(str(fname)+'*.dat'))
    for dat in datfiles:

        realfft_command = 'realfft -disk '+str(dat)
        run_realfft_cmd = subprocess.Popen([realfft_command],shell=True)
        run_realfft_cmd.wait()

    fftfiles = sorted(glob.glob(str(fname)+'*.fft'))

    if zaplist:

        for fft in fftfiles:

            baryv = subprocess.check_output(['prepdata -start 0.99 %s -o tmp |grep Average' %fil],shell=True).split()[-1].rstrip('\n')
            zaplist_command = 'zapbirds -zap -zapfile %s -baryv %f'%(zaplist,baryv) 

    if rednoise:

        for fft in fftfiles:

            rednoise_command = 'rednoise '+str(fft)
            run_rednoise_cmd = subprocess.Popen([rednoise_command],shell=True)
            run_rednoise_cmd.wait()
            os.rename(str(fft).rstrip('.fft')+'_red.fft',str(fft))

def run_accelsearch(fname,zmax,wmax,binary=True):

    fftfiles= sorted(glob.glob(str(fname)+'*.fft'))

    for fft in fftfiles:

        if binary == True:

            accelsearch_command = 'accelsearch -zmax %d -wmax %d %s' %(zmax,wmax,fft)
            run_accelsearch_cmd = subprocess.Popen([accelsearch_command],shell=True)
            run_accelsearch_cmd.wait()

        else:

            accelsearch_command = 'accelsearch -zmax 0 -numharm 32 %s' %(fft)
            run_accelsearch_cmd = subprocess.Popen([accelsearch_command],shell=True)
            run_accelsearch_cmd.wait()

def run_ffa(fname):

    for dm in pipeline_config.ffa_dm_set:

        if dm == 0.0:

            datfiles = sorted(glob.glob(fname+'*DM*'+str(dm)+'0.dat'))

        else:

            datfiles.extend(sorted(glob.glob(fname+'*DM*'+str(dm)+'0.dat')))

    for dat in datfiles:

        ffa_command = '/psr_scratch/common_utils/ffaGo/ffa.py '+str(dat)
        run_ffa_cmd = subprocess.Popen([ffa_command],shell=True)
        run_ffa_cmd.wait()

def run_sp(fname):
    sp_command = 'single_pulse_search.py -f %s*.dat' %(fname)
    run_sp_cmd = subprocess.Popen([sp_command],shell=True)
    run_sp_cmd.wait()

def run_accelsift(fname):

    accelsift_command = 'python /usr/local/src/presto/python/ACCEL_sift.py'
    with open(fname+'_ACCEL_sift_cands.lis','w+') as outfile:
        run_accelsift_cmd = subprocess.Popen([accelsift_command],stdout=outfile,shell=True)
        run_accelsift_cmd.wait()

def run_ffa_sift(fname):

    if not os.path.isfile(fname+'_rfifind.inf'):

        rfifind_empty_command = 'rfifind -blocks '+str(pipeline_config.rfiblocks)+' -ignorechan 0:1024 -o '+str(fname)+' '+str(fname)+'.fil'
        run_rfifind_empty_cmd = subprocess.Popen([rfifind_empty_command],shell=True)
        run_rfifind_empty_cmd.wait()

    ffa_sift_command = 'python /psr_scratch/common_utils/ffaGo/ffa_final.py '+str(fname)+'_rfifind.inf'
    run_ffa_sift_cmd = subprocess.Popen([ffa_sift_command],shell=True)
    run_ffa_sift_cmd.wait()

def run_prepfold(fname,accelfile,accelcand,candDM,candperiod,sk_mad=False):

    #check period for folding option
    prepfoldcmd = pipeline_config.foldplan(fname,accelfile,accelcand,candDM,candperiod,sk_mad)
    run_prepfold_cmd = subprocess.Popen([prepfoldcmd],shell=True)
    run_prepfold_cmd.wait()

def fold_candidates(fname,source_dm,coherent=True):

    if fft:

        with open(fname+'_ACCEL_sift_cands.lis','r') as candfile:

            for line in candfile:

                if '_ACCEL_' in line:

                    accelfile = line.split(':')[0]+'.cand'
                    accelcand = line.split(':')[1].split()[0]
                    candDM = accelfile.split('DM')[1].split('_')[0]
                    candperiod = float(line.split()[7]) / 1000.0
                    if coherent and np.abs(float(candDM)-float(source_dm)) > 20.0:
                        continue
                    else:
                        run_prepfold(fname,accelfile,accelcand,candDM,candperiod,sk_mad)
    
    if ffa:
    
        with open(fname+'_cands.ffa','r') as candfile:
   
            for line in candfile:

                if '.dat' in line:
 
                    candDM = line.split('DM')[1].split('.dat')[0]
                    candperiod = line.split()[2]
                    accelfile = line.split(':')[0].rstrip('.dat')
                    accelcand = None
                    if coherent and np.abs(float(candDM)-float(source_dm)) > 20.0:
                        continue
                    else:
                        run_prepfold(fname,accelfile,accelcand,candDM,candperiod,sk_mad)

def run_spp_regular(fil):
    dm=7
    fil = args.fil
    fname = fil.rstrip('.fil')

    #get some header details from the filterbank file
    filfile = FilterbankFile(fil)
    tsamp = filfile.dt
    nsamp = filfile.nspec
    nchan = filfile.nchan

    sk_mad=False
    coherent=False

    run_rfifind(fname)
    dm_list = [dm,dm+20,dm+40]
    for dm in dmlist:
        #don't need coherent dm value.
        run_prepsubband(fname,tsamp,nsamp,dm,123,sk_mad,coherent)
    run_sp(fname)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--fil', type=str, help='Input filterbank file')
    parser.add_argument('--dedisp',action='store_true',help='Run prepsubband and dedisperse the data')
    parser.add_argument('--sk_mad',action='store_true',help='Run sk_mad RFI excision instead of rfifind')
    parser.add_argument('--dm', type=float,help='DM of the candidate. This will determine the max DM to search for pulsar')
    parser.add_argument('--coherent',action='store_true',help='Use the coherent ddplan for searching')
    parser.add_argument('--no_fft',action='store_false',help='Do not run fft search (Default is to run)')
    parser.add_argument('--no_rednoise',action='store_false',help='Do not run rednoise removal on Fourier series (Default is to run)')
    parser.add_argument('--zaplist',nargs='?',default=None,help='zaplist to remove from the fft files')
    parser.add_argument('--binary',action='store_true',help='Run binary search (This only works if fft is run)')
    parser.add_argument('--zmax',nargs='?',default=100,help='zmax value for binary search')
    parser.add_argument('--wmax',nargs='?',default=0,help='wmax value for binary search')
    parser.add_argument('--ffa',action='store_true',help='Run Fast Folding Algorithm')
    parser.add_argument('--sp',action='store_true',help='Run single pulse search')
    parser.add_argument('--fold',action='store_true',help='Fold the candidates')
    parser.add_argument('--speg',action='store_true',help='creates the SPEGID files')
    parser.add_argument('--fetch',action='store_true',help='creates the FETCH files')
    parser.add_argument('--rfifind',action='store_true',help='Runs rfifind using the configuration in pipeline config')
    args = parser.parse_args()

    fil = args.fil
    fname = fil.rstrip('.fil')

    if os.path.islink(fil):
        fil = os.readlink(fil)
        if not os.path.isfile(fil):

            print('File does not exist')
            sys.exit()

    #get some header details from the filterbank file
    filfile = FilterbankFile(fil)
    tsamp = filfile.dt
    nsamp = filfile.nspec
    nchan = filfile.nchan

    sk_mad = args.sk_mad
    source_dm = args.dm
    coherent = args.coherent
    fft = args.no_fft
    rednoise = args.no_rednoise
    zaplist = args.zaplist
    binary = args.binary
    zmax = args.zmax
    wmax = args.wmax
    ffa = args.ffa
    sp = args.sp
    fold = args.fold
    speg = args.speg
    fetch =args.fetch 
    dedisp = args.dedisp
    rfifind = args.rfifind
    print('Running RFI mitigation')
    if sk_mad:
        if os.path.exists('{}_sk_mad.fil'.format(fname)):
            print('sk mad cleaned data already exists')
            fname = '{}_sk_mad'.format(fname)
        else:
            fname = run_sk_mad(fname,fil)
    if rfifind:
        run_rfifind(fname)
    if dedisp:
        #dedispersion
        if coherent:
            dmlist = [source_dm-i for i in pipeline_config.coherent_dm_set if source_dm-i > 0]
            dmlist.append(0)
        else:
            dmlist = [i for i in pipeline_config.dm_set if i < source_dm+20]
        for dm in dmlist:
            run_prepsubband(fname,tsamp,nsamp,dm,source_dm,coherent)

    #run fft
    if fft:
        run_realfft(fname,rednoise,zaplist)
        run_accelsearch(fname,zmax,wmax,binary)
        run_accelsift(fname)

    #run ffa
    if ffa:
        run_ffa(fname)
        run_ffa_sift(fname)

    if sp:
        run_sp(fname)

    if fold:
        fold_candidates(fname, source_dm, coherent=coherent)
    #run SPEGID on candidates
    if speg:
        from prep_speg import prep_speg
        #prep_speg
        #run SPEGID
        prep_speg(fname+'_rfifind.inf')
    if fetch:
        from prep_fetch import prep_fetch_csv
        prep_fetch_csv(fname+'.fil',rank=2)
