import numpy as np
import subprocess
import os, glob
import argparse
import re

from presto.filterbank import FilterbankFile
import pipeline_config_pilot as pipeline_config
from sk_mad_rficlean import sk_mad_rfi_excision
import sys

from ddplan_rejig import (frankenstein_ddplan, ddplan)

#original GWG pipeline written by Chiamin
def run_rfifind(fname):
    options = ""
    if pipeline_config.rfizerodm:
        options += "-zerodm"
    rfifind_command = 'rfifind -blocks %d -zapchan %s %s -o %s %s.fil' %(pipeline_config.rfiblocks,pipeline_config.zaplist,options,fname,fname)
    try:
        run_rfifind_cmd = subprocess.check_call([rfifind_command], shell=True)
    except subprocess.CalledProcessError:
        import traceback
        traceback.print_exc()
        [print(f) for f in os.listdir('.')]
        sys.exit(1)


def run_sk_mad(fname,fil):
    sk_mad_rfi_excision(fname,fil)
    fnamenew = str(fname)+'_sk_mad'
    return fnamenew

def run_sp_ddplan(fname,dm):
    if dm>20.1:
        dml=dm-20
    else:
        dml=0
    dmh=dm+20
    #run the ddplan in my current directory, it's got the rfi masking included
    import pathlib
    #run the ddplan that lies within the directory of this file because the default presto one can't do masks
    path=pathlib.Path(__file__).parent.absolute()
    ddplan_command = "python %s/DDplan.py -l %.2f -d %.2f -s 256 -o %s_ddplan -w %s.fil" %(path,dml,dmh,fname,fname)
    try:
        run_ddplan = subprocess.check_call([ddplan_command],shell=True)
        #run_ddplan.wait()
    except subprocess.CalledProcessError:
        import traceback
        traceback.print_exc()
        [print(f) for f in os.listdir('.')]
        sys.exit(1)

def run_make_dedisp_from_template_ddplan(filfname, ddplanfname, multiple_cDMs=False, mask=True):
    """Make a dedisp.py for filfname based on a pre-computed ddplan stored in a .npz (ddplanfname)
    If the ddplan is for multiple coherently dedispersed observations set multiple_cDMs=True"""
    if ddplanfname is not None:
        if multiple_cDMs:
            ddp = frankenstein_ddplan.read_from_npz(ddplanfname)
            ddp.write_dedisp_for(filfname, to_file=True, mask=mask)
        else:
            ddp = ddplan.read_from_npz(ddplanfname)
            ddp.write_dedisp_py(filfname, to_file=True, mask=mask)
    else:
        print("No template ddplan file given")
        sys.exit(1)


def run_dedisp_from_ddplan(fname):
    """Execute a dedisp_<fname>.py prepsubband script output by DDplan"""
    try:
        prepsubband_command = "python dedisp_%s.py" %(fname)
        run_ddplan_python = subprocess.check_call([prepsubband_command],shell=True)
    except subprocess.CalledProcessError:
        import traceback
        traceback.print_exc()
        [print(f) for f in os.listdir('.')]
        sys.exit(1)

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
    datfiles = sorted(glob.glob(str(fname)+'*DM*.dat'))
    # think the below is there because there's a different set of dms for ffa vs fft - do not know why
#    for dm in pipeline_config.ffa_dm_set:
#        if dm == 0.0:
#            datfiles = sorted(glob.glob(fname+'*DM*'+str(dm)+'0.dat'))
#
#        else:
#            datfiles.extend(sorted(glob.glob(fname+'*DM*'+str(dm)+'0.dat')))

    for dat in datfiles:
        ffa_command = '/psr_scratch/common_utils/ffaGo/ffa.py '+str(dat)
        run_ffa_cmd = subprocess.Popen([ffa_command],shell=True)
        run_ffa_cmd.wait()

def run_sp(fname):
    sp_command = 'single_pulse_search.py %s*.dat' %(fname)
    failed=True
    try:
        run_sp_cmd = subprocess.check_call([sp_command],shell=True)
        #run_sp_cmd.wait()
    except subprocess.CalledProcessError:
        [print(f) for f in os.listdir('.')]
        import traceback
        traceback.print_exc()
        sys.exit(1)


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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--fil', type=str, help='Input filterbank file')
    parser.add_argument('--dm', type=float,help='DM of the candidate. This will determine the max DM to search for pulsar')
    parser.add_argument('--slurm',type=str,help='specifies the root folder to output to, this can be useful on computecanada to reduce IO of files, we use the ${SLURM_TMPDIR} on CC')

    args = parser.parse_args()

    fil = args.fil
    source_dm = args.dm
    slurm = args.slurm
    if slurm:
        os.chdir(slurm)

    sk_mad = pipeline_config.run_sk_mad
    coherent = pipeline_config.use_coherent_ddplan
    fft = pipeline_config.run_fft
    rednoise = pipeline_config.rednoise
    zaplist = pipeline_config.fftzaplist
    binary = pipeline_config.run_binary
    zmax = pipeline_config.zmax
    wmax = pipeline_config.wmax
    ffa = pipeline_config.run_ffa
    sp = pipeline_config.run_sp
    fold = pipeline_config.fold_candidates
    speg = pipeline_config.run_prep_speg
    fetch = pipeline_config.run_prep_fetch
    ddplan_sp = pipeline_config.run_sp_ddplan
    dedisp_from_ddplan = pipeline_config.run_dedisp_from_ddplan
    prepsub = pipeline_config.run_prepsubband
    rfifind = pipeline_config.run_rfifind
    make_dedisp_from_template = pipeline_config.run_make_dedisp_from_template

    #get only the file name
    fname = fil.rstrip('.fil')
    fname = fname.split('/')
    fname = fname[-1]
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

    if sk_mad:
        if os.path.exists('{}_sk_mad.fil'.format(fname)):
            print('sk mad cleaned data already exists')
            fname = '{}_sk_mad'.format(fname)
        else:
            print('Running RFI mitigation - sk_mad')
            fname = run_sk_mad(fname,fil)

    if rfifind:
        print('Running RFI mitigation - rfifind')
        run_rfifind(fname)

    if ddplan_sp:
        # run ddplan
        print('Running DDplan: +-20 around target DM')
        run_sp_ddplan(fname,source_dm)

    if make_dedisp_from_template:
        print(f'Making dedisp.py for {fname} from {pipeline_config.template_ddplan_fname}')
        run_make_dedisp_from_template_ddplan(
            fname,
            pipeline_config.template_ddplan_fname,
            multiple_cDMs=pipeline_config.template_ddplan_has_multiple_cdms
            )

    if dedisp_from_ddplan:
        #run prepsubband plan output by ddplan
        print('Dedispersing from a dedisp_<filename>.py script')
        run_dedisp_from_ddplan(fname)

    if prepsub:
        # run prepsubband based on ddplan/coherent_ddplan defined in config
        print('Running prepsubband based on ddplan/coherent_ddplan defined in config')
        if coherent:
            dmlist = [source_dm-i for i in pipeline_config.coherent_dm_set if source_dm-i > 0]
            dmlist.append(0)
        else:
            dmlist = [i for i in pipeline_config.dm_set if (i < source_dm+20)&(i > source_dm-20)]
        for dm in dmlist:
            run_prepsubband(fname,tsamp,dm,source_dm,coherent)

    #run fft
    if fft:
        print('Running FFT search')
        run_realfft(fname,rednoise,zaplist)
        run_accelsearch(fname,zmax,wmax,binary)
        run_accelsift(fname)

    #run ffa
    if ffa:
        print('Running FFA search')
        run_ffa(fname)
        run_ffa_sift(fname)

    if sp:
        print('Running single pulse search')
        run_sp(fname)

    if fold:
        print('Folding candidates')
        fold_candidates(fname, source_dm, coherent=coherent)

    if speg:
        print('Running SPEGID on single pulse candidates')
        from prep_speg import prep_speg
        #prep_speg
        #run SPEGID
        prep_speg(fname+'_rfifind.inf')

    if fetch:
        print('Running FETCH')
        from prep_fetch import prep_fetch_csv
        prep_fetch_csv(fname+'.fil',rank=2)
