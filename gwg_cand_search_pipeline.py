import numpy as np
import subprocess
import os, glob
import argparse
import re

from presto.filterbank import FilterbankFile
import pipeline_config
from sk_mad_rficlean import sk_mad_rfi_excision
import sys
def write_mask_file(filename, extra_zap_chans, header):
    """Write a mask numpy array as a rfifind .mask file
    filename: filename to write mask to (will add a .mask extension if not present)
    maskarr: 2D mask array of shape (nint, nchan)
             where the channels follow rfifind convention (index 0 corresponds to the lowest frequency channel)
             and 1/True in the array means mask, 0/False means don't mask
    header: dictionary which must contain keys:
        'time_sig' (float) - from rfifind options, PRESTO default is 10
        'freq_sig' (float) - from rfifind options, PRESTO default is 4
        'MJD' (float) - starting MJD of the observation
        'dtint' (float) - length of one time interval in seconds
        'lofreq' (float) - center frequency of lowest channel
        'df' (float) - width of one frequency channel
        'nchan' (int) - number of frequency channels in the mask
        'nint' (int) - number of time intervals in the mask
        'ptsperint' (int) - number of time samples per interval
    """
    # Massage inputs
    if filename[-5:] != ".mask":
        filename += ".mask"
    header_params = [
        np.array(header.time_sig, dtype=np.float64),
        np.array(header.freq_sig, dtype=np.float64),
        np.array(header.MJD, dtype=np.float64),
        np.array(header.dtint, dtype=np.float64),
        np.array(header.lofreq, dtype=np.float64),
        np.array(header.df, dtype=np.float64),
        np.array(header.nchan, dtype=np.int32),
        np.array(header.nint, dtype=np.int32),
        np.array(header.ptsperint, dtype=np.int32),
    ]

    # Check maskarr shape matches nint and nchan in header

    zap_chans = list(header.mask_zap_chans)
    zap_ints = list(header.mask_zap_ints)
    zap_chans_per_int = header.mask_zap_chans_per_int
    # Write to file
    with open(filename, "wb") as fout:
        for variable in header_params:
            variable.tofile(fout)

        for zap in extra_zap_chans:
            if not (zap in zap_chans):
                zap_chans.append(zap)
        zap_chans.sort()
        zap_chans = np.array(zap_chans)
        zap_ints = np.array(zap_ints)
        for i,zap in enumerate(zap_chans_per_int):
            for z in zap_chans:
                if not(z in zap):
                    zap = np.append(zap,z)
                    zap_chans_per_int[i] = zap

        nzap_chan = np.asarray(zap_chans.size, dtype=np.int32)
        nzap_chan.tofile(fout)
        if nzap_chan:
            zap_chans.tofile(fout)

        nzap_int = np.asarray(zap_ints.size, dtype=np.int32)
        nzap_int.tofile(fout)
        if nzap_int:
            zap_ints.tofile(fout)

        nzap_per_int = []
        for an_arr in zap_chans_per_int:
            nzap_per_int.append(an_arr.size)
        if len(nzap_per_int) != header.nint:
            raise AttributeError("BUG: nzap_per_int should be of length nint!")
        nzpi = np.asarray(nzap_per_int, dtype=np.int32)
        nzpi.tofile(fout)

        # rfifind.py only calls fromfile if nzap != 0 and nzap != nchan
        for i in range(header.nint):
            if nzap_per_int[i]:
                if nzap_per_int[i] != header.nchan:
                    zap_chans_per_int[i].tofile(fout)



#original GWG pipeline written by Chiamin
def run_rfifind(fname,dead_gpus=''):
    pipeline_config_mask = pipeline_config.ignorelist.split(',')
    if dead_gpus!='':
        dead_gpu_mask = dead_gpus.split(',')
        #combine the two masks
        final_mask = []
        pipeline_config_mask = list(int(pgm) for pgm in pipeline_config_mask)
        dead_gpu_mask = list(int(dgm) for dgm in dead_gpu_mask)
        for dgm in dead_gpu_mask:
            if dgm in pipeline_config_mask:
                #do nothing
                pass
            else:
                #if something in the dead gpu mask isn't in the pipe config mask
                pipeline_config_mask.append(dgm)
                print('ignoring ',dgm)

    #conver pipeline config mask back into string
    ignore_chan_string = ''
    for i,chan in enumerate(pipeline_config_mask):
        if i==0:
            ignore_chan_string = str(chan)
        else:
            ignore_chan_string = ignore_chan_string+','+str(chan)
    print('ignoring these channels', ignore_chan_string)
    rfifind_command = 'rfifind -time 1 -timesig 10 -freqsig 10 -ignorechan %s -o %s %s.fil' %(ignore_chan_string,fname,fname)
    # rfifind_command = 'rfifind -time 1 -ignorechan %s -o %s %s.fil' %(ignore_chan_string,fname,fname)

    print(rfifind_command)
    try:
        run_rfifind_cmd = subprocess.check_call([rfifind_command], shell=True)
    except subprocess.CalledProcessError:
        import traceback
        traceback.print_exc()
        [print(f) for f in os.listdir('.')]
        sys.exit(1)


def run_ddplan(fname,dm):
    if dm>20.1:
        dml=dm-20
    else:
        dml=0
    dmh=dm+20
    #run the ddplan in my current directory, it's got the rfi masking included
    import pathlib
    #run the ddplan that lies within the directory of this file because the default presto one can't do masks
    path=pathlib.Path(__file__).parent.absolute()
    # ignorechan= pipeline_config.ignorechan
    fb_file = FilterbankFile(f"{fname}.fil")
    dt = fb_file.dt
    ddplan_command = f"python {path}/DDplan.py -r 0.1 -c {dm} -l {dml} -d {dmh} -s 256 -o {fname}_ddplan -w {fname}.fil"
    print(ddplan_command)
    # ddplan_command = "python %s/DDplan.py -l %.2f -d %.2f -s 256 -o %s_ddplan -w %s.fil" %(path,dml,dmh,fname,fname)
    try:
        run_ddplan = subprocess.check_call([ddplan_command],shell=True)
        #run_ddplan.wait()
        prepsubband_command = "python dedisp_%s.py" %(fname)
        run_ddplan_python = subprocess.check_call([prepsubband_command],shell=True)
        #run_ddplan_python.wait()
    except subprocess.CalledProcessError:
        import traceback
        traceback.print_exc()
        [print(f) for f in os.listdir('.')]
        sys.exit(1)


def run_sp(fname):
    #I set -m to 300, but I don't think I need 300 because it's in bins
    # sp_command = 'single_pulse_search.py -b -m 300 %s*.dat' %(fname)
    sp_command = 'single_pulse_search.py -t 8 %s*.dat' %(fname)
    print(sp_command)
    failed=True
    try:
        run_sp_cmd = subprocess.check_call([sp_command],shell=True)
        #run_sp_cmd.wait()
    except subprocess.CalledProcessError:
        [print(f) for f in os.listdir('.')]
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--fil', type=str, help='Input filterbank file')
    parser.add_argument('--dedisp',action='store_true',help='Run prepsubband and dedisperse the data')
    parser.add_argument('--dm', type=float,help='DM of the candidate. This will determine the max DM to search for pulsar')
    parser.add_argument('--sp',action='store_true',help='Run single pulse search')
    parser.add_argument('--speg',action='store_true',help='creates the SPEGID files')
    parser.add_argument('--fetch',action='store_true',help='creates the FETCH files')
    parser.add_argument('--rfifind',action='store_true',help='Runs rfifind using the configuration in pipeline config')
    parser.add_argument('--dead_gpu',type=str,default='',help='use this option if you want to input a mask for dead GPUs')
    parser.add_argument('--slurm',type=str,help='specifies the root folder to output to, this can be useful on computecanada to reduce IO of files, we use the ${SLURM_TMPDIR} on CC')

    args = parser.parse_args()

    fil = args.fil
    source_dm = args.dm
    sp = args.sp
    speg = args.speg
    fetch =args.fetch 
    dedisp = args.dedisp
    rfifind = args.rfifind
    dead_gpu=args.dead_gpu

    slurm=args.slurm
    if slurm:
        os.chdir(slurm)
    print('Running RFI mitigation')
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


    if rfifind:
        run_rfifind(fname,dead_gpu)
    if dedisp:
        #run ddplan
        run_ddplan(fname,source_dm) 
    if sp:
        run_sp(fname)
    #run SPEGID on candidates
    if speg:
        from prep_speg import prep_speg
        # prep_speg
        #run SPEGID
        prep_speg(fname+'_rfifind.inf')
    #prep the file needed for fetch
    if fetch:
        from prep_fetch import prep_fetch_csv
        prep_fetch_csv(fname+'.fil',rank=5,dm=source_dm)
