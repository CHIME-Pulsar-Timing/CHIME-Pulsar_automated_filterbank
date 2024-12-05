import numpy as np
import subprocess
import os, glob
import argparse
import re
import pipeline_config
import sys
import your_rfi_sk
import logging

#create a fake stream logger to redirect stdout and stderr
class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, level):
       self.logger = logger
       self.level = level
       self.linebuf = ''

    def write(self, buf):
       for line in buf.rstrip().splitlines():
          self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass


#handling exceptions
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.critical(
        "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
    )

#original GWG pipeline written by Chiamin
def run_rfifind(fname,ext,dead_gpus=''):
    pipeline_config_mask = pipeline_config.ignorelist.split(',')
    if dead_gpus!='':
        #legacy code to combine deagpu and ignore list, we don't use this anymore
        dead_gpu_mask = dead_gpus.split(',')
        if len(dead_gpu_mask)<400:
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
                    logging.info(f'ignoring dead gpus {dgm}')
        else:
            logging.info("dead GPU mask is huge, it's probably wrong, ignore")
    else:
        #find our own deadgpu mask
        from sigpyproc import readers as r
        #load 10s of data
        filf = r.FilReader(fname+ext)
        tsamp = filf.header.tsamp
        _ = filf.read_block(0,int(10/tsamp))
        #find the data that has std of 0
        stds = np.std(_,axis=1)
        pipeline_config_mask_sigpyproc = np.where(stds==0)[0]
        pipeline_config_mask = 1023-pipeline_config_mask_sigpyproc


    #conver pipeline config mask back into string
    ignore_chan_string = ''
    for i,chan in enumerate(pipeline_config_mask):
        if i==0:
            ignore_chan_string = str(chan)
        else:
            ignore_chan_string = ignore_chan_string+','+str(chan)
    logging.info(f"RFIFIND out fname {fname}")
    if ignore_chan_string=='':
        #don't have any ignore channels if there's no zaplist
        #print warning
        print("WARNING: no zaplist, so no channels will be ignored")
        rfifind_command = f"rfifind -blocks {pipeline_config.rfiblocks} -o {fname} {fname}{ext}"
    else:
        rfifind_command = f"rfifind -blocks {pipeline_config.rfiblocks} -ignorechan {ignore_chan_string} -o {fname} {fname}{ext}"
    logging.info(rfifind_command)
    try:
        run_rfifind_cmd = subprocess.check_call([rfifind_command], shell=True)
        return ignore_chan_string
    except subprocess.CalledProcessError:
        import traceback
        traceback.print_exc()
        [logging.info(f) for f in os.listdir('.')]
        sys.exit(1)

def run_gsk(fname,ext):
    from gsk import load_rfi_mask
    from gsk import perform_SK
    from rfifind_numpy_tools import write_new_mask_from
    rfimask_fn = f"{fname}_rfifind.mask"
    rfimask,nints = load_rfi_mask(rfimask_fn,"initial_mask.png")
    import os
    #copy the original mask to a new file
    os.system(f"cp {rfimask_fn} {fname}_original_rfifind.mask")
    mask = perform_SK(fname+ext,nints)
    write_new_mask_from(rfimask_fn, mask, rfimask, include_old=True, infstats_too=True)
    rfimask,nints = load_rfi_mask(rfimask_fn,"final_mask.png")

def run_kc_iqrm(fname,ext,ignorechans):
    #for now hardcode the link
    ##this doesn't work and I'm going to give up on it
    from gsk import load_rfi_mask
    import pathlib
    path=pathlib.Path(__file__).parent.absolute()
    rfimask_fn = fname+"_rfifind.mask"
    #copy the original mask to a new file
    os.system(f"cp {rfimask_fn} {fname}_before_kc_rfifind.mask")
    new_rfimask_name = "_01234568_rfifind.mask"
    command = f"python {path}/rfi_pipeline.py --option 0,1,2,3,4,5,6,8 --outfilename {fname}{new_rfimask_name} --rfac 16 --include_rfifind --ignorechans {ignorechans} {rfimask_fn}"
    run_kc_iqrm_cmd = subprocess.check_call([command], shell=True)
    rfimask,nints = load_rfi_mask(rfimask_fn,"kc_mask.png")
    return new_rfimask_name

def run_ddplan(fname,ext,dm,mask_name,ignorelist):
    from presto import infodata
    #run the ddplan in my current directory, it's got the rfi masking included
    import pathlib

    info = infodata.infodata(fname+'_rfifind.inf')
    dt = info.dt
    numchan = info.numchan
    if numchan>=1024:
        subband = int(numchan/4)
    else:
        subband = numchan
    if dm>26.1:
        dml=dm-20
    else:
        dml=6.1
    dmh=dm+20
    #run the ddplan that lies within the directory of this file because the default presto one can't do masks
    path=pathlib.Path(__file__).parent.absolute()
    # ignorechan= pipeline_config.ignorechan
    fname = fname.split('/')[-1]
    if ignorelist=='':
        ddplan_command = f"python {path}/DDplan.py -t {dt} -y {mask_name} -c {dm} -l {dml} -d {dmh} -s {subband} -o {fname}_ddplan -w {fname}{ext}"
    else:
        ddplan_command = f"python {path}/DDplan.py --ignore {ignorelist} -t {dt} -y {mask_name} -c {dm} -l {dml} -d {dmh} -s {subband} -o {fname}_ddplan -w {fname}{ext}"

    logging.info(ddplan_command)
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
        [logging.info(f) for f in os.listdir('.')]
        sys.exit(1)


def run_sp(fname):
    #I set -m to 300, but I don't think I need 300 because it's in bins
    sp_command = 'single_pulse_search.py %s*.dat' %(fname)
    logging.info(sp_command)
    failed=True
    try:
        run_sp_cmd = subprocess.check_call([sp_command],shell=True,executable="/bin/bash")
        #run_sp_cmd.wait()
    except subprocess.CalledProcessError:
        [logging.info(f) for f in os.listdir('.')]
        import traceback
        traceback.print_exc()
        sys.exit(1)

def edit_mask(fname,ext,mask_name):
    #finally edit the mask to include the pipeline things
    ignorelist = pipeline_config.ignorelist
    if ignorelist=='':
        print("WARNING: no zaplist, so no channels will be ignored")
        return None
    command = f"rfifind -mask {fname}{mask_name} -zapchan {ignorelist} -nocompute -o {fname} {fname}.{ext}"
    run_rfifind_cmd = subprocess.check_call([command], shell=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--sk_mask', action='store_true', help="apply a SK mask from 'your' onto the rfifind mask")
    parser.add_argument('--kc_iqrm', action='store_true', help="apply a Kathryn's IQRM mask")
    parser.add_argument('--fil', type=str, help='Input filterbank file')
    parser.add_argument('--dedisp',action='store_true',help='Run prepsubband and dedisperse the data')
    parser.add_argument('--dm', type=float,help='DM of the candidate. This will determine the max DM to search for pulsar')
    parser.add_argument('--sp',action='store_true',help='Run single pulse search')
    parser.add_argument('--speg',action='store_true',help='creates the SPEGID files')
    parser.add_argument('--fetch',action='store_true',help='creates the FETCH files')
    parser.add_argument('--rfifind',action='store_true',help='Runs rfifind using the configuration in pipeline config')
    parser.add_argument('--dead_gpu',type=str,default='',help='use this option if you want to input a mask for dead GPUs')
    parser.add_argument('--slurm',type=str,help='specifies the root folder to output to, this can be useful on computecanada to reduce IO of files, we use the ${SLURM_TMPDIR} on CC')
    parser.add_argument("--log",type=str,help="name of file to write log to",default="automated_filterbank_batch_")


    args = parser.parse_args()

    fil = args.fil
    source_dm = args.dm
    sp = args.sp
    speg = args.speg
    fetch =args.fetch 
    dedisp = args.dedisp
    rfifind = args.rfifind
    dead_gpu = args.dead_gpu
    slurm=args.slurm
    sk_mask = args.sk_mask
    kc_iqrm = args.kc_iqrm
    current_dir = os.getcwd()
    #get only the file name
    if fil.endswith(".fits"):
        fname = fil.rstrip('.fits')
        ext = '.fits'
    elif fil.endswith(".fil"):
        fname = fil.rstrip('.fil')
        ext = '.fil'
    fname = fname.split('/')
    fname = fname[-1]

    logging_fn = os.path.join(current_dir,args.log+fname+".log")
    print(f"logging file path {logging_fn}")
    logging.basicConfig(
        filename=logging_fn,
        filemode="a",
        format="%(asctime)s %(levelname)s:%(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.DEBUG,
        force=True
    )
    # log = logging.getLogger('stdlogger')
    # sys.stdout = StreamToLogger(log,logging.INFO)
    # sys.stderr = StreamToLogger(log,logging.ERROR)
    # logging.info("test logging info")
    # print('Test to standard out')

    if slurm:
        #this is a change to the SLURM tmpdir directory
        os.chdir(slurm)
        slurm_dir = os.getcwd()
        logging.info(f"current cwd: {slurm_dir}")
        # fname = os.path.join(slurm_dir,fname)
        logging.info(f"fname:{fname}")

    if os.path.islink(fil):
        fil = os.readlink(fil)
        if not os.path.isfile(fil):
            logging.info('File does not exist')
            sys.exit()

    logging.info('Running RFI mitigation')
    if rfifind:
        logging.info("Running rfifind")
        ignore_chan_string = run_rfifind(fname,ext,dead_gpu)
        mask_name = "_rfifind.mask"
    if sk_mask:
        run_gsk(fname,ext)
    if kc_iqrm:
        mask_name = run_kc_iqrm(fname,ext,ignore_chan_string)
    if dedisp:
        #run ddplan
        logging.info("Running DDplan")
        run_ddplan(fname,ext,source_dm,mask_name,ignorelist =ignore_chan_string)
    if sp:
        logging.info("Running single pulse search")
        run_sp(fname)
    #run SPEGID on candidates
    if speg:
        from prep_speg import prep_speg
        # prep_speg
        #run SPEGID
        logging.info("running prep speg")
        prep_speg(fname+'_rfifind.inf')
    #prep the file needed for fetch
    if fetch:
        logging.info("running prep fetch")
        from prep_fetch import prep_fetch_csv
        #this script needs the actual file name
        fname = fname.split('/')[-1]
        prep_fetch_csv(fname+ext,float(source_dm),rank=5)

    edit_mask(fname,ext,mask_name)
