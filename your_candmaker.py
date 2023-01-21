#!/usr/bin/env python3

import argparse
import glob
import logging
import os

from rich.logging import RichHandler

from your.utils.math import normalise

os.environ[
    "OPENBLAS_NUM_THREADS"
] = "1"  # stop numpy multithreading regardless of the backend
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

from multiprocessing import Pool
from datetime import datetime

import numpy as np
import pandas as pd

from your.candidate import Candidate, crop
from your.utils.gpu import gpu_dedisp_and_dmt_crop
from your.utils.misc import YourArgparseFormatter
from plot_h5 import plot_h5
import textwrap

logger = logging.getLogger()


def cpu_dedisp_dmt(cand, args):

    cand.dmtime(target="GPU")
    print("Made DMT")

    cand.dedisperse()
    print("Made Dedispersed profile")
    #resize
    #crop out1 the center data
    tsamp = cand.your_header.tsamp
    #gotta crop to the widths we want
    width = args.ws/1000 #convert to seconds
    extra = cand.dedispersed.shape[0]-int(width/tsamp)
    start = int(extra//2)
    length = int(width//tsamp)
    cand.dedispersed = crop(cand.dedispersed, start, length, 0)
    cand.dmt = crop(cand.dmt, start, length, 1)

    #this might not be 100% accurate... but whatever
    downsample_factor = cand.dedispersed.shape[0]//args.time_size+4
    cand.your_header.time_decimation_factor = downsample_factor

    cand.resize(
        key="ft", size=args.time_size+4, axis=0, anti_aliasing=True, mode="constant"
    )
    cand.resize(
        key="dmt", size=args.time_size+4, axis=1, anti_aliasing=True, mode="constant"
    )
    cand.dedispersed = crop(cand.dedispersed, 2, args.time_size, 0)
    cand.dmt = crop(cand.dmt, 2, args.time_size, 1)
    return cand


def cand2h5(cand_val):
    """
    TODO: Add option to use cand.resize for reshaping FT and DMT
    Generates h5 file of candidate with resized frequency-time and DM-time arrays
    :param cand_val: List of candidate parameters (filename, snr, width, dm, label, tcand(s))
    :type cand_val: Candidate
    :return: None
    """
    (
        filename,
        snr,
        width,
        dm,
        label,
        tcand,
        kill_mask_path,
        num_files,
        args,
        gpu_id,
    ) = cand_val
    if os.path.exists(str(kill_mask_path)):
        logger.info(f"Using mask {kill_mask_path}")
        kill_chans = np.loadtxt(kill_mask_path, dtype=np.int)
    else:
        logger.debug("No Kill Mask")

    fname, ext = os.path.splitext(filename)
    if ext == ".fits" or ext == ".sf":
        if num_files == 1:
            files = [filename]
        else:
            files = glob.glob(fname[:-5] + "*fits")
            if len(files) != num_files:
                raise ValueError(
                    "Number of fits files found was not equal to num_files in cand csv."
                )
    elif ext == ".fil":
        files = [filename]
    else:
        raise TypeError("Can only work with list of fits file or filterbanks")

    logger.debug(f"Source file list: {files}")
    print("loading candidate")
    cand = Candidate(
        files,
        snr=snr,
        width=width,
        dm=dm,
        label=label,
        tcand=tcand,
        device=gpu_id,
        spectral_kurtosis_sigma=args.spectral_kurtosis_sigma,
        savgol_frequency_window=args.savgol_frequency_window,
        savgol_sigma=args.savgol_sigma,
        flag_rfi=args.flag_rfi,
        min_samp=25000,
    )
    # get the tsamp
    tsamp = cand.your_header.tsamp
    cand.min_samp = (args.ws/1000)//tsamp
    print(f"min width: {cand.min_samp} samples")
    if os.path.exists(str(kill_mask_path)):
        kill_mask = np.zeros(cand.nchans, dtype=np.bool_)
        kill_mask[kill_chans] = True
        cand.kill_mask = kill_mask
    print("Getting chunk")
    cand.get_chunk(for_preprocessing=True)
    if cand.format == "fil":
        cand.fp.close()
    logger.info("Got Chunk")
    print("dedispersing")
    cand = cpu_dedisp_dmt(cand, args)
    print("Resizing")
    cand.resize(
        key="ft", size=args.frequency_size, axis=1, anti_aliasing=True, mode="constant"
    )
    logger.info(f"Resized Frequency axis of FT to fsize: {cand.dedispersed.shape[1]}")
    print("normalising dmt")
    cand.dmt = normalise(cand.dmt)
    print("normalising ft")
    cand.dedispersed = normalise(cand.dedispersed)
    fout = cand.save_h5(out_dir=args.fout)
    fout2 = plot_h5(fout,save=True,detrend_ft=True)
    logger.debug(f"Filesize of {fout} is {os.path.getsize(fout)}")
    if not os.path.isfile(fout):
        raise IOError(f"File with {cand.id} not written")
    if os.path.getsize(fout) < 100 * 1024:
        raise ValueError(f"File with id: {cand.id} has issues! Its size is too less.")
    logger.info(fout)
    del cand
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="your_candmaker.py",
        description="Your candmaker! Make h5 candidates from the candidate csv files",
        formatter_class=YourArgparseFormatter,
        epilog=textwrap.dedent(
            """\
        `your_candmaker.py` can be used to make candidate cutout files. Some additional notes for this script: 
        - These files are generated in HDF file format. 
        - The output candidates have been preprocessed and consists of Dedispersed Frequency-Time and DM-Time information of the candidate. 
        - The input should be a csv file containing the parameters of the candidates. The input csv file should contain the following fields: 
                - file: Filterbank or PSRFITs file containing the data. In case of multiple files, this should contain the name of first file. 
                - snr: Signal to Noise of the candidate.
                - width: Width of candidate as log2(number of samples). 
                - dm: DM of candidate
                - label: Label of candidate (can be just set to 0, if not known)
                - stime: Start time (seconds) of the candidate.
                - chan_mask_path: Path of the channel mask file. 
                - num_files: Number of files. 
            """
        ),
    )
    parser.add_argument("-v", "--verbose", help="Be verbose", action="store_true")
    parser.add_argument(
        "-fs",
        "--frequency_size",
        type=int,
        help="Frequency size after rebinning",
        default=256,
    )
    parser.add_argument(
        "-g",
        "--gpu_id",
        help="GPU ID (use -1 for CPU). To use multiple GPUs (say with id 2 and 3 use -g 2 3",
        nargs="+",
        required=False,
        default=[-1],
        type=int,
    )
    parser.add_argument(
        "-ts", "--time_size", type=int, help="Time length after rebinning", default=256
    )
    parser.add_argument(
        "-c",
        "--cand_param_file",
        help="csv file with candidate parameters",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-n",
        "--nproc",
        type=int,
        help="number of processors to use in parallel (default: 2)",
        default=2,
    )
    parser.add_argument(
        "-o",
        "--fout",
        help="Output file directory for candidate h5",
        type=str,
        default=".",
    )
    parser.add_argument(
        "-r",
        "--flag_rfi",
        help="Turn on RFI flagging",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-sksig",
        "--spectral_kurtosis_sigma",
        help="Sigma for spectral kurtosis filter",
        type=float,
        default=4,
        required=False,
    )
    parser.add_argument(
        "-sgsig",
        "--savgol_sigma",
        help="Sigma for savgol filter",
        type=float,
        default=4,
        required=False,
    )
    parser.add_argument(
        "-sgfw",
        "--savgol_frequency_window",
        help="Filter window for savgol filter (MHz)",
        type=float,
        default=15,
        required=False,
    )
    parser.add_argument(
        "-opt",
        "--opt_dm",
        dest="opt_dm",
        help="Optimise DM",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-ws",
        "--window_size",
        dest="ws",
        help="window size of the outcoming plot (in ms)",
        type=float,
        default=500,
    )

    parser.add_argument(
        "--no_log_file", help="Do not write a log file", action="store_true"
    )
    values = parser.parse_args()

    logging_format = (
        "%(asctime)s - %(funcName)s -%(name)s - %(levelname)s - %(message)s"
    )
    log_filename = (
        values.fout
        + "/"
        + datetime.utcnow().strftime("your_candmaker_%Y_%m_%d_%H_%M_%S_%f.log")
    )
    print(log_filename)

    if not values.no_log_file:
        if values.verbose:
            logging.basicConfig(
                filename=log_filename,
                level=logging.DEBUG,
                format=logging_format,
            )
        else:
            logging.basicConfig(
                filename=log_filename, level=logging.INFO, format=logging_format
            )
    else:
        if values.verbose:
            logging.basicConfig(
                level=logging.DEBUG,
                format=logging_format,
                handlers=[RichHandler(rich_tracebacks=True)],
            )
        else:
            logging.basicConfig(
                level=logging.INFO,
                format=logging_format,
                handlers=[RichHandler(rich_tracebacks=True)],
            )

    logging.info("Input Arguments:-")
    for arg, value in sorted(vars(values).items()):
        logging.info("%s: %r", arg, value)

    if -1 not in values.gpu_id:
        from numba.cuda.cudadrv.driver import CudaAPIError

        for gpu_ids in values.gpu_id:
            logger.info(f"Using the GPU {gpu_ids}")
        if len(values.gpu_id) > 1:
            from itertools import cycle

            gpu_id_cycler = cycle(range(len(values.gpu_id)))
    else:
        logger.info(f"Using CPUs only")
    cand_pars = pd.read_csv(values.cand_param_file)
    # Randomly shuffle the candidates, this is so that the high DM candidates are spread through out
    # Else they will clog the GPU memory at once
    cand_pars.sample(frac=1).reset_index(drop=True)
    process_list = []
    for index, row in cand_pars.iterrows():
        if len(values.gpu_id) > 1:
            # If there are more than one GPUs cycle the candidates between them.
            gpu_id = next(gpu_id_cycler)
        else:
            gpu_id = values.gpu_id[0]
        process_list.append(
            [
                row["file"],
                row["snr"],
                2 ** row["width"],
                row["dm"],
                row["label"],
                row["stime"],
                row["chan_mask_path"],
                row["num_files"],
                values,
                gpu_id,
            ]
        )

    with Pool(processes=values.nproc) as pool:
        pool.map(cand2h5, process_list, chunksize=1)
    # for p in process_list:
        # cand2h5(p)
