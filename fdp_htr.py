#!/usr/bin/env python3
from presto import sigproc
import sys
import numpy as np
import logging


def get_dtype(nbits):
    """
    Returns:
        dtype of the data
    """
    if nbits == 8:
        return np.uint8
    elif nbits == 16:
        return np.uint16
    elif nbits == 32:
        return np.float32
    else:
        raise RuntimeError(f"nbits={nbits} not supported")

def get_nbits(dtype):
    """
    Returns:
        number of bits of the data
    """
    if dtype == np.uint8:
        return 8
    elif dtype == np.uint16:
        return 16
    elif dtype == np.float32:
        return 32
    else:
        raise RuntimeError(f"dtype={dtype} not supported")

def write_header(header, outfile):
    header_list = list(header.keys())
    manual_head_start_end = False
    if header_list[0] != "HEADER_START" or header_list[-1] != "HEADER_END":
        logging.debug(
            f"HEADER_START not first and/or HEADER_END not last in header_list, removing them from header_list (if present) and writing them manually"
        )
        try_remove("HEADER_START", header_list)
        try_remove("HEADER_END", header_list)
        manual_head_start_end = True

    if manual_head_start_end:
        outfile.write(sigproc.addto_hdr("HEADER_START", None))
    for paramname in header_list:
        if paramname not in sigproc.header_params:
            # Only add recognized parameters
            continue
        logging.debug("Writing header param (%s)" % paramname)
        value = header[paramname]
        outfile.write(sigproc.addto_hdr(paramname, value))
    if manual_head_start_end:
        outfile.write(sigproc.addto_hdr("HEADER_END", None))


log = logging.getLogger(__name__)

fn = sys.argv[1]
chunk_size = 30674 # for 40.96us data this is ~1.25s, may want to change
thresh_sig = 4.5

header, hdrlen = sigproc.read_header(fn)
nspecs = int(sigproc.samples_per_file(fn, header, hdrlen))
nchans = header["nchans"]
arr_dtype = get_dtype(header["nbits"])

if header["nifs"] != 1:
    log.error(f"Code not written to deal with unsummed polarization data")


#loop through chunks
loop_iters = int(nspecs//chunk_size)
if nspecs % chunk_size:
    loop_iters += 1
fn_clean = fn.strip('.fil')
fdp_fn = f"{fn_clean}_fdp.fil"
new_fil = open(fdp_fn, "wb")
write_header(header, new_fil)
fil = open(fn, "rb")
fil.seek(hdrlen)

for i in range(loop_iters):
    print(f"{i+1}/{loop_iters}", end='\r', flush=True)
    spec = (
        np.fromfile(fil, count=chunk_size*nchans, dtype=arr_dtype)
        .reshape(-1, nchans)
    )
    med = np.median(spec,axis=0)
    std = np.std(spec,axis=0)
    #set the thresold
    threshold = med - thresh_sig*std
    #find values below threshold and replace with the median
    mask = np.where(spec < threshold)
    spec[mask] = med[mask[1]]
    new_fil.write(spec.ravel().astype(arr_dtype))

fil.close()
new_fil.close()
