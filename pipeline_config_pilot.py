# options for what to run in the pipeline:
run_sk_mad = False  # run sk_mad RFI excision instead of rfifind (NB don't think this works atm)
run_rfifind = True  # run rfifind, using configuration found below
run_sp_ddplan = False  # run DDplan using range +-20 around target DM. Suitable for single pulse searches (SPEGID & FETCH pipeline)
# = False  # make the dedsip.py for the file based on a ddplan saved as an npz based on configurations below
run_dedisp_from_ddplan = True  # run a dedisp_<filename>.py script output by DDplan, aka dedisperse the data using prepsubband
run_fft = False  # run FFT search
run_ffa = False # run  FFA (fast folding algorithm) search
fold_candidates = False  # fold FFA and/or FFT candidates
run_sp = False  # run single pulse search
run_prep_speg = False  # create the files needed to run SPEGID
run_prep_fetch = False  # create the files needed to run FETCH

# rfifind configuration:
#as of 22/4/2021,merged kill file from brad and adam
zaplist = "1023,989,988,987,986,985,984,983,982,981,980,979,978,977,976,910,909,908,907,906,905,904,903,902,901,900,899,898,897,896,895,894,893,892,891,890,889,888,887,886,885,884,883,882,881,880,879,878,877,876,875,874,873,872,871,870,869,868,867,866,865,864,863,862,861,860,859,858,857,856,855,854,853,852,851,850,849,848,847,846,845,838,837,836,835,834,833,832,831,830,829,828,827,826,825,824,823,822,821,820,819,818,817,816,815,814,813,805,804,803,802,801,800,799,798,797,796,795,794,793,792,791,790,789,788,787,786,785,784,783,782,781,780,779,778,777,776,775,774,773,772,771,770,769,618,617,616,615,614,613,612,611,610,609,608,607,606,605,604,603,602,601,600,599,598,597,596,595,594,593,592,591,590,589,588,587,586,585,584,583,582,581,580,579,578,577,576,575,574,573,572,571,570,569,568,567,566,565,564,563,562,561,560,559,558,557,556,555,554,471,470,469,468,467,466,465,464,463,462,461,460,459,458,457,456,455,451,448,439,438,437,436,435,434,433,432,431,430,429,428,427,426,392,391,390,389,388,387,386,385,384,383,382,381,380,379,346,345,344,343,342,341,340,339,338,337,336,335,334,333,332,331,330,269,268,267,266,265,264,263,262,261,235,234,233,232,177,169,168,167,166,165,164,163,150,149,148,140,136"
rfiblocks = 4
rfizerodm = True  # run rfifind with option -zerodm
# NOT QUITE SURE IF THIS PLAYS NICE WITH COHERENTLY DEDISPERSED DATA

# dedispersion configuration:
# for make_dedisp_from_template
#template_ddplan_fname = None  # "/home/kcrowter/survey/amalgamated_ddplan_cdms_30_91_152_213_273_334_395_456.npz"
#template_ddplan_has_multiple_cdms = True  # set True if template_ddplan is or multiple coherently dedispersed passes


# fft configuration:
fftzaplist = None  # zaplist to remove from fft files (NB do not understand how this is currently implemented. need birdies in there too)
rednoise = True  # remove rednoise
run_binary = False  # run binary (acceleration) search (only works if run_fft is True)
zmax = 100  # zmax for accelsearch
wmax = 0  # wmax for accelsearch

def foldplan(fname, accelfile, accelcand, canddm, candperiod, sk_mad=False):
    if candperiod < 0.02:
        n, npfact, ndmfact = 32, 2, 2
        npart = 90
        otheropts = ""
    elif candperiod < 0.1:
        n, npfact, ndmfact = 64, 2, 1
        npart = 90
        otheropts = "-pstep 1 -pdstep 2 -dmstep 2"
    elif candperiod < 1:
        n, npfact, ndmfact = 128, 1, 1
        npart = 90
        otheropts = "-pstep 1 -pdstep 2 -dmstep 1"
    elif candperiod < 10:
        n, npfact, ndmfact = 256, 1, 1
        npart = 45
        otheropts = "-nosearch -slow"
    elif candperiod < 60:
        n, npfact, ndmfact = 256, 1, 1
        npart = 30
        otheropts = "-nosearch -slow"
    else:
        n, npfact, ndmfact = 256, 1, 1
        npart = 15
        otheropts = "-nosearch -slow"
    nsub = 128
    if not sk_mad:
        otheropts += " -mask {0}_rfifind.mask {0}.fil".format(fname)
    elif sk_mad:
        otheropts += " {0}.fil".format(fname)
    if accelcand:
        outname = accelfile.rstrip('.cand')
        prepfold_opt = "prepfold -noxwin -accelfile {0} -accelcand {1} -dm {2} \
                        -nsub {3} -npart {4} -n {5} -npfact {6} -ndmfact {7} \
                        {8} -o {9}".format(accelfile, accelcand, canddm, nsub, \
                        npart, n, npfact, ndmfact, otheropts, outname)
    else:
        outname = accelfile
        prepfold_opt = "prepfold -noxwin -p {0} -dm {1} -nsub {2} -npart {3} \
                        -n {4} -npfact {5} -ndmfact {6} {7} -o {8}".format( \
                        candperiod, canddm, nsub, npart, n, npfact, ndmfact, \
                        otheropts, outname)

    return prepfold_opt
