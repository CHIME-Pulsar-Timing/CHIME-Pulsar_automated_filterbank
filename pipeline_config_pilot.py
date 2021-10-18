# options for what to run in the pipeline:
run_sk_mad = False  # run sk_mad RFI excision instead of rfifind (NB don't think this works atm)
run_rfifind = True  # run rfifind, using configuration found below
run_ddplan = False  # run DDplan - uses range +-20 of target DM
run_dedisp = False  # run a dedisp_<filename>.py script output by DDplan, aka dedisperse the data using prepsubband
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

# dedispersion configuration:
use_coherent_ddplan = True  # use coherent ddplan defined below
# . . . from *_cand_search_pipeline looks like this doesn't actually get used

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

def coherent_ddplan(tsamp, dm, coherent_dm):

    if tsamp < 0.000327:
        tsamp_fac = 2
    else:
        tsamp_fac = 1

    if coherent_dm - 80.0 < dm:
        dms = 0.05
        ds = 1
        sb = 32
    elif coherent_dm - 140.0 < dm:
        dms = 0.10
        ds = 2
        sb = 32
    elif coherent_dm - 300.0 < dm:
        dms = 0.20
        ds = 4
        sb = 32
    elif coherent_dm - 540.0 < dm:
        dms = 0.40
        ds = 8
        sb = 32
    elif coherent_dm - 780.0 < dm:
        dms = 0.40
        ds = 8
        sb = 64
    else:
        dms = 0.40
        ds = 16
        sb = 64
    if dm >= 260.0:
        sb = sb * tsamp_fac

    return dms, ds, sb

def ddplan(tsamp, dm):

    if tsamp < 0.000327:
        tsamp_fac = 2

    if dm < 80.0:
        dms = 0.05
        ds = 1
        sb = 32
    elif dm < 140.0:
        dms = 0.10
        ds = 2
        sb = 32
    elif dm < 300.0:
        dms = 0.20
        ds = 4
        sb = 32
    elif dm < 540.0:
        dms = 0.40
        ds = 8
        sb = 32
    elif dm < 780.0:
        dms = 0.40
        ds = 8
        sb = 64
    else:
        dms = 0.40
        ds = 16
        sb = 64
    if dm >= 260.0:
        sb = sb * tsamp_fac

    return dms, ds, sb

coherent_dm_set = [-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,100,120,140,180,220,260,300,380,460,540,620,700,780,860,940]
dm_set = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,160,180,200,220,240,260,280,300,380,460,540,620,700,780,860,940]
ffa_dm_set = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.0,9.2,9.4,9.6,9.8] #[0.0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.0,6.4,6.8,7.2,7.6,8.0,8.4,8.8,9.2,9.6]
