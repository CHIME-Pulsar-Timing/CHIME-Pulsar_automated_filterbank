#!/bin/bash
#SBATCH --account=rrg-istairs-ad
#SBATCH --export=NONE
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=4096M
#SBATCH --cpus-per-task=1
#SBATCH --job-name=automated_filterbank
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
set -euo pipefail
LOCAL=false
while getopts "ld:a:p:" flag
do
    case "${flag}" in
        l) LOCAL=true;;
        d) DM=${OPTARG};;
        a) AFP=${OPTARG};;
        p) p=${OPTARG};;
    esac
done

PULSAR=$(echo "$p" | rev | cut -f2- -d '.' | rev)
EXT="${p##*.}"
if [ "$LOCAL" != true ]; then
    module use /project/6004902/chimepsr-software/v2/environment-modules
    module load presto
    module load chime-psr
    source ~/projects/rrg-istairs-ad/Your_060524/bin/activate
else
    #set slurm tmpdir to current directory
    SLURM_TMPDIR='./'
    SLURM_JOB_ID=1

fi
# module load cuda

#make sure that $p is a file
if test -f "$p"; then
    #rename it FIL
    FIL=$p
    if [ "$LOCAL" != true ]; then
        #copy it to the scratch directory
        cp -d $FIL ${SLURM_TMPDIR}
        echo $FIL
        if [ ! -d ${SLURM_TMPDIR} ]; then
            mkdir ${SLURM_TMPDIR}
        fi
    fi
    n=0
    #basically try catch
    until [ "$n" -ge 1 ]
    do
        echo $EXT
        if [ $EXT == "fits" ]; then
            #set dead gpu string to empty if using fits
            echo "python $AFP/gwg_cand_search_pipeline.py --dm $DM --speg --fetch --rfifind --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}" && break"
            python $AFP/gwg_cand_search_pipeline.py --dm $DM --speg --fetch --rfifind --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}" && break
        else
            if [ "$LOCAL" != true ]; then
                DEAD_GPU=$(get_bad_channel_list.py --fmt presto --type filterbank $FIL)
                echo "python $AFP/gwg_cand_search_pipeline.py --dm $DM --speg --fetch --rfifind --sk_mask --dead_gpu $DEAD_GPU --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}" && break"
                python $AFP/gwg_cand_search_pipeline.py --dm $DM --speg --fetch --rfifind --dead_gpu $DEAD_GPU --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}" && break
            else
                echo "python $AFP/gwg_cand_search_pipeline.py --dm $DM --speg --fetch --rfifind --sk_mask --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}" && break"
                python $AFP/gwg_cand_search_pipeline.py --dm $DM --speg --fetch --rfifind --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}" && break
            fi
        fi

        #ANYTHING BEYOND HERE IS DEBUGGING
        #FETCH and SPEGID are slow so lets like ignore that for now
        #python $AFP/gwg_cand_search_pipeline.py --dm $DM --rfifind --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}" && break
        #for rapid tests, only do rfifind
        # python $AFP/gwg_cand_search_pipeline.py --dm $DM --rfifind --fil $FIL --slurm "${SLURM_TMPDIR}" && break
        # Should never get here. if it does exit with an error
        n=$((n+1))
        # sleep 15
        #if it fails, lets copy all the things to my scratch directory then exit with error code
        # PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
        # ERRORS=~/"scratch/errors/${PULSAR}_${SLURM_JOB_ID}"
        # echo "copying error files to ${ERRORS}"
        # df -h
        #mkdir -p $ERRORS
        #cp -r -d ${SLURM_TMPDIR}/* $ERRORS
        echo "ERRORS IN PROCESSING CHECK LOGS OR ENABLE DEBUGGING"
        exit 1
        #remove the extra fil files
    done

    PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
    # rm "${SLURM_TMPDIR}"/$FIL
    #remove the .dat files
    rm "${SLURM_TMPDIR}"/*.dat
    #tarball the infs and singlepulse files
    tar -cf "${SLURM_TMPDIR}/${PULSAR}_singlepulse.tar" "${SLURM_TMPDIR}/"*.singlepulse
    tar -cf "${SLURM_TMPDIR}/${PULSAR}_inf.tar" "${SLURM_TMPDIR}/"*DM*.inf
    rm "${SLURM_TMPDIR}"/*DM*.inf
    rm "${SLURM_TMPDIR}"/*DM*.singlepulse
    #chown -R adamdong:rrg-istairs-ad ${SLURM_TMPDIR}/*
    if [ "$LOCAL" == true ]; then
    # clean up - not needed on compute canada, but nice to run clean up when on my own computer
        #gotta do this for some weird reason on my local machine
        touch "${PULSAR}"_singlepulse.ps
    else
        cp -r ${SLURM_TMPDIR}/* .
    fi
    exit 0
    echo "ALL DONE"
fi
#didn't find file, throw error
echo "filterbank file not found"
exit 1
