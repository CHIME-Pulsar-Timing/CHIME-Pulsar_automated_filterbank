#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=8:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=automated_filterbank
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#1 is number of splits
#2 is DM
#3 is filterbank file
#4 is the location of the scripts
#if we are not splitting then just set out as the fil file

#path to automated filterbank file script locations
#this is set when you run a batch script by default
#load the module needed
module use /project/6004902/modulefiles
module load presto
AFP=$4
#check that the filterbank file exists this prevents accidental deletion of files with the later rm command
#SLURM_TMPDIR='new'
if test -f "$3"; then
    if [ $1 -gt 1 ]
    then
        OUT=`python $AFP/split_filterbank.py $1 $3 ${SLURM_TMPDIR}`
    else
        OUT=$3
        cp -d $3 ${SLURM_TMPDIR}
    fi
    i=0
    for FIL in $OUT;
    do
        echo $FIL
        SPFILES="${SLURM_TMPDIR}/$i"
        if [ ! -d $SPFILES ]; then
            mkdir $SPFILES
        fi
        FILFILE="${SLURM_TMPDIR}/$FIL"
        mv $FILFILE $SPFILES
        #copy the killfile into the folder
        #run pipeline and prep_fetch prep spegID
        #python $AFP/gwg_cand_search_pipeline.py --dm $2 --speg --fetch --no_fft --rfifind --sk_mad --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}/$i"
        #don't run sk_mad
        python $AFP/gwg_cand_search_pipeline.py --dm $2 --speg --fetch --no_fft --rfifind --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}/$i"

        #remove the extra fil files
        #rm "$SPFILES/$FIL"
        #rm "$SPFILES/"*sk_mad.fil
        #remove the .dat files
        rm "$SPFILES"/*.dat
        #tarball the infs and singlepulse files
        tar cf "$SPFILES/${FIL}_singlepulse.tar.gz" "$SPFILES/"*.singlepulse
        tar cf "$SPFILES/${FIL}_inf.tar.gz" "$SPFILES/"*DM*.inf
        rm "$SPFILES"/*DM*.inf
        rm "$SPFILES"/*DM*.singlepulse
        ((i=i+1))
    done
    #uncomment this code if you want to make a folder and shove everything in there, if you're using process_all_fil.sh, it already makes folder for you.
    #now copy all the files back
    #if [ ! -d $3 ]; then
    	#FN=$(echo "$3" | cut -f 1 -d '.')
	#mkdir $FN
    #fi
    #cp -r ${SLURM_TMPDIR}/* $FN
    cp -r ${SLURM_TMPDIR}/* .
fi
