#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=1:00:00
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=automated_filterbank
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#1 is number of splits
#2 is filterbank file
#3 is DM
#if we are not splitting then just set out as the fil file

#path to automated filterbank file script locations
#this is set when you run a batch script by default
#load the module needed
module use /project/6004902/modulefiles
module load presto
#unload scipy-stack because CC has an outdates version of scipy
#module unload scipy-stack
#load my own scipy stack
#source ~/SPEGID/bin/activate
AFP=$4
#check that the filterbank file exists this prevents accidental deletion of files with the later rm command
if test -f "$2"; then
    if [ $1 -gt 1 ]
    then
        OUT=`python $AFP/split_filterbank.py $1 $2 ${SLURM_TMPDIR}`
    else
        OUT=$2
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
        cp $FILFILE $SPFILES
        #run pipeline and prep_fetch prep spegID
        python $AFP/gwg_cand_search_pipeline.py --dm $3 --speg --fetch --no_fft --rfifind --sk_mad --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}/$i"
        #remove the extra fil files
        #rm "$SPFILES/$FIL"
        #remove the .dat files
        rm "$SPFILES"/*.dat
        ((i=i+1))
    done
    #uncomment this code if you want to make a folder and shove everything in there, if you're using process_all_fil.sh, it already makes folder for you.
    #now copy all the files back
    #if [ ! -d $2 ]; then
    	#FN=$(echo "$2" | cut -f 1 -d '.')
	#mkdir $FN
    #fi
    #cp -r ${SLURM_TMPDIR}/* $FN
    cp -r ${SLURM_TMPDIR}/* .
fi
