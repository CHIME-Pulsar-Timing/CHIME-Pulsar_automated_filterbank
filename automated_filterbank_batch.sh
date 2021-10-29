#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=00:10:00
#SBATCH --mem=256M
#SBATCH --cpus-per-task=1
#SBATCH --job-name=automated_filterbank
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-user=kcrowter@phas.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
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
#********************THIS IS THE LAZY WAY OUT!!!
#PULSAR=$(echo "$3" | cut -f 1 -d '.')
# KC will need to figure out the file movement better but for now:
FN=${3%.fil}
#SLURM_TMPDIR="/home/kcrowter/scratch/survey/$FN"
#SLURM_TMPDIR='/home/adamdong/scratch/tmpdir/'$PULSAR
#SLURM_TMPDIR='/media/adam/1c126a4b-fb16-4471-909f-4b0fda74a5d2/tmpdir/'$PULSAR
#mkdir -p $SLURM_TMPDIR
#SLURM_JOB_ID=1
if test -f "$3"; then
  FIL=$3
  TMP_OUTDIR=${SLURM_TMPDIR}/0
  mkdir ${TMP_OUTDIR}
  cp -d ${FIL} ${TMP_OUTDIR}
  cp 0/${FN}_rfifind* ${TMP_OUTDIR}
  cp 0/dedisp_${FN}.py ${TMP_OUTDIR}

	n=0
		#basically try catch
		until [ "$n" -ge 1 ]
		do
			#python $AFP/gwg_cand_search_pipeline.py --dm $2 --speg --fetch --no_fft --rfifind --sk_mad --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}/$i" && break
			python $AFP/pilot_survey_cand_search_pipeline.py --dm $2 --fil $FIL --slurm ${TMP_OUTDIR} && break
			n=$((n+1))
			sleep 15
			#if it fails, lets copy all the things to my scratch directory then exit with error code
			PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
			ERRORS=~/"scratch/errors/${PULSAR}_${SLURM_JOB_ID}"
			echo "copying error files to ${ERRORS}"
			df -h
			mkdir -p $ERRORS
			   cp -r -d ${SLURM_TMPDIR}/* $ERRORS
			exit 1
		done
        # KC commented out for testing, uncomment later
        # Need to remove everything you don't want copied back
        #remove the extra fil files
        #rm ${TMP_OUTDIR}/$FIL
        #rm ${TMP_OUTDIR}/*sk_mad.fil
        #tarball the .dat files
        tar -czf ${TMP_OUTDIR}/${FIL}_dat.tar -C ${TMP_OUTDIR} *.dat
        rm ${TMP_OUTDIR}/*.dat
        #tarball the infs and singlepulse files
        #tar -czf "$SPFILES/${FIL}_singlepulse.tar" -C "$SPFILES" *.singlepulse
        tar -cff ${TMP_OUTDIR}/${FIL}_inf.tar -C ${TMP_OUTDIR} *DM*.inf
        rm ${TMP_OUTDIR}/*DM*.inf
        #rm "$SPFILES"/*DM*.singlepulse

    #uncomment this code if you want to make a folder and shove everything in there, if you're using process_all_fil.sh, it already makes folder for you.
    #now copy all the files back
    #if [ ! -d $3 ]; then
    	#FN=$(echo "$3" | cut -f 1 -d '.')
	#mkdir $FN
    #fi
    #cp -r ${SLURM_TMPDIR}/* $FN
    cp -r ${SLURM_TMPDIR}/* .  # KC commented out while reset SLURM_TMPDIR
    #clean up
    #rm -r ${SLURM_TMPDIR}  # KC commented out while reset SLURM_TMPDIR
fi
