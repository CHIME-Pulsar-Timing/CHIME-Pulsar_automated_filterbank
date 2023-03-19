#!/bin/sh

#Run this on the directory where the filterbank file is (softlinked .fil file should work)

#Use filterbank file name without .fil extension as first input
PSRFILE=$1

#The number of trial DM to search (in steps of 0.1 pc/cc) as second input
NUMDMS=$2

#run rfifind on the data
rfifind -blocks 2 -o ${PSRFILE} ${PSRFILE}.fil
wait

#dedisperse the data (assuming data are 15 minutes long)
prepsubband -lodm 0 -numdms ${NUMDMS} -dmstep 0.1 -numout 3072000 -mask ${PSRFILE}_rfifind.mask -o ${PSRFILE} ${PSRFILE}.fil
wait

#run realfft and search for periodicity candidates
for f in ${PSRFILE}*.dat; do
    realfft -disk $f;
    rednoise ${f%.dat}.fft;
    mv ${f%.dat}_red.fft ${f%.dat}.fft;
    accelsearch -zmax 0 -numharm 32 ${f%.dat}.fft;
done
wait

#sift the candidates and put into a file
python /usr/local/src/presto/python/ACCEL_sift.py > ${PSRFILE}_candidates.lis
