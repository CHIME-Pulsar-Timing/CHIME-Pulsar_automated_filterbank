#!/bin/bash
#if we are not splitting then just set out as the fil file
#argument:
#1) segments
#2) fil name
#3) DM
if [ $1 -gt 1 ]
then
    out=`python $(dirname "$BASH_SOURCE")/split_filterbank.py $1 $2`
else
    out=$2
fi
i=0
echo $out
for fil in $out;
do
    echo $fil
    mkdir $i
    mv $fil $i
    cd $i
    #run pipeline and prep_fetch prep spegID
    python $(dirname "$BASH_SOURCE")/gwg_cand_search_pipeline.py --dm $3 --speg --fetch --no_fft --rfifind --sk_mad --dedisp --sp --fil $fil
    mkdir data
    #remove all the large files that we no longer need
    rm *.dat
    if [ $1 -gt 1 ]
    then
        rm $fil
    fi
    #run FETCH
    #the following code is only valid for Adam's personal computer
    #source ~/anaconda3/etc/profile.d/conda.sh    
    #conda activate fetch
    candmaker.py --frequency_size 256 --time_size 256 --cand_param_file cands.csv --plot --fout data/
    predict.py --data_dir data/ --model a
    #conda deactivate
    cd ..
    ((i=i+1))
done
