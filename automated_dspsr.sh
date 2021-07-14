#!/bin/bash
#
echo $1
basename=$1
DM=$2
destination=fluence_archives
mkdir -p $destination


while IFS=, read -r mjddec obs_length pulse_toa
do
    #find the filterbank and rfi files
    mjd=${mjddec%.*}
    filbase="$basename"_"$mjd"_pow
    filfile="$filbase".fil
    filfolder="$filbase"/0/
    mkdir -p "$destination"/"$mjddec"
    cp $filfile/*rfifind* $destination/$mjddec
    cp -d $filfile $destination/$mjddec
    #this should have set up the folder, now to run dspsr
    #pulse_toa gives the timestamps
    echo entering fil folder
    cd $filfolder
    dspsr -c 4 -N $basename -s -S $pulse_toa -T 8 -k chime -D $DM $filfile
    #run rfistats
    rfifind_stats.py "$filbase"_rfifind.mask
    #run paz to remove rfi
    paz -e zap -z $( python ~/CHIME-Pulsar_automated_filterbank/read_zapchan.py "$basename"_rfifind.zapchans ) *.ar
    #psrspa
    psrspa -a above:threshold=6 -N 127 *.zap
    PSRSPA_OUTPUT=($(python ~/CHIME-Pulsar_automated_filterbank/read_psrspa_output.py psrspa_pulses.dat | tr -d '[],'))
    mkdir -p PSRSPA_PULSES
    for OUT in "${PSRSPA_OUTPUT[@]}"
    do
        mv "$OUT" PSRSPA_PULSES
    done
    cd ..
done < extracted_bursts.csv
