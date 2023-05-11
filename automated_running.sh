#!/bin/bash
while true; do
    input="pulsar"
    filterbank="test/"
    AFP="$(dirname $(readlink -f $0))"
    while IFS= read -r my_pulsar
    do   
        echo "Process $my_pulsar"
        #check if that pulsar has a folder
        ln -s $filterbank/"$my_pulsar"* .
        if [ ! -d $my_pulsar]; then
            mkdir $my_pulsar
        fi
        $AFP/check_single_pulse.sh -b -i "$my_pulsar"/
    done < "$input"
done
