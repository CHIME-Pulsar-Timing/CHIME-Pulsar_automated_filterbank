#!/usr/bin/env bash
#THE FIRST ARGUMENT SETS THE PERMISSIONS OF THE TARBALL
set -euo pipefail
FILFILES=*.fil
for FIL in $FILFILES;
do
    #strip the extension
    PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
    if [ -d $PULSAR ]; then
        #lists all the splits in the pulsar's directory
        cd $PULSAR
        nsub="nsub_128_0.tar"
        nsub_1="nsub_128_1.tar"
        nsub_short="nsub_128_short.tar"
        if [ -f $nsub ]; then
            tar -xf $nsub
            if [ $? -eq 0 ]; then
                mv ${PULSAR}/nsub_128_0 .
            else
                echo "untar failed for ${PULSAR}"
            fi
        fi
        if [ -f $nsub_1 ]; then
            tar -xf $nsub_1
            if [ $? -eq 0 ]; then
                mv ${PULSAR}/nsub_128_1 .
            else
                echo "untar failed for ${PULSAR}"
            fi
        fi
        if [ -f $nsub_short ]; then
            tar -xf $nsub_short
            if [ $? -eq 0 ]; then
                mv ${PULSAR}/nsub_128_0_short .
            else
                echo "untar failed for ${PULSAR}"
            fi
        fi
        cd ..
    fi
done
