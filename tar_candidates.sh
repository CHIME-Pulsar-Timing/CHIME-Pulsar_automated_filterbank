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
        nsub="${PULSAR}/nsub_128_0"
        nsub_1="${PULSAR}/nsub_128_1"
        nsub_short="${PULSAR}/nsub_128_0_short"
        if [ -d $nsub ]; then
            tar -cf "${PULSAR}/nsub_128_0.tar" $nsub
            if [ $? -eq 0 ]; then
                rm -r $nsub
                chown -R $1 "${PULSAR}/nsub_128_0.tar"
            else
                echo "tar failed for ${PULSAR}"
            fi
        fi
        if [ -d $nsub_1 ]; then
            tar -cf "${PULSAR}/nsub_128_1.tar" $nsub_1
            if [ $? -eq 0 ]; then
                rm -r $nsub_1
                chown -R $1 "${PULSAR}/nsub_128_1.tar"
            else
                echo "tar failed for ${PULSAR}"
            fi
        fi
        if [ -d $nsub_short ]; then
            tar -cf "${PULSAR}/nsub_128_short.tar" $nsub_short
            if [ $? -eq 0 ]; then
                rm -r $nsub_short
                chown -R $1 "${PULSAR}/nsub_128_short.tar"
            else
                echo "tar failed for ${PULSAR}"
            fi
        fi
    fi
done
