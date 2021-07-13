#!/usr/bin/env bash
set -euo pipefail
#this serves as the automated dspsr program
dspsr -c 0.4 -N J0012+54 -s -S 627.6 -T 0.4 -D 130.4 -k chime J0012+54_59358_pow.fil
#need to read the zapchan file
paz -e zap -z $( python ~/CHIME-Pulsar_automated_filterbank/read_zapchan.py J0012+54_59358_pow_rfifind.zapchans ) start_59358.6671434013635.ar
