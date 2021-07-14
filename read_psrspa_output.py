#!/usr/bin/env python

import csv
import sys

with open(sys.argv[1],'r') as psrspa_file:
    reader = csv.reader(psrspa_file,delimiter=' ')
    obs = []
    for row in reader:
        obs.append(row[0])
    #get unique observations
    obs=set(obs)
    #print them
    (print(ob) for ob in obs)
