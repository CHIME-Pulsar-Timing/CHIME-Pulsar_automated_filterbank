import numpy as np
import matplotlib.pyplot as plt
import glob
import os
#get all the .ar files
ar_files = glob.glob('*.ar')
masked_channel_arr = []
mjds = []
for my_timer in ar_files:
    command = f"get_bad_channel_list.py $my_timer"
    mjds.append(float(my_timer.split('_')[4]))
    output = os.popen(command).read()
    channels = output.split('_')
    masked_channel_arr.append(channels)
import pdb; pdb.set_trace()
