# heretofore known as CHIme/Pulsar Single-pulse PIPEline (CHIPSPIPE)
Video tutorials of how to run CHIPSPIPE 

Introduction and on Cedar/SLURM based super computers: https://youtu.be/pZlYtUMod4A

On local machines: https://youtu.be/qDWkxYvxOns

Analysing the output: https://youtu.be/tbM_24_3Tbo


# Preparing the data
As of March 2022 we found that there are sometimes dropped packets in filterbank data (see https://bao.chimenet.ca/doc/documents/1624)
Therefore I have developed a script to correct the dropped packets. _This is reccommended for all observations untill it is fixed up stream_

`fdp_submit_jobs.sh *.fil` will submit a job for each filterbank file. The script expects you to give it all the filterbank files you're interesting in. For example I will do the following to correct observations of J0012+54 (if you're using high time resolution data, you should try to use fdp_htr.py instead)


`fdp_submit_jobs.sh J0012+54*.fil` 

# Setting up your folder structure
Place all the filterbank files in some folder (default to filterbank directory in project space)
symlink it to your scratch sapce

`ln -s $FULL_PATH_TO_FILTERBANK $DESTINATION`

# CHIME-Pulsar_automated_filterbank
*make sure your have adjusted your channel mask in pipeline config.py*

the best way to run is to run the following command, it will submit two jobs per `observation.fil`, the first one runs the presto pipeline, the second will run the FETCH pipeline, which is dependent on the success of the first job.

`check_single_pulse.sh -b -d $DM $Filterbank/Fits_files`

to run locally use the -l flag (on your own computer i.e. not on a SLURM controlled computer)

`check_single_pulse.sh -l -b -d $DM $Filterbank/Fits_files`

If for some reason FETCH failed you can run it again with the following command (`check_single_pulse.sh` will check if FETCH ran)
MAKE SURE YOU HAVE NO MODULES LOADED ON CC. THE SCRIPT WILL LOAD THEM FOR YOU.


`check_single_pulses.sh -f *.fil` 

similarly

`check_single_pulses.sh -l -f *.fil` to run locally

will run FETCH on any unprocessed `.fil` files in your current directory.

This script will only run the pipeline incrementally on `.fil` files, meaning that it will not rerun the same `.fil` file twice. If you want to do that copy your `.fil` file to an empty directory (symbolic link works).

`get_bright_bursts.sh -i .`

will grab all the astrophysical bursts and put them into `positive_bursts/` and create a `positive_burst.csv` file
