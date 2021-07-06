# Setting up your folder structure
Place all the filterbank files in some folder (default to filterbank directory in project space)
symlink it to your scratch sapce

`ln -s $FULL_PATH_TO_FILTERBANK $DESTINATION`

# CHIME-Pulsar_automated_filterbank
the best way to run is in two parts
`check_single_pulse.sh -b -d $DM`

to run the initial prep for fetch

MAKE SURE YOU HAVE NO MODULES LOADED ON CC. THE SCRIPT WILL LOAD THEM FOR YOU.

`check_single_pulses.sh -f` 

will run FETCH on any unprocessed `.fil` files in your current directory.

This script will only run the pipeline incrementally on `.fil` files, meaning that it will not rerun the same `.fil` file twice. If you want to do that either use the method below or copy your `.fil` file to an empty directory (symbolic link works).

`get_bright_bursts.sh -i .`

will grab all the astrophysical bursts and put them into `positive_bursts/` and create a `positive_burst.csv` file


The method outlined below still works but is deprecated.

To run this on a set of filterbank files

On compute canada the command is 

`.process_all_fil $split_size $DM *.fil`

wait for the submitted jobs to finish
then to create the candidates and grade

MAKE SURE YOU HAVE NO MODULES LOADED ON CC. THE SCRIPT WILL LOAD THEM FOR YOU.

`sbatch automated_filterbank_FETCH.batch .`

The split size should be determined by the size of the filterbank file you're using, a general rule of thumb is to try get the split filterbank file to ~800MB for every 40GB of RAM you use. If you are not sure about this, just use a split size of 1 and then disable the SK portion of the code in `automated_filterbank_batch.sh` 
todo: make the split size an automated feature
i.e. How I've configured it 4GB -> split_size=5

To analyse the results you should see the following folder structure

`$Pulsar_Name/*`

Within this folder is everything you will need to look at

`$Pulsar_Name/{0..split_size-1}/*`

This contains all the `.singlepulse` and rfifind files, along with the candidate images

`$Pulsar_Name/combined_results.csv`

This file contains the output from the FETCH pipeline

So I would expect the workflow to look something like this
1) look at `combined_results.csv`
2) If any candidate looks like it has a astrophysical probability >0.5 then we go look at the plots for it


