# CHIME-Pulsar_automated_filterbank
To run this on a set of filterbank files

On compute canada the command is 

`.process_all_fil $split_size $DM *.fil`

wait for the submitted jobs to finish
then to create the candidates and grade

`sbatch automated_filterbank_FETCH.batch .`

The split size should be determined by the size of the filterbank file you're using, a general rule of thumb is to try get the split filterbank file to ~800MB for every 40GB of RAM you use. 
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


