# CHIME-Pulsar_automated_filterbank
To run this on a set of filterbank files

On compute canada the command is 

`.process_all_fil $split_size $DM *.fil`

wait for the submitted jobs to finish
then to create the candidates and grade
`sbatch automated_filterbank_FETCH.batch .`

The split size should be determined by the size of the filterbank file you're using, a general rule of thumb is to try get the split filterbank file to ~800MB for every 40GB of RAM you use. 
i.e. How I've configured it 4GB -> split_size=5

