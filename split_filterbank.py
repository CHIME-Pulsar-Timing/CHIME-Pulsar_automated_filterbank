import sys
from presto import filterbank as fb
from presto.filterbank import FilterbankFile as fbf

def split_fb(segments,filfile,slurm=''):
    '''
    This function splits a filterbank file up for easier processing
    input
    segments: integer number of segments
    fil: fil file

    output
    string array of filenames
    '''
    fil = fbf(filfile,mode='read')
    total_samples = fil.nspec
    segment_length = int(fil.nspec/segments)
    fnames=[]
    for i in range(segments):
        if (i*segment_length+segment_length)>total_samples:
            #prevent us from running off the end
            segment_length = total_samples-i*segment_length
        my_specs = fil.get_spectra(i*segment_length,segment_length)
        if slurm:
            import os
            fname=os.path.join(slurm,filfile.rstrip('.fil')+'_'+str(i)+'.fil')
        
        out_fname=filfile.rstrip('.fil')+'_'+str(i)+'.fil'

        fil.header['tstart'] = fil.header['tstart']+((segment_length*float(fil.header['tsamp']))/(60*60*24))
        fb.create_filterbank_file(fname,fil.header,spectra=my_specs.data.T,nbits=fil.header['nbits'])
        fnames.append(out_fname)
    return fnames

filnames=split_fb(int(sys.argv[1]),sys.argv[2],sys.argv[3])
returnstr=''
for n in filnames:
    returnstr+=n+' '
print(returnstr)
