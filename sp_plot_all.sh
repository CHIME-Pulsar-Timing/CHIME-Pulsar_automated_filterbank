sh /psr_scratch/chiamin/presto_pipeline/sp_plot.sh 0.00 80.00
sh /psr_scratch/chiamin/presto_pipeline/sp_plot.sh 80.00 140.00
sh /psr_scratch/chiamin/presto_pipeline/sp_plot.sh 140.00 300.00
sh /psr_scratch/chiamin/presto_pipeline/sp_plot.sh 300.00 540.00
sh /psr_scratch/chiamin/presto_pipeline/sp_plot.sh 540.00 780.00
sh /psr_scratch/chiamin/presto_pipeline/sp_plot.sh 780.00 1020.00

for f in *single*.ps; do convert $f ${f%.ps}.png; done
