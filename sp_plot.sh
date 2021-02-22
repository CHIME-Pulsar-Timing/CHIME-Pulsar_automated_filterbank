lodm=$1
hidm=$2


rm to_process
for f in *.singlepulse; do dm=`echo $f |awk -F"DM" '{print $2}' | awk -F".s" '{print $1}'`; if (( $(echo "$dm > $lodm && $dm < $hidm" | bc -l) )); then echo $f >> to_process; fi; done

rm too_process
echo `cat to_process` > too_process

single_pulse_search.py -t 6 `cat too_process`

spname=`ls *_singlepulse.ps`
mv $spname ${spname%.ps}_DM${lodm}-${hidm}.ps
