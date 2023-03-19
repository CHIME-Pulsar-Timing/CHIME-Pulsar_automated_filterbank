while read p; do PER=`cut -d' ' -f1 <<< $p`; DM=`cut -d' ' -f2 <<< $p`; PERIOD=`echo 'scale=4;'${PER}'*0.001' | bc`; prepfold -p ${PERIOD} -dm ${DM} -nosearch -slow $1 -noxwin; done < fold
