#!/bin/bash

rm -r out_job*
rm -r jobsub*.sh

for i in {0..2}
do
	echo "start with $i th"
for j in {0..10}
do
#if [[  $j -ge 9 ]];then
	echo '#!/bin/bash' >> jobsub_$i$j.sh
	echo 'cd /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/' >> jobsub_$i$j.sh
	echo 'eval `scramv1 runtime -sh`' >> jobsub_$i$j.sh
	echo 'cd /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/PlusEta4' >> jobsub_$i$j.sh
	echo "export GTH=$j" >> jobsub_$i$j.sh
	echo "export STH=$i" >> jobsub_$i$j.sh
	echo 'root -l -b <<EOF' >> jobsub_$i$j.sh
	echo '.x store.C' >> jobsub_$i$j.sh
	echo 'EOF' >> jobsub_$i$j.sh
	chmod +x jobsub_$i$j.sh
	echo -n "$j	"
	./jobsub_$i$j.sh
	#bsub -o out_job$i$j -R "pool>30000" -q 8nh -J  jobsub_$i$j < jobsub_$i$j.sh
	rm jobsub_$i$j.sh
#fi
done
	echo ""
done
