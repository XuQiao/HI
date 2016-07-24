#!/bin/sh
nfiles=300
nfileperjob=15
njobs=`echo "$nfiles/$nfileperjob" | bc`
echo "split into $njobs jobs, $nfileperjob files per job"
for i in $( seq 1 $njobs )
do
#if [[ $i -ge 11 ]];then
echo $i
start=`echo "($i-1)*$nfileperjob" | bc`
end=`echo "$i*$nfileperjob" | bc`
qsub -v START=$start,END=$end jobsub.sh
#fi
done
