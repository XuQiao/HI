#!/bin/sh
SumorProd="Prod"
Vorv="V"

nfiles=300
nfileperjob=20
njobs=`echo "$nfiles/$nfileperjob" | bc`
echo "split into $njobs jobs, $nfileperjob files per job"
for i in $( seq 1 $njobs )
do
#if [[ $i == 2 ]];then
start=`echo "($i-1)*$nfileperjob" | bc`
end=`echo "$i*$nfileperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
qsub -v I=$i,START=$start,END=$end,SUMORPROD=$SumorProd,VORV=$Vorv -N jobsub$i -z jobsub.sh
#fi
done
