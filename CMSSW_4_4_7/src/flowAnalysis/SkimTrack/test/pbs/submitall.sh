#!/bin/bash
nfiles=`cat filelistAll.dat | wc -l`
nfilesperjob=20
njobs=`echo "$nfiles/$nfilesperjob" | bc`
echo "split into $(($njobs+1)) jobs, $nfilesperjob files per job"

for i in $( seq 0 $njobs );do
#if [[ $i == 0 ]];then
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
qsub -v I=$i,BEGIN=$begin,END=$end -N job$i -z jobsub.pbs
#fi
done
