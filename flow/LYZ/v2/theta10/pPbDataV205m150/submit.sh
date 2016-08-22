#!/bin/sh
nfiles=300
nfileperjob=10
njobs=`echo "$nfiles/$nfileperjob" | bc`
echo "split into $njobs jobs, $nfileperjob files per job"
for i in `seq 1 $njobs`
do
#if [[ $i == 11 || $i == 14 || $i == 17 ||  $i == 24 ]];then
#echo $i
start=`echo "($i-1)*$nfileperjob" | bc`
end=`echo "$i*$nfileperjob" | bc`
#qsub -v I=$i,START=$start,END=$end -N jobsub$i -z jobsub.sh
sbatch -J job$i jobsub.slurm $i $start $end
#fi
done
