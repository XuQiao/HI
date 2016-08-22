#!/bin/sh
nfiles=427
nfileperjob=30
njobs=`echo "$nfiles/$nfileperjob+1" | bc`
echo "split into $njobs jobs, $nfileperjob files per job"
for i in `seq 2 $njobs`
do
#if [[ $i == 1 ]];then
echo $i
start=`echo "($i-1)*$nfileperjob" | bc`
end=`echo "$i*$nfileperjob" | bc`
if [[ $i == $njobs ]];then 
    end=$nfiles
fi
export I=$i
export START=$start
export END=$end
#qsub -v I=$i,START=$start,END=$end -N jobsub$i -z jobsub.sh
sbatch -J PbPb$i -o job$i.out jobsub.slurm
#fi
done
