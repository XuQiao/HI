#!/bin/sh
nfiles=427
nfileperjob=20
njobs=`echo "$nfiles/$nfileperjob" | bc`
echo "split into $njobs jobs, $nfileperjob files per job"
for i in `seq 0 $njobs`
do
#if [[ $i != 27 ]];then
echo $i
start=`echo "$i*$nfileperjob" | bc`
end=`echo "($i+1)*$nfileperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
export I=$i
export START=$start
export END=$end
#qsub -v I=$i,START=$start,END=$end -N jobsub$i -z jobsub.sh
sbatch -J PbPbMC$i -o job$i.out jobsub.slurm
#fi
done
