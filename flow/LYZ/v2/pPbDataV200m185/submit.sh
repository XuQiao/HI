#!/bin/sh
nfiles=300
nfileperjob=10
njobs=`echo "$nfiles/$nfileperjob" | bc`
echo "split into $njobs jobs, $nfileperjob files per job"
for i in `seq 1 $njobs`
do
#if [[ $i == 5 || $i == 8 || $i == 9 || $i == 12 || $i == 15 || $i == 17 || $i == 18 || $i == 21 || $i == 20  ]];then
echo $i
start=`echo "($i-1)*$nfileperjob" | bc`
end=`echo "$i*$nfileperjob" | bc`
export I=$i
export START=$start
export END=$end
#qsub -v I=$i,START=$start,END=$end -N jobsub$i -z jobsub.sh
sbatch -J ToyMC$i -o job$i.out jobsub.slurm
#fi
done
