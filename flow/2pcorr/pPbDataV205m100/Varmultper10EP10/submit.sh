#!/bin/bash
nfiles=300
nfilesperjob=1
njobs=`echo "$nfiles/$nfilesperjob" | bc`
echo "split into $(($njobs+1)) jobs, $nfilesperjob files per job"

for i in $( seq 0 $njobs );do
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
echo -e $begin "to" $end '\t'
export I=$i
export BEGIN=$begin
export END=$end
sbatch -o job$i.out jobsub.slurm
done
