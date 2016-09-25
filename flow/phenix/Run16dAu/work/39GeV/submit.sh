#!/bin/sh
njobs=108
jobs=(0 101 106)
#for i in $( seq 0 $njobs );do
for i in ${jobs[@]};do
    export I=$i
    sbatch -o job$i.out -J 39GeVjob$i jobsub.slurm
done
