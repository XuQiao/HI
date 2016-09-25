#!/bin/sh
njobs=73
jobs=(13 34 35 39 40 49 50 56 65 70 2 3 5)
#for i in $( seq 0 $njobs );do
for i in ${jobs[@]};do
    export I=$i
    sbatch -o job$i.out -J 62GeVjob$i jobsub.slurm
done
