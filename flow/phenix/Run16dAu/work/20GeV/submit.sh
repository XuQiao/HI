#!/bin/sh
njobs=79
for i in $( seq 0 $njobs );do
    export I=$i
    sbatch -o job$i.out -J 20GeVjob$i jobsub.slurm
done
