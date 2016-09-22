#!/bin/sh
njobs=79
for i in $( seq 0 $njobs );do
    export I=$i
    sbatch -o job$i.out -J job$i jobsub.slurm
done
