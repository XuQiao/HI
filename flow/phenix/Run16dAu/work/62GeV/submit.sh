#!/bin/sh
njobs=73
for i in $( seq 0 $njobs );do
    export I=$i
    sbatch -o job$i.out -J job$i jobsub.slurm
done
