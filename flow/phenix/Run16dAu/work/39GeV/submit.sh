#!/bin/sh
njobs=108
for i in $( seq 0 $njobs );do
    export I=$i
    sbatch -o job$i.out -J 39GeVjob$i jobsub.slurm
done
