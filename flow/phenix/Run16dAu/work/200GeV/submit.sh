#!/bin/sh
njobs=106
for i in $( seq 0 $njobs );do
    export I=$i
    sbatch -o job$i.out -J 200GeVjob$i jobsub.slurm
done
