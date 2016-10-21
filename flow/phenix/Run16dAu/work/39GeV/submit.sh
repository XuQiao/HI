#!/bin/sh
njobs=108
#for i in $( seq 0 $njobs );do
for i in `echo 16 27 38 47 56 58 64 69 71 73 74 75 76 78 79 88`;do
    export I=$i
    if [[ ${1} == "EP" ]]; then
        sbatch -o job$i.out -J 39${1}job$i jobsub.slurm
    elif [[ ${1} == "Ridge" ]]; then
        sbatch -o job$i.out -J 39${1}job$i jobsubRidge.slurm
    elif [[ ${1} == "Perf" ]]; then
        sbatch -o job$i.out -J 39${1}job$i jobsubPerf.slurm
    else
        echo "Wrong parameter!"
    fi
done
