#!/bin/sh
njobs=73
for i in $( seq 0 $njobs );do
#for i in `echo 10 22 24 27 35 36 39 46 48 55 56 59 6 9`; do
    export I=$i
    if [[ ${1} == "EP" ]]; then
        sbatch -o job$i.out -J 62${1}job$i jobsub.slurm
    elif [[ ${1} == "Ridge" ]]; then
        sbatch -o job$i.out -J 62${1}job$i jobsubRidge.slurm
    elif [[ ${1} == "Perf" ]]; then
        sbatch -o job$i.out -J 62${1}job$i jobsubPerf.slurm
    else
        echo "Wrong parameter!"
    fi
done
