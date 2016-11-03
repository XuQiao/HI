#!/bin/sh
njobs=79
for i in $( seq 0 $njobs );do
#for i in `echo 20 33 50 61 66 71`;do
    export I=$i
    if [[ ${1} == "EP" ]]; then
        sbatch -o job$i.out -J 20${1}job$i jobsub.slurm
    elif [[ ${1} == "Ridge" ]]; then
        sbatch -o job$i.out -J 20${1}job$i jobsubRidge.slurm
    elif [[ ${1} == "Perf" ]]; then
        sbatch -o job$i.out -J 20${1}job$i jobsubPerf.slurm
    else
        echo "Wrong parameter!"
    fi
done
