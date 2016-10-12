#!/bin/sh
njobs=106
for i in $( seq 0 $njobs );do
    export I=$i
    if [[ ${1} == "EP" ]]; then
        sbatch -o job$i.out -J 200${1}job$i jobsub.slurm
    elif [[ ${1} == "Ridge" ]]; then
        sbatch -o job$i.out -J 200${1}job$i jobsubRidge.slurm
    elif [[ ${1} == "Perf" ]]; then
        sbatch -o job$i.out -J 200${1}job$i jobsubPerf.slurm
    else
        echo "Wrong parameter!"
    fi
done
