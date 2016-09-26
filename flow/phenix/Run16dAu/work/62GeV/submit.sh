#!/bin/sh
njobs=73
for i in $( seq 0 $njobs );do
    export I=$i
    if [[ ${1} == "EP" ]]; then
        sbatch -o job$i.out -J 62${1}job$i jobsub.slurm
    elif [[ ${1} == "Ridge" ]]; then
        sbatch -o job$i.out -J 62${1}job$i jobsubRidge.slurm
    else
        echo "Wrong parameter!"
    fi
done
