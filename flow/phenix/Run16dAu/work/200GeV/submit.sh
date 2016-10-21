#!/bin/sh
njobs=106
#for i in $( seq 0 $njobs );do
for i in `echo 101 104 105 2 21 24 34 40 42 46 50 6 60 64 65 67 7 82 83 87 88 9`; do
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
