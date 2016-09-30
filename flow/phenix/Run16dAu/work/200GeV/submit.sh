#!/bin/sh
njobs=106
for i in $( seq 0 $njobs );do
    export I=$i
    if [[ ${1} == "EP" ]]; then
        sbatch -o job$i.out -J 200${1}job$i jobsub.slurm
        #echo 200${1}job$i
    elif [[ ${1} == "Ridge" ]]; then
        sbatch -o job$i.out -J 200${1}job$i jobsubRidge.slurm
    else
        echo "Wrong parameter!"
    fi
done
