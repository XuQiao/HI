#!/bin/sh

for dir in `ls`
do
if [[ -d $dir && $dir == M* ]];then
#qsub -v DIR=$dir -N save$dir -z shisto.pbs
export DIR=$dir
sbatch -J save$dir -o $dir/save.out shisto.slurm
fi
done
