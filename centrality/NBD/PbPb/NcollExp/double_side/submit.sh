#!/bin/bash

for i in {0..0}
do
	echo "start with $i th"
for j in {5..9}
do
#if [[  $i == 1 || $i == 2 ]]; then
#if [[  $j == 10 ]]; then
	export GTH=$j;	export STH=$i;
        #./jobsub.sh	
	#qsub -v GTH=$j,STH=$i jobsub.sh
	sbatch -J job$j -o out$STH$GTH jobsub.slurm
#fi
#fi
done
	echo ""
done
