#!/bin/bash

for i in {0..0}
do
	echo "start with $i th"
for j in {0..10}
do
#if [[  $i == 1 || $i == 2 ]]; then
#if [[  $j == 0 ]]; then
	export GTH=$j;	export STH=$i;
        ./jobsub.sh	
	#qsub -v GTH=$j,STH=$i jobsub.sh
	#sbatch -o out$STH$GTH jobsub.slurm
#fi
#fi
done
	echo ""
done
