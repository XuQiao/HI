#!/bin/sh

mult=100
for i in {0..300}
do
    export MULT=$mult
    export IFILE=$i
#qsub -v MULT=$mult,IFILE=$i proSTEG.pbs
#sbatch -o /dev/null -J pro0_${i}_${mult} proSTEG0.slurm
#sbatch -o job5_${i}_${mult}.out -J pro_${i}_${mult} proSTEG4.slurm
#sbatch -o jobwithnf_${i}_${mult}.out -J pro_${i}_${mult} proSTEGwithnfv3.slurm
sbatch -o job_${i}_${mult}.out -J pro_${i}_${mult} proSTEG.slurm
done
