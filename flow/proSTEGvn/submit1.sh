#!/bin/sh

mult=185
for i in {0..299}
do
    export MULT=$mult
    export IFILE=$i
#qsub -v MULT=$mult,IFILE=$i proSTEG.pbs
#sbatch -o /dev/null proSTEG1.slurm
sbatch -o /dev/null -J pro1_${i}_${mult} proSTEG1.slurm
#root -l -b -q proSTEGvn1.C
done
