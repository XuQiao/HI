#!/bin/sh

export I=$1
sbatch --account=cms_stage2 -o job$I.out -J 200checkjob$I jobsub.slurm

#root -b << EOF
#int nhar = 3
#.L Getvn2D.C+
#int ihar = $1 % nhar
#int UseCNTEP = $1 / nhar
##Getvn2D(0,ihar,UseCNTEP)
#EOF
