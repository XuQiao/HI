#!/bin/bash
#PBS -l nodes=1
#PBS -l mem=4000mb
#PBS -l pmem=4000mb
#PBS -l walltime=5:00:00
#PBS -o out$STH$GTH
#PBS -j oe

cd /home/xuq7/HI/CMSSW_5_3_20/src/
eval `scramv1 runtime -sh`
cd /home/xuq7/HI/centrality/Bias
root -l -b <<EOF
.x store.C
EOF
