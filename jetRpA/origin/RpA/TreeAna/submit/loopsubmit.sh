#!/bin/bash
cd /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/TreeAna/submit

array=(Jet20 Jet40 Jet60 Jet80 Jet100 MB)

for i in `seq 0 5` 
  do
 #   qsub submit${array[i]}.pbs
    ./submit${array[i]}.pbs
done
