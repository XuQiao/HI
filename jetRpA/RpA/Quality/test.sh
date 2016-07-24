#!/bin/sh

rm ./pic/test/*
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/CMSSW_6_2_3_patch1/src/
eval `scramv1 runtime -sh`

cd /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality

for ptfile in `ls | grep pteta`; do
root -l -q $ptfile
done
 
exit 0
