#!/bin/bash
SumorProd="Prod"
Vorv="V"

cd /home/xuq7/HI/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`
cd /home/xuq7/HI/flow/LYZ/v2/PFcandpt03to6

for dir in `ls`;do
if [[ -d $dir && $dir == M150120 ]];then
cd $dir
export SUMORPROD=$SumorProd
export DIR=$dir
if [[ $Vorv == "V" ]];then
root -l -q getResVsub.C
else
root -l -q getResv.C
fi
cd ..
fi
done