#!/bin/bash
SumorProd="Prod"
Vorv="V"

cd /home/xuq7/HI/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`
cd /home/xuq7/HI/flow/LYZ/v2/PFcandpt01to10

for dir in `ls`;do
if [[ -d $dir && $dir == M260220 ]];then
cd $dir
export SUMORPROD=$SumorProd
export DIR=$dir
if [[ $Vorv == "V" ]];then
root -l -q getResV.C
else
root -l -q getResv.C
fi
cd ..
fi
done
