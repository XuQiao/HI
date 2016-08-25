#!/bin/bash
SumorProd="Prod"
Vorv="v"

cd /home/xuq7/HI/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`
cd /home/xuq7/HI/flow/LYZ/v2/tracktc/tracknormcpt03to3tracknormcpt03to6/

for dir in `ls`;do
if [[ -d $dir && $dir == M* ]];then
cd $dir
export SUMORPROD=$SumorProd
export DIR=$dir
if [[ $Vorv == "V" ]];then
root -l -q getResV.C
#root -l -q getResVsub.C
#root -l -q nsubvsV2.C
elif [[ $Vorv == "v" ]];then
#root -l -q getResv.C
root -l -q "getResvsub.C(1)"
elif [[ $Vorv == "veta" ]];then
root -l -q "getResvsub.C(0)"
fi
cd ..
fi
done
