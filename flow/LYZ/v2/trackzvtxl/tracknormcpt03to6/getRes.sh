#!/bin/bash
SumorProd="Prod"
Vorv="vetac"

cd /home/xuq7/HI/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`
cd /home/xuq7/HI/flow/LYZ/v2/trackzvtxl/tracknormcpt03to6/

for dir in `ls`;do
if [[ -d $dir && $dir == M* ]];then
cd $dir
export SUMORPROD=$SumorProd
export DIR=$dir
if [[ $Vorv == "V" ]];then
root -l -q getResV.C
#root -l getResVsub.C
#root -l -q getResVsub.C
#root -l -q nsubvsV2.C
fi
if [[ $Vorv == "v" ]];then
root -l -q "getResvsub.C(1)"
fi
#root -l -q getResvsub.C
if [[ $Vorv == "veta" ]];then
root -l -q "getResvsub.C(0)"
fi
if [[ $Vorv == "vetac" ]];then
root -l -q "getResvsub.C(-1)"
fi
cd ..
fi
done
