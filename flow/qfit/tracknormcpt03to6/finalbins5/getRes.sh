#!/bin/bash
Vorv="V"

cd /home/xuq7/HI/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`
cd /home/xuq7/HI/flow/qfit/tracknormcpt03to6/finalbins5

for dir in `ls`;do
if [[ -d $dir && $dir == M* ]];then
cd $dir
export DIR=$dir
if [[ $Vorv == "V" ]];then
root -l <<EOF
.L getResV.C
getResV()
.L plotqfit.C
plotqfit(1)
EOF
#root -l -q getResVsub.C
#root -l -q nsubvsV2.C
fi
if [[ $Vorv == "v" ]];then
root -l <<EOF
.L getResv.C
getResv(1)
EOF
fi
#root -l -q getResvsub.C
if [[ $Vorv == "veta" ]];then
root -l <<EOF
.L getResv.C
getResv(0)
EOF
fi
cd ..
fi
done
