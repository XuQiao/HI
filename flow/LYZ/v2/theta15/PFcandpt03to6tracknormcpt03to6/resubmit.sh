#!/bin/bash

for idir in M150120 M185150 M220185 M260220 M300260 ;do
arr=()
arrindex=-1
for i in `seq 0 889`;do
    if [[ ! -f /scratch/xuq7/flow/pbsjoboutput/theta15/PFcandpt03to6tracknormcpt03to6/${idir}/AnaV_Sum_$i.root ]];then
        j=$((i/25))
        if [[ $arrindex == -1 || ${arr[$arrindex]} != $j ]];then
            arr=(${arr[@]} $j)
            arrindex=$((arrindex+1))
        fi
    fi
done
echo  $idir ${arr[@]}
#./submitall.sh $idir ${arr[@]}
done
