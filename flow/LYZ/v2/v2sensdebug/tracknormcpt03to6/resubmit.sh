#!/bin/bash

for iloop in 0 1 2 3 4 5; do
    for idir in M150120 M185150 M220185 M260220 M300260 ;do
arr=()
arrindex=-1
for i in `seq 0 461`;do
    if [[ ! -f /scratch/xuq7/flow/pbsjoboutput/v2sensdebug/tracknormcpt03to6/loop${iloop}/${idir}/Anav_Prod2_$i.root ]];then
        j=$((i/25))
        if [[ $arrindex == -1 || ${arr[$arrindex]} != $j ]];then
            arr=(${arr[@]} $j)
            arrindex=$((arrindex+1))
        fi
    fi
done
echo $iloop $idir ${arr[@]}
#./submitall.sh $iloop $idir ${arr[@]}
done
done
