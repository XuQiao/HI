#!/bin/bash

arr=()
arrindex=-1
for i in `seq 1377 1512`;do
#for i in `seq 0 338`;do
    if [[ ! -f /scratch/xuq7/phenix/pAu/Ridgecntbbc/AnapAumbcentral_$i.root ]];then
        j=$((i/3))
        if [[ $arrindex == -1 || ${arr[$arrindex]} != $j ]];then
            arr=(${arr[@]} $j)
            arrindex=$((arrindex+1))
        fi
    fi
done
echo  ${arr[@]}
#./submit.sh ${arr[@]}
