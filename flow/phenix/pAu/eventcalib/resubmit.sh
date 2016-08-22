#!/bin/bash

arr=()
arrindex=-1
for i in `seq 0 1375`;do
    if [[ ! -f output/AnapAumbst_$i.root ]];then
        j=$((i/10))
        if [[ $arrindex == -1 || ${arr[$arrindex]} != $j ]];then
            arr=(${arr[@]} $j)
            arrindex=$((arrindex+1))
        fi
    fi
done
echo  ${arr[@]}
./submit.sh ${arr[@]}
