#!/bin/bash

arr=()
arrindex=-1
for i in `seq 0 442`;do
    if [[ ! -f output/AnapAlmbst_$i.root ]];then
        j=$((i/3))
        if [[ $arrindex == -1 || ${arr[$arrindex]} != $j ]];then
            arr=(${arr[@]} $j)
            arrindex=$((arrindex+1))
        fi
    fi
done
echo  ${arr[@]}
./submit.sh ${arr[@]}
