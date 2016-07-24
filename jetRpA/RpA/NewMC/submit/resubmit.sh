#!/bin/bash

arr=()
arrindex=-1
for i in `seq 0 360`;do
    if [[ `grep Break AnaTree$i.out || grep error AnaTree$i.out || grep Error AnaTree$i.out` ]];then
        if [[ $arrindex == -1 || ${arr[$arrindex]} != $j ]];then
            arr=(${arr[@]} $i)
            arrindex=$((arrindex+1))
        fi
    fi
done
echo  ${arr[@]}
echo $((arrindex+1))
#./loopbatchsubAnaTree.sh ${arr[@]}
