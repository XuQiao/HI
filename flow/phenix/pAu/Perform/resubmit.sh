#!/bin/bash

arr=()
arrindex=-1
count=0
for i in `seq 0 1512`;do
    if [[ -f /scratch/xuq7/phenix/pAu/Perform/AnapAumbst_$i.root ]]; then
        if test `find "/scratch/xuq7/phenix/pAu/Perform/AnapAumbst_$i.root" -mmin +2400`; then
        j=$((i/1))
        if [[ $arrindex == -1 || ${arr[$arrindex]} != $j ]];then
            arr=(${arr[@]} $j)
            arrindex=$((arrindex+1))
            count=$((count+1))
        fi
    fi
fi
done
echo  ${arr[@]}
echo $count
#./submit.sh ${arr[@]}
