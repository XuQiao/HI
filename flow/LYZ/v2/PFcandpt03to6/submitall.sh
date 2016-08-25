#!/bin/bash
SumorProd="Prod"
Vorv="V"

nfilesperjob=40
nfiles=890	#1480
njobs=`echo "$nfiles/$nfilesperjob" | bc`
echo "split into $(($njobs+1)) jobs, $nfilesperjob files per job"

echo $SumorProd $Vorv

for dir in `ls`;do
if [[ -d $dir && $dir == M* &&$dir != M150120 ]];then
echo $dir

for i in $( seq 0 $njobs );do
#if [[ $i == 7 || $i == 8 ]];then
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
echo -e $begin "to" $end '\t'
qsub -v I=$i,BEGIN=$begin,END=$end,DIR=$dir,SUMORPROD=$SumorProd,VORV=$Vorv -N $dir$i -z jobsub.sh
#fi
done

fi
done
