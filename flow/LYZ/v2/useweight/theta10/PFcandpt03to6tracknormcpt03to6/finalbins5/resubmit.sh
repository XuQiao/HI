#!/bin/sh
SumorProd="Prod"
Vorv="V"

echo $SumorProd $Vorv

for dir in `ls`;do
if [[ -d $dir ]];then
echo $dir
for i in `seq 0 461`;do
if [[ !(-a $dir/jobout/Ana${Vorv}_${SumorProd}_$i.root) ]];then
echo -ne $i '\t'

begin=`echo "$i" | bc`
end=`echo "($i+1)" | bc`


#qsub -v I=$i,BEGIN=$begin,END=$end,DIR=$dir,SUMORPROD=$SumorProd,VORV=$Vorv -N re$dir$i -z jobsub.sh

fi
done
echo 

fi
done
cd ..
