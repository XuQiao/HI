#!/bin/sh
SumorProd="Prod"
Vorv="v"
for i in `seq 0 299`;do
if [[ !(-a jobout/Ana${Vorv}_${SumorProd}_$i.root) ]];then
echo -ne $i '\t'

begin=`echo "$i" | bc`
end=`echo "($i+1)" | bc`
export START=$begin
export END=$end

root -l -q dosecond.C

fi
done

