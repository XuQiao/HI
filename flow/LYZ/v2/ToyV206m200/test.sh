#!/bin/sh

for file in `ls jobout`
do
if [[ $file == Anav_Prod2_*.root ]];then
file1=${file/Prod2/Prod}
mv jobout/$file jobout/$file1
fi
done
