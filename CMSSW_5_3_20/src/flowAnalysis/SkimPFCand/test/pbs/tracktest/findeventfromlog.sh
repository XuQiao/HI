#!/bin/bash

totalevents=0
for file in `ls joboutput`;do
if [[ $file == *pPbReReco_*.out ]];then
x=`grep "passed =" joboutput/$file`
y=${x%%"failed"*}
total=${y##* passed =}
echo -n " "$total
totalevents=$((totalevents+total))
fi
done
echo
echo $totalevents
