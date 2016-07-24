#!/bin/bash
HLTTrig="M300260"
#dataset="Pbp"
#dataset="pPbReReco"
dataset="pPb"
trigname=${HLTTrig#$"M"}

nfiles=`cat filelist${dataset}.dat | wc -l`
nfilesperjob=10
njobs=`echo "$nfiles/$nfilesperjob" | bc`
echo "split into $(($njobs+1)) jobs, $nfilesperjob files per job"

echo "Now skim" $HLTTrig for data $dataset

#for i in $( seq 0 $njobs );do
for i in {496,499,502,503,504,505,507,510,512,513,516,520,524,525,526,527,528,529,530,532,533,534,535,536,537,538,539,540,541,543,544,545,546,547,548,549,554};do
#if [[ $i == 7 ]];then
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
qsub -v I=$i,BEGIN=$begin,END=$end,DATASET=$dataset,TRIGNAME=$trigname -N ${HLTTrig}${dataset}$i -z jobsub.pbs
#fi
done
