##!/bin/bash
#trig="pAumbst"
trig="pAumbcentral"
#trig="pAufvtxand"
#trig="pAufvtxor"
#trig="pAufvtxsouth"
#trig="pAufvtxnorth"
if [[ $trig == *mb* ]];then
nfiles=1513
fi
if [[ $trig == *fvtx* ]];then
nfiles=339
fi
nfilesperjob=3
njobs=`echo "$nfiles/$nfilesperjob" | bc`
echo "split into $(($njobs+1)) jobs, $nfilesperjob files per job"

arr=("$@")
for i in $( seq 0 $njobs );do
    for j in ${arr[@]};do
if [[ $i == $j ]];then
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
echo -e $begin "to" $end '\t'
export I=$i
export BEGIN=$begin
export END=$end
export TRIG=$trig
sbatch -J  ${trig}out$i -o ${trig}job$i.out jobsub.slurm
fi
done
done
