#!/bin/bash
SumorProd="Prod"
Vorv="v"

nfilesperjob=25
if [[ $Vorv == "V" ]]; then
nfiles=890
else
nfiles=462
fi
njobs=`echo "$nfiles/$nfilesperjob" | bc`
echo "split into $(($njobs+1)) jobs, $nfilesperjob files per job"

echo $SumorProd $Vorv

for dir in `ls`;do
if [[ -d $dir && $dir == M300260 ]];then
echo $dir
cd $dir
rm LYZPF_C*
rm LYZtrack_C*
root -l -b<<EOF
.L LYZPF.C+
EOF
root -l -b<<EOF
.L LYZtrack.C+
EOF
cd ..

for i in $( seq 0 $njobs );do
if [[ $i == 16 ]];then
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
echo -e $begin "to" $end '\t'
qsub -v I=$i,BEGIN=$begin,END=$end,DIR=$dir,SUMORPROD=$SumorProd,VORV=$Vorv -N $dir$i -z jobsub.pbs
fi
done

fi
done
