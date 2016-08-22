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
if [[ -d $dir && $dir == M220185 ]];then
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
if [[ $i == 0 || $i == 1 || $i == 7 || $i == 8 || $i == 11 || $i == 12 || $i == 14 || $i == 15 ]];then
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
echo -e $begin "to" $end '\t'
#qsub -v I=$i,BEGIN=$begin,END=$end,DIR=$dir,SUMORPROD=$SumorProd,VORV=$Vorv -N $dir$i -z jobsub.pbs
sbatch -J $dir$i jobsub.slurm $i $begin $end $dir $SumorProd $Vorv
fi
done

fi
done
