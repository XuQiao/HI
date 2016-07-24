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
arr=("$@")
par1=$1
unset arr[0]

for dir in `ls`;do
if [[ -d $dir && $dir == $par1 ]];then
echo $dir
cd $dir
rm LYZPF_C*
rm LYZtrack_C*
root -l <<EOF
.L LYZPF.C+
EOF
root -l <<EOF
.L LYZtrack.C+
EOF
cd ..
for i in $( seq 0 $njobs );do
#if [[ $i != 0 ]];then
    for j in ${arr[@]};do
if [[ $i == $j ]];then
begin=`echo "$i*$nfilesperjob" | bc`
end=`echo "($i+1)*$nfilesperjob" | bc`
if [[ $i == $njobs ]];then
end=$nfiles
fi
echo -e $begin "to" $end '\t'
#qsub -v I=$i,BEGIN=$begin,END=$end,DIR=$dir,SUMORPROD=$SumorProd,VORV=$Vorv -N $dir$i -z jobsub.pbs
export I=$i
export BEGIN=$begin
export END=$end
export DIR=$dir
export SUMORPROD=$SumorProd
export VORV=$Vorv
sbatch -J $dir$i -o $DIR/job$I$SUMORPROD.out jobsub.slurm
#fi
fi
done
done
fi
done
