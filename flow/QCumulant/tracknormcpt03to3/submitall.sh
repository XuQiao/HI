#!/bin/bash
nfilesperjob=25
nfiles=462
njobs=`echo "$nfiles/$nfilesperjob" | bc`
echo "split into $(($njobs+1)) jobs, $nfilesperjob files per job"
arr=("$@")
par1=$1
unset arr[0]
for dir in `ls`;do
if [[ -d $dir && $dir == $par1 ]];then
echo $dir
cd $dir
root -l -b<<EOF
gInterpreter->AddIncludePath("/home/xuq7/HI/flow/QCumulant/correlations-0.9_test");
.L Analyze.C+
EOF
cd ..
for i in $( seq 0 $njobs );do
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
export START=$begin
export END=$end
export DIR=$dir
sbatch -J $dir$i -o $DIR/job$I.out jobsub.slurm
#./jobsub.slurm
fi
done
done
fi
done
