#!/bin/bash
SumorProd="Prod"
Vorv="v"

nfilesperjob=25
nfiles=462
njobs=`echo "$nfiles/$nfilesperjob" | bc`
echo "split into $(($njobs+1)) jobs, $nfilesperjob files per job"

echo $SumorProd $Vorv
arr=("$@")
par1=$1
par2=$2
echo loop = $par1
unset arr[0]
unset arr[1]
for dir in `ls`;do
if [[ -d $dir && $dir == $par2 ]];then
echo $dir
cd $dir
rm LYZ_C*
root -l -b<<EOF
.L LYZ.C+
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
export BEGIN=$begin
export END=$end
export DIR=$dir
export SUMORPROD=$SumorProd
export VORV=$Vorv
sbatch -J $par1$dir$i -o $DIR/job$I$SUMORPROD.out jobsub$par1.slurm
#sbatch -J 0$dir$i -o $DIR/job$I$SUMORPROD.out jobsub.slurm
#sbatch -J 2$dir$i -o $DIR/job$I$SUMORPROD.out jobsub2.slurm
#sbatch -J 3$dir$i -o $DIR/job$I$SUMORPROD.out jobsub3.slurm
#sbatch -J 4$dir$i -o $DIR/job$I$SUMORPROD.out jobsub4.slurm
#sbatch -J 5$dir$i -o $DIR/job$I$SUMORPROD.out jobsub5.slurm
fi
done
done
fi
done
