#!/bin/sh


k0=8;
Nfiles=`cat filelistAll.dat | wc -l`
echo "Totally $Nfiles submit"
Njobs=`echo "$Nfiles/$k0" | bc`
echo "split into $(($Njobs+1)) jobs, $k0 files per job"
dir="/home/xuq7/tmp/flow/CMSSW_5_3_8_HI_patch2/src/Analysis/LYZsecond/test/pbs/"

cd $dir/jobs/Prod
############
rm *.py 
rm *.pbs 
rm *.output 
rm *.dat
############

for i in $( seq 0 $Njobs )
#for i in $( seq 28 28)
do
	temp1=`expr $i \* ${k0} + 1`
	temp2=`expr \( $i + 1 \) \* ${k0}`
	sed -n "$temp1,${temp2}p" $dir/filelistAll.dat > filelist_$i.dat
	cp $dir/lyzsecond_cfi.py lyzsecond_cfi_$i.py
	cp $dir/jobsub.pbs jobsub_$i.pbs
	
	sed -i "s/filelistAll.dat/filelist_$i.dat/g" lyzsecond_cfi_$i.py
	sed -i "s/pPbAnaTree.root/pPbAnaTree_$i.root/g" lyzsecond_cfi_$i.py
	sed -i "s/lyzsecond_cfi.py/lyzsecond_cfi_$i.py/g" jobsub_$i.pbs
	sed -i "s/jobsub.output/jobsub_$i.output/g" jobsub_$i.pbs
	qsub jobsub_$i.pbs
done
