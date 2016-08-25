#!/bin/bash
#PBS -M qiao.xu@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=1
#PBS -l mem=4000mb
#PBS -l pmem=4000mb
#PBS -l walltime=12:00:00
#PBS -o submitTrackTrig.output
#PBS -j oe

cd $PBS_O_WORKDIR
pwd
echo $PBS_NODEFILE
echo $HOSTNAME

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`
#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/CMSSW_6_2_3_patch1/src/

#cmsenv
eval `scramv1 runtime -sh`

cd /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA

ptmin=(15 30 50 80 120 170 220 280)
ptmax=(30 50 80 120 170 220 280 9999)
#for i in `seq 1 7`;
for ((i=0; i<2 ; i++)) #${#ptmin[@]}; i++))
     do
        export PTHAT=${ptmin[$i]}
        export PTMAX=${ptmax[$i]}
        echo "i = $i and pthat = $PTHAT"
echo "Processing..."
echo "root -l -b -q anaTrackAsy.C+"

#Do the analysis
#root -b > runjob.log <<EOF
#.L anaInclusiveJS.C+
#anaInclusiveJS()
#.q
#EOF

root -b > submit/runjobTrack.log <<EOF
.L anaTrackAsy.C+
anaTrackAsy()
.q
EOF
     echo "Done! $i"
     done

echo "Done all jobs!"

#echo "Copying output files to " $destination
