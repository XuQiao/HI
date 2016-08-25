#!/bin/bash
#PBS -l nodes=1
#PBS -l mem=8000mb
#PBS -l pmem=8000mb
#PBS -l walltime=24:00:00
#PBS -o 1.out
#PBS -j oe

cd $PBS_O_WORKDIR
pwd
echo $PBS_NODEFILE
echo $HOSTNAME

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`

#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/HI/CMSSW_5_3_20/src
#cmsenv
eval `scramv1 runtime -sh`

cd /home/xuq7/HI/flow/proSTEGvn

root -l -b proSTEGvn.C
echo "Done!"

exit 0
