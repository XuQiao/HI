#!/bin/bash
#PBS -l nodes=1
#PBS -l mem=4000mb
#PBS -l pmem=4000mb
#PBS -l walltime=2:00:00
#PBS -o job$END.out
#PBS -j oe

cd $PBS_O_WORKDIR
pwd
echo $PBS_NODEFILE
echo $HOSTNAME

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`

#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc434

cd /home/xuq7/CMSSW_5_3_8_HI/src
#cmsenv
eval `scramv1 runtime -sh`

cd /home/xuq7/CMSSW_5_3_8_HI/src/UserCode/flow/LYZ/v2/ToyV206m150

#root -l -b dofirst.C
root -l -b dosecond.C
echo "Done!"

exit 0
