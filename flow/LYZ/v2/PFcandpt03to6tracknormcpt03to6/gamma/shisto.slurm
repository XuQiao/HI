#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=0-06:00:00

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST

#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/HI/CMSSW_5_3_20/src
#cmsenv
eval `scramv1 runtime -sh`
cd /home/xuq7/HI/flow/LYZ/v2/PFcandpt03to6tracknormcpt03to6/$DIR

root -l -b savehisto.C

echo "Done!"

exit 0
