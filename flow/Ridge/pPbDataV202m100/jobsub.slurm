#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --mem-per-cpu=32G
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
cd /home/xuq7/HI/flow/Ridge/pPbDataV202m100

root -l -b <<EOF
.x dofirst_v3.C
EOF

echo "Done!"

exit 0
