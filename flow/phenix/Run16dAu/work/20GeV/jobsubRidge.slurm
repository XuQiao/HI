#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1-00:00:00

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST

#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/HI/flow/phenix/Run16dAu/work/20GeV

root -l <<EOF 
.L ../RidgedAuRun16.C+
.L WritehistoRidge.C
Writehisto()
EOF
echo "Done!"

exit 0
