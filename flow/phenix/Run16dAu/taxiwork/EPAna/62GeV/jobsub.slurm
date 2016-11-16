#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-01:00:00

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST

#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/HI/flow/phenix/Run16dAu/taxiwork/EPAna/62GeV

root -b <<EOF 
int i=atoi(getenv("I"));
.L checkflat.C
.L checkrec.C
checkflat(0,i)
checkrec(0,i)
EOF
echo "Done!"

exit 0
