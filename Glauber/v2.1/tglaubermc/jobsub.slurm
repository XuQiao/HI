#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --mem-per-cpu=4G
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
cd /home/xuq7/HI/Glauber/v2.1/tglaubermc

root -l -b <<EOF
gSystem->Load("libMathMore")
.L runglauber_v2.1_test.C+
runAndSaveNtuple(1000000,"Au","Au",42,-1,0.4,"glau_AuAu_ntuple_1M.root")
EOF
echo "Done!"

exit 0
