#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST

#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc491

#cd /home/xuq7/HI/CMSSW_5_3_20/src
cd /home/xuq7/HI/CMSSW_7_5_0_pre2/src
#cmsenv
eval `scramv1 runtime -sh`
cd /home/xuq7/HI/centrality/NBD/PbPb/NcollExp/double_side


root -l -b <<EOF
.L /home/xuq7/HI/centrality/NBD/PbPb/NcollExp/NBDclassNcollExp.C+
.L store.C
store(0,$STH,$GTH)
EOF
