#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=0-06:00:00

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`

#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/HI/CMSSW_5_3_20/src/
#cd /home/xuq7/HI/CMSSW_7_4_0_pre6/src/
#cmsenv
eval `scramv1 runtime -sh`

cd /home/xuq7/HI/jetRpA/RpA/NewMC/produceandcheck
if [[ ! -d /tmp/xuq7 ]];then
mkdir /tmp/xuq7
fi
#  echo "First = $FIRST and last file = $LAST"   
#.L analyzeTreesPA.C+
#analyzeTreesPA(1, 0, 0, 1, 1, 1, 0, -1, false, 1)
root -b > ../submit/AnaTree.log <<EOF
.L skimTreesPAsimple.C+
skimTreesPAsimple(1)
.q
EOF

echo "Done all jobs!"

#echo "Copying output files to " $destination
