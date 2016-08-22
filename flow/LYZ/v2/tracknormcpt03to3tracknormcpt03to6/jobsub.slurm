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
cd /home/xuq7/HI/flow/LYZ/v2/tracknormcpt03to3tracknormcpt03to6/$DIR

root -l -b <<EOF
.L ../Do.C
string Vorv=getenv("VORV")
if(Vorv=="V"){
gROOT->ProcessLine(".L LYZ.C+");
dofirst();
}
if(Vorv=="v"){
gROOT->ProcessLine(".L LYZ.C+");
dosecond();
}
if(Vorv=="veta"){
gROOT->ProcessLine(".L LYZ.C+");
dothird();
}
EOF

echo "Done!"

exit 0
