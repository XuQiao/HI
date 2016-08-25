#!/bin/bash
#PBS -l nodes=1
#PBS -l mem=4000mb
#PBS -l pmem=4000mb
#PBS -l walltime=4:00:00
#PBS -o $DIR/job$I$SUMORPROD.out
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
cd /home/xuq7/HI/flow/LYZ/v2/PFcandpt01to10/OnlyCharged/$DIR

#gSystem->AddIncludePath("-D__USE_XOPEN2K8")
root -l -b <<EOF
.L LYZ.C+
.L ../Do.C
string Vorv=getenv("VORV")
if(Vorv=="V") dofirst()
if(Vorv=="v") dosecond()
EOF

echo "Done!"

exit 0
