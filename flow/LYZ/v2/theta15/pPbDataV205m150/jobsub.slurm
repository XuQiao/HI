#!/bin/bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-01:00:00
#SBATCH --output /dev/null
args=("$@")
echo ${args[0]} ${args[1]} ${args[2]}
export I=${args[0]}
export START=${args[1]}
export END=${args[2]}

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

cd /home/xuq7/HI/flow/LYZ/v2/theta8/pPbDataV205m150

#root -l -b dofirst.C
root -l -b dosecond.C
#root -l -b dosecondsub.C
echo "Done!"

exit 0
