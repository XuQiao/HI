#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --mem-per-cpu=32G
#SBATCH --time=0-02:00:00

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
cd /home/xuq7/HI/flow/phenix/PP/Perform
#hadd merged_Anappfvtxor.root output/Anappfvtxor_*.root
#hadd merged_AnaMWGppfvtxand.root /store/user/qixu/.phenix/PP/Perform/output/AnaMWGppfvtxand_*.root
#hadd merged_AnaMWGppfvtxor.root /store/user/qixu/.phenix/PP/Perform/output/AnaMWGppfvtxor_*.root
#hadd merged_AnaMWGppmb.root /store/user/qixu/.phenix/PP/Perform/output/AnaMWGppmb_*.root
#hadd merged_AnaMWGppfvtxand_NODROP.root /store/user/qixu/.phenix/PP/Perform/physicsoutput/AnaMWGppfvtxand_*.root
#hadd merged_AnaMWGppfvtxor_NODROP.root /store/user/qixu/.phenix/PP/Perform/physicsoutput/AnaMWGppfvtxor_*.root
hadd merged_AnaMWGppmb_NODROP.root /store/user/qixu/.phenix/PP/Perform/physicsoutput/AnaMWGppmb_*.root
#./haddwithchoice.sh
echo "Done!"
exit 0
