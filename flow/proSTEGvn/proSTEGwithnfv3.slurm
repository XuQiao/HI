#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --account=cms_stage2

echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NODELIST

#source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`

#source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/HI/CMSSW_5_3_20/src
#cmsenv
eval `scramv1 runtime -sh`

TMPDIR=/tmp/xuq7/pPbDataV205m$MULT
if [[ ! -d $TMPDIR ]]; then
mkdir /tmp/xuq7
mkdir $TMPDIR
fi
cd /home/xuq7/HI/flow/proSTEGvn

root -l -b proSTEGvnwithnfv3.C
echo "Done!"
if [[ ! -d /cms/store/user/qixu/flow/STEGwithnf/pPbDataV205m$MULT ]]; then
    mkdir /cms/store/user/qixu/flow/STEGwithnf/pPbDataV205m$MULT
fi
if [[ ! -d /cms/store/user/qixu/flow/STEGwithnf/pPbDataV205m$MULT/Varmultper10EP10 ]]; then
    mkdir /cms/store/user/qixu/flow/STEGwithnf/pPbDataV205m$MULT/Varmultper10EP10
fi
export X509_USER_CERT=/home/xuq7/.globus/usercert.pem
/usr/local/cms-stageout/vandyCp.sh $TMPDIR/vndata_50k_$IFILE.root /cms/store/user/qixu/flow/STEGwithnf/pPbDataV205m$MULT/Varmultper10EP10
rm -r $TMPDIR/vndata_50k_$IFILE.root
if [[ "$(ls -A $TMPDIR)" ]];then
    echo "NOT EMPTY!"
else
    echo "EMPTY!"
    rm -rf $TMPDIR
fi
exit 0
