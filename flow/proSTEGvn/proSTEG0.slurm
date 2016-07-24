#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --reservation=xuq7_1

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

TMPDIR=/tmp/${USER}-${SLURM_JOB_ID}
if [[ ! -d $TMPDIR ]]; then
mkdir $TMPDIR
fi
cd /home/xuq7/HI/flow/proSTEGvn

root -l -b "proSTEGvn0.C(\"${TMPDIR}/vndata_50k_${IFILE}.root\")"
echo "Done!"
if [[ ! -d /cms/store/user/qixu/flow/STEG/pPbDataV200m$MULT ]]; then
    mkdir /cms/store/user/qixu/flow/STEG/pPbDataV200m$MULT
fi
export X509_USER_CERT=/home/xuq7/.globus/usercert.pem
/usr/local/cms-stageout/vandyCp.sh $TMPDIR/vndata_50k_$IFILE.root /cms/store/user/qixu/flow/STEG/pPbDataV200m$MULT/vndata_50k_$IFILE.root
rm -r $TMPDIR/vndata_50k_$IFILE.root
if [[ "$(ls -A $TMPDIR)" ]];then
    echo "NOT EMPTY!"
else
    echo "EMPTY!"
    rm -r $TMPDIR
fi
exit 0
