#!/bin/bash
#PBS -l nodes=1
#PBS -l mem=4000mb
#PBS -l pmem=4000mb
#PBS -l walltime=5:00:00
#PBS -o job0.out
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
filename=(
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2918_1_mNZ.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2919_1_4c6.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2920_1_sx9.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2921_1_VXm.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2922_1_vYh.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2923_1_MWR.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2924_1_22d.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2925_1_Sd3.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2926_1_AD6.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2927_1_w7r.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2928_1_UFa.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2930_1_YZU.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2931_1_AIw.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2932_1_50J.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2933_1_MmS.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2934_1_uma.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2935_1_b5b.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2937_1_HoM.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2938_1_Ily.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2939_1_jOk.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2942_1_oBV.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2943_1_ago.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2944_1_1OI.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2945_1_Led.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2946_1_JAb.root'
       '/store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_2947_1_5oM.root'
)

for ifile in $(seq 0 ${#filename[@]})
#for ifile in $(seq 564 564)
do
#xrdcp root://xrootd.unl.edu/$file /lio/lfs/cms/store/user/qixu/Glauber
xrdcp  root://cmsxrootd.fnal.gov/${filename[$ifile]} /lio/lfs/cms/store/user/qixu/flow/PACorrHM/PbPb2011MBData
done

echo "Done!"

exit 0
