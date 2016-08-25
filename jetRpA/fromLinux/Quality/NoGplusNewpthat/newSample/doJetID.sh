#!/bin/sh
export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/xuq7/CMSSW_6_2_3_patch1/src/
#cmsenv
eval `scramv1 runtime -sh`

#rm /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/pic/JetIDcut/merged/*.png

cd /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/NoGplusNewpthat/newSample

for (( i=13;i<14;i++ )) 
#for (( i=0;i<1;i++ )) 

do
export LIST=$i

root -l -b -q ptJetID.C
#root -l -q overlayDataMC.C
#root -l -b -q overlayProfile.C
#root -l -q overlayJetIDcut.C

done
#exit 0

