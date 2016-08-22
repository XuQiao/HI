#!/bin/sh
#cd /home/xuq7/HI/jetRpA/RpA/NewMC/JetID

for (( i=0;i<22;i++ )) 
do
export LIST=$i

root -l -b -q ptJetIDOnlyPPb_Etabin.C
#root -l -b -q ptJetIDOnlyPPb.C
#root -l -b -q overlayDataMC.C
#root -l -b -q overlayProfile.C
#root -l -q overlayJetIDcut.C

done
#exit 0

