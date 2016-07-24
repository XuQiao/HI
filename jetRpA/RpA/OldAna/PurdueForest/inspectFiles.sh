#!/bin/bash
trigname=MB
rm -f /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/FailFile$trigname.txt
cd /scratch/xuq7/RpA/$trigname
list=$(find . -name 'DATAPPb*.root'| sed -e "s|./DATAPPb_$trigname||g" | sed -e  's/akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PurdueForest_//g' |  sed -e 's/.root//g'| sort -n)
max=$(find . -name 'DATAPPb*.root'| sed -e "s|./DATAPPb_$trigname||g" | sed -e  's/akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PurdueForest_//g' |  sed -e 's/.root//g' | sort -n | tail -n 1)
#diff <(echo $list) <(seq 0 $max) | grep '>' | cut -c 3-
diff <(printf "%s\n" $list) <(printf "%s\n" `seq 0 $max`) >> /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/FailFile$trigname.txt
