#!/bin/bash
trigname="MB"
rm -f /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/checkFile$trigname.txt
cd /scratch/xuq7/RpA/$trigname
max=`find . -name "DATAPPb_*.root" | sed -e "s|./DATAPPb_$trigname||g" | sed -e  's/akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PurdueForest_//g' |  sed -e 's/.root//g' | sort -n | tail -n 1`
echo $max
for i in `seq 0 $max`
do
flag=0
for filename in `ls`
do
if [[ -f $filename ]]; then
#for filename in `ls /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/NewForest`
filename=`echo $filename | sed -e "s|DATAPPb_$trigname||g" | sed -e  's/akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PurdueForest_//g' |  sed -e 's/.root//g'`
#echo $filename
if [ $i -eq $filename ]; then
flag=1
break
fi
fi
done

if [ $flag -eq 0 ]; then
echo "NOT found $i!" >> /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/FailFile$trigname.txt
#else
#echo "found $i OK!"
fi
done
