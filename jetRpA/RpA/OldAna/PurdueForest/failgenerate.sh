#!/bin/bash
cd /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/
trigname="MB"
cd $trigname/script
for pbsname in `ls`
do
if [[ -f $pbsname && $pbsname =~ "fail" ]]; then
rm $pbsname
fi
done
sed -i '/a/d' "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/FailFile$trigname.txt"
while read fail
do
   fail=`echo $fail | sed -e 's/> //g'`
   filename="$trigname"_fail"$fail"
echo $filename
   cp ../../submit.pbs ../temp.pbs
   sed -i "s/.output/..\/output\/$filename.output/g" ../temp.pbs
   cp  ../temp.pbs $filename.pbs
   rm ../temp.pbs
       echo "fail="$fail >> $filename.pbs
       echo 'trigname="'$trigname'"' >> $filename.pbs
       echo "rm -f ~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/$trigname/Log/$filename.log" >> $filename.pbs
       echo 'export NFILE=$fail' >> $filename.pbs
       echo 'export TRIG=$trigname' >> $filename.pbs
       echo "root -b >>~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/$trigname/Log/$filename.log<<EOF" >> $filename.pbs
       echo ".L anaDijetMassRpA_test.C+">>$filename.pbs
       echo "anaDijetMassRpA_test()">>$filename.pbs
       echo ".q">>$filename.pbs
       echo "EOF" >>$filename.pbs
       echo "rm /tmp/x509up_u171345" >> $filename.pbs
       echo "exit 0" >> $filename.pbs
done < "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/FailFile$trigname.txt"
