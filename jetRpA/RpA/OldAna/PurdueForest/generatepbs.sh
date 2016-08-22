trigarr=(MB Jet20 Jet40 Jet60 Jet80 Jet100)
for trigname in ${trigarr[*]}
do
cd /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest
cd $trigname
rm -f script/*
cd script
#cd /home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/script
for i in `seq 0 1030`
     do	
	filename="$trigname"_"$i"
	cp ../../submit.pbs ../temp.pbs
	sed -i "s/.output/..\/output\/$filename.output/g" ../temp.pbs
	cp  ../temp.pbs $filename.pbs
	rm ../temp.pbs
       echo "i="$i >> $filename.pbs
       echo 'trigname="'$trigname'"' >> $filename.pbs
       echo "rm -f ~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/$trigname/Log/$filename.log" >> $filename.pbs
       echo 'for j in `seq 0 9`' >> $filename.pbs
       echo "do" >> $filename.pbs
       echo 'ifile=$(( $i*10+$j ))' >> $filename.pbs
       echo 'export NFILE=$ifile' >> $filename.pbs
       echo 'export TRIG=$trigname' >> $filename.pbs
       echo "root -b >>~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/PurdueForest/$trigname/Log/$filename.log<<EOF" >> $filename.pbs
       echo ".L anaDijetMassRpA_test.C+">>$filename.pbs
       echo "anaDijetMassRpA_test()">>$filename.pbs
       echo ".q">>$filename.pbs
       echo "EOF" >>$filename.pbs 
       echo "done" >> $filename.pbs
       echo "rm /tmp/x509up_u171345" >> $filename.pbs
       echo "exit 0" >> $filename.pbs
done
done 
