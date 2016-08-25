#!/bin/bash
CWD=`pwd`
cd /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/
eval `scramv1 runtime -sh`
cd $CWD
for i in {0..0}
do
	echo "start with $i th"
for j in {0..8}
do
	export STH=$i
	export GTH=$j

root -l -b -q getpar.C	>>run.log <<EOF
EOF
	echo -n ""
done
done
