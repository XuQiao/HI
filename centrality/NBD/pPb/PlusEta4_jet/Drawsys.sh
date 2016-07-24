#!/bin/bash

rm -rf graph.root
rm -rf sys.txt

CWD=`pwd`
cd /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/
eval `scramv1 runtime -sh`
cd ${CWD}
root -l -b <<EOF
.x Drawsys.C
EOF
