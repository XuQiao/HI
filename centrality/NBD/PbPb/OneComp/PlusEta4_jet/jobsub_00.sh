#!/bin/bash
cd /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/PlusEta4_jet
export GTH=0
export STH=0
root -l -b <<EOF
.x store.C
EOF
