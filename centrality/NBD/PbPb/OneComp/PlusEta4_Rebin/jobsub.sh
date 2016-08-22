#!/bin/bash
cd /home/xuq7/CMSSW_6_2_3_patch1/src/
eval `scramv1 runtime -sh`
cd /home/xuq7/CMSSW_6_2_3_patch1/src/Centrality/NBD/PlusEta4_Rebin
root -l -b <<EOF
.x store.C
EOF
