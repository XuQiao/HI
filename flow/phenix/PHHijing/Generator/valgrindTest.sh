#!/bin/sh

# Here is a wrapper script for doing the valgrind test in preparation for the 
# train running.  That is why it is specific to CNT/PWG running.

#run=${1:-119919}
nevents=${2:-10}
logfile=valgrindTest.log
#rootfile=valgrind_nDST_$run.root

# This just makes sure that the pi04all lib will be loaded.  This should
# NOT be here -- it couples the script to my area!!!
#export LD_LIBRARY_PATH=/phenix/u/workarea/winter/install-pi04all/lib:$LD_LIBRARY_PATH

exec 4>$logfile

echo Running valgrind test >&4 
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH >&4

valgrind -v --num-callers=20 --tool=memcheck --leak-check=full --error-limit=no --log-file=valgrind.log --suppressions=$ROOTSYS/root.supp --leak-resolution=high root.exe -b >&4 2>&1 <<EOF
.L phhijing.C
phhijing($nevents)
EOF

