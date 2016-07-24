#!/bin/bash
source /gpfs21/grid/grid-app/cmssoft/cms/cmsset_default.sh
eval `scramv1 runtime -sh`
arr=(Jet20 Jet40 Jet60 Jet80 Jet100)
for var in "${arr[@]}"
do
  export TRIG=${var}
root -l -b <<EOF
.L anaEv.C+
anaEv()
EOF
done
