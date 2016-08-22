#!/bin/bash
for i in `seq 0 140`; do
    start=$((i*50))
    end=$(((i+1)*50))
    root -l -b /cms/store/user/qixu/jetRpA/skimTree/DATAPPbakPu3PFskimfile${start}_${end}.root <<EOF
    
done
