#/bin/bash
for i in 30 40 50 65 80 100 120 150 300;do
    dir=pPbDataV202m$i
    cp pPbDataV202m30/merge.C pPbDataV202m$i
    cd $dir
    root -l -q merge.C
    cd ..
done


