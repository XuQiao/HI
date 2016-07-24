#cd /home/nih/CMSSW_5_3_8_HI_patch2/src
#eval `scramv1 runtime -sh`

grid-proxy-init

#cd /home/nih/CMSSW_5_3_8_HI_patch2/src/RpA
for i in `seq 0 1`;
     do 
        export NFILE=$i
     qsub submitNewForest.pbs
done 
