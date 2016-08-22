#!/bin/sh
ptmin=(15 30 50 80 120 170 220 280)
ptmax=(30 50 80 120 170 220 280 9999)

for (( i=0; i<${#ptmin[@]}; i++))
     do
        echo "i = $i and pthat = ${ptmin[$i]} "
     qsub -v PTHAT=${ptmin[$i]},PTMAX=${ptmax[$i]} -N MC$i -z submitJetTrig.pbs
done 
