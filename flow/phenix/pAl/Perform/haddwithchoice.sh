#!/bin/bash
lowlumirun=(436763 436764 436765 436778 436781 436782 436783 436935 436936 436939 437029 437075 437076 437110 437111 437112 437440 438160 438162 438163 438263 438264 438404 438405 438421)
#mkdir lowoutput
prefix=AnaPUsub
trig=pAlmbst_
suffix=.root
for ifile in `ls output`;do
        number=${ifile#$prefix}
        number=${number#$trig}
        number=${number%$suffix}
        number=$((number+1))
        filename=`sed "${number}q;d" ../filelistmb1.dat` 
        for irun in ${lowlumirun[@]};do
        if [[ $filename == *${irun}*.root ]];then
            mv output/${ifile} lowoutput/
            echo $filename
        fi
    done
done
