#!/bin/bash
dir=/store/user/qixu/.phenix/Run15pp200MuonMinBias/6794/data/
newdir=/store/user/qixu/.phenix/Run15pp200MuonMinBias/6794/new/
cd $dir
for ifile in `ls`;do
    if [[ $ifile == *_0.root ]];then
        irun=${ifile%"_0.root"}
        if [[ ! -f ${irun}_1.root ]];then
            cp $ifile ${newdir}/${irun}.root
            echo cp $ifile ${newdir}/${irun}.root
        else
            hadd ${newdir}/${irun}.root ${irun}*.root
            echo hadd ${newdir}/${irun}.root ${irun}*.root
        fi
    fi
done
