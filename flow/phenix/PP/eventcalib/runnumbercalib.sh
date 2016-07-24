#!/bin/bash
#dir1=/store/user/qixu/.phenix/Run15pp200CAMinBias/6700/new/
#dir2=/store/user/qixu/.phenix/Run15pp200MuonMinBias/6794/new/
dir2=/store/user/qixu/.phenix/Run15pp200CAFVTX/6699/new/
dir1=/store/user/qixu/.phenix/Run15pp200MuonFVTX/6793/new/
cd $dir1
for ifile in `ls`;do
    if [[ -f $dir2/$ifile ]];then
        echo "$dir2/$ifile" >> /home/xuq7/HI/flow/phenix/PP/eventcalib/filelistfvtx2.dat
    fi
    #if [[ ! -f $dir2/$ifile ]];then
    #    echo "$dir1/$ifile" >> /home/xuq7/HI/flow/phenix/PP/eventcalib/filelistDIffmb1.dat
    #fi
done

