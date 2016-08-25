#!/bin/sh
rm inputcheck.txt

sAsyjet="/afs/cern.ch/user/d/dgulhan/public/Corrections/Casym_pPb_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root"

echo "Asyjet_input:">>inputcheck.txt
echo $sAsyjet>>inputcheck.txt
echo `edmDumpEventContent $sAsyjet` >>inputcheck.txt

