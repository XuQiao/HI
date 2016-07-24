#!/bin/sh

rm pPbHist_Epos*.root
rm pPbHist_Hijing*.root
rm pPbHistPlusEta4.root

root -l -q Trig_PlusEta4.C
root -l -q Trig_eff.C
root -l -q Eff_Cor.C
exit 0
