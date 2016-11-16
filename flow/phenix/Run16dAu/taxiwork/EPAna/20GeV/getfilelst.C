#include "RpPar.h"
void getfilelst(){
    TString pro="";
    system(Form("dire=`pwd` && 
            echo $dire &&
            cd /home/xuq7/HI/flow/phenix/Run16dAu/work/20GeV/output && ls EP*.root > $dire/%s%s.Lst &&
            cd $dire",dataset.Data(),pro.Data(),dataset.Data(),pro.Data()));
}
    

