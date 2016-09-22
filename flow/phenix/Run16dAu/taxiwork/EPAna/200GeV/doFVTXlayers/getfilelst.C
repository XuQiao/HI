#include "RpPar.h"
void getfilelst(){
    TString pro="";
    system(Form("dire=`pwd` && 
            echo $dire &&
            cd /gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/work/200GeV/output && ls output_EP*.root > $dire/%s%s.Lst &&
            cd $dire",dataset.Data(),pro.Data(),dataset.Data(),pro.Data()));
}
    

