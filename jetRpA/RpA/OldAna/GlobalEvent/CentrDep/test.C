#include "TF1.h"
void test(){
TF1 *fitEff_HFEnergyPlus4= new TF1("fitEff_HFEnergyPlus4","[0]-[1]*exp([2]*x)-[3]*exp(-[4]*x)",0,300);
fitEff_HFEnergyPlus4->SetParameters(1.,0.367894,-3.36629,0.104605,0.690272);
cout<<fitEff_HFEnergyPlus4->Eval(350)<<endl;
}

