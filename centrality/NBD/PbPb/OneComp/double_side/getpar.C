#include "../parameter.h"
#include "par.h"
void getpar(){ 
	int sth=0;
        TFile *f = TFile::Open(outG);
        TString dirname;
	if(sth==0){ dirname = "std";}
        else if(sth==1){dirname ="Gri055";}
        else {dirname ="Gri101";}
        TString name;
	for(int iGlau=0; iGlau<nGlau+2; iGlau++){
	if(iGlau==0)
		name = "G0";
	else if(iGlau<nGlau)
		name = Form("Glau_%d",iGlau);
	else
		name = Form("bin_%d",iGlau-nGlau+1);
		TVectorD* mubest = (TVectorD*)f->Get(Form("%s/%s/mubest",dirname.Data(),name.Data()));
		TVectorD* kbest = (TVectorD*)f->Get(Form("%s/%s/kbest",dirname.Data(),name.Data()));
		TVectorD* chis = (TVectorD*)f->Get(Form("%s/%s/chis",dirname.Data(),name.Data()));
	        TVectorD* Ndf = (TVectorD*)f->Get(Form("%s/%s/Ndf",dirname.Data(),name.Data()));
	//	TVectorD* NcollAver = (TVectorD*)f->Get(Form("%s/Glau_%d/NcollAver",dirname.Data(),iGlau));
//		TVectorD* BAver = (TVectorD*)f->Get(Form("%s/Glau_%d/BAver",dirname.Data(),iGlau));
		cout<<(*mubest)[0]<<","<<(*kbest)[0]<<"},{";
		cout<<(*chis)[0]<<","<<(*Ndf)[0]<<endl;
	//	NcollAver->Print();
	//	BAver->Print();
	}
}
