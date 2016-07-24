#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
void getpar(){ 
	int sth=0;
        TFile *f = TFile::Open(Form("G.root"));
	if(sth==0){TString dirname = "std";}
        else if(sth==1){TString dirname ="Gri055";}
        else {TString dirname ="Gri101";}
	for(int iGlau=0; iGlau<nGlau; iGlau++){
		TVectorD* mubest = (TVectorD*)f->Get(Form("%s/Glau_%d/mubest",dirname.Data(),iGlau));
		TVectorD* kbest = (TVectorD*)f->Get(Form("%s/Glau_%d/kbest",dirname.Data(),iGlau));
		TVectorD* chis = (TVectorD*)f->Get(Form("%s/Glau_%d/chis",dirname.Data(),iGlau));
		TVectorD* Ndf = (TVectorD*)f->Get(Form("%s/Glau_%d/Ndf",dirname.Data(),iGlau));
		TVectorD* NcollAver = (TVectorD*)f->Get(Form("%s/Glau_%d/NcollAver",dirname.Data(),iGlau));
		TVectorD* BAver = (TVectorD*)f->Get(Form("%s/Glau_%d/BAver",dirname.Data(),iGlau));
		cout<<(*mubest)[0]<<","<<(*kbest)[0]<<"},{";
		cout<<(*chis)[0]<<","<<(*Ndf)[0]<<endl;
	//	NcollAver->Print();
	//	BAver->Print();
	}
}
