#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
#include <vector>
void Drawgrid(){

        if (!TClass::GetDict("NBD")) {
                gROOT->ProcessLine(".L /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C+");
        }
	
	int sth=0;
	int Gth=0;
	TFile *fs = TFile::Open(Form("G%d%d.root",sth,Gth),"Read");
        tl = fs->GetListOfKeys();
        key = ((TKey*) tl->At(0));
        NBD *n0 = (NBD*)key->ReadObj();
	cout<<n0->mubest<<"\t"<<n0->kbest<<endl;
	cout<<n0->chis<<"\t"<<n0->Ndf<<endl;
	TGraph2D *gr2;
	TCanvas *c1 = new TCanvas();
//	c1->SetLogz();
	gr2 = n0->gr2;
	int N = gr2->GetN();
	double *z = gr2->GetZ();
	double *x = gr2->GetX();
	double *y = gr2->GetY();
	for(int i=0;i<N;i++){
		z[i]=-z[i];
		y[i]=-y[i];
	}
	gr2 = new TGraph2D(N,x,y,z);
	gr2->Draw("surf1");
}
