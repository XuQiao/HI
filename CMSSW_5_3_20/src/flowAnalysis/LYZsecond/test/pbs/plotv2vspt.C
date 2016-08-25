#include "par.h"
void plotv2vspt(){

	V2vsPt->SetParameters(4.81159,1.80783,3.69272,3.11889,0.931485);        //Real data V~0.05
int ibin=0;
TFile *f = TFile::Open("mergedv_Prod.root");
TVectorD *vecDpt = (TVectorD*)f->Get(Form("D_%d/avgpt",ibin));
TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
double *pt = vecDpt->GetMatrixArray();
double *v2 = vecDv2->GetMatrixArray();
double *v2err = vecDv2err->GetMatrixArray();
TGraphErrors *gr=new TGraphErrors(npt,pt,v2,0,0);
gr->SetTitle("v_{2} vs momentum");
gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr->GetYaxis()->SetTitle("v_{2}");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);
V2vsPt->SetLineColor(2);
V2vsPt->Draw();
gr->Draw("Psame");

f->Close();

}

