#include "par.h"
void plotV2vsNtrkvs(){

TFile *f = TFile::Open("mergedV_Sum.root");
TVectorD *vecDV2 = (TVectorD*)f->Get(Form("Vmean"));
TVectorD *vecDV2err = (TVectorD*)f->Get(Form("deltaVmean"));
TVectorD *vecDNtrk = (TVectorD*)f->Get(Form("avgmult"));
double *V2 = vecDV2->GetMatrixArray();
double *V2err = vecDV2err->GetMatrixArray();
cout<<V2err[0]<<endl;
double *Ntrk = vecDNtrk->GetMatrixArray();
TGraphErrors *gr=new TGraphErrors(nbin,Ntrk,V2,0,V2err);
gr->SetTitle("V_{2} vs Ntrkoffline");
gr->GetXaxis()->SetTitle("Ntrkoffline");
gr->GetYaxis()->SetTitle("V_{2}");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);
gr->Draw("AP");

f->Close();

}

