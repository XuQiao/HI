#include "par.h"
void plotV2vsNtrk(){

TFile *f = TFile::Open("mergedV_Sum.root");
TVectorD *vecDV2 = (TVectorD*)f->Get(Form("Vmean"));
TVectorD *vecDV2err = (TVectorD*)f->Get(Form("deltaVmean"));
//TVectorD *vecDNtrk = (TVectorD*)f->Get(Form("avgmult"));
double Ntrk[npt];
for(int iptbin=0;iptbin<npt;iptbin++){ 
Ntrk[iptbin]=(double)(trkbin[iptbin]+trkbin[iptbin+1])/2;
}
double *V2 = vecDV2->GetMatrixArray();
double *V2err = vecDV2err->GetMatrixArray();
//double *Ntrk = vecDNtrk->GetMatrixArray();
TCanvas *c1 = new TCanvas();
TGraphErrors *gr=new TGraphErrors(nbin,Ntrk,V2,0,0);
gr->GetYaxis()->SetRangeUser(0.04,0.06);
gr->SetTitle("V_{2} vs Ntrkoffline");
gr->GetXaxis()->SetTitle("Ntrkoffline");
gr->GetYaxis()->SetTitle("V_{2}");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);
gr->Draw("AP");
c1->Print("V2vsNtrk.png");
f->Close();

}

