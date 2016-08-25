#include "par.h"
void plotv2vspt(){

int ibin=0;
TFile *f = TFile::Open("mergedv_Prod.root");
TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("D_%d/avgpt",ibin));
double *avgpt = vecDavgpt->GetMatrixArray();
double *v2 = vecDv2->GetMatrixArray();
double *v2err = vecDv2err->GetMatrixArray();
TCanvas *c1 = new TCanvas();
V2vsPt->SetLineColor(2);
V2vsPt->SetTitle("v_{2} vs momentum");
V2vsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
V2vsPt->GetYaxis()->SetTitle("v_{2}");
V2vsPt->Draw();
TGraphErrors *gr=new TGraphErrors(npt,avgpt,v2,0,v2err);
gr->SetTitle("v_{2} vs momentum");
gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr->GetYaxis()->SetTitle("v_{2}");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);
gr->Draw("Psame");
TLegend *tl = new TLegend(0.6,0.2,0.8,0.35);

tl->SetFillColor(0);
tl->SetBorderSize(0);
tl->AddEntry(V2vsPt,"input v_{2}","lp");
tl->AddEntry(gr,"LYZ method","lp");
tl->Draw("same");
c1->Print("v2vspt.png");
f->Close();

}

