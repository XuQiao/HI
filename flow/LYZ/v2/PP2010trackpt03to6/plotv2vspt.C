#include "par.h"
void plotv2vspt(){
gStyle->SetOptStat(kFALSE);
int ibin=0;
TFile *f = TFile::Open("mergedv_Prod.root");
TFile *fProd = TFile::Open("mergedv_Prod2.root");
TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("D_%d/avgpt",ibin));

TVectorD *vecDv2_ = (TVectorD*)fProd->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err_ = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",ibin));
TVectorD *vecDavgpt_ = (TVectorD*)fProd->Get(Form("D_%d/avgpt",ibin));

double *avgpt = vecDavgpt->GetMatrixArray();
double *v2 = vecDv2->GetMatrixArray();
double *v2err = vecDv2err->GetMatrixArray();

double *avgpt_ = vecDavgpt_->GetMatrixArray();
double *v2_ = vecDv2_->GetMatrixArray();
double *v2err_ = vecDv2err_->GetMatrixArray();

TCanvas *c1 = new TCanvas();
TGraphErrors *gr=new TGraphErrors(nptv,avgpt,v2,0,v2err);
TH1D* hFrame = new TH1D("hFrame","",1000,0,10);
hFrame->SetTitle("v_{2} vs momentum");
hFrame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
hFrame->GetYaxis()->SetTitle("v_{2}");
hFrame->GetXaxis()->SetRangeUser(0,6.0);
hFrame->GetYaxis()->SetRangeUser(0,0.2);
TGraphErrors *grProd=new TGraphErrors(nptv,avgpt_,v2_,0,v2err_);
gr->SetMarkerSize(1.3);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);

grProd->SetMarkerStyle(29);
grProd->SetMarkerColor(4);
hFrame->Draw();
//agr->Draw("Psame");
grProd->Draw("Psame");
TLegend *tl = new TLegend(0.6,0.2,0.8,0.35);

tl->SetFillColor(0);
tl->SetBorderSize(0);
tl->SetTextSize(0.05);
//tl->AddEntry(gr,"LYZ Sum method","lp");
tl->AddEntry(grProd,"LYZ Prod method","lp");
tl->Draw("same");
TLatex *tlx = new TLatex(0.2,0.3,Form("PP 2010 HM Data, %d<mult<%d",trkbin[1],trkbin[0]));
tlx->SetNDC();
tlx->SetTextSize(0.045);
tlx->Draw("same");
c1->Print("v2vspt.png");
f->Close();

}

