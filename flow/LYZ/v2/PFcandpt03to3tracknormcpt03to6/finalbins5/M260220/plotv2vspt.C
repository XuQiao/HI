#include "par.h"
void plotv2vspt(){
const int npt24 = 9;
const double pt24[npt24]={0.4,0.7,1.2,1.7,2.2,2.7,3.4,4.4,5.4};
const double v24[npt24]={0.02507,0.04214,0.06461,0.08294,0.09241,0.09917,0.09662,0.08479,0.06271};
const double v24err[npt24]={0.00079,0.00057,0.00084,0.0013,0.002,0.0032,0.004,0.0075,0.012};

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

double v2corr[nptv],v2corr_[nptv];
for(int i=0;i<nptv;i++){
v2corr[i]=v2[i];
v2corr_[i]=v2_[i];
}

TCanvas *c1 = new TCanvas();
TGraphErrors *gr=new TGraphErrors(nptv,avgpt,v2corr,0,v2err);

TGraphErrors *grProd=new TGraphErrors(nptv,avgpt_,v2corr_,0,v2err_);
TGraphErrors *grv24 = new TGraphErrors(npt24,pt24,v24,v24err);
gr->SetTitle("v_{2} vs momentum");
gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr->GetYaxis()->SetTitle("v_{2}");
gr->GetYaxis()->SetRangeUser(-0.1,0.3);
gr->SetMarkerSize(1.3);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);

grv24->SetMarkerSize(1.3);
grv24->SetMarkerColor(2);
grv24->SetMarkerStyle(20);

grProd->SetMarkerStyle(29);
grProd->SetMarkerColor(4);
gr->Draw("AP");
grProd->Draw("Psame");
grv24->Draw("Psame");
TLegend *tl = new TLegend(0.6,0.2,0.8,0.35);

tl->SetFillColor(0);
tl->SetBorderSize(0);
tl->SetTextSize(0.05);
tl->AddEntry(gr,"LYZ Sum method","lp");
tl->AddEntry(grProd,"LYZ Prod method","lp");
tl->AddEntry(grv24,"v_{2} 4-particle cumulant","lp");
tl->Draw("same");
TLatex *tlx = new TLatex(0.2,0.3,"");
tlx->SetNDC();
tlx->SetTextSize(0.045);
tlx->Draw("same");
c1->Print("v2vspt.png");
f->Close();

}

