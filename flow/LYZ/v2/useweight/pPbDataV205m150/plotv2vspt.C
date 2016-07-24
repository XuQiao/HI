#include "par.h"
void plotv2vspt(){

int ibin=0;
TFile *fProd = TFile::Open("mergedv_Prod2.root");
TVectorD *vecDv2_ = (TVectorD*)fProd->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err_ = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",ibin));
TVectorD *vecDavgpt_ = (TVectorD*)fProd->Get(Form("D_%d/avgpt",ibin));

double *avgpt_ = vecDavgpt_->GetMatrixArray();
double *v2_ = vecDv2_->GetMatrixArray();
double *v2err_ = vecDv2err_->GetMatrixArray();

double v2corr[nptv],v2corr_[nptv];
for(int i=0;i<nptv;i++){
v2corr_[i]=v2_[i];
}

TCanvas *c1 = new TCanvas();
V2vsPt->SetLineColor(2);
V2vsPt->SetTitle("v_{2} vs momentum");
V2vsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
V2vsPt->GetYaxis()->SetTitle("v_{2}");
V2vsPt->GetYaxis()->SetRangeUser(-0.1,0.13);
V2vsPt->Draw();
TGraphErrors *grProd=new TGraphErrors(nptv,avgpt_,v2corr_,0,v2err_);
grProd->SetTitle("v_{2} vs momentum");
grProd->GetXaxis()->SetTitle("p_{T} (GeV/c)");
grProd->GetYaxis()->SetTitle("v_{2}");

grProd->SetMarkerStyle(29);
grProd->SetMarkerColor(4);
grProd->Draw("Psame");
TLegend *tl = new TLegend(0.6,0.2,0.8,0.35);

tl->SetFillColor(0);
tl->SetBorderSize(0);
tl->SetTextSize(0.05);
tl->AddEntry(V2vsPt,"input v_{2}","lp");
tl->AddEntry(grProd,"LYZ Prod method","lp");
tl->Draw("same");
TLatex *tlx = new TLatex(0.2,0.3,Form("STEG, 15M events, mult=%d",(trkbin[0]+trkbin[1])/2));
tlx->SetNDC();
tlx->SetTextSize(0.045);
tlx->Draw("same");
c1->Print("v2vspt.png");

}

