#include "par.h"
void plotv2vseta(){

int ibin=0;
TFile *f = TFile::Open("mergedv_Prod2_eta.root");
TFile *fProd = TFile::Open("mergedv_Prod2_eta.root");
TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
TVectorD *vecDavgeta = (TVectorD*)f->Get(Form("D_%d/avgeta",ibin));

TVectorD *vecDv2_ = (TVectorD*)fProd->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err_ = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",ibin));
TVectorD *vecDavgeta_ = (TVectorD*)fProd->Get(Form("D_%d/avgeta",ibin));

double *avgeta = vecDavgeta->GetMatrixArray();
double *v2 = vecDv2->GetMatrixArray();
double *v2err = vecDv2err->GetMatrixArray();

double *avgeta_ = vecDavgeta_->GetMatrixArray();
double *v2_ = vecDv2_->GetMatrixArray();
double *v2err_ = vecDv2err_->GetMatrixArray();

double v2corr[netav],v2corr_[netav];
for(int i=0;i<netav;i++){
v2corr[i]=v2[i];
v2corr_[i]=v2_[i];
}

TCanvas *c1 = new TCanvas();
/*V2vsPt->SetLineColor(2);
V2vsPt->SetTitle("v_{2} vs momentum");
V2vsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
V2vsPt->GetYaxis()->SetTitle("v_{2}");
V2vsPt->GetYaxis()->SetRangeUser(0,0.1);
V2vsPt->Draw();
*/
TGraphErrors *gr=new TGraphErrors(netav,avgeta,v2corr,0,v2err);

TGraphErrors *grProd=new TGraphErrors(netav,avgeta_,v2corr_,0,v2err_);
gr->SetTitle("v_{2} vs #eta");
gr->GetXaxis()->SetTitle("#eta");
gr->GetYaxis()->SetTitle("v_{2}");
gr->SetMarkerSize(1.3);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);

grProd->SetMarkerStyle(29);
grProd->SetMarkerColor(4);
gr->Draw("AP");
grProd->Draw("Psame");
TLegend *tl = new TLegend(0.6,0.2,0.8,0.35);

tl->SetFillColor(0);
tl->SetBorderSize(0);
tl->SetTextSize(0.05);
//tl->AddEntry(V2vsPt,"input v_{2}","lp");
tl->AddEntry(gr,"LYZ Sum method","lp");
tl->AddEntry(grProd,"LYZ Prod method","lp");
tl->Draw("same");
TLatex *tlx = new TLatex(0.2,0.3,"260<multi<300");
tlx->SetNDC();
tlx->SetTextSize(0.045);
tlx->Draw("same");
c1->Print("v2vseta.png");
f->Close();

}

