#include "par.h"
void plotv2vseta(){

int ibin=0;
//TFile *f = TFile::Open("mergedv_Prod2_inV2.root");
TFile *fProd = TFile::Open("mergedv_Prod2_eta.root");
//TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/vmean",ibin));
//TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
//TVectorD *vecDavgeta = (TVectorD*)f->Get(Form("D_%d/avgeta",ibin));

TVectorD *vecDv2_ = (TVectorD*)fProd->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err_ = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",ibin));
TVectorD *vecDavgeta_ = (TVectorD*)fProd->Get(Form("D_%d/avgeta",ibin));

//double *avgeta = vecDavgeta->GetMatrixArray();
//double *v2 = vecDv2->GetMatrixArray();
//double *v2err = vecDv2err->GetMatrixArray();

double *avgeta_ = vecDavgeta_->GetMatrixArray();
double *v2_ = vecDv2_->GetMatrixArray();
double *v2err_ = vecDv2err_->GetMatrixArray();

double v2corr[netav],v2corr_[netav];
for(int i=0;i<netav;i++){
//v2corr[i]=v2[i];
v2corr_[i]=v2_[i];
}

TCanvas *c1 = new TCanvas();
//TGraphErrors *gr=new TGraphErrors(netav,avgeta,v2corr,0,v2err);
TF1 *V2vseta = new TF1("V2vseta",Form("%f",inV2),-2.4,2.4);
V2vseta->SetTitle("");
V2vseta->GetXaxis()->SetTitle("#eta");
V2vseta->GetYaxis()->SetTitle("v_{2}");
V2vseta->GetXaxis()->CenterTitle();
V2vseta->GetYaxis()->CenterTitle();
V2vseta->GetXaxis()->SetTitleSize(0.05);
V2vseta->GetYaxis()->SetTitleSize(0.05);
V2vseta->GetXaxis()->SetTitleOffset(0.9);
V2vseta->GetYaxis()->SetTitleOffset(0.9);
V2vseta->Draw();

TGraphErrors *grProd=new TGraphErrors(netav,avgeta_,v2corr_,0,v2err_);
grProd->SetTitle("");
grProd->GetXaxis()->SetTitle("#eta");
grProd->GetYaxis()->SetTitle("v_{2}");
//gr->SetMarkerSize(1.3);
//gr->SetMarkerColor(1);
//gr->SetMarkerStyle(20);

grProd->SetMarkerStyle(20);
grProd->SetMarkerColor(4);
grProd->SetMarkerSize(1.3);
//gr->Draw("Psame");
grProd->Draw("Psame");
TLegend *tl = new TLegend(0.6,0.75,0.8,0.85);

tl->SetFillColor(0);
tl->SetBorderSize(0);
tl->SetTextSize(0.05);
//tl->AddEntry(gr,"LYZ Prod method inV2","lp");
tl->AddEntry(V2vseta,"input v_{2}","lp");
tl->AddEntry(grProd,"v_{2}{LYZ}","lp");
tl->Draw("same");
TLatex *tlx = new TLatex(0.12,0.8,Form("STEG, 15M events, mult=%d",(trkbin[0]+trkbin[1])/2));
tlx->SetNDC();
tlx->SetTextSize(0.045);
tlx->Draw("same");
c1->Print("v2vseta.png");
c1->Print("v2vseta.pdf");
//f->Close();

}

