#include <iomanip>
#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
void DrawVz(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
ofstream fstr("EventSel.txt");
const int N=6; const int NSel=10;
static const int Color[N]={
   1, 2, 4, 46,6, 7,8
};
static const int Style[N] = {
    20, 34, 33, 25,27, 28//,24
};

TFile *fpPb = TFile::Open("/home/xuq/Documents/HI/RpA/output/TreeAna/JetTrigvzandcent.root");
//TFile *fPbp = TFile::Open("/home/xuq/Documents/HI/RpA/output/TreeAna/DATAPPbakPu3PFskimUniv.root");
TFile *fPbp = TFile::Open("/home/xuq/Documents/HI/RpA/output/TreeAna/PbPJetTrigvzandcent.root");
hVzpPb = (TH1F*)fpPb->Get("vz");
hVzPbp = (TH1F*)fPbp->Get("vz");
hVzpPb->Rebin(5);
hVzPbp->Rebin(5);
hVzpPb->Scale(1./hVzpPb->Integral());
hVzPbp->Scale(1./hVzPbp->Integral());

TH1F* hFrame=new TH1F("","",1000,-50,50);
hFrame->GetXaxis()->SetRangeUser(-15,15);
hFrame->GetYaxis()->SetTitle("Event Fraction");
fixedFontHist(hFrame,1.2,1.7);
TLegend *leg=new TLegend(0.30,0.75,0.70,0.90);
TLine *l=new TLine(-15,1,15,1);
l->SetLineStyle(2);
l->SetLineColor(1);

leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextFont(42);
leg->SetTextSize(0.04);
leg->SetHeader("After Event Selection");
c1 = new TCanvas("c1"," ",600,600);
c2 = new TCanvas("c2"," ",600,600);
makeMultiPanelCanvas(c1,1,1,-0.08,0,0.12,0.1,0.03);
makeMultiPanelCanvas(c2,1,1,-0.08,0,0.12,0.1,0.03);

c1->cd(1);
hFrame->GetXaxis()->SetTitle("Vz (cm)");
hFrame->GetYaxis()->SetRangeUser(1e-4,6e-2);
hFrame->DrawCopy();
hVzpPb->SetMarkerSize(1.2);
hVzpPb->SetMarkerStyle(20);
hVzpPb->SetMarkerColor(1);
hVzpPb->SetLineColor(1);

hVzPbp->SetMarkerSize(1.2);
hVzPbp->SetMarkerStyle(24);
hVzPbp->SetMarkerColor(2);
hVzPbp->SetLineColor(2);
leg->AddEntry(hVzpPb,"pPb: Pb going positive Data","lp");
leg->AddEntry(hVzPbp,"Pbp: Proton going positive Data","lp");
hVzpPb->Draw("same");
hVzPbp->Draw("same");
leg->Draw("same");
c2->cd(1);
//hFrame->GetYaxis()->SetRangeUser(0.9,1.25);
hFrame->GetYaxis()->SetRangeUser(0.75,1.09);
hFrame->GetYaxis()->SetTitle("Ratio Pbp/pPb");
hFrame->DrawCopy();
hRatio = (TH1F*)hVzPbp->Clone();
hRatio->Divide(hVzpPb);
hRatio->SetMarkerColor(1);
hRatio->SetLineColor(1);
hRatio->SetMarkerStyle(20);
hRatio->Draw("same");
l->Draw("same");
c1->Print("AftEvSelVz_PbppPb.pdf");
c1->Print("AftEvSelVz_PbppPb.png");
c2->Print("AftEvSelVz_PbppPb_ratio.pdf");
c2->Print("AftEvSelVz_PbppPb_ratio.png");
}
