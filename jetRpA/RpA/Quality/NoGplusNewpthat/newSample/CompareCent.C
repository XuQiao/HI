#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"

void CompareCent(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

TString histoname1="pfjet_Cent";
TString histoname2="pfjet_CentWeighted";

TH1F* histo1=(TH1F*)fPPb->Get(histoname1);
TH1F* histo2=(TH1F*)fPPb->Get(histoname2);
TH1F* histodata=(TH1F*)fdata->Get("cent");

//histo1=(TH1F*)histo1->Rebin(Nbin_vz,"histo1",binbound_vz);
//histo2=(TH1F*)histo2->Rebin(Nbin_vz,"histo2",binbound_vz);
histo1->Scale(1/histo1->Integral());
histo2->Scale(1/histo2->Integral());
histodata->Scale(1/histodata->Integral());
histo1->SetMarkerStyle(24);
histo1->SetMarkerSize(1.2);
histo1->SetMarkerColor(2);
histo1->SetLineColor(2);
histo2->SetMarkerStyle(0);
histo2->SetMarkerSize(0);
histo2->SetFillStyle(3004);
histo2->SetFillColor(2);
histo2->SetLineColor(2);
histo2->SetMarkerColor(0);
//histo2->SetMarkerColor(2);
histodata->SetMarkerStyle(20);
histodata->SetMarkerSize(1.2);
histodata->SetMarkerColor(1);
histodata->SetLineColor(1);

c1 = new TCanvas("c1"," ",500,500);
makeMultiPanelCanvas(c1,1,1,-0.1,0,0.16,0.14,0.03);

TH1F* hFrame=new TH1F("","",20000,-1000,1000);
fixedFontHist(hFrame,1.4,1.8);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("Centrality");
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetXaxis()->SetLimits(0,100);
//hFrame->GetXaxis()->SetLimits(-3,3);
hFrame->GetYaxis()->SetRangeUser(0,5e-2);
c1->cd(1);
hFrame->DrawCopy();
histodata->Draw("same");
histo1->Draw("same");
histo2->Draw("HIST same");
TLegend *leg1=new TLegend(0.40,0.80,0.85,0.92);
TLegend *leg2=new TLegend(0.50,0.9,0.80,0.98);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetTextSize(0.04);
leg1->AddEntry(histo1,"Before Centrality weighting","lp");
leg1->AddEntry(histo2,"After Centrality weighting","lfp");
leg1->AddEntry(histodata,"Data","lp");
leg1->Draw("same");
TLatex *T1=new TLatex(0.25,0.92,"");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
c1->Print("CentCompare_DataPPb.png");
c1->Print("CentCompare_DataPPb.pdf");
}
