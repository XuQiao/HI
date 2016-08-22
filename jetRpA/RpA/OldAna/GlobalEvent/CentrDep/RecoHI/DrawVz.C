#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void DrawVz(){

gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

TFile *f=TFile::Open("NoVRw/Closure_v1.root","ReadOnly");
TH1F* hVz_MC=(TH1F*)f->Get("demo/vz");
hVz_MC->GetXaxis()->SetRangeUser(-15,15);
TH1F* hFrame=new TH1F("","",300,-30,30);
fixedFontHist(hFrame,1.9,2.8);
TString Data_file="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/PPbCombinedJetTrigakPu3PFJetTrkEffHIN12017v5TrkCorr2DCutAllHistHFsumEta4Bin1.root";
TFile *fdata=TFile::Open(Data_file);
c1 = new TCanvas("c1","",400,600);
makeMultiPanelCanvas(c1,1,2,-0.14,0,0.2,0.14,0.03);
c1->cd(1);
hVz_MC->SetFillStyle(3004);
hVz_MC->SetFillColor(2);
hFrame->SetTitle("");
hFrame->SetMaximum(2e-2);
hFrame->SetMinimum(5e-5);
hFrame->GetXaxis()->SetRangeUser(-20,20);
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->DrawCopy();
hVz_Data=(TH1F*)fdata->Get("pfjet_Vz");
hVz_Data->Rebin(2);
hVz_Data->GetXaxis()->SetRangeUser(-15,15);
hVz_Data->SetMarkerStyle(20);
hVz_Data->SetMarkerSize(1.1);
hVz_Data->SetMarkerColor(1);
normalizeByBinWidth(hVz_MC);
normalizeByBinWidth(hVz_Data);
hVz_MC->Scale(1/hVz_MC->Integral());
hVz_Data->Scale(1/hVz_Data->Integral());
hVz_Data->Draw("same");
hVz_MC->Draw("HIST same");
TLegend *leg=new TLegend(0.75,0.75,0.9,0.9);
leg->AddEntry(hVz_Data,"Data","lp");
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.06);
leg->AddEntry(hVz_MC,"MC","lp");
leg->Draw("same");
c1->cd(2);
hFrame->SetMaximum(2.29);
hFrame->SetMinimum(0.1);
hFrame->GetXaxis()->SetTitle("Vz (cm)");
hFrame->GetYaxis()->SetTitle("Data/MC");
hFrame->DrawCopy();
TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
fVz->SetLineColor(2);
hVz_MC_cp=(TH1F*)hVz_Data->Clone("hVz_MC_cp");
for(int i=0;i<hVz_MC_cp->GetNbinsX();i++)
hVz_MC_cp->SetBinContent(i,hVz_MC->GetBinContent(50+i));
hVz_ratio=(TH1F*)hVz_Data->Clone("hVz_ratio");
hVz_ratio->GetYaxis()->SetTitle("Data/MC");
hVz_ratio->GetXaxis()->SetTitle("Vz (cm)");
hVz_ratio->Divide(hVz_MC_cp);
hVz_ratio->GetXaxis()->SetRangeUser(-20,20);
hVz_ratio->SetMinimum(0);
hVz_ratio->Draw("same");
fVz->Draw("same");
TLegend *leg1=new TLegend(0.45,0.75,0.65,0.95);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.05);
leg1->AddEntry(hVz_ratio,"ratio","lp");
//leg1->AddEntry(fVz,"#splitline{reweighting function}{from jetRpA analysis}","lp");
leg1->AddEntry(fVz,"reweighting function","lp");
leg1->Draw("same");
c1->Print("Hijing_VzWeighting.png");

}

