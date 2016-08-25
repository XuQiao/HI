#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void DrawVz(){

gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

TFile *f=TFile::Open("output/Hijing_vz.root","ReadOnly");
TH1F* hVz_MC=(TH1F*)f->Get("hVz_MC");
TH2F* hVzHF4=(TH2F*)f->Get("hVzHF4");
TH2F* hVzHFPlus4=(TH2F*)f->Get("hVzHFPlus4");
TH2F* hVzHFPlus=(TH2F*)f->Get("hVzHFPlus");
fixedFontHist(hVz_MC);
fixedFontHist(hVzHF4);
fixedFontHist(hVzHFPlus4);
fixedFontHist(hVzHFPlus);
TString Data_file="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/latest/PPbCombinedJetTrigakPu3PFJetTrkEffHIN12017v5TrkCorr2DCutAllHistHFsumEta4Bin1.root";
TFile *fdata=TFile::Open(Data_file);
c1 = new TCanvas("c1","",800,800);
c2 = new TCanvas("c2","",600,600);
c3 = new TCanvas("c3","",600,600);
c4 = new TCanvas("c4","",600,600);
makeMultiPanelCanvas(c1,1,2,0.03,0.03,0.12,0.12,0.03);
makeMultiPanelCanvas(c2,1,1,0.03,0.04,0.02,0.04,0.03);
makeMultiPanelCanvas(c3,1,1,0.03,0.04,0.02,0.04,0.03);
makeMultiPanelCanvas(c4,1,1,0.03,0.04,0.02,0.04,0.03);
c1->cd(1);
hVz_MC->SetFillStyle(3004);
hVz_MC->SetFillColor(2);
hVz_MC->SetTitle("");
hVz_MC->SetMaximum(1);
hVz_MC->GetXaxis()->SetTitle("Vz (cm)");
hVz_MC->GetXaxis()->SetTitleOffset(1.1);
hVz_MC->GetYaxis()->SetTitle("Event Fraction");
hVz_Data=(TH1F*)fdata->Get("pfjet_Vz");
hVz_Data->GetYaxis()->SetTitle("Event Fraction");
hVz_Data->SetMarkerStyle(20);
hVz_Data->SetMarkerSize(1.2);
hVz_Data->SetMarkerColor(1);
normalizeByBinWidth(hVz_MC);
normalizeByBinWidth(hVz_Data);
hVz_MC->Scale(1/hVz_MC->Integral());
hVz_Data->Scale(1/hVz_Data->Integral());
hVz_Data->Draw();
hVz_MC->Draw("HIST same");
TLegend *leg=new TLegend(0.8,0.8,0.95,0.9);
leg->AddEntry(hVz_Data,"Data","lp");
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.05);
leg->AddEntry(hVz_MC,"MC","lp");
leg->Draw("same");
c1->cd(2);

TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
fVz->SetLineColor(2);

hVz_ratio=(TH1F*)hVz_Data->Clone("hVz_ratio");
hVz_ratio->GetYaxis()->SetTitle("Data/MC");
hVz_ratio->GetXaxis()->SetTitle("Vz (cm)");
hVz_ratio->Divide(hVz_MC);
fixedFontHist(hVz_ratio);
hVz_ratio->GetXaxis()->SetRangeUser(-20,20);
hVz_ratio->SetMaximum(2.49);
hVz_ratio->SetMinimum(0);
hVz_ratio->Draw();
fVz->Draw("same");
TLegend *leg1=new TLegend(0.7,0.75,0.9,0.9);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.05);
leg1->AddEntry(hVz_ratio,"ratio","lp");
leg1->AddEntry(fVz,"#splitline{reweighting}{function}","lp");
leg1->Draw("same");
c1->Print("pic/Hijing_VzWeighting.png");

c2->cd();
hVzHF4->SetTitle("");
hVzHF4->GetYaxis()->SetTitleOffset(0.9);
hVzHF4->GetXaxis()->SetTitle("Vz (cm)");
hVzHF4->GetYaxis()->SetTitle("HF E_{T} |#eta|>4 (GeV/c)");
hVzHF4->Draw("colz");
c2->Print("pic/Hijing_VzHF4.png");

c3->cd();
hVzHFPlus4->SetTitle("");
hVzHFPlus4->GetYaxis()->SetTitleOffset(0.9);
hVzHFPlus4->GetXaxis()->SetTitle("Vz (cm)");
hVzHFPlus4->GetYaxis()->SetTitle("HF E_{T} #eta>4 (GeV/c)");
hVzHFPlus4->Draw("colz");

c4->cd();
hVzHFPlus->SetTitle("");
hVzHFPlus->GetYaxis()->SetTitleOffset(0.9);
hVzHFPlus->GetXaxis()->SetTitle("Vz (cm)");
hVzHFPlus->GetYaxis()->SetTitle("HF E_{T} #eta>0 (GeV/c)");
hVzHFPlus->Draw("colz");

}

