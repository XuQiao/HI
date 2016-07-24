#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void PbppPbCompareVz(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TString filename1="JetTrigvzandcent.root";

TString filename2="JetTrigvzandcentPbp.root";

TString dir="/scratch/xuq7/RpA/TreeAna";

TString histoname="vz";

TFile *file1=TFile::Open(Form("%s/%s",dir.Data(),filename1.Data()));
TFile *file2=TFile::Open(Form("%s/%s",dir.Data(),filename2.Data()));

TH1F* histo1=(TH2F*)file1->Get(histoname);
TH1F* histo2=(TH2F*)file2->Get(histoname);
histo1->Rebin(5);
histo2->Rebin(5);
normalizeByBinWidth(histo1);
normalizeByBinWidth(histo2);
histo1->Scale(1/histo1->Integral());
histo2->Scale(1/histo2->Integral());
histo1->SetMarkerStyle(24);
histo1->SetMarkerSize(1.2);
histo1->SetMarkerColor(1);
histo2->SetMarkerStyle(20);
histo2->SetMarkerSize(1.2);
histo2->SetMarkerColor(2);

c1 = new TCanvas("c1"," ",500,800);
makeMultiPanelCanvas(c1,1,2,-0.08,0,0.14,0.16,0.03);

TH1F* hFrame=new TH1F("","",20000,-1000,1000);
fixedFontHist(hFrame,2,2.2);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetXaxis()->SetLimits(-15,15);
//hFrame->GetXaxis()->SetLimits(-3,3);
hFrame->GetYaxis()->SetRangeUser(1e-4,1e-1);
c1->cd(1);
hFrame->DrawCopy();
histo1->Draw("E1same");
histo2->Draw("E1same");
TLegend *leg1=new TLegend(0.60,0.75,0.85,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(histo1,"pPb","lp");
leg1->AddEntry(histo2,"Pbp","lp");
leg1->Draw("same");
c1->cd(2);
TLatex *T1=new TLatex(0.5,0.92,"");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
TH1D* ratio=(TH1D*)histo2->Clone("ratio");
hFrame->GetXaxis()->SetTitle("vz (cm)");
hFrame->GetYaxis()->SetTitle("Ratio: Pbp/pPb");
hFrame->GetYaxis()->SetRangeUser(0.6,1.39);
hFrame->DrawCopy();
ratio->Divide(histo1);
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.2);
ratio->SetMarkerColor(1);
ratio->SetTitle("");
ratio->Draw("E1same");
TLine *l =new TLine(-15,1,15,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T1->Draw("same");
c1->Print("pic/VzCompare_PbppPb.png");
}
