#include "root_setting.h"
void checkBump(){
const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/DATAPPb_Jet60akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";

//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet100akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
TString PPbfile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/mergedCSdiff_MCPPb_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";
TFile *f = TFile::Open(PPbfile);
TH2F *hSumSumpt=(TH2F*)f->Get("jetptSumSumpt_0-100%");
TH2F *hSumSumrawpt=(TH2F*)f->Get("jetptSumSumrawpt_0-100%");
TH2F *hptJES=(TH2F*)f->Get("jetptJES_0-100%");
bool isRebin=kFALSE;
TProfile* hProSumSumpt=(TProfile*)hSumSumpt->ProfileX();
TProfile* hProSumSumrawpt=(TProfile*)hSumSumrawpt->ProfileX();
TProfile* hProptJES=(TProfile*)hptJES->ProfileX();
if(isRebin){
hProSumSumpt=(TProfile*)hProSumSumpt->Rebin(Nbin_pt,hProSumSumpt->GetName(),binbound_pt);
hProSumSumrawpt=(TProfile*)hProSumSumrawpt->Rebin(Nbin_pt,hProSumSumrawpt->GetName(),binbound_pt);
hProptJES=(TProfile*)hProptJES->Rebin(Nbin_pt,hProptJES->GetName(),binbound_pt);
}
double SumSumpt[1000],SumSumrawpt[1000], ptJES[1000], pt[1000], x[1000];
for(int nbin=1;nbin<=hProSumSumpt->GetNbinsX();nbin++){
if(isRebin)
x[nbin-1]=(binbound_pt[nbin-1]+binbound_pt[nbin])/2;
else
x[nbin-1]=nbin-0.5;
SumSumpt[nbin-1]=hProSumSumpt->GetBinContent(nbin);
SumSumrawpt[nbin-1]=hProSumSumrawpt->GetBinContent(nbin);
ptJES[nbin-1]=hProptJES->GetBinContent(nbin);
if(ptJES[nbin-1])
pt[nbin-1]=SumSumpt[nbin-1]/ptJES[nbin-1];
}

TGraph *gpt = new TGraph(hProSumSumpt->GetNbinsX(),x,pt);

c1 = new TCanvas("c1"," ",800,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.05,0.07,0.03);
gpt->GetXaxis()->SetRangeUser(binbound_pt[5],binbound_pt[Nbin_pt]);
gpt->GetXaxis()->SetTitle("p_{T}^{jet}");
gpt->GetXaxis()->CenterTitle();
gpt->GetYaxis()->SetRangeUser(0.75,1.25);
gpt->SetMarkerStyle(20);
gpt->SetMarkerColor(1);
gpt->SetMarkerSize(1.2);
gpt->Draw("AP");

hProptJES->SetMarkerStyle(30);
hProptJES->SetMarkerColor(4);
hProptJES->SetMarkerSize(1.2);
hProptJES->Draw("same");

hProSumSumpt->SetMarkerStyle(29);
hProSumSumpt->SetMarkerColor(4);
hProSumSumpt->SetMarkerSize(1.2);
hProSumSumpt->Draw("same");

hProSumSumrawpt->SetMarkerStyle(24);
hProSumSumrawpt->SetMarkerColor(2);
hProSumSumrawpt->SetMarkerSize(1.2);
hProSumSumrawpt->Draw("same");

TLegend *leg=new TLegend(0.7,0.7,0.9,0.9);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(hProSumSumpt,"SumSum/jetpt","lp");
leg->AddEntry(hProSumSumrawpt,"SumSum/rawpt","lp");
leg->AddEntry(hProptJES,"JES rawpt/jetpt","lp");
leg->AddEntry(gpt,"after JES corr","lp");
leg->Draw("same");
//c1->Print("~/Rebin_JES100.png");
}
