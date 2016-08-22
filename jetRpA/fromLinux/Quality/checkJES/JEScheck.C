#include "root_setting.h"
void JEScheck(){
gStyle->SetOptStat(kFALSE);
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

TString datafile20="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
TString datafile40="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet40akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
TString datafile80="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
TString datafile100="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet100akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";

TFile *f20=TFile::Open(datafile20);
TFile *f40=TFile::Open(datafile40);
TFile *f80=TFile::Open(datafile80);
TFile *f100=TFile::Open(datafile100);

TH2F* h20JES=(TH2F*)f20->Get("jetptJES_0-100%");
TH2F* h40JES=(TH2F*)f40->Get("jetptJES_0-100%");
TH2F* h80JES=(TH2F*)f80->Get("jetptJES_0-100%");
TH2F* h100JES=(TH2F*)f100->Get("jetptJES_0-100%");

TProfile *Pro_h20=(TProfile*)h20JES->ProfileX("h20JES");
TProfile *Pro_h40=(TProfile*)h40JES->ProfileX("h40JES");
TProfile *Pro_h80=(TProfile*)h80JES->ProfileX("h80JES");
TProfile *Pro_h100=(TProfile*)h100JES->ProfileX("h100JES");

Pro_h20=(TProfile*)Pro_h20->Rebin(Nbin_pt,Pro_h20->GetName(),binbound_pt);
Pro_h40=(TProfile*)Pro_h40->Rebin(Nbin_pt,Pro_h40->GetName(),binbound_pt);
Pro_h80=(TProfile*)Pro_h80->Rebin(Nbin_pt,Pro_h80->GetName(),binbound_pt);
Pro_h100=(TProfile*)Pro_h100->Rebin(Nbin_pt,Pro_h100->GetName(),binbound_pt);

c1 = new TCanvas("c1"," ",800,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.05,0.07,0.03);
Pro_h20->GetXaxis()->SetTitle("p_{T} (GeV/c)");
Pro_h20->GetYaxis()->SetRangeUser(0.8,1);
Pro_h20->SetMarkerColor(1);
Pro_h20->SetLineColor(1);
Pro_h20->SetMarkerStyle(20);

Pro_h40->SetMarkerColor(6);
Pro_h40->SetLineColor(6);
Pro_h40->SetMarkerStyle(22);

Pro_h80->SetMarkerColor(2);
Pro_h80->SetLineColor(2);
Pro_h80->SetMarkerStyle(24);

Pro_h100->SetMarkerColor(4);
Pro_h100->SetLineColor(4);
Pro_h100->SetMarkerStyle(29);
TLegend *leg=new TLegend(0.7,0.8,0.9,0.9);
leg->AddEntry(Pro_h20,"JetTrig20","lp");
leg->AddEntry(Pro_h40,"JetTrig40","lp");
leg->AddEntry(Pro_h80,"JetTrig80","lp");
leg->AddEntry(Pro_h100,"JetTrig100","lp");

Pro_h20->Draw();
Pro_h40->Draw("same");
Pro_h80->Draw("same");
Pro_h100->Draw("same");
leg->Draw("same");
}

