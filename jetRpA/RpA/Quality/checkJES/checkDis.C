#include "root_setting.h"
void checkDis(){
gStyle->SetOptStat(kFALSE);
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet40akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet100akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/CombinedJetTrigakPu3PFJetAllHist.root";
TFile *f=TFile::Open(datafile);

TH2F* h=(TH2F*)f->Get("jetptSumSumpt_0-100%");
TH1D* h1= (TH1D*)h->ProjectionY("h1",80,90);
TH1D* h2= (TH1D*)h->ProjectionY("h2",90,95);
TH1D* h3= (TH1D*)h->ProjectionY("h3",96,100);
TH1D* h4= (TH1D*)h->ProjectionY("h4",100,110);

h1->Scale(1./h1->Integral());
h2->Scale(1./h2->Integral());
h3->Scale(1./h3->Integral());
h4->Scale(1./h4->Integral());


c1 = new TCanvas("c1","",800,800);
makeMultiPanelCanvas(c1,1,1,0,0,0.06,0.08,0.03);
c1->cd();
c1->SetLogy();
h1->GetXaxis()->SetRangeUser(0,2);
h1->GetXaxis()->SetTitle("SumSum/p_{T}");
h1->SetTitle("");
h1->SetMarkerStyle(20);
h1->SetMarkerColor(2);
h1->SetLineColor(2);
h2->SetMarkerStyle(24);
h2->SetMarkerColor(1);
h2->SetLineColor(1);
h3->SetMarkerStyle(29);
h3->SetMarkerColor(4);
h3->SetLineColor(4);
h4->SetMarkerStyle(30);
h4->SetMarkerColor(6);
h4->SetLineColor(6);

TLegend *leg=new TLegend(0.7,0.7,0.9,0.8);
leg->AddEntry(h1,"p_{T}=92-93","lp");
leg->AddEntry(h2,"p_{T}=94-95","lp");
leg->AddEntry(h3,"p_{T}=95-96","lp");
leg->AddEntry(h4,"p_{T}=97-98","lp");
h1->Draw();
leg->Draw("same");
h2->Draw("same");
h3->Draw("same");
h4->Draw("same");
cout<<h1->GetMean()<<endl;
cout<<h2->GetMean()<<endl;
cout<<h3->GetMean()<<endl;
cout<<h4->GetMean()<<endl;
}
