#include "root_setting.h"
void checkDisSingle(){
gStyle->SetOptStat(kFALSE);
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
TString datafile20="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
TString datafile40="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet40akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
TString datafile80="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
TString datafile100="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet100akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";

TFile *f20=TFile::Open(datafile20);
TFile *f40=TFile::Open(datafile40);
TFile *f80=TFile::Open(datafile80);
TFile *f100=TFile::Open(datafile100);
TString histoname="jetptSumSumpt_0-100%";

TH2F* h20=(TH2F*)f20->Get(histoname);
TH2F* h40=(TH2F*)f40->Get(histoname);
TH2F* h80=(TH2F*)f80->Get(histoname);
TH2F* h100=(TH2F*)f100->Get(histoname);
int bin1, bin2;
bin1=70;
bin2=95;
TH1D *h1=(TH1D*)h20->ProjectionY("h1",bin1,bin2);
TH1D *h2=(TH1D*)h40->ProjectionY("h2",bin1,bin2);
TH1D *h3=(TH1D*)h80->ProjectionY("h3",bin1,bin2);
TH1D *h4=(TH1D*)h100->ProjectionY("h4",bin1,bin2);
/*
h1=(TH1D*)h1->Rebin(Nbin_JetID,"h1",binbound_JetID);
h2=(TH1D*)h2->Rebin(Nbin_JetID,"h2",binbound_JetID);
h3=(TH1D*)h3->Rebin(Nbin_JetID,"h3",binbound_JetID);
h4=(TH1D*)h4->Rebin(Nbin_JetID,"h4",binbound_JetID);
*/
h1->Scale(1./h1->Integral());
h2->Scale(1./h2->Integral());
h3->Scale(1./h3->Integral());
h4->Scale(1./h4->Integral());

c1 = new TCanvas("c1","",800,800);
makeMultiPanelCanvas(c1,1,1,0,0,0.08,0.08,0.03);
c1->cd();
c1->SetLogy();
h1->GetXaxis()->SetRangeUser(0.6,1.4);
h1->GetXaxis()->SetTitle("SumSum/p_{T}");
h1->GetYaxis()->SetTitle("Event Fraction");
h1->GetXaxis()->CenterTitle();
h1->GetYaxis()->CenterTitle();
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

TLegend *leg=new TLegend(0.7,0.7,0.88,0.88);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(h1,"Jet20","lp");
leg->AddEntry(h2,"Jet40","lp");
leg->AddEntry(h3,"Jet80","lp");
leg->AddEntry(h4,"Jet100","lp");
h1->Draw();
leg->Draw("same");
h2->Draw("same");
h3->Draw("same");
h4->Draw("same");
TLatex *T1=new TLatex(0.4,0.85,Form("p_{T}: %d -- %d",bin1,bin2));
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");

c1->Print("pic/JetIDcut/SumSumptDis_JetTrig.png");
cout<<h1->GetMean()<<endl;
cout<<h2->GetMean()<<endl;
cout<<h3->GetMean()<<endl;
cout<<h4->GetMean()<<endl;
}
