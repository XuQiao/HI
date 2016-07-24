#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void Bumpcheck(){
gStyle->SetOptStat(kFALSE);
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
static const int nColor = 7;
static const int colorCode[nColor] = {
    1, 2, 4, 46,6, 7,8
};
static const int markerCode[nColor] = {
    20, 34, 33, 25,27, 28,24
};
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/DATAPPb_Jet60akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
//TString datafile="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
bool IsMC=kTRUE; TString algo="akPu3PF";
const int Ndata=4, NMC=10;
TString datafile[Ndata]; TString MCfile[NMC];
TString JetTrig[Ndata]={"JetTrig20","JetTrig40","JetTrig80","JetTrig100"};
TString pthat[NMC]={"pthat15","pthat30","pthat50","pthat80","pthat120","pthat170","pthat220","pthat280","pthat370","merged"};
int Marker[NMC]={20,34,33,25,27,28,24,33,34,20};
int Color[NMC]={1,2,4,46,6,7,8,8,1,2};

if(algo=="akPu3PF"){
datafile[0]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
datafile[1]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet40akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
datafile[2]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
datafile[3]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet100akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
}
else if(algo=="ak3PF"){
datafile[0]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet20ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
datafile[1]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet40ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
datafile[2]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet80ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
datafile[3]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/newsubmitted/noleadingcut/DATAPPb_Jet100ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
}

MCfile[0]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt15_HP04_prod16_v77_merged_forest_0.root";
MCfile[1]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt30_HP04_prod16_v77_merged_forest_0.root";
MCfile[2]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt50_HP04_prod16_v77_merged_forest_0.root";
MCfile[3]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt80_HP04_prod16_v77_merged_forest_0.root";
MCfile[4]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt120_HP04_prod16_v77_merged_forest_0.root";
MCfile[5]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt170_HP04_prod16_v77_merged_forest_0.root";
MCfile[6]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt220_HP04_prod16_v77_merged_forest_0.root";
MCfile[7]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt280_HP04_prod16_v77_merged_forest_0.root";
MCfile[8]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/MCPPb_akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_pt370_HP04_prod16_v77_merged_forest_0.root";
MCfile[9]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/NoResidual/mergedCSdiff_MCPPb_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";
//MCfile[9]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/mergedCSdiff_MCPPb_ak3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";

if(!IsMC){
TFile *f[Ndata];
TH2F* hSumSumpt[Ndata]; TH2F* hptJES[Ndata];	TProfile* hProSumSumpt[Ndata]; TProfile* hProptJES[4]; int N=Ndata;}
else{
TFile *f[NMC];
TH2F* hSumSumpt[NMC]; TH2F* hptJES[NMC];	TProfile* hProSumSumpt[NMC]; TProfile* hProptJES[NMC]; int N=NMC;}

TLegend *leg1=new TLegend(0.6,0.2,0.9,0.5);
TLegend *leg2=new TLegend(0.6,0.2,0.9,0.5);
c1 = new TCanvas("c1"," ",800,600);
c2 = new TCanvas("c2"," ",800,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.05,0.07,0.03);
makeMultiPanelCanvas(c2,1,1,0,0,0.05,0.07,0.03);
TLatex *T1=new TLatex(0.45,0.85,algo);
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
c1->cd();
for(int i=0;i<N;i++){
if(!IsMC)
f[i] = TFile::Open(datafile[i]);
else
f[i] = TFile::Open(MCfile[i]);
hSumSumpt[i]=(TH2F*)f[i]->Get("jetptSumSumpt_0-100%");
hptJES[i]=(TH2F*)f[i]->Get("jetptJES_0-100%");
hProSumSumpt[i]=(TProfile*)hSumSumpt[i]->ProfileX();
hProptJES[i]=(TProfile*)hptJES[i]->ProfileX();
hProSumSumpt[i]=(TProfile*)hProSumSumpt[i]->Rebin(Nbin_pt,hProSumSumpt[i]->GetName(),binbound_pt);
hProptJES[i]=(TProfile*)hProptJES[i]->Rebin(Nbin_pt,hProptJES[i]->GetName(),binbound_pt);
}
for(i=0;i<N;i++){
hProSumSumpt[i]->GetXaxis()->SetRangeUser(binbound_pt[12],binbound_pt[Nbin_pt-1]);
hProSumSumpt[i]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hProSumSumpt[i]->GetYaxis()->SetTitle("SumSum/p_{T}^{jet}");
hProSumSumpt[i]->SetTitle("");
hProSumSumpt[i]->GetXaxis()->CenterTitle();
hProSumSumpt[i]->GetYaxis()->CenterTitle();
hProSumSumpt[i]->GetYaxis()->SetTitleOffset(1.2);
hProSumSumpt[i]->GetYaxis()->SetRangeUser(0.9,1.1);
if(IsMC){
hProSumSumpt[i]->SetMarkerStyle(Marker[i]);
hProSumSumpt[i]->SetMarkerColor(Color[i]);
}
else{
hProSumSumpt[i]->SetMarkerStyle(markerCode[i]);
hProSumSumpt[i]->SetMarkerColor(colorCode[i]);
}
hProSumSumpt[i]->SetMarkerSize(1.2);
if(!IsMC)
leg1->AddEntry(hProSumSumpt[i],JetTrig[i].Data(),"lp");
else
leg1->AddEntry(hProSumSumpt[i],pthat[i].Data(),"lp");
if(i==0)
hProSumSumpt[i]->Draw();
else
hProSumSumpt[i]->Draw("same");
}
T1->Draw("same");

leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->Draw("same");
if(!IsMC)
c1->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/pic/JetIDcut/JES/SumSumpt_%s.png",algo.Data()));
else
c1->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/pic/JetIDcut/JES/MC_SumSumpt_%s.png",algo.Data()));

c2->cd();
for(int i=0;i<N;i++){
hProptJES[i]->GetXaxis()->SetRangeUser(binbound_pt[12],binbound_pt[Nbin_pt-1]);
hProptJES[i]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hProptJES[i]->GetYaxis()->SetTitle("p_{T}^{raw}/p_{T}^{jet}");
hProptJES[i]->SetTitle("");
hProptJES[i]->GetXaxis()->CenterTitle();
hProptJES[i]->GetYaxis()->CenterTitle();
hProptJES[i]->GetYaxis()->SetTitleOffset(1.2);
hProptJES[i]->GetYaxis()->SetRangeUser(0.88,0.95);
if(IsMC){
hProptJES[i]->SetMarkerStyle(Marker[i]);
hProptJES[i]->SetMarkerColor(Color[i]);
}
else{
hProptJES[i]->SetMarkerStyle(markerCode[i]);
hProptJES[i]->SetMarkerColor(colorCode[i]);
}

hProptJES[i]->SetMarkerSize(1.2);
if(!IsMC)
leg2->AddEntry(hProptJES[i],JetTrig[i].Data(),"lp");
else
leg2->AddEntry(hProptJES[i],pthat[i].Data(),"lp");
if(i==0)
hProptJES[i]->Draw();
else
hProptJES[i]->Draw("same");
}
T1->Draw("same");
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->Draw("same");
if(!IsMC)
c2->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/pic/JetIDcut/JES/JES_%s.png",algo.Data()));
else
c2->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/pic/JetIDcut/JES/MC_JES_%s_pthat.png",algo.Data()));
}
