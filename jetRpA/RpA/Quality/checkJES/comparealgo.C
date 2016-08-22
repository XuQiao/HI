#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void comparealgo(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
static const int nColor = 8;
static const int colorCode[nColor] = {
    1, 2, 4, 46,6, 7,8,4
};
static const int markerCode[nColor] = {
    20, 34, 33, 25,27, 28,24,33
};

const int Ndata=4;
const int Nalgo=2;
TString datafile[Ndata][Nalgo];

TString algoName[Nalgo]={"akPu3PF","ak3PF"};
TString TrigName[Ndata]={"JetTrig80","JetTrig40","JetTrig20","PYTHIA+HIJING"};

datafile[0][0]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/NoResidual/noleadingcut/DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
datafile[0][1]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/noleadingcut/DATAPPb_Jet80ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
datafile[1][0]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/NoResidual/noleadingcut/DATAPPb_Jet40akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
datafile[1][1]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/noleadingcut/DATAPPb_Jet40ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
datafile[2][0]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/NoResidual/noleadingcut/DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
datafile[2][1]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/noleadingcut/DATAPPb_Jet20ak3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
datafile[3][0]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/akPu3PF/NoResidual/mergedCSdiff_MCPPb_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";
datafile[3][1]="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/mergedCSdiff_MCPPb_ak3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";
TFile *f[Ndata][Nalgo];
TH2F* hptJES[Ndata][Nalgo];	 TProfile* hProptJES[Ndata][Nalgo]; 

TLegend *leg1=new TLegend(0.45,0.2,0.6,0.5);
TLegend *leg2=new TLegend(0.45,0.2,0.6,0.5);
leg1->SetTextSize(0.04);
leg2->SetTextSize(0.04);
c1 = new TCanvas("c1"," ",1000,520);
makeMultiPanelCanvas(c1,2,1,0,0,0.2,0.18,0.03);
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(30.1,599.9);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
hFrame->GetYaxis()->SetTitle("p_{T}^{jet}/p_{T}^{raw}");
hFrame->GetYaxis()->SetRangeUser(0.9,1.2);
fixedFontHist(hFrame,1.2,1.5);
TLatex *T1=new TLatex(0.55,0.83,"Data: No leadingjet cut");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);
for(int i=0;i<Ndata;i++){
for(int j=0;j<Nalgo;j++){
f[i][j] = TFile::Open(datafile[i][j]);
hptJES[i][j]=(TH2F*)f[i][j]->Get("jetrawptJES_0-100%");
hptJES[i][j]->SetName(Form("hrawptJES_%d_%d",i,j));
hProptJES[i][j]=(TProfile*)hptJES[i][j]->ProfileX();
hProptJES[i][j]=(TProfile*)hProptJES[i][j]->Rebin(Nbin_pt,hProptJES[i][j]->GetName(),binbound_pt);
fixedFontHist(hProptJES[i][j],1.2,1.5);
}
}

c1->cd(1);
hFrame->GetYaxis()->SetTitle("");
hFrame->DrawCopy();
for(int i=0;i<Ndata;i++){
if(i==1 || i==2) continue;
for(int j=0;j<Nalgo;j++){
hProptJES[i][j]->SetMarkerStyle(markerCode[i+4*j]);
hProptJES[i][j]->SetMarkerColor(colorCode[i+4*j]);
hProptJES[i][j]->SetMarkerSize(1.4);
leg1->AddEntry(hProptJES[i][j],Form("%s %s",algoName[j].Data(),TrigName[i].Data()),"lp");
hProptJES[i][j]->DrawCopy("same");
}
}

T1->Draw("same");
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->Draw("same");
drawCMS(0.4,0.9);

c1->cd(2);
hFrame->GetYaxis()->SetTitle("");
hFrame->DrawCopy();
for(int i=0;i<Ndata;i++){
if(i==3) continue;
for(int j=0;j<Nalgo;j++){
hProptJES[i][j]->GetYaxis()->SetRangeUser(1.0,1.2);
hProptJES[i][j]->SetMarkerStyle(markerCode[i+4*j]);
hProptJES[i][j]->SetMarkerColor(colorCode[i+4*j]);
hProptJES[i][j]->SetMarkerSize(1.4);
leg2->AddEntry(hProptJES[i][j],Form("%s %s",algoName[j].Data(),TrigName[i].Data()),"lp");
hProptJES[i][j]->Draw("same");
}
}

T1->Draw("same");
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->Draw("same");
drawCMS(0.4,0.9);
c1->Print("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/pic/JetIDcut/JES/JES_compalgo.png");
}
