#include "root_setting.h"

void overlay(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;


//TString filename1="mergedCSdiff_MCPP_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_hiforest77_hiSignal.root";
//TString filename1="mergedCSdiff_MCPPb_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";
//TString filename1="DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
//TString filename1="DATAPPb_Jet40akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
//TString filename1="DATAPPb_Jet60akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
//TString filename1="DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
//TString filename1="DATAPPb_MBakPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedMB_pPb_SingleTrack_Full_UsingKKForest_v1.root";
TString filename1="CombinedJetTrigak3PFJetAllHist.root";

//TString filename2="mergedCSdiff_MCPP_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_hiforest77_hiSignal.root";
//TString filename2="mergedCSdiff_MCPPb_akPu3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";
//TString filename2="DATAPPb_Jet20akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root";
//TString filename2="DATAPPb_Jet40akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
//TString filename2="DATAPPb_Jet80akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
//TString filename2="DATAPPb_Jet60akPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root";
//TString filename2="DATAPPb_MBakPu3PFDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedMB_pPb_SingleTrack_Full_UsingKKForest_v1.root";
//TString filename2="PPbCombinedJetTrigakPu3PFJetTrkEffHIN12017v5TrkCorr2DCutAllHistHFsumEta4Bin1.root";
TString filename2="PPbCombinedJetTrigak3PFJetTrkEffHIN12017v5TrkCorr2DCutAllHistHFsumEta4Bin1.root";

TString dir1="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual";
TString dir2="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX";	//From Yaxian

TString histoname="jetpt_0-100%";
TString histoname1="jetptJES_0-100%";
//TString histoname1=Form("jetpt%s_0-100%%",JetIDName.Data());

TFile *file1=TFile::Open(Form("%s/%s",dir1.Data(),filename1.Data()));
TFile *file2=TFile::Open(Form("%s/%s",dir2.Data(),filename2.Data()));

TH2F* histo12F=(TH2F*)file1->Get(histoname1);
TH2F* histo22F=(TH2F*)file2->Get(histoname1);
//TH1D* histo1=(TH1D*)histo12F->ProjectionY("histo1",histo12F->GetXaxis()->FindBin(binbound_pt[5]),histo12F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
//TH1D* histo2=(TH1D*)histo22F->ProjectionY("histo2",histo22F->GetXaxis()->FindBin(binbound_pt[5]),histo22F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
TH1D* histo1=(TH1D*)histo12F->ProjectionX("histo1");
TH1D* histo2=(TH1D*)histo22F->ProjectionX("histo2");

//TH1D* histo1=(TH1D*)file1->Get(histoname);
//TH1D* histo2=(TH1D*)file2->Get(histoname);

rehisto1=(TH1D*)histo1->Clone("rehisto1");
rehisto2=(TH1D*)histo2->Clone("rehisto2");
rehisto1=(TH1D*)rehisto1->Rebin(Nbin_pt,"rehisto1",binbound_pt);
rehisto2=(TH1D*)rehisto2->Rebin(Nbin_pt,"rehisto2",binbound_pt);
rehisto1->GetXaxis()->SetRangeUser(binbound_pt[5],binbound_pt[Nbin_pt-1]);
rehisto2->GetXaxis()->SetRangeUser(binbound_pt[5],binbound_pt[Nbin_pt-1]);
normalizeByBinWidth(rehisto1);
normalizeByBinWidth(rehisto2);
//rehisto1->Scale(1/rehisto1->Integral());
//rehisto2->Scale(1/rehisto2->Integral());
rehisto1->SetMarkerStyle(24);
rehisto1->SetMarkerSize(1.2);
rehisto1->SetMarkerColor(1);
rehisto1->GetYaxis()->SetTitle("Number of Events");
rehisto2->SetMarkerStyle(20);
rehisto2->SetMarkerSize(1.2);
rehisto2->SetMarkerColor(2);
fixedFontHist(rehisto1,2,2);
fixedFontHist(rehisto2,2,2);

c1 = new TCanvas("c1","",800,800);
makeMultiPanelCanvas(c1,1,2,0.1,0.12,0.12,0.1,0.03);

c1->cd(1)->SetLogy();
//rehisto1->GetXaxis()->SetTitle(histo12F->GetYaxis()->GetTitle());
//rehisto1->GetXaxis()->SetRangeUser(0.8,1.2);
//rehisto2->GetXaxis()->SetRangeUser(0.8,1.2);
//rehisto1->SetMinimum(1e-8);
//rehisto1->SetMaximum(1);
rehisto1->SetTitle("");
rehisto1->Draw("E1");
rehisto2->Draw("E1same");
TLegend *leg1=new TLegend(0.60,0.75,0.85,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(rehisto1,"Qiao","lp");
leg1->AddEntry(rehisto2,"Yaxian","lp");
leg1->Draw("same");
c1->cd(2);
TH1D* ratio=(TH1D*)rehisto2->Clone("ratio");
ratio->Divide(rehisto1);
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.2);
ratio->SetMarkerColor(1);
ratio->SetMinimum(0.6);
ratio->SetMaximum(1.39);
ratio->SetTitle("");
ratio->GetYaxis()->SetTitle("Ratio: Qiao/Yaxian");
ratio->Draw("E1");
TLine *l =new TLine(binbound_pt[0],1,binbound_pt[Nbin_pt-1],1);
//TLine *l =new TLine(-3,1,3,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
//T->Draw("same");

}
