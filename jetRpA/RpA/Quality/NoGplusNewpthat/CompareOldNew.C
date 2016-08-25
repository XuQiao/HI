#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void CompareOldNew(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TString gen="PPb";
//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

TString filename1=Form("MC%sNoGplus_ak3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root",gen.Data());
if(gen=="PP")
TString filename2="mergedCSdiff_MCPP_ak3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_hiforest77_hiSignal.root";
else if(gen=="PPb")
TString filename2="mergedCSdiff_MCPPb_ak3PFIncJet_TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_HP04_prod16_v77_merged_forest_0.root";

TString dir1="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/NoGplusNewpthat";
TString dir2="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/";
//TString dir2="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX";

TString histoname="jetpt_0-100%";
TString histoname1="jetptEta_0-100%";
//TString histoname1=Form("jetpt%s_0-100%%",JetIDName.Data());

TFile *file1=TFile::Open(Form("%s/%s",dir1.Data(),filename1.Data()));
TFile *file2=TFile::Open(Form("%s/%s",dir2.Data(),filename2.Data()));

TH2F* histo12F=(TH2F*)file1->Get(histoname1);
TH2F* histo22F=(TH2F*)file2->Get(histoname1);
//TH1D* histo1=(TH1D*)histo12F->ProjectionY("histo1",histo12F->GetXaxis()->FindBin(binbound_pt[5]),histo12F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
//TH1D* histo2=(TH1D*)histo22F->ProjectionY("histo2",histo22F->GetXaxis()->FindBin(binbound_pt[5]),histo22F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
TH1D* histo1=(TH1D*)histo12F->ProjectionX("histo1");
TH1D* histo2=(TH1D*)histo22F->ProjectionX("histo2");

TH1D* histo1=(TH1D*)file1->Get(histoname);
TH1D* histo2=(TH1D*)file2->Get(histoname);

rehisto1=(TH1D*)histo1->Clone("rehisto1");
rehisto2=(TH1D*)histo2->Clone("rehisto2");
rehisto1=(TH1D*)rehisto1->Rebin(Nbin_pt,"rehisto1",binbound_pt);
rehisto2=(TH1D*)rehisto2->Rebin(Nbin_pt,"rehisto2",binbound_pt);
normalizeByBinWidth(rehisto1);
normalizeByBinWidth(rehisto2);
//rehisto1->Scale(1/rehisto1->Integral());
//rehisto2->Scale(1/rehisto2->Integral());
rehisto1->SetMarkerStyle(24);
rehisto1->SetMarkerSize(1.2);
rehisto1->SetMarkerColor(1);
rehisto2->SetMarkerStyle(20);
rehisto2->SetMarkerSize(1.2);
rehisto2->SetMarkerColor(2);

c1 = new TCanvas("c1"," ",500,800);
makeMultiPanelCanvas(c1,1,2,-0.14,0,0.2,0.14,0.03);
if(gen=="PP")
TLatex *T1=new TLatex(0.5,0.92,"PYTHIA");
else if(gen=="PPb")
TLatex *T1=new TLatex(0.5,0.92,"PYTHIA+HIJING");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,2,3.1);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T}d#eta}");
hFrame->GetXaxis()->SetLimits(30,600);
hFrame->GetYaxis()->SetRangeUser(1.01e-10,1e-1);
c1->cd(1)->SetLogy();
hFrame->DrawCopy();
rehisto1->Draw("E1same");
rehisto2->Draw("E1same");
TLegend *leg1=new TLegend(0.70,0.75,0.9,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(rehisto1,"New","lp");
leg1->AddEntry(rehisto2,"Old","lp");
leg1->Draw("same");
T1->Draw("same");
c1->cd(2);
TLatex *T1=new TLatex(0.5,0.92,"|#eta|<1");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
TH1D* ratio=(TH1D*)rehisto2->Clone("ratio");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Ratio: New/Old");
hFrame->GetYaxis()->SetRangeUser(0.6,1.39);
hFrame->DrawCopy();
ratio->Divide(rehisto1);
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.2);
ratio->SetMarkerColor(1);
ratio->SetTitle("");
ratio->Draw("E1same");
TLine *l =new TLine(30,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T1->Draw("same");
if(gen=="PP")
c1->Print("OldNewComparePP.png");
else if(gen=="PPb")
c1->Print("OldNewComparePPb.png");
}

