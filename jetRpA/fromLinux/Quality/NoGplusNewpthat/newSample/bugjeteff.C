#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void bugjeteff(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

TString filename1="MCPPbNoGpluscW_akPu3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root";
//TString filename2="MCPPbakPu3PFjetana_pPb_Full_BTagForest30_Fix3output.root";
//TString filename2="MCPPbakPu3PFetacutjetana_pPb_Full_BTagForest30_Fix3output.root";
TString filename2="MCPPbakPu3PFjetana_merged.root";
//TString filename2="MCPPbakPu3PFjeteff_pPb_Full_BTagForest30_Fix3output.root";
TString dir1="/scratch/xuq7/RpA/akPu3PF/NoResidual/NoGplusNewpthat/newSample/";
TString dir2="/scratch/xuq7/RpA/akPu3PF/NoResidual/NoGplusNewpthat/newSample/";

TString histoname="pfjet_Cent";
//TString histoname="refjetpt_0-100%";
TString histoname2="rawptJES_0-100%";

TFile *file1=TFile::Open(Form("%s/%s",dir1.Data(),filename1.Data()));
TFile *file2=TFile::Open(Form("%s/%s",dir2.Data(),filename2.Data()));

TH1D* histo1=(TH1D*)file1->Get(histoname);
TH1D* histo2=(TH1D*)file2->Get(histoname);
/*
TH2F* histo12=(TH2F*)file1->Get(histoname2);
TH2F* histo22=(TH2F*)file2->Get(histoname2);
TH1D* histo1=(TH1D*)histo12->ProjectionX();
TH1D* histo2=(TH1D*)histo22->ProjectionX();
histo1=(TH1D*)histo1->Rebin(Nbin_pt,"histo1",binbound_pt);
histo2=(TH1D*)histo2->Rebin(Nbin_pt,"histo2",binbound_pt);
normalizeByBinWidth(histo1);
normalizeByBinWidth(histo2);
*/
//histo1->Scale(1/histo1->Integral());
//histo2->Scale(1/histo2->Integral());
histo1->SetMarkerStyle(24);
histo1->SetMarkerSize(1.2);
histo1->SetMarkerColor(1);
histo1->SetLineColor(1);
histo2->SetMarkerStyle(20);
histo2->SetMarkerSize(1.2);
histo2->SetMarkerColor(2);
histo2->SetLineColor(2);

c1 = new TCanvas("c1"," ",500,800);
makeMultiPanelCanvas(c1,1,2,-0.14,0,0.2,0.14,0.03);
TLatex *T1=new TLatex(0.5,0.92,"PYTHIA+HIJING");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
TH1F* hFrame=new TH1F("","",20000,-1000,1000);
fixedFontHist(hFrame,2,3.1);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T}d#eta}");
hFrame->GetXaxis()->SetRangeUser(0,100);
//hFrame->GetXaxis()->SetRangeUser(-15,15);
//hFrame->GetXaxis()->SetRangeUser(1,1.3);
hFrame->GetYaxis()->SetRangeUser(1e-5,1e-2);
c1->cd(1);//->SetLogy();
hFrame->DrawCopy();
histo1->Draw("E1same");
histo2->Draw("E1same");
TLegend *leg1=new TLegend(0.70,0.75,0.9,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(histo1,"Old","lp");
leg1->AddEntry(histo2,"New","lp");
leg1->Draw("same");
T1->Draw("same");
c1->cd(2);
TLatex *T1=new TLatex(0.5,0.92,"|#eta|<1");
TH1D* ratio=(TH1D*)histo1->Clone("ratio");
hFrame->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Ratio: Old/New");
hFrame->GetYaxis()->SetRangeUser(0.6,1.39);
hFrame->DrawCopy();
ratio->Divide(histo2);
//for(int ibin=0;ibin<ratio->GetNbinsX();ibin++){
//cout<<ratio->GetBinLowEdge(ibin)<<'\t';
//cout<<ratio->GetBinContent(ibin)<<endl;
//}
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.2);
ratio->SetMarkerColor(1);
ratio->SetLineColor(1);
ratio->SetTitle("");
ratio->Draw("E1same");
TLine *l =new TLine(30,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T1->Draw("same");
//c1->Print("OldNewComparePPb.png");
}

