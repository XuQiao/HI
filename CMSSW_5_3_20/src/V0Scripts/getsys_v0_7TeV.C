{
#include <tdrstyle.C>
#include <sstream>
  //#include <iostream>
#include <string>

gROOT->Reset();
setTDRStyle();
//gStyle->SetLabelOffset(0.017, "X");
gStyle->SetStatFormat("5.5g");
gStyle->SetStatY(0.94);
gStyle->SetOptFit(111);
gStyle->SetFitFormat("5.5g");
gStyle->SetStatFontSize(0.03);
gStyle->SetTitleSize(0.05, "XYZ");
gStyle->SetLabelSize(0.04, "XYZ");
gStyle->SetTitleYOffset(1.4);
gStyle->SetStatH(0.16);
gStyle->SetStatW(0.3);
gStyle->SetTextSize(0.038);
gStyle->SetTextFont(42);
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);
//gStyle->SetStatFormat("NELU",AutoPrecision(3));


bool isMC = false;
bool kshort = false;

TString ksSysNames[6] = {"tuneD6T_7TeV_weight101", "tuneP0_7TeV_weight103", "pythia8_7TeV_weight102_fullwgt", 
			 "pythia8_7TeV_weight2", "pythia8_7TeV_weight102_trigBSC", "pythia8_7TeV_weight102_SDp50" };
TString lamSysNames[7] = {"tuneD6T_7TeV_weight101_subgeant2", "tuneP0_7TeV_weight103_subgeant2", "pythia8_7TeV_weight102_subgeant2_fullwgt", 
			  "pythia8_7TeV_weight2_subgeant2", "pythia8_7TeV_weight102_subgeant2_trigBSC", "pythia8_7TeV_weight102_subgeant2_SDp50", 
			  "pythia8_7TeV_weight102_subgeant"};
TString inHist_pt;
TString inFileBase_pt;
TString inFileSys_pt[10];
TString outFileSysPre_pt[10];

TString inHist_y;
TString inFileBase_y;
TString inFileSys_y[10];
TString outFileSysPre_y[10];

int numSys = 0;

if ( kshort ) {
  numSys = 6;
  inHist_pt = "ks_dNdpT_scale";
  inFileBase_pt = "plots/kspt_cor_data_pythia8_7TeV_weight102.root";
  for ( int i = 0; i < numSys; ++i ) {
    inFileSys_pt[i] = "plots/kspt_cor_data_" + ksSysNames[i] + ".root";
    outFileSysPre_pt[i] = "final_plots/kspt_sys_" + ksSysNames[i];
  }
  inHist_y = "ks_dNdy";
  inFileBase_y = "plots/ksy_cor_data_pythia8_7TeV_weight102.root";
  for ( int i = 0; i < numSys; ++i ) {
    inFileSys_y[i] = "plots/ksy_cor_data_" + ksSysNames[i] + ".root";
    outFileSysPre_y[i] = "final_plots/ksy_sys_" + ksSysNames[i];
  }
} else {
  numSys = 7;
  inHist_pt = "lam_dNdpT_scale";
  inFileBase_pt = "plots/lampt_cor_data_pythia8_7TeV_weight102_subgeant2.root";
  for ( int i = 0; i < numSys; ++i ) {
    inFileSys_pt[i] = "plots/lampt_cor_data_" + lamSysNames[i] + ".root";
    outFileSysPre_pt[i] = "final_plots/lampt_sys_" + lamSysNames[i];
  }
  inHist_y = "lam_dNdy";
  inFileBase_y = "plots/lamy_cor_data_pythia8_7TeV_weight102_subgeant2.root";
  for ( int i = 0; i < numSys; ++i ) {
    inFileSys_y[i] = "plots/lamy_cor_data_" + lamSysNames[i] + ".root";
    outFileSysPre_y[i] = "final_plots/lamy_sys_" + lamSysNames[i];
  }
}


TFile* fin = TFile::Open(inFileBase_pt, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* v0_dNdpT_cor_data_base = (TH1F*)c1->GetPrimitive(inHist_pt);
fin->Close();
delete c1;
delete fin;

TFile* fin = TFile::Open(inFileBase_y, "READ");
TCanvas *c1 = (TCanvas*)fin.Get("c1");
TH1F* v0_dNdy_cor_data_base = (TH1F*)c1->GetPrimitive(inHist_y);
fin->Close();
delete c1;
delete fin;


for (int i = 0; i < numSys; ++i ) {
  TFile* fin = TFile::Open(inFileSys_pt[i], "READ");
  TCanvas *c1 = (TCanvas*)fin.Get("c1");
  TH1F* v0_dNdpT_cor_data_sys = (TH1F*)c1->GetPrimitive(inHist_pt);
  fin->Close();
  delete c1;
  delete fin;
  v0_dNdpT_cor_data_sys->Divide(v0_dNdpT_cor_data_base);
  v0_dNdpT_cor_data_sys->GetXaxis()->SetTitle("V^{0}  p_{T}  [GeV/c]");
  v0_dNdpT_cor_data_sys->GetYaxis()->SetTitle("Systematic / Default");
  v0_dNdpT_cor_data_sys->SetAxisRange(0.8, 1.2, "Y");
  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");
  v0_dNdpT_cor_data_sys->Draw();
  c1->SaveAs(outFileSysPre_pt[i] + ".root");
  c1->SaveAs(outFileSysPre_pt[i] + ".eps");
  c1->SaveAs(outFileSysPre_pt[i] + ".png");

  TFile* fin = TFile::Open(inFileSys_y[i], "READ");
  TCanvas *c1 = (TCanvas*)fin.Get("c1");
  TH1F* v0_dNdy_cor_data_sys = (TH1F*)c1->GetPrimitive(inHist_y);
  fin->Close();
  delete c1;
  delete fin;
  v0_dNdy_cor_data_sys->Divide(v0_dNdy_cor_data_base);
  v0_dNdy_cor_data_sys->GetXaxis()->SetTitle("V^{0}  rapidity");
  v0_dNdy_cor_data_sys->GetYaxis()->SetTitle("Systematic / Default");
  v0_dNdy_cor_data_sys->SetAxisRange(0.8, 1.2, "Y");
  TCanvas* c1 = new TCanvas("c1", "V0 y distribution");
  v0_dNdy_cor_data_sys->Draw();
  c1->SaveAs(outFileSysPre_y[i] + ".root");
  c1->SaveAs(outFileSysPre_y[i] + ".eps");
  c1->SaveAs(outFileSysPre_y[i] + ".png");
}


}



