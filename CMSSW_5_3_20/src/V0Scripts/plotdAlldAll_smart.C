{
#include <tdrstyle.C>
#include <sstream>
//#include <iostream>
#include <iomanip.h>
#include <string>
#include <vector>

#include "ksfit.C"
#include "lamfit.C"
  //#include "kshortrec.C"
#include "constants.h"


  gSystem->Load("kshortrec_C.so");

  vector<double> ksptBinMinima;
  vector<double> lamptBinMinima;
  vector<double> kspBinMinima;
  vector<double> lampBinMinima;
  vector<double> kslifeBinMinima;
  vector<double> lamlifeBinMinima;
  vector<double> ntrkBinMinima;
  
  for ( int ndx = 0; ndx <= ksptNbins; ndx++ ) {
    ksptBinMinima.push_back( ksptarray[ndx] );
  }
  
  for ( int ndx = 0; ndx <= lamptNbins; ndx++ ) {
    lamptBinMinima.push_back( lamptarray[ndx] );
  }

  for ( int ndx = 0; ndx <= kspNbins; ndx++ ) {
    kspBinMinima.push_back( ksparray[ndx] );
  }

  for ( int ndx = 0; ndx <= lampNbins; ndx++ ) {
    lampBinMinima.push_back( lamparray[ndx] );
  }

  for ( int ndx = 0; ndx <= kslifeNbins; ndx++ ) {
    kslifeBinMinima.push_back( kslifearray[ndx] );
  }

  for ( int ndx = 0; ndx <= lamlifeNbins; ndx++ ) {
    lamlifeBinMinima.push_back( lamlifearray[ndx] );
  }

  for ( int ndx = 0; ndx <= ntrkNbins; ndx++ ) {
    ntrkBinMinima.push_back( ntrkarray[ndx] );
  }

  // Switches
  bool dodNdy = true;
  bool dodNdpT = true;
  bool dodNdlife = true;
  bool dodNdp = true;
  bool dodNdy_ks;
  bool dodNdpT_ks;
  bool dodNdlife_ks;
  bool dodNdy_lam;
  bool dodNdpT_lam;
  bool dodNdlife_lam;
  bool isMC = false;
  bool debug = false;
  bool kshort = true;

  bool dodNdy_ks = kshort && dodNdy;
  bool dodNdpT_ks = kshort && dodNdpT;
  bool dodNdp_ks = kshort && dodNdp;
  bool dodNdlife_ks = kshort && dodNdlife;

  bool dodNdy_lam = !kshort && dodNdy;
  bool dodNdpT_lam = !kshort && dodNdpT;
  bool dodNdp_lam = !kshort && dodNdp;
  bool dodNdlife_lam = !kshort && dodNdlife;

  TString mcname = "pythia8";
  //TString mcname = "tuneD6T";

  // Style and other ROOT commands
  setTDRStyle();
  gStyle->SetOptStat(1);
  TGaxis::SetMaxDigits(3);
  gROOT->SetBatch(kTRUE);

  TString finstr;
  TString outFileName;

  if( !isMC ) {
    finstr = "/nfs/data37/cms/v0ntuple/data_7TeV_june4fixesGood.root";
    //    fin = "../data_7TeV_allNtuples_400Jobs.root";
    if ( kshort ) {
      outFileName = "data_dNdX_kshort.root";
    } 
    else {
      outFileName = "data_dNdX_lambda.root";
    }
  }
  else {
    //fin = "/nfs/data37/cms/v0ntuple/mc_7TeV_" + mcname + "_50M_START3X_V26B_june3fixes_PARTIAL.root";
    finstr = "/nfs/data37/cms/v0ntuple/mc_7TeV_" + mcname + "_START3X_V26B_june4fixesGood.root";
    //fin = "/nfs/data35/cms/v0ntuple/mc_7TeV_tuneD6T_may22fixes_partial_2.root";
    //fin = "/nfs/data35/cms/v0ntuple/mc_7TeV_tuneP0_may22fixes.root";
    //    fin = "/nfs/data35/cms/v0ntuple/mc_500kEvents.root";
    //   fin = "../mc_7TeV_pythia8.root";
    //    fin = "../mc_500kEvents.root";
    if (kshort) {
      outFileName = "sim_dNdX_kshort_" + mcname + ".root";
    } 
    else {
      outFileName = "sim_dNdX_lambda_" + mcname + ".root";
    }
  }
  TFile* fout = TFile::Open(outFileName, "RECREATE");
  TFile* fin = TFile::Open(finstr, "READ");


  //cout << "Files are open." << endl;
  // Constants
  double piMass = 0.13957018;
  double protonMass = 0.93827203;
  double ksMassConst = 0.497614;
  double lamMassConst = 1.115683;

  int fit_command = 0;

  // TPaves for labels
  TText* labltxt1_1;
  TText* labltxt1_2;
  TText* labltxt2_1;
  TText* labltxt2_2;

  TPaveText* labl = new TPaveText(0.25, 0.75, 0.55, 0.88, "brNDC");
  labl->SetBorderSize(0);
  labl->SetFillColor(0);
  labl->SetFillStyle(0);
  labl->SetTextSizePixels(20);
  TPaveText* labl2 = new TPaveText(0.25, 0.79, 0.55, 0.92, "brNDC");
  labl2->SetBorderSize(0);
  labl2->SetFillColor(0);
  labl2->SetFillStyle(0);
  labl2->SetTextSizePixels(20);

  if( isMC ) {
    labltxt1_1 = labl->AddText("CMS Simulation " + mcname + "");
    labltxt1_2 = labl->AddText("#sqrt{s} = 7 TeV");
    labltxt2_1 = labl2->AddText("CMS Simulation " + mcname + "");
    labltxt2_2 = labl2->AddText("#sqrt{s} = 7 TeV");
  }
  else {
    labltxt1_1 = labl->AddText("CMS Preliminary");
    labltxt1_2 = labl->AddText("#sqrt{s} = 7 TeV");
    labltxt2_1 = labl2->AddText("CMS Preliminary");
    labltxt2_2 = labl2->AddText("#sqrt{s} = 7 TeV");
  }
  /*
  // Histogram quantities
  Double_t ksMassXmin = 0.340;
  Double_t ksMassXmax = 0.740;
  int ksMassNbins = 400;
  Double_t ksMassBinWidth = (ksMassXmax - ksMassXmin) / ksMassNbins;
  Double_t ksMassBinWidth_DG = ksMassBinWidth / sqrt(2*TMath::Pi());

  Double_t ksMassXmin_narrow = 0.440;
  Double_t ksMassXmax_narrow = 0.640;
  int ksMassNbins_narrow = 200;

  Double_t lamMassXmin = 1.088;
  Double_t lamMassXmax = 1.188;
  int lamMassNbins = 200;
  Double_t lamMassBinWidth = (lamMassXmax - lamMassXmin) / lamMassNbins;
  Double_t lamMassBinWidth_DG = lamMassBinWidth / sqrt(2*TMath::Pi());

  Double_t lamMassXmin_narrow = 1.09;
  Double_t lamMassXmax_narrow = 1.15;
  int lamMassNbins_narrow = 120;

  Double_t ksyXmin = 0.;//-2.5;
  Double_t ksyXmax = 2.4;
  int ksyNbins = 23;
  Double_t ksyBinWidth = (ksyXmax - ksyXmin) / ksyNbins;

  Double_t ksyXmin_signed = -2.4;//-2.5;
  Double_t ksyXmax_signed = 2.4;
  int ksyNbins_signed = 46;
  Double_t ksyBinWidth_signed = (ksyXmax_signed - ksyXmin_signed) / ksyNbins_signed;

  Double_t lamyXmin = 0.;
  Double_t lamyXmax = 2.4;
  int lamyNbins = 23;
  Double_t lamyBinWidth = (lamyXmax - lamyXmin) / lamyNbins;

  Double_t lamyXmin_signed = -2.4;
  Double_t lamyXmax_signed = 2.4;
  int lamyNbins_signed = 46;
  Double_t lamyBinWidth_signed = (lamyXmax_signed - lamyXmin_signed) / lamyNbins;

  double ksptarray[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0,10.0};
  int ksptNbins = sizeof(ksptarray)/sizeof(double) -1;
  vector<Double_t> ksptBinMinima;
  for ( int ndx = 0; ndx <= ksptNbins; ndx++ ) {
    ksptBinMinima.push_back( ksptarray[ndx] );
  }
  Double_t kspt_default_binWidth = 0.1;

  double lamptarray[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0,10.0};
  int lamptNbins = sizeof(lamptarray)/sizeof(double) -1;
  vector<Double_t> lamptBinMinima;
  for ( int ndx = 0; ndx <= lamptNbins; ndx++ ) {
    lamptBinMinima.push_back( lamptarray[ndx] );
  }
  Double_t lampt_default_binWidth = 0.1;

  double ksparray[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.4,5.8,6.2,6.6,7.0,7.4,7.8,8.2,8.6,9.0,10.0,11.0,12.0,13.0,14.0,16.0,18.0,20.0};
  int kspNbins = sizeof(ksparray)/sizeof(double) -1;
  vector<Double_t> kspBinMinima;
  for ( int ndx = 0; ndx <= kspNbins; ndx++ ) {
    kspBinMinima.push_back( ksparray[ndx] );
  }
  Double_t ksp_default_binWidth = 0.2;

  double lamparray[] = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.4,5.8,6.2,6.6,7.0,7.4,7.8,8.2,8.6,9.0,10.0,11.0,12.0,13.0,14.0,16.0,18.0,20.0};
  int lampNbins = sizeof(lamparray)/sizeof(double) -1;
  vector<Double_t> lampBinMinima;
  for ( int ndx = 0; ndx <= lampNbins; ndx++ ) {
    lampBinMinima.push_back( lamparray[ndx] );
  }
  Double_t lamp_default_binWidth = 0.2;

  double kslifearray[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.4,6.8,7.2,7.6,8.0,8.8,9.6,10.4,11.2,12.0};
  int kslifeNbins = sizeof(kslifearray)/sizeof(double) -1;
  vector<Double_t> kslifeBinMinima;
  for ( int ndx = 0; ndx <= kslifeNbins; ndx++ ) {
    kslifeBinMinima.push_back( kslifearray[ndx] );
  }
  Double_t kslife_default_binWidth = 0.1;

  double lamlifearray[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,6.4,6.8,7.2,7.6,8.0,8.8,9.6,10.4,11.2,12.0,13.0,14.0,15.0,16.0,18.0,20.0};
  int lamlifeNbins = sizeof(lamlifearray)/sizeof(double) -1;
  vector<Double_t> lamlifeBinMinima;
  for ( int ndx = 0; ndx <= lamlifeNbins; ndx++ ) {
    lamlifeBinMinima.push_back( lamlifearray[ndx] );
  }
  Double_t lamlife_default_binWidth = 0.1;

  double ntrkarray[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,84,88,92,96,100,110,120,140,160,200};
  int ntrkNbins = sizeof(ntrkarray)/sizeof(double) -1;
  vector<double> ntrkBinMinima;
  for ( int ndx = 0; ndx <= ntrkNbins; ndx++ ) {
    ntrkBinMinima.push_back( ntrkarray[ndx] );
  }
  Double_t ntrk_default_binWidth = 0.1;
  */
  TH1F* ntrk = new TH1F("ntrk", "Track multiplicity", ntrkNbins, &ntrkBinMinima[0]);
  TH1F* ntrk_scale = new TH1F("ntrk_scale", "Track multiplicity", ntrkNbins, &ntrkBinMinima[0]);
  TH1F* trkscale = new TH1F("trkscale", "ntrk scale factor", ntrkNbins, &ntrkBinMinima[0]);
  for( int ndx = 0; ndx < ntrkNbins; ndx++ ) {
    double binWidth = ntrkBinMinima[ndx+1] - ntrkBinMinima[ndx];
    double scalingFactor = binWidth / ntrk_default_binWidth;
    trkscale->SetBinContent( ndx + 1, scalingFactor );
    trkscale->SetBinError( ndx + 1, 0. );
  }
  
  if (kshort) {
    TH1F* vee_nocuts = new TH1F("vee_nocuts", "vee_nocuts", ksMassNbins, ksMassXmin, ksMassXmax);
    TH1F* vee_withcuts = new TH1F("vee_withcuts", "vee_withcuts", ksMassNbins, ksMassXmin, ksMassXmax);
    TH1F* veeReg_nocuts = 0;
    TH1F* veeBar_nocuts = 0;
    TH1F* veeReg_withcuts = 0;
    TH1F* veeBar_withcuts = 0;
  } 
  else {
    TH1F* vee_nocuts = new TH1F("vee_nocuts", "vee_nocuts", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* vee_withcuts = new TH1F("vee_withcuts", "vee_withcuts", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* veeReg_nocuts = new TH1F("veeReg_nocuts", "veeReg_nocuts", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* veeReg_withcuts = new TH1F("veeReg_withcuts", "veeReg_withcuts", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* veeBar_nocuts = new TH1F("veeBar_nocuts", "veeBar_nocuts", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* veeBar_withcuts = new TH1F("veeBar_withcuts", "veeBar_withcuts", lamMassNbins, lamMassXmin, lamMassXmax);
  }

  //cout << "Creating histos..." << endl;
  // Create all of the Kshort y mass histograms
  TH1F** ksPtrs;
  ksPtrs = new TH1F*[ksyNbins];
  for( int ndx = 0; ndx < ksyNbins; ndx++ ) {
    double min = ksyXmin + ndx*ksyBinWidth;
    double max = ksyXmin + (ndx+1)*ksyBinWidth;
    ostringstream tempHistoName;
    tempHistoName << "ksM_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(), tempHistoName.str().c_str(), ksMassNbins, ksMassXmin, ksMassXmax);
    ksPtrs[ndx] = tempHist;
    ksPtrs[ndx]->SetXTitle("m_{#pi^{+}#pi^{-}} [GeV/c^{2}]");
    ksPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    ksPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** ksPtrs_signed;
  ksPtrs_signed = new TH1F*[ksyNbins_signed];
  for( int ndx = 0; ndx < ksyNbins_signed; ndx++ ) {
    double min = ksyXmin_signed + ndx*ksyBinWidth_signed;
    ostringstream tempHistoName;
    tempHistoName << "ksMsgn_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(), tempHistoName.str().c_str(), ksMassNbins, ksMassXmin, ksMassXmax);
    ksPtrs_signed[ndx] = tempHist;
    ksPtrs_signed[ndx]->SetXTitle("m_{#pi^{+}#pi^{-}} [GeV/c^{2}]");
    ksPtrs_signed[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    ksPtrs_signed[ndx]->SetNdivisions(506);
  }

  TH1F** lamPtrs;
  lamPtrs = new TH1F*[lamyNbins];
  for( int ndx = 0; ndx < lamyNbins; ndx++ ) {
    double min = lamyXmin + ndx*lamyBinWidth;
    double max = lamyXmin + (ndx+1)*lamyBinWidth;
    ostringstream tempHistoName;
    tempHistoName << "lamM_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(), tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    lamPtrs[ndx] = tempHist;
    lamPtrs[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** lamPtrs_signed;
  lamPtrs_signed = new TH1F*[lamyNbins_signed];
  for( int ndx = 0; ndx < lamyNbins_signed; ndx++ ) {
    double min = lamyXmin_signed + ndx*lamyBinWidth_signed;
    double max = lamyXmin_signed + (ndx+1)*lamyBinWidth_signed;
    ostringstream tempHistoName;
    tempHistoName << "lamMsgn_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(), tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    lamPtrs_signed[ndx] = tempHist;
    lamPtrs_signed[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamPtrs_signed[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamPtrs_signed[ndx]->SetNdivisions(506);
  }

  TH1F** lamRegPtrs;
  lamRegPtrs = new TH1F*[lamyNbins];
  for( int ndx = 0; ndx < lamyNbins; ndx++ ) {
    double min = lamyXmin + ndx*lamyBinWidth;
    double max = lamyXmin + (ndx+1)*lamyBinWidth;
    ostringstream tempHistoName;
    tempHistoName << "lamRegM_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(), tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    lamRegPtrs[ndx] = tempHist;
    lamRegPtrs[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamRegPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamRegPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** lamRegPtrs_signed;
  lamRegPtrs_signed = new TH1F*[lamyNbins_signed];
  for( int ndx = 0; ndx < lamyNbins_signed; ndx++ ) {
    double min = lamyXmin_signed + ndx*lamyBinWidth_signed;
    double max = lamyXmin_signed + (ndx+1)*lamyBinWidth_signed;
    ostringstream tempHistoName;
    tempHistoName << "lamRegMsgn_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(), tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    lamRegPtrs_signed[ndx] = tempHist;
    lamRegPtrs_signed[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamRegPtrs_signed[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamRegPtrs_signed[ndx]->SetNdivisions(506);
  }

  TH1F** lamBarPtrs;
  lamBarPtrs = new TH1F*[lamyNbins];
  for( int ndx = 0; ndx < lamyNbins; ndx++ ) {
    double min = lamyXmin + ndx*lamyBinWidth;
    double max = lamyXmin + (ndx+1)*lamyBinWidth;
    ostringstream tempHistoName;
    tempHistoName << "lamBarM_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(), tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    lamBarPtrs[ndx] = tempHist;
    lamBarPtrs[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamBarPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamBarPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** lamBarPtrs_signed;
  lamBarPtrs_signed = new TH1F*[lamyNbins_signed];
  for( int ndx = 0; ndx < lamyNbins_signed; ndx++ ) {
    double min = lamyXmin_signed + ndx*lamyBinWidth_signed;
    double max = lamyXmin_signed + (ndx+1)*lamyBinWidth_signed;
    ostringstream tempHistoName;
    tempHistoName << "lamBarMsgn_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(), tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    lamBarPtrs_signed[ndx] = tempHist;
    lamBarPtrs_signed[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamBarPtrs_signed[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamBarPtrs_signed[ndx]->SetNdivisions(506);
  }

  TH1F** ksptPtrs;
  ksptPtrs = new TH1F*[ksptNbins];
  for( int ndx = 0; ndx < ksptNbins; ndx++ ) {
    double min = ksptBinMinima[ndx];
    double max = ksptBinMinima[ndx+1];
    double binWidth = max - min;
    double scalingFactor = binWidth / kspt_default_binWidth;
    ostringstream tempHistoName;
    tempHistoName << "ksMpT_" << min << "_" << max;
    if ( ndx == 0 ) {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), 
				ksMassNbins_narrow, ksMassXmin_narrow, ksMassXmax_narrow);
    } else {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), ksMassNbins, ksMassXmin, ksMassXmax);
    }
    ksptPtrs[ndx] = tempHist;
    ksptPtrs[ndx]->SetXTitle("m_{#pi^{+}#pi^{-}} [GeV/c^{2}]");
    ksptPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    ksptPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** lamptPtrs;
  lamptPtrs = new TH1F*[lamptNbins];
  for( int ndx = 0; ndx < lamptNbins; ndx++ ) {
    ostringstream lamSelection;
    double min = lamptBinMinima[ndx];
    double max = lamptBinMinima[ndx+1];
    double binWidth = max - min;
    double scalingFactor = binWidth / lampt_default_binWidth;
    ostringstream tempHistoName;
    tempHistoName << "lamMpT_" << min << "_" << max;
    if ( ndx == 0 ) {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), 
				lamMassNbins_narrow, lamMassXmin_narrow, lamMassXmax_narrow);
    } else {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    }
    lamptPtrs[ndx] = tempHist;
    lamptPtrs[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamptPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamptPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** lamRegptPtrs;
  lamRegptPtrs = new TH1F*[lamptNbins];
  for( int ndx = 0; ndx < lamptNbins; ndx++ ) {
    ostringstream lamSelection;
    double min = lamptBinMinima[ndx];
    double max = lamptBinMinima[ndx+1];
    double binWidth = max - min;
    double scalingFactor = binWidth / lampt_default_binWidth;
    ostringstream tempHistoName;
    tempHistoName << "lamRegMpT_" << min << "_" << max;
    if ( ndx == 0 ) {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), 
				lamMassNbins_narrow, lamMassXmin_narrow, lamMassXmax_narrow);
    } else {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    }
    lamRegptPtrs[ndx] = tempHist;
    lamRegptPtrs[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamRegptPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamRegptPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** lamBarptPtrs;
  lamBarptPtrs = new TH1F*[lamptNbins];
  for( int ndx = 0; ndx < lamptNbins; ndx++ ) {
    ostringstream lamSelection;
    double min = lamptBinMinima[ndx];
    double max = lamptBinMinima[ndx+1];
    double binWidth = max - min;
    double scalingFactor = binWidth / lampt_default_binWidth;
    ostringstream tempHistoName;
    tempHistoName << "lamBarMpT_" << min << "_" << max;
    if ( ndx == 0 ) {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), 
				lamMassNbins_narrow, lamMassXmin_narrow, lamMassXmax_narrow);
    } else {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    }
    lamBarptPtrs[ndx] = tempHist;
    lamBarptPtrs[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lamBarptPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamBarptPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** kspPtrs;
  kspPtrs = new TH1F*[kspNbins];
  for( int ndx = 0; ndx < kspNbins; ndx++ ) {
    double min = kspBinMinima[ndx];
    double max = kspBinMinima[ndx+1];
    double binWidth = max - min;
    double scalingFactor = binWidth / ksp_default_binWidth;
    ostringstream tempHistoName;
    tempHistoName << "ksMp_" << min << "_" << max;
    if ( ndx == 0 ) {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), 
				ksMassNbins_narrow, ksMassXmin_narrow, ksMassXmax_narrow);
    } else {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), ksMassNbins, ksMassXmin, ksMassXmax);
    }
    kspPtrs[ndx] = tempHist;
    kspPtrs[ndx]->SetXTitle("m_{#pi^{+}#pi^{-}} [GeV/c^{2}]");
    kspPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    kspPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** lampPtrs;
  lampPtrs = new TH1F*[lampNbins];
  for( int ndx = 0; ndx < lampNbins; ndx++ ) {
    ostringstream lamSelection;
    double min = lampBinMinima[ndx];
    double max = lampBinMinima[ndx+1];
    double binWidth = max - min;
    double scalingFactor = binWidth / lamp_default_binWidth;
    ostringstream tempHistoName;
    tempHistoName << "lamMp_" << min << "_" << max;
    if ( ndx == 0 ) {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), 
				lamMassNbins_narrow, lamMassXmin_narrow, lamMassXmax_narrow);
    } else {
      TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    }
    lampPtrs[ndx] = tempHist;
    lampPtrs[ndx]->SetXTitle("m_{p#pi^{-}}(+c.c) [GeV/c^{2}]");
    lampPtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lampPtrs[ndx]->SetNdivisions(506);
  }

  TH1F** kslifePtrs;
  kslifePtrs = new TH1F*[kslifeNbins];
  for( int ndx = 0; ndx < kslifeNbins; ndx++ ) {
    double min = kslifeBinMinima[ndx];
    double max = kslifeBinMinima[ndx+1];
    double binWidth = max - min;
    double scalingFactor = binWidth / kslife_default_binWidth;
    ostringstream tempHistoName;
    tempHistoName << "ksMlife_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), ksMassNbins, ksMassXmin, ksMassXmax);
    kslifePtrs[ndx] = tempHist;
    kslifePtrs[ndx]->SetXTitle("m_{#pi^{+}#pi^{-}} [GeV/c^{2}]");
    kslifePtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    kslifePtrs[ndx]->SetNdivisions(506);
  }

  TH1F** lamlifePtrs;
  lamlifePtrs = new TH1F*[lamlifeNbins];
  for( int ndx = 0; ndx < lamlifeNbins; ndx++ ) {
    double min = lamlifeBinMinima[ndx];
    double max = lamlifeBinMinima[ndx+1];
    double binWidth = max - min;
    double scalingFactor = binWidth / lamlife_default_binWidth;
    ostringstream tempHistoName;
    tempHistoName << "lamMlife_" << min << "_" << max;
    TH1F* tempHist = new TH1F(tempHistoName.str().c_str(),tempHistoName.str().c_str(), lamMassNbins, lamMassXmin, lamMassXmax);
    lamlifePtrs[ndx] = tempHist;
    lamlifePtrs[ndx]->SetXTitle("m_{p^{+}#pi^{-}} [GeV/c^{2}]");
    lamlifePtrs[ndx]->SetYTitle("Candidates / 0.001 GeV/c^{2}");
    lamlifePtrs[ndx]->SetNdivisions(506);
  }


  // Find number of events.
  TH1F* nV0sH = new TH1F("nV0sH", "nV0sH", 100, 0., 100);
  TH1F* numvee = new TH1F("numvee", "numvee", 100, 0., 100);

  // kshortData t;
  // if ( isMC ) {
  //   kshortrec t;
  // } else {
  //   kshortrec t;
  // }
  //cout << "Kshortrec" << endl;
  TTree* theTree;
  if( kshort ) {
    fin->cd("analyzeKshort");
    theTree = (TTree*)gDirectory->Get("ntuple");
  }
  else {
    fin->cd("analyzeLambda");
    theTree = (TTree*)gDirectory->Get("ntuple");
  }
  kshortrec *t = new kshortrec(theTree);
  t->Loop(kshort, nV0sH, numvee, ntrk, vee_nocuts, veeReg_nocuts, veeBar_nocuts,
	 vee_withcuts, veeReg_withcuts, veeBar_withcuts,
	 ksPtrs, ksPtrs_signed, ksptPtrs, kspPtrs, kslifePtrs,
	 lamPtrs, lamRegPtrs, lamBarPtrs, lamptPtrs, lamRegptPtrs, lamBarptPtrs,
	 lampPtrs, lamlifePtrs);

  fout->cd();

  cout << "Finished loop" << endl;
  ntrk->Sumw2();
  ntrk_scale->Divide(ntrk,trkscale);
  TCanvas* c1 = new TCanvas("c1", "Track distribution");
  gStyle->SetOptTitle(1);
  fout->cd();
  ntrk_scale->SetXTitle("Number of tracks");
  ntrk_scale->SetYTitle("Events / track");
  ntrk_scale->SetNdivisions(506);
  ntrk_scale->SetMinimum(0);
  ntrk_scale->Draw();
  labl->Draw("same");
  if( isMC ) {
    c1->SaveAs("ntrk_MC_" + mcname + ".png");
    c1->SaveAs("ntrk_MC_" + mcname + ".eps");
    c1->SaveAs("ntrk_MC_" + mcname + ".root");
    c1->SaveAs("ntrk_MC_" + mcname + ".C");
  }
  else {
    c1->SaveAs("ntrk.png");
    c1->SaveAs("ntrk.eps");
    c1->SaveAs("ntrk.root");
    c1->SaveAs("ntrk.C");
  }
  delete c1;
 
  Double_t nEvts = nV0sH->GetEntries();
  cout << "Precision: " << cout.precision() << endl;
  cout.precision(0);
  cout.setf(ios::fixed);
  cout << "Running over " << nEvts << " events..." << endl;
  cout.unsetf(ios::fixed);
  cout.precision(6);


  double avgsigma;

  gStyle->SetOptTitle(1);

  // Now we only need to create this object once, and we can refit with it
  //   as many times as we like.  Just pass the TH1F* to the Fit method.
  cout << "Creating..." << endl;
  ksfit* ksfitobj = new ksfit();
  lamfit* lamfitobj = new lamfit();
  cout << "Created." << endl;
  if (kshort) {
    cout << "fitting full Kshort histogram\n";
    int gfitresult = ksfitobj->Fit(vee_withcuts);
    cout << "Total yield = " << ksfitobj->sig() << " ± " << ksfitobj->sig_err() << endl;
  } 
  else {
    cout << "fitting full Lambda histogram\n";
    int gfitresult = lamfitobj->Fit(vee_withcuts);
    cout << "Total yield = " << lamfitobj->sig() << " ± " << lamfitobj->sig_err() << endl;
  }

  if( dodNdy_ks ) {
    cout << "Doing dN/dy..." << endl;

    int ksCanvasX = 8;
    int ksCanvasY = ksyNbins / ksCanvasX;
    if( (ksyNbins)%ksCanvasX != 0 ) ksCanvasY += 1;
    TCanvas* ksyCanvas = new TCanvas("ksyFits", "K0S mass fits in y bins", 
				       300*ksCanvasX, 300*ksCanvasY);
    ksyCanvas->Divide( ksCanvasX, ksCanvasY );
    
    TH1F* ks_dNdy = new TH1F("ks_dNdy", "dN/dy", ksyNbins, ksyXmin, ksyXmax);
    TH1F* ks_dNdy_neg = new TH1F("ks_dNdy_neg", "dN/dy", 
				 ksyNbins, 0., ksyXmax_signed);
    TH1F* ks_dNdy_pos = new TH1F("ks_dNdy_pos", "dN/dy",
				 ksyNbins, 0., ksyXmax_signed);
    TH1F* ks_dStNdy = new TH1F("ks_dStNdy", "(S/N) / dy", ksyNbins, ksyXmin, ksyXmax);
    TH1F* ks_sigdy = new TH1F("ks_sigdy", "#sigma / dy", ksyNbins, ksyXmin, ksyXmax);
    TH1F* ks_massdy = new TH1F("ks_massdy", "Mass / dy", ksyNbins, ksyXmin, ksyXmax);
    TH1F* ks_fitdy = new TH1F("ks_fitdy", "Fit / dy", ksyNbins, ksyXmin, ksyXmax);
    TH1F* ks_chi2dy = new TH1F("ks_chi2dy", "#chi^{2}/dof / dy", ksyNbins, ksyXmin, ksyXmax);

    for( int ndx = 0; ndx < ksyNbins; ndx++ ) {
      ksyCanvas->cd(ndx + 1);
      double min = ksyXmin + ndx*ksyBinWidth;
      double max = ksyXmin + (ndx+1)*ksyBinWidth;
      cout << "Doing fit for K0S with " << min << " < y < " << max << endl;

      int fitresult = ksfitobj->Fit(ksPtrs[ndx]);
      cout << "Fit result = " << fitresult << endl;

      ks_dNdy->SetBinContent( ndx + 1, ksfitobj->sig() );
      ks_dNdy->SetBinError( ndx + 1, ksfitobj->sig_err() );
      double ston =  0.95*ksfitobj->sig()/ksfitobj->bg();
      double stone = ston*sqrt( pow(ksfitobj->sig_err() / ksfitobj->sig(),2)
				+ pow(ksfitobj->bg_err() / ksfitobj->bg(),2) );
      ks_dStNdy->SetBinContent( ndx + 1, ston );
      ks_dStNdy->SetBinError( ndx + 1, stone );
      ks_sigdy->SetBinContent( ndx + 1, ksfitobj->res() );
      ks_sigdy->SetBinError( ndx + 1, ksfitobj->res_err() );
      ks_massdy->SetBinContent( ndx + 1, ksfitobj->mass() );
      ks_massdy->SetBinError( ndx + 1, ksfitobj->mass_err() );
      ks_fitdy->SetBinContent( ndx + 1, ksfitobj->fitbg() );
      ks_chi2dy->SetBinContent( ndx + 1, ksfitobj->nChi2() );

      ksPtrs[ndx]->Draw();
      cout << "---------------" << endl;
    }
    if( isMC ) {
      ksyCanvas->SaveAs("plots/KsMassFitsInyBins_MC_" + mcname + ".png");
      ksyCanvas->SaveAs("plots/KsMassFitsInyBins_MC_" + mcname + ".eps");
      ksyCanvas->SaveAs("plots/KsMassFitsInyBins_MC_" + mcname + ".root");
      ksyCanvas->SaveAs("plots/KsMassFitsInyBins_MC_" + mcname + ".pdf");
    }
    else {
      ksyCanvas->SaveAs("plots/KsMassFitsInyBins.png");
      ksyCanvas->SaveAs("plots/KsMassFitsInyBins.eps");
      ksyCanvas->SaveAs("plots/KsMassFitsInyBins.root");
      ksyCanvas->SaveAs("plots/KsMassFitsInyBins.pdf");
    }
    delete ksyCanvas;

    ksCanvasY = ksyNbins_signed / ksCanvasX;
    if( (ksyNbins_signed)%ksCanvasX != 0 ) ksCanvasY += 1;
    ksyCanvas = new TCanvas("ksyFits_signed", "K0S mass fits in y bins, signed",
			    300*ksCanvasX, 300*ksCanvasY);
    ksyCanvas->Divide( ksCanvasX, ksCanvasY );

    for( int ndx = 0; ndx < ksyNbins_signed; ndx++ ) {
      ksyCanvas->cd(ndx + 1);
      double min = ksyXmin_signed + ndx*ksyBinWidth_signed;
      double max = ksyXmin + (ndx+1)*ksyBinWidth_signed;
      cout << "Doing for for K0S with " << min << " < y < " << max << endl;

      int fitresult = ksfitobj->Fit(ksPtrs_signed[ndx]);
      if( ndx < ksyNbins ) {
	ks_dNdy_neg->SetBinContent( (ksyNbins-ndx), fabs(ksfitobj->sig()) );
	ks_dNdy_neg->SetBinError( (ksyNbins-ndx), ksfitobj->sig_err() );
      }
      else {
	ks_dNdy_pos->SetBinContent( (ndx-ksyNbins)+1, ksfitobj->sig() );
	ks_dNdy_pos->SetBinError( (ndx-ksyNbins)+1, ksfitobj->sig_err() );
      }
      ksPtrs_signed[ndx]->Draw();
    }
    if( isMC ) {
      ksyCanvas->SaveAs("plots/KsMassFitsInSignedyBins_MC_" + mcname + ".png");
      ksyCanvas->SaveAs("plots/KsMassFitsInSignedyBins_MC_" + mcname + ".eps");
      ksyCanvas->SaveAs("plots/KsMassFitsInSignedyBins_MC_" + mcname + ".root");
      ksyCanvas->SaveAs("plots/KsMassFitsInSignedyBins_MC_" + mcname + ".pdf");
    }
    else {
      ksyCanvas->SaveAs("plots/KsMassFitsInSignedyBins.png");
      ksyCanvas->SaveAs("plots/KsMassFitsInSignedyBins.eps");
      ksyCanvas->SaveAs("plots/KsMassFitsInSignedyBins.root");
      ksyCanvas->SaveAs("plots/KsMassFitsInSignedyBins.pdf");
    }
    delete ksyCanvas;

    TCanvas* c1 = new TCanvas("c1", "K0S yields");

    gStyle->SetOptTitle(0);

    fout->cd();

    ks_dNdy_neg->SetXTitle("K^{0}_{S} rapidity");
    ks_dNdy_neg->SetYTitle("Candidates / 0.1");
    ks_dNdy_neg->SetNdivisions(506);
    ks_dNdy_neg->SetMinimum(0);
    ks_dNdy_neg->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dNdy_neg_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dNdy_neg_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dNdy_neg_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dNdy_neg_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dNdy_neg.png");
      c1->SaveAs("plots/ks_dNdy_neg.eps");
      c1->SaveAs("plots/ks_dNdy_neg.root");
      c1->SaveAs("plots/ks_dNdy_neg.C");
    }

    ks_dNdy_pos->SetXTitle("K^{0}_{S} rapidity");
    ks_dNdy_pos->SetYTitle("Candidates / 0.1");
    ks_dNdy_pos->SetNdivisions(506);
    ks_dNdy_pos->SetMinimum(0);
    ks_dNdy_pos->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dNdy_pos_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dNdy_pos_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dNdy_pos_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dNdy_pos_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dNdy_pos.png");
      c1->SaveAs("plots/ks_dNdy_pos.eps");
      c1->SaveAs("plots/ks_dNdy_pos.root");
      c1->SaveAs("plots/ks_dNdy_pos.C");
    }

    ks_dNdy->SetXTitle("K^{0}_{S} rapidity");
    ks_dNdy->SetYTitle("Candidates / 0.1");
    ks_dNdy->SetNdivisions(506);
    ks_dNdy->SetMinimum(0);
    ks_dNdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dNdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dNdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dNdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dNdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dNdy.png");
      c1->SaveAs("plots/ks_dNdy.eps");
      c1->SaveAs("plots/ks_dNdy.root");
      c1->SaveAs("plots/ks_dNdy.C");
    }

    ks_dStNdy->SetXTitle("K^{0}_{S} rapidity");
    ks_dStNdy->SetYTitle("Signal-to-noise");
    ks_dStNdy->SetNdivisions(506);
    ks_dStNdy->SetMinimum(0);
    ks_dStNdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dStNdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dStNdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dStNdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dStNdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dStNdy.png");
      c1->SaveAs("plots/ks_dStNdy.eps");
      c1->SaveAs("plots/ks_dStNdy.root");
      c1->SaveAs("plots/ks_dStNdy.C");
    }

    ks_sigdy->SetXTitle("K^{0}_{S} rapidity");
    ks_sigdy->SetYTitle("#sigma");
    ks_sigdy->SetNdivisions(506);
    ks_sigdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_sigdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_sigdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_sigdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_sigdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_sigdy.png");
      c1->SaveAs("plots/ks_sigdy.eps");
      c1->SaveAs("plots/ks_sigdy.root");
      c1->SaveAs("plots/ks_sigdy.C");
    }

    ks_chi2dy->SetXTitle("K^{0}_{S} rapidity");
    ks_chi2dy->SetYTitle("#chi^2/dof");
    ks_chi2dy->SetNdivisions(506);
    ks_chi2dy->SetMinimum(0);
    ks_chi2dy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_chi2dy_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_chi2dy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_chi2dy_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_chi2dy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_chi2dy.png");
      c1->SaveAs("plots/ks_chi2dy.eps");
      c1->SaveAs("plots/ks_chi2dy.root");
      c1->SaveAs("plots/ks_chi2dy.C");
    }

    ks_massdy->SetXTitle("K^{0}_{S} rapidity");
    ks_massdy->SetYTitle("Mass  [GeV/c^2]");
    ks_massdy->SetNdivisions(506);
    ks_massdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_massdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_massdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_massdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_massdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_massdy.png");
      c1->SaveAs("plots/ks_massdy.eps");
      c1->SaveAs("plots/ks_massdy.root");
      c1->SaveAs("plots/ks_massdy.C");
    }

    ks_fitdy->SetXTitle("K^{0}_{S} rapidity");
    ks_fitdy->SetYTitle("Fit type");
    ks_fitdy->SetNdivisions(506);
    ks_fitdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_fitdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_fitdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_fitdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_fitdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_fitdy.png");
      c1->SaveAs("plots/ks_fitdy.eps");
      c1->SaveAs("plots/ks_fitdy.root");
      c1->SaveAs("plots/ks_fitdy.C");
    }

    delete c1;

    //    for( int ndx = 0; ndx < ksyNbins; ndx++ ) {
    //      delete ksPtrs[ndx];
    //    }
    //    delete ksPtrs[];
  }


  //  Lambda y distribution
  if( dodNdy_lam ) {
    cout << "Doing lambda dN/dy..." << endl;

    int lamCanvasX = 8;
    int lamCanvasY = lamyNbins / lamCanvasX;
    if( (lamyNbins)%lamCanvasX != 0 ) lamCanvasY += 1;
    TCanvas* lamyCanvas = new TCanvas("lamyFits", "Lambda mass fits in y bins", 
				       300*lamCanvasX, 300*lamCanvasY);
    lamyCanvas->Divide( lamCanvasX, lamCanvasY );
    
    TH1F* lam_dNdy = new TH1F("lam_dNdy", "dN/dy", lamyNbins, lamyXmin, lamyXmax);
    TH1F* lamReg_dNdy = new TH1F("lamReg_dNdy", "dN/dy", lamyNbins, lamyXmin, lamyXmax);
    TH1F* lamBar_dNdy = new TH1F("lamBar_dNdy", "dN/dy", lamyNbins, lamyXmin, lamyXmax);
    TH1F* lam_dStNdy = new TH1F("lam_dStNdy", "(S/N) / dy", lamyNbins, lamyXmin, lamyXmax);
    TH1F* lam_sigdy = new TH1F("lam_sigdy", "#sigma / dy", lamyNbins, lamyXmin, lamyXmax);
    TH1F* lam_massdy = new TH1F("lam_massdy", "Mass / dy", lamyNbins, lamyXmin, lamyXmax);
    TH1F* lam_fitdy = new TH1F("lam_fitdy", "Fit / dy", lamyNbins, lamyXmin, lamyXmax);
    TH1F* lam_chi2dy = new TH1F("lam_chi2dy", "#chi^{2}/dof / dy", lamyNbins, lamyXmin, lamyXmax);

    for( int ndx = 0; ndx < lamyNbins; ndx++ ) {
      lamyCanvas->cd(ndx + 1);
      double min = lamyXmin + ndx*lamyBinWidth;
      double max = lamyXmin + (ndx+1)*lamyBinWidth;
      cout << "Doing fit for Lambda with " << min << " < y < " << max << endl;

      int fitresult = lamfitobj->Fit(lamPtrs[ndx]);
      cout << "Fit result = " << fitresult << endl;

      lam_dNdy->SetBinContent( ndx + 1, lamfitobj->sig() );
      lam_dNdy->SetBinError( ndx + 1, lamfitobj->sig_err() );
      double ston =  0.95*lamfitobj->sig()/lamfitobj->bg();
      double stone = ston*sqrt( pow(lamfitobj->sig_err() / lamfitobj->sig(),2)
				+ pow(lamfitobj->bg_err() / lamfitobj->bg(),2) );
      lam_dStNdy->SetBinContent( ndx + 1, ston );
      lam_dStNdy->SetBinError( ndx + 1, stone );
      lam_sigdy->SetBinContent( ndx + 1, lamfitobj->res() );
      lam_sigdy->SetBinError( ndx + 1, lamfitobj->res_err() );
      lam_massdy->SetBinContent( ndx + 1, lamfitobj->mass() );
      lam_massdy->SetBinError( ndx + 1, lamfitobj->mass_err() );
      lam_fitdy->SetBinContent( ndx + 1, lamfitobj->fitbg() );
      lam_chi2dy->SetBinContent( ndx + 1, lamfitobj->nChi2() );

      lamPtrs[ndx]->Draw();
      cout << "---------------" << endl;
    }
    if( isMC ) {
      lamyCanvas->SaveAs("plots/LamMassFitsInyBins_MC_" + mcname + ".png");
      lamyCanvas->SaveAs("plots/LamMassFitsInyBins_MC_" + mcname + ".eps");
      lamyCanvas->SaveAs("plots/LamMassFitsInyBins_MC_" + mcname + ".root");
      lamyCanvas->SaveAs("plots/LamMassFitsInyBins_MC_" + mcname + ".pdf");
    }
    else {
      lamyCanvas->SaveAs("plots/LamMassFitsInyBins.png");
      lamyCanvas->SaveAs("plots/LamMassFitsInyBins.eps");
      lamyCanvas->SaveAs("plots/LamMassFitsInyBins.root");
      lamyCanvas->SaveAs("plots/LamMassFitsInyBins.pdf");
    }
    delete lamyCanvas;

    TCanvas* c1 = new TCanvas("c1", "Lambda yields");

    gStyle->SetOptTitle(0);

    fout->cd();

    lam_dNdy->SetXTitle("#Lambda^{0} rapidity");
    lam_dNdy->SetYTitle("Candidates / 0.1");
    lam_dNdy->SetNdivisions(506);
    lam_dNdy->SetMinimum(0);
    lam_dNdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_dNdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_dNdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_dNdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_dNdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_dNdy.png");
      c1->SaveAs("plots/lam_dNdy.eps");
      c1->SaveAs("plots/lam_dNdy.root");
      c1->SaveAs("plots/lam_dNdy.C");
    }

    lam_dStNdy->SetXTitle("#Lambda^{0} rapidity");
    lam_dStNdy->SetYTitle("Signal-to-noise");
    lam_dStNdy->SetNdivisions(506);
    lam_dStNdy->SetMinimum(0);
    lam_dStNdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_dStNdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_dStNdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_dStNdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_dStNdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_dStNdy.png");
      c1->SaveAs("plots/lam_dStNdy.eps");
      c1->SaveAs("plots/lam_dStNdy.root");
      c1->SaveAs("plots/lam_dStNdy.C");
    }

    lam_sigdy->SetXTitle("#Lambda^{0} rapidity");
    lam_sigdy->SetYTitle("#sigma");
    lam_sigdy->SetNdivisions(506);
    lam_sigdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_sigdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_sigdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_sigdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_sigdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_sigdy.png");
      c1->SaveAs("plots/lam_sigdy.eps");
      c1->SaveAs("plots/lam_sigdy.root");
      c1->SaveAs("plots/lam_sigdy.C");
    }

    lam_chi2dy->SetXTitle("#Lambda^{0} rapidity");
    lam_chi2dy->SetYTitle("#chi^2/dof");
    lam_chi2dy->SetNdivisions(506);
    lam_chi2dy->SetMinimum(0);
    lam_chi2dy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_chi2dy_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_chi2dy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_chi2dy_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_chi2dy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_chi2dy.png");
      c1->SaveAs("plots/lam_chi2dy.eps");
      c1->SaveAs("plots/lam_chi2dy.root");
      c1->SaveAs("plots/lam_chi2dy.C");
    }

    lam_massdy->SetXTitle("#Lambda^{0} rapidity");
    lam_massdy->SetYTitle("Mass  [GeV/c^2]");
    lam_massdy->SetNdivisions(506);
    lam_massdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_massdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_massdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_massdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_massdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_massdy.png");
      c1->SaveAs("plots/lam_massdy.eps");
      c1->SaveAs("plots/lam_massdy.root");
      c1->SaveAs("plots/lam_massdy.C");
    }

    lam_fitdy->SetXTitle("#Lambda^{0} rapidity");
    lam_fitdy->SetYTitle("Fit type");
    lam_fitdy->SetNdivisions(506);
    lam_fitdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_fitdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_fitdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_fitdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_fitdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_fitdy.png");
      c1->SaveAs("plots/lam_fitdy.eps");
      c1->SaveAs("plots/lam_fitdy.root");
      c1->SaveAs("plots/lam_fitdy.C");
    }

    delete c1;

    //    for( int ndx = 0; ndx < lamyNbins; ndx++ ) {
    //      delete lamPtrs[ndx];
    //    }
    //    delete lamPtrs[];
    lamyCanvas = new TCanvas("lamyfits", "Lambda mass fits in y bins",
			     300*lamCanvasX, 300*lamCanvasY);
    lamyCanvas->Divide(lamCanvasX, lamCanvasY);
    for( int ndx = 0; ndx < lamyNbins; ndx++ ) {
      lamyCanvas->cd(ndx + 1);
      double min = lamyXmin + ndx*lamyBinWidth;
      double max = lamyXmin + (ndx+1)*lamyBinWidth;
      cout << "Doing fit for Lambda (REG) with " << min << " < y < " << max << endl;

      int fitresult = lamfitobj->Fit(lamRegPtrs[ndx]);
      cout << "Fit result = " << fitresult << endl;

      lamReg_dNdy->SetBinContent( ndx + 1, lamfitobj->sig() );
      lamReg_dNdy->SetBinError( ndx + 1, lamfitobj->sig_err() );
      /*double ston =  0.95*lamfitobj->sig()/lamfitobj->bg();
      double stone = ston*sqrt( pow(lamfitobj->sig_err() / lamfitobj->sig(),2)
				+ pow(lamfitobj->bg_err() / lamfitobj->bg(),2) );
      lam_dStNdy->SetBinContent( ndx + 1, ston );
      lam_dStNdy->SetBinError( ndx + 1, stone );
      lam_sigdy->SetBinContent( ndx + 1, lamfitobj->res() );
      lam_sigdy->SetBinError( ndx + 1, lamfitobj->res_err() );
      lam_massdy->SetBinContent( ndx + 1, lamfitobj->mass() );
      lam_massdy->SetBinError( ndx + 1, lamfitobj->mass_err() );
      lam_fitdy->SetBinContent( ndx + 1, lamfitobj->fitbg() );
      lam_chi2dy->SetBinContent( ndx + 1, lamfitobj->nChi2() );*/

      lamRegPtrs[ndx]->Draw();
      cout << "---------------" << endl;
    }
    if( isMC ) {
      lamyCanvas->SaveAs("plots/LamRegMassFitsInyBins_MC_" + mcname + ".png");
      lamyCanvas->SaveAs("plots/LamRegMassFitsInyBins_MC_" + mcname + ".eps");
      lamyCanvas->SaveAs("plots/LamRegMassFitsInyBins_MC_" + mcname + ".root");
      lamyCanvas->SaveAs("plots/LamRegMassFitsInyBins_MC_" + mcname + ".pdf");
    }
    else {
      lamyCanvas->SaveAs("plots/LamRegMassFitsInyBins.png");
      lamyCanvas->SaveAs("plots/LamRegMassFitsInyBins.eps");
      lamyCanvas->SaveAs("plots/LamRegMassFitsInyBins.root");
      lamyCanvas->SaveAs("plots/LamRegMassFitsInyBins.pdf");
    }
    delete lamyCanvas;

    c1 = new TCanvas("c1", "Lambda yields");

    gStyle->SetOptTitle(0);

    fout->cd();

    lamReg_dNdy->SetXTitle("#Lambda^{0} rapidity");
    lamReg_dNdy->SetYTitle("Candidates / 0.1");
    lamReg_dNdy->SetNdivisions(506);
    lamReg_dNdy->SetMinimum(0);
    lamReg_dNdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lamReg_dNdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/lamReg_dNdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lamReg_dNdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/lamReg_dNdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lamReg_dNdy.png");
      c1->SaveAs("plots/lamReg_dNdy.eps");
      c1->SaveAs("plots/lamReg_dNdy.root");
      c1->SaveAs("plots/lamReg_dNdy.C");
    }

    delete c1;


    lamyCanvas = new TCanvas("lamyfits", "Lambda mass fits in y bins",
			     300*lamCanvasX, 300*lamCanvasY);
    lamyCanvas->Divide(lamCanvasX, lamCanvasY);
    for( int ndx = 0; ndx < lamyNbins; ndx++ ) {
      lamyCanvas->cd(ndx + 1);
      double min = lamyXmin + ndx*lamyBinWidth;
      double max = lamyXmin + (ndx+1)*lamyBinWidth;
      cout << "Doing fit for Lambda (REG) with " << min << " < y < " << max << endl;

      int fitresult = lamfitobj->Fit(lamBarPtrs[ndx]);
      cout << "Fit result = " << fitresult << endl;

      lamBar_dNdy->SetBinContent( ndx + 1, lamfitobj->sig() );
      lamBar_dNdy->SetBinError( ndx + 1, lamfitobj->sig_err() );
      /*double ston =  0.95*lamfitobj->sig()/lamfitobj->bg();
      double stone = ston*sqrt( pow(lamfitobj->sig_err() / lamfitobj->sig(),2)
				+ pow(lamfitobj->bg_err() / lamfitobj->bg(),2) );
      lam_dStNdy->SetBinContent( ndx + 1, ston );
      lam_dStNdy->SetBinError( ndx + 1, stone );
      lam_sigdy->SetBinContent( ndx + 1, lamfitobj->res() );
      lam_sigdy->SetBinError( ndx + 1, lamfitobj->res_err() );
      lam_massdy->SetBinContent( ndx + 1, lamfitobj->mass() );
      lam_massdy->SetBinError( ndx + 1, lamfitobj->mass_err() );
      lam_fitdy->SetBinContent( ndx + 1, lamfitobj->fitbg() );
      lam_chi2dy->SetBinContent( ndx + 1, lamfitobj->nChi2() );*/

      lamBarPtrs[ndx]->Draw();
      cout << "---------------" << endl;
    }
    if( isMC ) {
      lamyCanvas->SaveAs("plots/LamBarMassFitsInyBins_MC_" + mcname + ".png");
      lamyCanvas->SaveAs("plots/LamBarMassFitsInyBins_MC_" + mcname + ".eps");
      lamyCanvas->SaveAs("plots/LamBarMassFitsInyBins_MC_" + mcname + ".root");
      lamyCanvas->SaveAs("plots/LamBarMassFitsInyBins_MC_" + mcname + ".pdf");
    }
    else {
      lamyCanvas->SaveAs("plots/LamBarMassFitsInyBins.png");
      lamyCanvas->SaveAs("plots/LamBarMassFitsInyBins.eps");
      lamyCanvas->SaveAs("plots/LamBarMassFitsInyBins.root");
      lamyCanvas->SaveAs("plots/LamBarMassFitsInyBins.pdf");
    }
    delete lamyCanvas;

    c1 = new TCanvas("c1", "Lambda yields");

    gStyle->SetOptTitle(0);

    fout->cd();

    lamBar_dNdy->SetXTitle("#Lambda^{0} rapidity");
    lamBar_dNdy->SetYTitle("Candidates / 0.1");
    lamBar_dNdy->SetNdivisions(506);
    lamBar_dNdy->SetMinimum(0);
    lamBar_dNdy->Draw();
    labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lamBar_dNdy_MC_" + mcname + ".png");
      c1->SaveAs("plots/lamBar_dNdy_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lamBar_dNdy_MC_" + mcname + ".root");
      c1->SaveAs("plots/lamBar_dNdy_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lamBar_dNdy.png");
      c1->SaveAs("plots/lamBar_dNdy.eps");
      c1->SaveAs("plots/lamBar_dNdy.root");
      c1->SaveAs("plots/lamBar_dNdy.C");
    }

    delete c1;
  }


  // Kshort pt distribution

  if( dodNdpT_ks ) {
    cout << "Doing dN/dpT..." << endl;

    int ksCanvasX = 8;
    int ksCanvasY = ksptNbins / ksCanvasX;
    if( ksptNbins%ksCanvasX != 0 ) ksCanvasY += 1;
    TCanvas* kspTCanvas = new TCanvas("kspTFits", "K0S mass fits in p_{T} bins",
				      300*ksCanvasX, 300*ksCanvasY);
    kspTCanvas->Divide( ksCanvasX, ksCanvasY );

    TH1F* ks_ptscale = new TH1F("ks_ptscale", "pT scaling factor", ksptNbins, &ksptBinMinima[0]);
    TH1F* ks_dNdpT = new TH1F("ks_dNdpT", "dN / dp_{T}", ksptNbins, &ksptBinMinima[0]);
    TH1F* ks_dStNdpT = new TH1F("ks_dStNdpT", "(S/N) / dp_{T}", ksptNbins, &ksptBinMinima[0]);
    TH1F* ks_sigdpT = new TH1F("ks_sigdpT", "#sigma / dp_{T}", ksptNbins, &ksptBinMinima[0]);
    TH1F* ks_massdpT = new TH1F("ks_massdpT", "Mass / dp_{T}", ksptNbins, &ksptBinMinima[0]);
    TH1F* ks_fitdpT = new TH1F("ks_fitdpT", "Fit / dp_{T}", ksptNbins, &ksptBinMinima[0]);
    TH1F* ks_chi2dpT = new TH1F("ks_chi2dpT", "#chi^{2}/dof / dp_{T}", ksptNbins, &ksptBinMinima[0]);

    for( int ndx = 0; ndx < ksptNbins; ndx++ ) {
      kspTCanvas->cd(ndx + 1);
      double min = ksptBinMinima[ndx];
      double max = ksptBinMinima[ndx+1];
      double binWidth = max - min;
      double scalingFactor = binWidth / kspt_default_binWidth;
      ks_ptscale->SetBinContent( ndx + 1, scalingFactor );
      ks_ptscale->SetBinError( ndx + 1, 0. );
      cout << "Doing fit for K0S with " << min << " < pT < " << max << endl;

      if ( ndx == 0 || ndx >= 31 ) {
	fit_command = 1;
      } else {
	fit_command = 0;
      }
      int fitresult = ksfitobj->Fit(ksptPtrs[ndx],fit_command);
      cout << "Fit result = " << fitresult << endl;

      ks_dNdpT->SetBinContent( ndx + 1, ksfitobj->sig() );
      ks_dNdpT->SetBinError( ndx + 1, ksfitobj->sig_err() );
      double ston = 0.95*ksfitobj->sig() / ksfitobj->bg();
      double stone = ston*sqrt( pow(ksfitobj->sig_err() / ksfitobj->sig(),2)
				+ pow(ksfitobj->bg_err() / ksfitobj->bg(),2) );
      ks_dStNdpT->SetBinContent( ndx + 1, ston );
      ks_dStNdpT->SetBinError( ndx + 1, stone );
      ks_sigdpT->SetBinContent( ndx + 1, ksfitobj->res() );
      ks_sigdpT->SetBinError( ndx + 1, ksfitobj->res_err() );
      ks_massdpT->SetBinContent( ndx + 1, ksfitobj->mass() );
      ks_massdpT->SetBinError( ndx + 1, ksfitobj->mass_err() );
      ks_fitdpT->SetBinContent( ndx + 1, ksfitobj->fitbg() );
      ks_chi2dpT->SetBinContent( ndx + 1, ksfitobj->nChi2() );

      ksptPtrs[ndx]->Draw();
      cout << "-----------------" << endl;
    }
    if( isMC ) {
      kspTCanvas->SaveAs("plots/KsMassFitsInpTBins_MC_" + mcname + ".png");
      kspTCanvas->SaveAs("plots/KsMassFitsInpTBins_MC_" + mcname + ".eps");
      kspTCanvas->SaveAs("plots/KsMassFitsInpTBins_MC_" + mcname + ".root");
      kspTCanvas->SaveAs("plots/KsMassFitsInpTBins_MC_" + mcname + ".pdf");
    }
    else {
      kspTCanvas->SaveAs("plots/KsMassFitsInpTBins.png");
      kspTCanvas->SaveAs("plots/KsMassFitsInpTBins.eps");
      kspTCanvas->SaveAs("plots/KsMassFitsInpTBins.root");
      kspTCanvas->SaveAs("plots/KsMassFitsInpTBins.pdf");
    }
    delete kspTCanvas;

    TCanvas* c1 = new TCanvas("c1", "K0S yields");

    gStyle->SetOptTitle(0);

    fout->cd();
    TH1F* ks_dNdpT_scale = ks_dNdpT->Clone("ks_dNdpT_scale");
    ks_dNdpT_scale->Divide(ks_ptscale);

    ks_dNdpT_scale->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
    ks_dNdpT_scale->SetYTitle("Candidates / 0.1 GeV/c");
    ks_dNdpT_scale->SetMinimum(0);
    ks_dNdpT_scale->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dNdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dNdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dNdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dNdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dNdpT.png");
      c1->SaveAs("plots/ks_dNdpT.eps");
      c1->SaveAs("plots/ks_dNdpT.root");
      c1->SaveAs("plots/ks_dNdpT.C");
    }

    ks_dStNdpT->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
    ks_dStNdpT->SetYTitle("Signal-to-noise");
    ks_dStNdpT->SetMinimum(0);
    ks_dStNdpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dStNdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dStNdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dStNdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dStNdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dStNdpT.png");
      c1->SaveAs("plots/ks_dStNdpT.eps");
      c1->SaveAs("plots/ks_dStNdpT.root");
      c1->SaveAs("plots/ks_dStNdpT.C");
    }

    ks_sigdpT->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
    ks_sigdpT->SetYTitle("#sigma");
    ks_sigdpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_sigdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_sigdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_sigdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_sigdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_sigdpT.png");
      c1->SaveAs("plots/ks_sigdpT.eps");
      c1->SaveAs("plots/ks_sigdpT.root");
      c1->SaveAs("plots/ks_sigdpT.C");
    }

    ks_chi2dpT->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
    ks_chi2dpT->SetYTitle("#chi^2/dof");
    ks_chi2dpT->SetMinimum(0);
    ks_chi2dpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_chi2dpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_chi2dpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_chi2dpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_chi2dpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_chi2dpT.png");
      c1->SaveAs("plots/ks_chi2dpT.eps");
      c1->SaveAs("plots/ks_chi2dpT.root");
      c1->SaveAs("plots/ks_chi2dpT.C");
    }

    ks_massdpT->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
    ks_massdpT->SetYTitle("Mass  [GeV/c^2]");
    ks_massdpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_massdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_massdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_massdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_massdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_massdpT.png");
      c1->SaveAs("plots/ks_massdpT.eps");
      c1->SaveAs("plots/ks_massdpT.root");
      c1->SaveAs("plots/ks_massdpT.C");
    }

    ks_fitdpT->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
    ks_fitdpT->SetYTitle("Fit type");
    ks_fitdpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_fitdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_fitdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_fitdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_fitdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_fitdpT.png");
      c1->SaveAs("plots/ks_fitdpT.eps");
      c1->SaveAs("plots/ks_fitdpT.root");
      c1->SaveAs("plots/ks_fitdpT.C");
    }

    delete c1;

    //    for( int ndx = 0; ndx < ksptNbins; ndx++ ) {
    //      delete ksptPtrs[ndx];
    //    }
    //    delete ksptPtrs[];
  }


  // Lambda pt distribution
  if( dodNdpT_lam ) {
    cout << "Doing dN/dpT..." << endl;

    int lamCanvasX = 8;
    int lamCanvasY = lamptNbins / lamCanvasX;
    if( lamptNbins%lamCanvasX != 0 ) lamCanvasY += 1;
    TCanvas* lampTCanvas = new TCanvas("lampTFits", "Lambda mass fits in p_{T} bins",
				      300*lamCanvasX, 300*lamCanvasY);
    lampTCanvas->Divide( lamCanvasX, lamCanvasY );

    TH1F* lam_ptscale = new TH1F("lam_ptscale", "pT scaling factor", lamptNbins, &lamptBinMinima[0]);
    TH1F* lam_dNdpT = new TH1F("lam_dNdpT", "dN / dp_{T}", lamptNbins, &lamptBinMinima[0]);
    TH1F* lamReg_dNdpT = new TH1F("lamReg_dNdpT", "dN / dp_{T}", lamptNbins, &lamptBinMinima[0]);
    TH1F* lamBar_dNdpT = new TH1F("lamBar_dNdpT", "dN / dp_{T}", lamptNbins, &lamptBinMinima[0]);
    TH1F* lam_dStNdpT = new TH1F("lam_dStNdpT", "(S/N) / dp_{T}", lamptNbins, &lamptBinMinima[0]);
    TH1F* lam_sigdpT = new TH1F("lam_sigdpT", "#sigma / dp_{T}", lamptNbins, &lamptBinMinima[0]);
    TH1F* lam_massdpT = new TH1F("lam_massdpT", "Mass / dp_{T}", lamptNbins, &lamptBinMinima[0]);
    TH1F* lam_fitdpT = new TH1F("lam_fitdpT", "Fit / dp_{T}", lamptNbins, &lamptBinMinima[0]);
    TH1F* lam_chi2dpT = new TH1F("lam_chi2dpT", "#chi^{2}/dof / dp_{T}", lamptNbins, &lamptBinMinima[0]);

    for( int ndx = 0; ndx < lamptNbins; ndx++ ) {
      lampTCanvas->cd(ndx + 1);
      double min = lamptBinMinima[ndx];
      double max = lamptBinMinima[ndx+1];
      double binWidth = max - min;
      double scalingFactor = binWidth / lampt_default_binWidth;
      lam_ptscale->SetBinContent( ndx + 1, scalingFactor );
      lam_ptscale->SetBinError( ndx + 1, 0. );
      cout << "Doing fit for Lambda with " << min << " < pT < " << max << endl;

      if ( ndx <= 1 || ndx >= 30 ) {
	fit_command = 1;
      } else {
	fit_command = 0;
      }
      int fitresult = lamfitobj->Fit(lamptPtrs[ndx],fit_command);
      cout << "Fit result = " << fitresult << endl;

      lam_dNdpT->SetBinContent( ndx + 1, lamfitobj->sig() );
      lam_dNdpT->SetBinError( ndx + 1, lamfitobj->sig_err() );
      double ston = 0.95*lamfitobj->sig() / lamfitobj->bg();
      double stone = ston*sqrt( pow(lamfitobj->sig_err() / lamfitobj->sig(),2)
				+ pow(lamfitobj->bg_err() / lamfitobj->bg(),2) );
      lam_dStNdpT->SetBinContent( ndx + 1, ston );
      lam_dStNdpT->SetBinError( ndx + 1, stone );
      lam_sigdpT->SetBinContent( ndx + 1, lamfitobj->res() );
      lam_sigdpT->SetBinError( ndx + 1, lamfitobj->res_err() );
      lam_massdpT->SetBinContent( ndx + 1, lamfitobj->mass() );
      lam_massdpT->SetBinError( ndx + 1, lamfitobj->mass_err() );
      lam_fitdpT->SetBinContent( ndx + 1, lamfitobj->fitbg() );
      lam_chi2dpT->SetBinContent( ndx + 1, lamfitobj->nChi2() );

      lamptPtrs[ndx]->Draw();
      cout << "-----------------" << endl;
    }
    if( isMC ) {
      lampTCanvas->SaveAs("plots/LamMassFitsInpTBins_MC_" + mcname + ".png");
      lampTCanvas->SaveAs("plots/LamMassFitsInpTBins_MC_" + mcname + ".eps");
      lampTCanvas->SaveAs("plots/LamMassFitsInpTBins_MC_" + mcname + ".root");
      lampTCanvas->SaveAs("plots/LamMassFitsInpTBins_MC_" + mcname + ".pdf");
    }
    else {
      lampTCanvas->SaveAs("plots/LamMassFitsInpTBins.png");
      lampTCanvas->SaveAs("plots/LamMassFitsInpTBins.eps");
      lampTCanvas->SaveAs("plots/LamMassFitsInpTBins.root");
      lampTCanvas->SaveAs("plots/LamMassFitsInpTBins.pdf");
    }
    delete lampTCanvas;

    TCanvas* c1 = new TCanvas("c1", "Lambda yields");

    gStyle->SetOptTitle(0);

    fout->cd();
    TH1F* lam_dNdpT_scale = lam_dNdpT->Clone("lam_dNdpT_scale");
    lam_dNdpT_scale->Divide(lam_ptscale);

    lam_dNdpT_scale->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
    lam_dNdpT_scale->SetYTitle("Candidates / 0.1 GeV/c");
    lam_dNdpT_scale->SetMinimum(0);
    lam_dNdpT_scale->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_dNdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_dNdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_dNdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_dNdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_dNdpT.png");
      c1->SaveAs("plots/lam_dNdpT.eps");
      c1->SaveAs("plots/lam_dNdpT.root");
      c1->SaveAs("plots/lam_dNdpT.C");
    }

    lam_dStNdpT->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
    lam_dStNdpT->SetYTitle("Signal-to-noise");
    lam_dStNdpT->SetMinimum(0);
    lam_dStNdpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_dStNdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_dStNdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_dStNdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_dStNdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_dStNdpT.png");
      c1->SaveAs("plots/lam_dStNdpT.eps");
      c1->SaveAs("plots/lam_dStNdpT.root");
      c1->SaveAs("plots/lam_dStNdpT.C");
    }

    lam_sigdpT->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
    lam_sigdpT->SetYTitle("#sigma");
    lam_sigdpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_sigdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_sigdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_sigdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_sigdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_sigdpT.png");
      c1->SaveAs("plots/lam_sigdpT.eps");
      c1->SaveAs("plots/lam_sigdpT.root");
      c1->SaveAs("plots/lam_sigdpT.C");
    }

    lam_chi2dpT->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
    lam_chi2dpT->SetYTitle("#chi^2/dof");
    lam_chi2dpT->SetMinimum(0);
    lam_chi2dpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_chi2dpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_chi2dpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_chi2dpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_chi2dpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_chi2dpT.png");
      c1->SaveAs("plots/lam_chi2dpT.eps");
      c1->SaveAs("plots/lam_chi2dpT.root");
      c1->SaveAs("plots/lam_chi2dpT.C");
    }

    lam_massdpT->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
    lam_massdpT->SetYTitle("Mass  [GeV/c^2]");
    lam_massdpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_massdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_massdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_massdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_massdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_massdpT.png");
      c1->SaveAs("plots/lam_massdpT.eps");
      c1->SaveAs("plots/lam_massdpT.root");
      c1->SaveAs("plots/lam_massdpT.C");
    }

    lam_fitdpT->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
    lam_fitdpT->SetYTitle("Fit type");
    lam_fitdpT->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_fitdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_fitdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_fitdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_fitdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_fitdpT.png");
      c1->SaveAs("plots/lam_fitdpT.eps");
      c1->SaveAs("plots/lam_fitdpT.root");
      c1->SaveAs("plots/lam_fitdpT.C");
    }

    delete c1;

    lampTCanvas = new TCanvas("lamptfits", "Lambda mass fits in pt bins",
			      300*lamCanvasX, 300*lamCanvasY);
    lampTCanvas->Divide(lamCanvasX, lamCanvasY);
    for( int ndx = 0; ndx < lamptNbins; ndx++ ) {
      lampTCanvas->cd(ndx + 1);
      double min = lamptBinMinima[ndx];
      double max = lamptBinMinima[ndx+1];
      double binWidth = max - min;
      double scalingFactor = binWidth / lampt_default_binWidth;
      //lam_ptscale->SetBinContent( ndx + 1, scalingFactor );
      //lam_ptscale->SetBinError( ndx + 1, 0. );
      cout << "Doing fit for Lambda with " << min << " < pT < " << max << endl;

      if ( ndx <= 1 || ndx >= 30 ) {
	fit_command = 1;
      } else {
	fit_command = 0;
      }
      int fitresult = lamfitobj->Fit(lamRegptPtrs[ndx],fit_command);
      cout << "Fit result = " << fitresult << endl;

      lamReg_dNdpT->SetBinContent( ndx + 1, lamfitobj->sig() );
      lamReg_dNdpT->SetBinError( ndx + 1, lamfitobj->sig_err() );
      /*double ston = 0.95*lamfitobj->sig() / lamfitobj->bg();
      double stone = ston*sqrt( pow(lamfitobj->sig_err() / lamfitobj->sig(),2)
				+ pow(lamfitobj->bg_err() / lamfitobj->bg(),2) );
      lam_dStNdpT->SetBinContent( ndx + 1, ston );
      lam_dStNdpT->SetBinError( ndx + 1, stone );
      lam_sigdpT->SetBinContent( ndx + 1, lamfitobj->res() );
      lam_sigdpT->SetBinError( ndx + 1, lamfitobj->res_err() );
      lam_massdpT->SetBinContent( ndx + 1, lamfitobj->mass() );
      lam_massdpT->SetBinError( ndx + 1, lamfitobj->mass_err() );
      lam_fitdpT->SetBinContent( ndx + 1, lamfitobj->fitbg() );
      lam_chi2dpT->SetBinContent( ndx + 1, lamfitobj->nChi2() );*/

      lamRegptPtrs[ndx]->Draw();
      cout << "-----------------" << endl;
    }
    if( isMC ) {
      lampTCanvas->SaveAs("plots/LamRegMassFitsInpTBins_MC_" + mcname + ".png");
      lampTCanvas->SaveAs("plots/LamRegMassFitsInpTBins_MC_" + mcname + ".eps");
      lampTCanvas->SaveAs("plots/LamRegMassFitsInpTBins_MC_" + mcname + ".root");
      lampTCanvas->SaveAs("plots/LamRegMassFitsInpTBins_MC_" + mcname + ".pdf");
    }
    else {
      lampTCanvas->SaveAs("plots/LamRegMassFitsInpTBins.png");
      lampTCanvas->SaveAs("plots/LamRegMassFitsInpTBins.eps");
      lampTCanvas->SaveAs("plots/LamRegMassFitsInpTBins.root");
      lampTCanvas->SaveAs("plots/LamRegMassFitsInpTBins.pdf");
    }
    delete lampTCanvas;

    c1 = new TCanvas("c1", "Lambda yields");

    gStyle->SetOptTitle(0);

    fout->cd();
    TH1F* lamReg_dNdpT_scale = lamReg_dNdpT->Clone("lamReg_dNdpT_scale");
    lamReg_dNdpT_scale->Divide(lam_ptscale);

    lamReg_dNdpT_scale->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
    lamReg_dNdpT_scale->SetYTitle("Candidates / 0.1 GeV/c");
    lamReg_dNdpT_scale->SetMinimum(0);
    lamReg_dNdpT_scale->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lamReg_dNdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/lamReg_dNdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lamReg_dNdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/lamReg_dNdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lamReg_dNdpT.png");
      c1->SaveAs("plots/lamReg_dNdpT.eps");
      c1->SaveAs("plots/lamReg_dNdpT.root");
      c1->SaveAs("plots/lamReg_dNdpT.C");
    }

    delete c1;

    lampTCanvas = new TCanvas("lamptfits", "Lambda mass fits in pt bins",
			      300*lamCanvasX, 300*lamCanvasY);
    lamyCanvas->Divide(lamCanvasX, lamCanvasY);
    for( int ndx = 0; ndx < lamptNbins; ndx++ ) {
      lampTCanvas->cd(ndx + 1);
      double min = lamptBinMinima[ndx];
      double max = lamptBinMinima[ndx+1];
      double binWidth = max - min;
      double scalingFactor = binWidth / lampt_default_binWidth;
      //lam_ptscale->SetBinContent( ndx + 1, scalingFactor );
      //lam_ptscale->SetBinError( ndx + 1, 0. );
      cout << "Doing fit for Lambda with " << min << " < pT < " << max << endl;

      if ( ndx <= 1 || ndx >= 30 ) {
	fit_command = 1;
      } else {
	fit_command = 0;
      }
      int fitresult = lamfitobj->Fit(lamBarptPtrs[ndx],fit_command);
      cout << "Fit result = " << fitresult << endl;

      lamBar_dNdpT->SetBinContent( ndx + 1, lamfitobj->sig() );
      lamBar_dNdpT->SetBinError( ndx + 1, lamfitobj->sig_err() );
      /*double ston = 0.95*lamfitobj->sig() / lamfitobj->bg();
      double stone = ston*sqrt( pow(lamfitobj->sig_err() / lamfitobj->sig(),2)
				+ pow(lamfitobj->bg_err() / lamfitobj->bg(),2) );
      lam_dStNdpT->SetBinContent( ndx + 1, ston );
      lam_dStNdpT->SetBinError( ndx + 1, stone );
      lam_sigdpT->SetBinContent( ndx + 1, lamfitobj->res() );
      lam_sigdpT->SetBinError( ndx + 1, lamfitobj->res_err() );
      lam_massdpT->SetBinContent( ndx + 1, lamfitobj->mass() );
      lam_massdpT->SetBinError( ndx + 1, lamfitobj->mass_err() );
      lam_fitdpT->SetBinContent( ndx + 1, lamfitobj->fitbg() );
      lam_chi2dpT->SetBinContent( ndx + 1, lamfitobj->nChi2() );*/

      lamBarptPtrs[ndx]->Draw();
      cout << "-----------------" << endl;
    }
    if( isMC ) {
      lampTCanvas->SaveAs("plots/LamBarMassFitsInpTBins_MC_" + mcname + ".png");
      lampTCanvas->SaveAs("plots/LamBarMassFitsInpTBins_MC_" + mcname + ".eps");
      lampTCanvas->SaveAs("plots/LamBarMassFitsInpTBins_MC_" + mcname + ".root");
      lampTCanvas->SaveAs("plots/LamBarMassFitsInpTBins_MC_" + mcname + ".pdf");
    }
    else {
      lampTCanvas->SaveAs("plots/LamBarMassFitsInpTBins.png");
      lampTCanvas->SaveAs("plots/LamBarMassFitsInpTBins.eps");
      lampTCanvas->SaveAs("plots/LamBarMassFitsInpTBins.root");
      lampTCanvas->SaveAs("plots/LamBarMassFitsInpTBins.pdf");
    }
    delete lampTCanvas;

    c1 = new TCanvas("c1", "Lambda yields");

    gStyle->SetOptTitle(0);

    fout->cd();
    TH1F* lamBar_dNdpT_scale = lamBar_dNdpT->Clone("lamBar_dNdpT_scale");
    lamBar_dNdpT_scale->Divide(lam_ptscale);

    lamBar_dNdpT_scale->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
    lamBar_dNdpT_scale->SetYTitle("Candidates / 0.1 GeV/c");
    lamBar_dNdpT_scale->SetMinimum(0);
    lamBar_dNdpT_scale->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lamBar_dNdpT_MC_" + mcname + ".png");
      c1->SaveAs("plots/lamBar_dNdpT_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lamBar_dNdpT_MC_" + mcname + ".root");
      c1->SaveAs("plots/lamBar_dNdpT_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lamBar_dNdpT.png");
      c1->SaveAs("plots/lamBar_dNdpT.eps");
      c1->SaveAs("plots/lamBar_dNdpT.root");
      c1->SaveAs("plots/lamBar_dNdpT.C");
    }

    delete c1;

    //    for( int ndx = 0; ndx < lamptNbins; ndx++ ) {
    //      delete lamptPtrs[ndx];
    //    }
    //    delete lamptPtrs[];
  }



  // Kshort p distribution

  if( dodNdp_ks ) {
    cout << "Doing dN/dp..." << endl;

    int ksCanvasX = 8;
    int ksCanvasY = kspNbins / ksCanvasX;
    if( kspNbins%ksCanvasX != 0 ) ksCanvasY += 1;
    TCanvas* kspCanvas = new TCanvas("kspFits", "K0S mass fits in p_{T} bins",
				     300*ksCanvasX, 300*ksCanvasY);
    kspCanvas->Divide( ksCanvasX, ksCanvasY );

    TH1F* ks_pscale = new TH1F("ks_pscale", "p scaling factor", kspNbins, &kspBinMinima[0]);
    TH1F* ks_dNdp = new TH1F("ks_dNdp", "dN/dp", kspNbins, &kspBinMinima[0]);
    TH1F* ks_dStNdp = new TH1F("ks_dStNdp", "(S/N) / dp", kspNbins, &kspBinMinima[0]);
    TH1F* ks_sigdp = new TH1F("ks_sigdp", "#sigma / dp", kspNbins, &kspBinMinima[0]);
    TH1F* ks_massdp = new TH1F("ks_massdp", "Mass / dp", kspNbins, &kspBinMinima[0]);
    TH1F* ks_fitdp = new TH1F("ks_fitdp", "Fit / dp", kspNbins, &kspBinMinima[0]);
    TH1F* ks_chi2dp = new TH1F("ks_chi2dp", "#chi^{2}/dof / dp", kspNbins, &kspBinMinima[0]);

    for( int ndx = 0; ndx < kspNbins; ndx++ ) {
      kspCanvas->cd(ndx + 1);
      double min = kspBinMinima[ndx];
      double max = kspBinMinima[ndx+1];
      double binWidth = max - min;
      double scalingFactor = binWidth / ksp_default_binWidth;
      ks_pscale->SetBinContent( ndx + 1, scalingFactor );
      ks_pscale->SetBinError( ndx + 1, 0. );
      cout << "Doing fit for K0S with " << min << " < p < " << max << endl;

      if ( ndx == 0 || ndx >= 40 ) {
	fit_command = 1;
      } else {
	fit_command = 0;
      }
      int fitresult = ksfitobj->Fit(kspPtrs[ndx],fit_command);
      cout << "Fit result = " << fitresult << endl;

      ks_dNdp->SetBinContent( ndx + 1, ksfitobj->sig() );
      ks_dNdp->SetBinError( ndx + 1, ksfitobj->sig_err() );
      double ston = 0.95*ksfitobj->sig() / ksfitobj->bg();
      double stone = ston*sqrt( pow(ksfitobj->sig_err() / ksfitobj->sig(),2)
				+ pow(ksfitobj->bg_err() / ksfitobj->bg(),2) );
      ks_dStNdp->SetBinContent( ndx + 1, ston );
      ks_dStNdp->SetBinError( ndx + 1, stone );
      ks_sigdp->SetBinContent( ndx + 1, ksfitobj->res() );
      ks_sigdp->SetBinError( ndx + 1, ksfitobj->res_err() );
      ks_massdp->SetBinContent( ndx + 1, ksfitobj->mass() );
      ks_massdp->SetBinError( ndx + 1, ksfitobj->mass_err() );
      ks_fitdp->SetBinContent( ndx + 1, ksfitobj->fitbg() );
      ks_chi2dp->SetBinContent( ndx + 1, ksfitobj->nChi2() );

      kspPtrs[ndx]->Draw();
      cout << "-----------------" << endl;
    }
    if( isMC ) {
      kspCanvas->SaveAs("plots/KsMassFitsInpBins_MC_" + mcname + ".png");
      kspCanvas->SaveAs("plots/KsMassFitsInpBins_MC_" + mcname + ".eps");
      kspCanvas->SaveAs("plots/KsMassFitsInpBins_MC_" + mcname + ".root");
      kspCanvas->SaveAs("plots/KsMassFitsInpBins_MC_" + mcname + ".pdf");
    }
    else {
      kspCanvas->SaveAs("plots/KsMassFitsInpBins.png");
      kspCanvas->SaveAs("plots/KsMassFitsInpBins.eps");
      kspCanvas->SaveAs("plots/KsMassFitsInpBins.root");
      kspCanvas->SaveAs("plots/KsMassFitsInpBins.pdf");
    }
    delete kspCanvas;

    TCanvas* c1 = new TCanvas("c1", "K0S yields");

    gStyle->SetOptTitle(0);

    fout->cd();
    TH1F* ks_dNdp_scale = ks_dNdp->Clone("ks_dNdp_scale");
    ks_dNdp_scale->Divide(ks_pscale);

    ks_dNdp_scale->SetXTitle("K^{0}_{S} momentum [GeV/c]");
    ks_dNdp_scale->SetYTitle("Candidates / 0.1 GeV/c");
    ks_dNdp_scale->SetMinimum(0);
    ks_dNdp_scale->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dNdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dNdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dNdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dNdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dNdp.png");
      c1->SaveAs("plots/ks_dNdp.eps");
      c1->SaveAs("plots/ks_dNdp.root");
      c1->SaveAs("plots/ks_dNdp.C");
    }

    ks_dStNdp->SetXTitle("K^{0}_{S} momentum [GeV/c]");
    ks_dStNdp->SetYTitle("Signal-to-noise");
    ks_dStNdp->SetMinimum(0);
    ks_dStNdp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dStNdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dStNdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dStNdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dStNdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dStNdp.png");
      c1->SaveAs("plots/ks_dStNdp.eps");
      c1->SaveAs("plots/ks_dStNdp.root");
      c1->SaveAs("plots/ks_dStNdp.C");
    }

    ks_sigdp->SetXTitle("K^{0}_{S} momentum [GeV/c]");
    ks_sigdp->SetYTitle("#sigma");
    ks_sigdp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_sigdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_sigdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_sigdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_sigdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_sigdp.png");
      c1->SaveAs("plots/ks_sigdp.eps");
      c1->SaveAs("plots/ks_sigdp.root");
      c1->SaveAs("plots/ks_sigdp.C");
    }

    ks_chi2dp->SetXTitle("K^{0}_{S} momentum [GeV/c]");
    ks_chi2dp->SetYTitle("#chi^2/dof");
    ks_chi2dp->SetMinimum(0);
    ks_chi2dp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_chi2dp_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_chi2dp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_chi2dp_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_chi2dp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_chi2dp.png");
      c1->SaveAs("plots/ks_chi2dp.eps");
      c1->SaveAs("plots/ks_chi2dp.root");
      c1->SaveAs("plots/ks_chi2dp.C");
    }

    ks_massdp->SetXTitle("K^{0}_{S} momentum [GeV/c]");
    ks_massdp->SetYTitle("Mass  [GeV/c^2]");
    ks_massdp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_massdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_massdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_massdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_massdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_massdp.png");
      c1->SaveAs("plots/ks_massdp.eps");
      c1->SaveAs("plots/ks_massdp.root");
      c1->SaveAs("plots/ks_massdp.C");
    }

    ks_fitdp->SetXTitle("K^{0}_{S} momentum [GeV/c]");
    ks_fitdp->SetYTitle("Fit type");
    ks_fitdp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_fitdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_fitdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_fitdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_fitdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_fitdp.png");
      c1->SaveAs("plots/ks_fitdp.eps");
      c1->SaveAs("plots/ks_fitdp.root");
      c1->SaveAs("plots/ks_fitdp.C");
    }

   delete c1;

    //    for( int ndx = 0; ndx < kspNbins; ndx++ ) {
    //      delete kspPtrs[ndx];
    //    }
    //    delete kspPtrs[];
  }


  // Lambda p distribution
  if( dodNdp_lam ) {
    cout << "Doing dN/dp..." << endl;

    int lamCanvasX = 8;
    int lamCanvasY = lampNbins / lamCanvasX;
    if( lampNbins%lamCanvasX != 0 ) lamCanvasY += 1;
    TCanvas* lampCanvas = new TCanvas("lampFits", "Lambda mass fits in p_{T} bins",
				     300*lamCanvasX, 300*lamCanvasY);
    lampCanvas->Divide( lamCanvasX, lamCanvasY );

    TH1F* lam_pscale = new TH1F("lam_pscale", "p scaling factor", lampNbins, &lampBinMinima[0]);
    TH1F* lam_dNdp = new TH1F("lam_dNdp", "dN/dp", lampNbins, &lampBinMinima[0]);
    TH1F* lam_dStNdp = new TH1F("lam_dStNdp", "(S/N) / dp", lampNbins, &lampBinMinima[0]);
    TH1F* lam_sigdp = new TH1F("lam_sigdp", "#sigma / dp", lampNbins, &lampBinMinima[0]);
    TH1F* lam_massdp = new TH1F("lam_massdp", "Mass / dp", lampNbins, &lampBinMinima[0]);
    TH1F* lam_fitdp = new TH1F("lam_fitdp", "Fit / dp", lampNbins, &lampBinMinima[0]);
    TH1F* lam_chi2dp = new TH1F("lam_chi2dp", "#chi^{2}/dof / dp", lampNbins, &lampBinMinima[0]);

    for( int ndx = 0; ndx < lampNbins; ndx++ ) {
      lampCanvas->cd(ndx + 1);
      double min = lampBinMinima[ndx];
      double max = lampBinMinima[ndx+1];
      double binWidth = max - min;
      double scalingFactor = binWidth / lamp_default_binWidth;
      lam_pscale->SetBinContent( ndx + 1, scalingFactor );
      lam_pscale->SetBinError( ndx + 1, 0. );
      cout << "Doing fit for Lambda with " << min << " < p < " << max << endl;

      if ( ndx <= 1 || ndx >= 28 ) {
	fit_command = 1;
      } else {
	fit_command = 0;
      }
      int fitresult = lamfitobj->Fit(lampPtrs[ndx],fit_command);
      cout << "Fit result = " << fitresult << endl;

      lam_dNdp->SetBinContent( ndx + 1, lamfitobj->sig() );
      lam_dNdp->SetBinError( ndx + 1, lamfitobj->sig_err() );
      double ston = 0.95*lamfitobj->sig() / lamfitobj->bg();
      double stone = ston*sqrt( pow(lamfitobj->sig_err() / lamfitobj->sig(),2)
				+ pow(lamfitobj->bg_err() / lamfitobj->bg(),2) );
      lam_dStNdp->SetBinContent( ndx + 1, ston );
      lam_dStNdp->SetBinError( ndx + 1, stone );
      lam_sigdp->SetBinContent( ndx + 1, lamfitobj->res() );
      lam_sigdp->SetBinError( ndx + 1, lamfitobj->res_err() );
      lam_massdp->SetBinContent( ndx + 1, lamfitobj->mass() );
      lam_massdp->SetBinError( ndx + 1, lamfitobj->mass_err() );
      lam_fitdp->SetBinContent( ndx + 1, lamfitobj->fitbg() );
      lam_chi2dp->SetBinContent( ndx + 1, lamfitobj->nChi2() );

      lampPtrs[ndx]->Draw();
      cout << "-----------------" << endl;
    }
    if( isMC ) {
      lampCanvas->SaveAs("plots/LamMassFitsInpBins_MC_" + mcname + ".png");
      lampCanvas->SaveAs("plots/LamMassFitsInpBins_MC_" + mcname + ".eps");
      lampCanvas->SaveAs("plots/LamMassFitsInpBins_MC_" + mcname + ".root");
      lampCanvas->SaveAs("plots/LamMassFitsInpBins_MC_" + mcname + ".pdf");
    }
    else {
      lampCanvas->SaveAs("plots/LamMassFitsInpBins.png");
      lampCanvas->SaveAs("plots/LamMassFitsInpBins.eps");
      lampCanvas->SaveAs("plots/LamMassFitsInpBins.root");
      lampCanvas->SaveAs("plots/LamMassFitsInpBins.pdf");
    }
    delete lampCanvas;

    TCanvas* c1 = new TCanvas("c1", "Lambda yields");

    gStyle->SetOptTitle(0);

    fout->cd();
    TH1F* lam_dNdp_scale = lam_dNdp->Clone("lam_dNdp_scale");
    lam_dNdp_scale->Divide(lam_pscale);

    lam_dNdp_scale->SetXTitle("#Lambda^{0} momentum [GeV/c]");
    lam_dNdp_scale->SetYTitle("Candidates / 0.1 GeV/c");
    lam_dNdp_scale->SetMinimum(0);
    lam_dNdp_scale->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_dNdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_dNdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_dNdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_dNdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_dNdp.png");
      c1->SaveAs("plots/lam_dNdp.eps");
      c1->SaveAs("plots/lam_dNdp.root");
      c1->SaveAs("plots/lam_dNdp.C");
    }

    lam_dStNdp->SetXTitle("#Lambda^{0} momentum [GeV/c]");
    lam_dStNdp->SetYTitle("Signal-to-noise");
    lam_dStNdp->SetMinimum(0);
    lam_dStNdp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_dStNdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_dStNdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_dStNdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_dStNdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_dStNdp.png");
      c1->SaveAs("plots/lam_dStNdp.eps");
      c1->SaveAs("plots/lam_dStNdp.root");
      c1->SaveAs("plots/lam_dStNdp.C");
    }

    lam_sigdp->SetXTitle("#Lambda^{0} momentum [GeV/c]");
    lam_sigdp->SetYTitle("#sigma");
    lam_sigdp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_sigdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_sigdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_sigdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_sigdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_sigdp.png");
      c1->SaveAs("plots/lam_sigdp.eps");
      c1->SaveAs("plots/lam_sigdp.root");
      c1->SaveAs("plots/lam_sigdp.C");
    }

    lam_chi2dp->SetXTitle("#Lambda^{0} momentum [GeV/c]");
    lam_chi2dp->SetYTitle("#chi^2/dof");
    lam_chi2dp->SetMinimum(0);
    lam_chi2dp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_chi2dp_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_chi2dp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_chi2dp_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_chi2dp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_chi2dp.png");
      c1->SaveAs("plots/lam_chi2dp.eps");
      c1->SaveAs("plots/lam_chi2dp.root");
      c1->SaveAs("plots/lam_chi2dp.C");
    }

    lam_massdp->SetXTitle("#Lambda^{0} momentum [GeV/c]");
    lam_massdp->SetYTitle("Mass  [GeV/c^2]");
    lam_massdp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_massdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_massdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_massdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_massdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_massdp.png");
      c1->SaveAs("plots/lam_massdp.eps");
      c1->SaveAs("plots/lam_massdp.root");
      c1->SaveAs("plots/lam_massdp.C");
    }

    lam_fitdp->SetXTitle("#Lambda^{0} momentum [GeV/c]");
    lam_fitdp->SetYTitle("Fit type");
    lam_fitdp->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_fitdp_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_fitdp_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_fitdp_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_fitdp_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_fitdp.png");
      c1->SaveAs("plots/lam_fitdp.eps");
      c1->SaveAs("plots/lam_fitdp.root");
      c1->SaveAs("plots/lam_fitdp.C");
    }

   delete c1;

    //    for( int ndx = 0; ndx < lampNbins; ndx++ ) {
    //      delete lampPtrs[ndx];
    //    }
    //    delete lampPtrs[];
  }


  // Kshort lifetime distribution

  if( dodNdlife_ks ) {
    cout << "Doing dN/dlife..." << endl;

    int ksCanvasX = 8;
    int ksCanvasY = kslifeNbins / ksCanvasX;
    if( kslifeNbins%ksCanvasX != 0 ) ksCanvasY += 1;
    TCanvas* kslifeCanvas = new TCanvas("kslifeFits", "K0S mass fits in ct bins",
					300*ksCanvasX, 300*ksCanvasY);
    kslifeCanvas->Divide( ksCanvasX, ksCanvasY );

    TH1F* ks_lifescale = new TH1F("ks_lifescale", "life scaling factor", kslifeNbins, &kslifeBinMinima[0]);
    TH1F* ks_dNdlife = new TH1F("ks_dNdlife", "dN/dct", kslifeNbins, &kslifeBinMinima[0]);
    TH1F* ks_dStNdlife = new TH1F("ks_dStNdlife", "(S/N) / dct", kslifeNbins, &kslifeBinMinima[0]);
    TH1F* ks_sigdlife = new TH1F("ks_sigdlife", "#sigma / dct", kslifeNbins, &kslifeBinMinima[0]);
    TH1F* ks_massdlife = new TH1F("ks_massdlife", "Mass / dct", kslifeNbins, &kslifeBinMinima[0]);
    TH1F* ks_fitdlife = new TH1F("ks_fitdlife", "Fit / dct", kslifeNbins, &kslifeBinMinima[0]);
    TH1F* ks_chi2dlife = new TH1F("ks_chi2dlife", "#chi^{2}/dof / dct", kslifeNbins, &kslifeBinMinima[0]);

    for( int ndx = 0; ndx < kslifeNbins; ndx++ ) {
      kslifeCanvas->cd(ndx + 1);
      double min = kslifeBinMinima[ndx];
      double max = kslifeBinMinima[ndx+1];
      double binWidth = max - min;
      double scalingFactor = binWidth / kslife_default_binWidth;
      ks_lifescale->SetBinContent( ndx + 1, scalingFactor );
      ks_lifescale->SetBinError( ndx + 1, 0. );
      cout << "Doing fit for K0S with " << min << " < life < " << max << endl;

      int fitresult = ksfitobj->Fit(kslifePtrs[ndx]);
      cout << "Fit result = " << fitresult << endl;

      ks_dNdlife->SetBinContent( ndx + 1, ksfitobj->sig() );
      ks_dNdlife->SetBinError( ndx + 1, ksfitobj->sig_err() );
      double ston =  0.95*ksfitobj->sig()/ksfitobj->bg();
      double stone = ston*sqrt( pow(ksfitobj->sig_err() / ksfitobj->sig(),2)
				+ pow(ksfitobj->bg_err() / ksfitobj->bg(),2) );
      ks_dStNdlife->SetBinContent( ndx + 1, ston );
      ks_dStNdlife->SetBinError( ndx + 1, stone );
      ks_sigdlife->SetBinContent( ndx + 1, ksfitobj->res() );
      ks_sigdlife->SetBinError( ndx + 1, ksfitobj->res_err() );
      ks_massdlife->SetBinContent( ndx + 1, ksfitobj->mass() );
      ks_massdlife->SetBinError( ndx + 1, ksfitobj->mass_err() );
      ks_fitdlife->SetBinContent( ndx + 1, ksfitobj->fitbg() );
      ks_chi2dlife->SetBinContent( ndx + 1, ksfitobj->nChi2() );

      kslifePtrs[ndx]->Draw();
      cout << "-----------------" << endl;
    }
    if( isMC ) {
      kslifeCanvas->SaveAs("plots/KsMassFitsInlifeBins_MC_" + mcname + ".png");
      kslifeCanvas->SaveAs("plots/KsMassFitsInlifeBins_MC_" + mcname + ".eps");
      kslifeCanvas->SaveAs("plots/KsMassFitsInlifeBins_MC_" + mcname + ".root");
      kslifeCanvas->SaveAs("plots/KsMassFitsInlifeBins_MC_" + mcname + ".pdf");
    }
    else {
      kslifeCanvas->SaveAs("plots/KsMassFitsInlifeBins.png");
      kslifeCanvas->SaveAs("plots/KsMassFitsInlifeBins.eps");
      kslifeCanvas->SaveAs("plots/KsMassFitsInlifeBins.root");
      kslifeCanvas->SaveAs("plots/KsMassFitsInlifeBins.pdf");
    }
    delete kslifeCanvas;

    TCanvas* c1 = new TCanvas("c1", "K0S yields");

    gStyle->SetOptTitle(0);

    fout->cd();
    TH1F* ks_dNdlife_scale = ks_dNdlife->Clone("ks_dNdlife_scale");
    ks_dNdlife_scale->Divide(ks_lifescale);

    ks_dNdlife_scale->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
    ks_dNdlife_scale->SetYTitle("Candidates / 0.1 cm");
    ks_dNdlife_scale->SetMinimum(0);
    ks_dNdlife_scale->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dNdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dNdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dNdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dNdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dNdlife.png");
      c1->SaveAs("plots/ks_dNdlife.eps");
      c1->SaveAs("plots/ks_dNdlife.root");
      c1->SaveAs("plots/ks_dNdlife.C");
    }

    ks_dStNdlife->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
    ks_dStNdlife->SetYTitle("Signal-to-noise");
    ks_dStNdlife->SetMinimum(0);
    ks_dStNdlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_dStNdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_dStNdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_dStNdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_dStNdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_dStNdlife.png");
      c1->SaveAs("plots/ks_dStNdlife.eps");
      c1->SaveAs("plots/ks_dStNdlife.root");
      c1->SaveAs("plots/ks_dStNdlife.C");
    }

    ks_sigdlife->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
    ks_sigdlife->SetYTitle("#sigma");
    ks_sigdlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_sigdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_sigdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_sigdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_sigdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_sigdlife.png");
      c1->SaveAs("plots/ks_sigdlife.eps");
      c1->SaveAs("plots/ks_sigdlife.root");
      c1->SaveAs("plots/ks_sigdlife.C");
    }

    ks_chi2dlife->SetXTitle("K^{0}_{S} lifetime  [ct]");
    ks_chi2dlife->SetYTitle("#chi^2/dof");
    ks_chi2dlife->SetMinimum(0);
    ks_chi2dlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_chi2dlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_chi2dlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_chi2dlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_chi2dlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_chi2dlife.png");
      c1->SaveAs("plots/ks_chi2dlife.eps");
      c1->SaveAs("plots/ks_chi2dlife.root");
      c1->SaveAs("plots/ks_chi2dlife.C");
    }

    ks_massdlife->SetXTitle("K^{0}_{S} lifetime  [ct]");
    ks_massdlife->SetYTitle("Mass  [GeV/c^2]");
    ks_massdlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_massdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_massdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_massdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_massdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_massdlife.png");
      c1->SaveAs("plots/ks_massdlife.eps");
      c1->SaveAs("plots/ks_massdlife.root");
      c1->SaveAs("plots/ks_massdlife.C");
    }

    ks_fitdlife->SetXTitle("K^{0}_{S} lifetime  [ct]");
    ks_fitdlife->SetYTitle("Fit type");
    ks_fitdlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/ks_fitdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/ks_fitdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/ks_fitdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/ks_fitdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/ks_fitdlife.png");
      c1->SaveAs("plots/ks_fitdlife.eps");
      c1->SaveAs("plots/ks_fitdlife.root");
      c1->SaveAs("plots/ks_fitdlife.C");
    }

    delete c1;

    //    for( int ndx = 0; ndx < kslifeNbins; ndx++ ) {
    //      delete kslifePtrs[ndx];
    //    }
    //    delete kslifePtrs[];
  }


  // Lambda lifetime distribution
  if( dodNdlife_lam ) {
    cout << "Doing dN/dlife..." << endl;

    int lamCanvasX = 8;
    int lamCanvasY = lamlifeNbins / lamCanvasX;
    if( lamlifeNbins%lamCanvasX != 0 ) lamCanvasY += 1;
    TCanvas* lamlifeCanvas = new TCanvas("lamlifeFits", "Lambda mass fits in ct bins",
					300*lamCanvasX, 300*lamCanvasY);
    lamlifeCanvas->Divide( lamCanvasX, lamCanvasY );

    TH1F* lam_lifescale = new TH1F("lam_lifescale", "life scaling factor", lamlifeNbins, &lamlifeBinMinima[0]);
    TH1F* lam_dNdlife = new TH1F("lam_dNdlife", "dN/dct", lamlifeNbins, &lamlifeBinMinima[0]);
    TH1F* lam_dStNdlife = new TH1F("lam_dStNdlife", "(S/N) / dct", lamlifeNbins, &lamlifeBinMinima[0]);
    TH1F* lam_sigdlife = new TH1F("lam_sigdlife", "#sigma / dct", lamlifeNbins, &lamlifeBinMinima[0]);
    TH1F* lam_massdlife = new TH1F("lam_massdlife", "Mass / dct", lamlifeNbins, &lamlifeBinMinima[0]);
    TH1F* lam_fitdlife = new TH1F("lam_fitdlife", "Fit / dct", lamlifeNbins, &lamlifeBinMinima[0]);
    TH1F* lam_chi2dlife = new TH1F("lam_chi2dlife", "#chi^{2}/dof / dct", lamlifeNbins, &lamlifeBinMinima[0]);

    for( int ndx = 0; ndx < lamlifeNbins; ndx++ ) {
      lamlifeCanvas->cd(ndx + 1);
      double min = lamlifeBinMinima[ndx];
      double max = lamlifeBinMinima[ndx+1];
      double binWidth = max - min;
      double scalingFactor = binWidth / lamlife_default_binWidth;
      lam_lifescale->SetBinContent( ndx + 1, scalingFactor );
      lam_lifescale->SetBinError( ndx + 1, 0. );
      cout << "Doing fit for Lambda with " << min << " < life < " << max << endl;

      int fitresult = lamfitobj->Fit(lamlifePtrs[ndx]);
      cout << "Fit result = " << fitresult << endl;

      lam_dNdlife->SetBinContent( ndx + 1, lamfitobj->sig() );
      lam_dNdlife->SetBinError( ndx + 1, lamfitobj->sig_err() );
      double ston =  0.95*lamfitobj->sig()/lamfitobj->bg();
      double stone = ston*sqrt( pow(lamfitobj->sig_err() / lamfitobj->sig(),2)
				+ pow(lamfitobj->bg_err() / lamfitobj->bg(),2) );
      lam_dStNdlife->SetBinContent( ndx + 1, ston );
      lam_dStNdlife->SetBinError( ndx + 1, stone );
      lam_sigdlife->SetBinContent( ndx + 1, lamfitobj->res() );
      lam_sigdlife->SetBinError( ndx + 1, lamfitobj->res_err() );
      lam_massdlife->SetBinContent( ndx + 1, lamfitobj->mass() );
      lam_massdlife->SetBinError( ndx + 1, lamfitobj->mass_err() );
      lam_fitdlife->SetBinContent( ndx + 1, lamfitobj->fitbg() );
      lam_chi2dlife->SetBinContent( ndx + 1, lamfitobj->nChi2() );

      lamlifePtrs[ndx]->Draw();
      cout << "-----------------" << endl;
    }
    if( isMC ) {
      lamlifeCanvas->SaveAs("plots/LamMassFitsInlifeBins_MC_" + mcname + ".png");
      lamlifeCanvas->SaveAs("plots/LamMassFitsInlifeBins_MC_" + mcname + ".eps");
      lamlifeCanvas->SaveAs("plots/LamMassFitsInlifeBins_MC_" + mcname + ".root");
      lamlifeCanvas->SaveAs("plots/LamMassFitsInlifeBins_MC_" + mcname + ".pdf");
    }
    else {
      lamlifeCanvas->SaveAs("plots/LamMassFitsInlifeBins.png");
      lamlifeCanvas->SaveAs("plots/LamMassFitsInlifeBins.eps");
      lamlifeCanvas->SaveAs("plots/LamMassFitsInlifeBins.root");
      lamlifeCanvas->SaveAs("plots/LamMassFitsInlifeBins.pdf");
    }
    delete lamlifeCanvas;

    TCanvas* c1 = new TCanvas("c1", "Lambda yields");

    gStyle->SetOptTitle(0);

    fout->cd();
    TH1F* lam_dNdlife_scale = lam_dNdlife->Clone("lam_dNdlife_scale");
    lam_dNdlife_scale->Divide(lam_lifescale);

    lam_dNdlife_scale->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
    lam_dNdlife_scale->SetYTitle("Candidates / 0.1 cm");
    lam_dNdlife_scale->SetMinimum(0);
    lam_dNdlife_scale->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_dNdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_dNdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_dNdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_dNdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_dNdlife.png");
      c1->SaveAs("plots/lam_dNdlife.eps");
      c1->SaveAs("plots/lam_dNdlife.root");
      c1->SaveAs("plots/lam_dNdlife.C");
    }

    lam_dStNdlife->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
    lam_dStNdlife->SetYTitle("Signal-to-noise");
    lam_dStNdlife->SetMinimum(0);
    lam_dStNdlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_dStNdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_dStNdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_dStNdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_dStNdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_dStNdlife.png");
      c1->SaveAs("plots/lam_dStNdlife.eps");
      c1->SaveAs("plots/lam_dStNdlife.root");
      c1->SaveAs("plots/lam_dStNdlife.C");
    }

    lam_sigdlife->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
    lam_sigdlife->SetYTitle("#sigma");
    lam_sigdlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_sigdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_sigdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_sigdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_sigdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_sigdlife.png");
      c1->SaveAs("plots/lam_sigdlife.eps");
      c1->SaveAs("plots/lam_sigdlife.root");
      c1->SaveAs("plots/lam_sigdlife.C");
    }

    lam_chi2dlife->SetXTitle("#Lambda^{0} lifetime  [ct]");
    lam_chi2dlife->SetYTitle("#chi^2/dof");
    lam_chi2dlife->SetMinimum(0);
    lam_chi2dlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_chi2dlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_chi2dlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_chi2dlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_chi2dlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_chi2dlife.png");
      c1->SaveAs("plots/lam_chi2dlife.eps");
      c1->SaveAs("plots/lam_chi2dlife.root");
      c1->SaveAs("plots/lam_chi2dlife.C");
    }

    lam_massdlife->SetXTitle("#Lambda^{0} lifetime [ct]");
    lam_massdlife->SetYTitle("Mass  [GeV/c^2]");
    lam_massdlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_massdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_massdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_massdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_massdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_massdlife.png");
      c1->SaveAs("plots/lam_massdlife.eps");
      c1->SaveAs("plots/lam_massdlife.root");
      c1->SaveAs("plots/lam_massdlife.C");
    }

    lam_fitdlife->SetXTitle("#Lambda^{0} lifetime  [ct]");
    lam_fitdlife->SetYTitle("Fit type");
    lam_fitdlife->Draw();
    //labl->Draw("same");
    if( isMC ) {
      c1->SaveAs("plots/lam_fitdlife_MC_" + mcname + ".png");
      c1->SaveAs("plots/lam_fitdlife_MC_" + mcname + ".eps");
      c1->SaveAs("plots/lam_fitdlife_MC_" + mcname + ".root");
      c1->SaveAs("plots/lam_fitdlife_MC_" + mcname + ".C");
    }
    else {
      c1->SaveAs("plots/lam_fitdlife.png");
      c1->SaveAs("plots/lam_fitdlife.eps");
      c1->SaveAs("plots/lam_fitdlife.root");
      c1->SaveAs("plots/lam_fitdlife.C");
    }

    delete c1;

    //    for( int ndx = 0; ndx < lamlifeNbins; ndx++ ) {
    //      delete lamlifePtrs[ndx];
    //    }
    //    delete lamlifePtrs[];
  }

  fout->Write();

}


