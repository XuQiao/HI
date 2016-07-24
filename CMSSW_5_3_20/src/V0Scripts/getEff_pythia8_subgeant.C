{
#include <tdrstyle.C>
#include <sstream>
//#include <iostream>
#include <iomanip.h>
#include <string>
#include <vector>
#include "ks_efficiency_pythia8.C"
#include "lam_efficiency_pythia8_subgeant.C"
  namespace std;

  // Switches
  bool debug = false;
  bool kshort = false;

  int mc_choice = 2; // 1=D6T, 2=pythia8, 3=P0
                     // add 10 for 900 GeV
                     // add 100 for pt,y weighting

  int weight_choice = 0;
  vector< TString > variables(3);
  variables[0] = "y";
  variables[1] = "pT";
  variables[2] = "life";

  vector< TString > formats(4);
  formats[0] = "png";
  formats[1] = "eps";
  formats[2] = "root";
  formats[3] = "C";

//  vector< vector<TString> > output_rec_plots(3, vector<TString>(4));
//  vector< vector<TString> > output_rec_plots(3, vector<TString>(4));
//  vector< vector<TString> > output_rec_plots(3, vector<TString>(4));
  TString output_rec_plots[3][4];
  TString output_gen_plots[3][4];
  TString output_eff_plots[3][4];

  TString mctype = "def";
  TString energy = "7TeV";
  if ( mc_choice >= 10 ) energy = "900GeV";
  TString v0type = "ks";
  if (!kshort) v0type = "lam";

  if ( mc_choice == 1 ) { // D6T tune
    mctype = "tuneD6T";
    TFile* fin = TFile::Open("/nfs/data37/cms/v0ntuple/mc_7TeV_tuneD6T_50M_START3X_V26B_june3fixes_PARTIAL.root");
    weight_choice = 101;
    TString weightchoice = "101";

  } else if ( mc_choice == 2 ) { // pythia8
    mctype = "pythia8";
    TFile* fin = TFile::Open("/nfs/data37/cms/v0ntuple/mc_7TeV_pythia8_START3X_V26B_june3fixes_PARTIAL.root");
    weight_choice = 102;
    TString weightchoice = "102";

  } else if ( mc_choice == 3 ) { // p0
    mctype = "tuneP0";
    TFile* fin = TFile::Open("/nfs/data37/cms/v0ntuple/mc_7TeV_tuneP0_START3X_V26B_june3fixes.root");
    weight_choice = 103;
    TString weightchoice = "103";

  } else if ( mc_choice == 11 ) { // 900 GeV D6T tune
    mctype = "tuneD6T";
    TFile* fin = TFile::Open("/nfs/data37/cms/v0ntuple/mc_900GeV_START3X_V22A_june4fixes.root");
    weight_choice = 111;
    TString weightchoice = "111";
  }

  for ( Int_t i = 0; i < variables.size(); ++i ) {
    for ( Int_t j = 0; j < formats.size(); ++j ) {
      output_rec_plots[i][j] = "plots/" + v0type + "_dNd" + variables[i] + "_rec_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant." + formats[j];
      output_gen_plots[i][j] = "plots/" + v0type + "_dNd" + variables[i] + "_gen_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant." + formats[j];
      output_eff_plots[i][j] = "plots/" + v0type + "_dNd" + variables[i] + "_eff_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant." + formats[j];
    }
  }

  TString outFileName = "eff_dNdX_" + v0type + "_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant.root";

  // Style and other ROOT commands
  setTDRStyle();
  gStyle->SetOptStat(1);
  TGaxis::SetMaxDigits(3);
  gROOT->SetBatch(kTRUE);


  TFile* fout = TFile::Open(outFileName, "RECREATE");

  // Constants
  double piMass = 0.13957018;
  double protonMass = 0.93827203;
  double ksMassConst = 0.497614;
  double lamMassConst = 1.115683;

  // TPaves for labels
  TText* labltxt1_1;
  TText* labltxt1_2;
  TText* labltxt2_1;
  TText* labltxt2_2;

  TPaveText* labl = new TPaveText(0.25, 0.25, 0.55, 0.38, "brNDC");
  labl->SetBorderSize(0);
  labl->SetFillColor(0);
  labl->SetFillStyle(0);
  labl->SetTextSizePixels(20);
  TPaveText* labl2 = new TPaveText(0.25, 0.79, 0.55, 0.92, "brNDC");
  labl2->SetBorderSize(0);
  labl2->SetFillColor(0);
  labl2->SetFillStyle(0);
  labl2->SetTextSizePixels(20);

  labltxt1_1 = labl->AddText("CMS Simulation");
  labltxt1_2 = labl->AddText("#sqrt{s} = 7 TeV");
  labltxt2_1 = labl2->AddText("CMS Simulation");
  labltxt2_2 = labl2->AddText("#sqrt{s} = 7 TeV");

  // Histogram quantities
  Double_t ksyXmin = 0.;//-2.5;
  Double_t ksyXmax = 2.3;
  int ksyNbins = 23;
  Double_t ksyBinWidth = (ksyXmax - ksyXmin) / ksyNbins;

  Double_t lamyXmin = 0.;
  Double_t lamyXmax = 2.3;
  int lamyNbins = 23;
  Double_t lamyBinWidth = (lamyXmax - lamyXmin) / lamyNbins;

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


 gStyle->SetOptTitle(0);

 // Find number of events.
 TH1F* veematch = new TH1F("veematch", "veematch", 5, 0., 5);
 if ( kshort ) {
   cout << "Doing Kshort..." << endl;

   TH1F* ks_dNdy_rec = new TH1F("ks_dNdy_rec", "dN/dy rec", ksyNbins, ksyXmin, ksyXmax);
   TH1F* ks_dNdy_gen = new TH1F("ks_dNdy_gen", "dN/dy gen", ksyNbins, ksyXmin, ksyXmax);
   TH1F* ks_dNdy_eff = new TH1F("ks_dNdy_eff", "dN/dy eff", ksyNbins, ksyXmin, ksyXmax);

   TH1F* ks_ptscale = new TH1F("ks_ptscale", "pT scaling factor", ksptNbins, &ksptBinMinima[0]);
   for( int ndx = 0; ndx < ksptNbins; ndx++ ) {
     double min = ksptBinMinima[ndx];
     double max = ksptBinMinima[ndx+1];
     double binWidth = max - min;
     double scalingFactor = binWidth / kspt_default_binWidth;
     ks_ptscale->SetBinContent( ndx + 1, scalingFactor );
     ks_ptscale->SetBinError( ndx + 1, 0. );
   }
   TH1F* ks_dNdpT_rec = new TH1F("ks_dNdpT_rec", "dN/dp_{T} rec", ksptNbins, &ksptBinMinima[0]);
   TH1F* ks_dNdpT_gen = new TH1F("ks_dNdpT_gen", "dN/dp_{T} gen", ksptNbins, &ksptBinMinima[0]);
   TH1F* ks_dNdpT_eff = new TH1F("ks_dNdpT_eff", "dN/dp_{T} eff", ksptNbins, &ksptBinMinima[0]);

   TH1F* ks_lifescale = new TH1F("ks_lifescale", "life scaling factor", kslifeNbins, &kslifeBinMinima[0]);
   for( int ndx = 0; ndx < kslifeNbins; ndx++ ) {
     double min = kslifeBinMinima[ndx];
     double max = kslifeBinMinima[ndx+1];
     double binWidth = max - min;
     double scalingFactor = binWidth / kslife_default_binWidth;
     ks_lifescale->SetBinContent( ndx + 1, scalingFactor );
     ks_lifescale->SetBinError( ndx + 1, 0. );
   }
   TH1F* ks_dNdlife_rec = new TH1F("ks_dNdlife_rec", "dN/dct rec", kslifeNbins, &kslifeBinMinima[0]);
   TH1F* ks_dNdlife_gen = new TH1F("ks_dNdlife_gen", "dN/dct gen", kslifeNbins, &kslifeBinMinima[0]);
   TH1F* ks_dNdlife_eff = new TH1F("ks_dNdlife_eff", "dN/dct eff", kslifeNbins, &kslifeBinMinima[0]);

 } 
 else {
   cout << "Doing Lambda..." << endl;
 
   TH1F* lam_dNdy_rec = new TH1F("lam_dNdy_rec", "dN/dy rec", lamyNbins, lamyXmin, lamyXmax);
   TH1F* lam_dNdy_gen = new TH1F("lam_dNdy_gen", "dN/dy gen", lamyNbins, lamyXmin, lamyXmax);
   TH1F* lam_dNdy_eff = new TH1F("lam_dNdy_eff", "dN/dy eff", lamyNbins, lamyXmin, lamyXmax);

   TH1F* lam_ptscale = new TH1F("lam_ptscale", "pT scaling factor", lamptNbins, &lamptBinMinima[0]);
   for( int ndx = 0; ndx < lamptNbins; ndx++ ) {
     double min = lamptBinMinima[ndx];
     double max = lamptBinMinima[ndx+1];
     double binWidth = max - min;
     double scalingFactor = binWidth / lampt_default_binWidth;
     lam_ptscale->SetBinContent( ndx + 1, scalingFactor );
     lam_ptscale->SetBinError( ndx + 1, 0. );
   }
   TH1F* lam_dNdpT_rec = new TH1F("lam_dNdpT_rec", "dN/dp_{T} rec", lamptNbins, &lamptBinMinima[0]);
   TH1F* lam_dNdpT_gen = new TH1F("lam_dNdpT_gen", "dN/dp_{T} gen", lamptNbins, &lamptBinMinima[0]);
   TH1F* lam_dNdpT_eff = new TH1F("lam_dNdpT_eff", "dN/dp_{T} eff", lamptNbins, &lamptBinMinima[0]);

   TH1F* lam_lifescale = new TH1F("lam_lifescale", "life scaling factor", lamlifeNbins, &lamlifeBinMinima[0]);
   for( int ndx = 0; ndx < lamlifeNbins; ndx++ ) {
     double min = lamlifeBinMinima[ndx];
     double max = lamlifeBinMinima[ndx+1];
     double binWidth = max - min;
     double scalingFactor = binWidth / lamlife_default_binWidth;
     lam_lifescale->SetBinContent( ndx + 1, scalingFactor );
     lam_lifescale->SetBinError( ndx + 1, 0. );
   }
   TH1F* lam_dNdlife_rec = new TH1F("lam_dNdlife_rec", "dN/dct rec", lamlifeNbins, &lamlifeBinMinima[0]);
   TH1F* lam_dNdlife_gen = new TH1F("lam_dNdlife_gen", "dN/dct gen", lamlifeNbins, &lamlifeBinMinima[0]);
   TH1F* lam_dNdlife_eff = new TH1F("lam_dNdlife_eff", "dN/dct eff", lamlifeNbins, &lamlifeBinMinima[0]);

 }



#include "ntrk_wgt.C"
#include "kspt_wgt.C"
#include "lampt_wgt.C"
#include "ksy_wgt.C"
#include "lamy_wgt.C"

 fin->cd("matchVeesByDeltaR");
 TTree* ksTree = (TTree*)gDirectory->Get("kseff");
 TTree* lamTree = (TTree*)gDirectory->Get("lameff");
 if( kshort ) {
   ks_efficiency keff(ksTree);
   keff.Loop(weight_choice);
 }
 else {
   lam_efficiency leff(lamTree);
   leff.Loop(weight_choice);
 }

 
 Double_t nVee = veematch->GetEntries();
 cout << "Precision: " << cout.precision() << endl;
 cout.precision(0);
 cout.setf(ios::fixed);
 cout << "Running over " << nVee << " generated vees..." << endl;
 cout.unsetf(ios::fixed);
 cout.precision(6);
  
 fout->cd();

 if ( kshort ) {
   TCanvas* c1 = new TCanvas("c1", "K0S results");
   
   ks_dNdy_rec->Sumw2();
   ks_dNdy_gen->Sumw2();
   ks_dNdy_eff->Divide(ks_dNdy_rec,ks_dNdy_gen,1,1,"B");

   TH1F* ks_dNdy_rec_copy = ks_dNdy_rec->Clone("ks_dNdy_rec_copy");
   TH1F* ks_dNdy_gen_copy = ks_dNdy_gen->Clone("ks_dNdy_gen_copy");

   ks_dNdy_rec->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_rec->SetYTitle("Candidates / 0.1");
   ks_dNdy_rec->SetNdivisions(506);
   ks_dNdy_rec->SetMinimum(0);
   ks_dNdy_rec->Draw();
   labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[0][i]);
   }

   ks_dNdy_gen->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_gen->SetYTitle("Candidates / 0.1");
   ks_dNdy_gen->SetNdivisions(506);
   ks_dNdy_gen->SetMinimum(0);
   ks_dNdy_gen->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[0][i]);
   }

   ks_dNdy_eff->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_eff->SetYTitle("Efficiency");
   ks_dNdy_eff->SetNdivisions(506);
   ks_dNdy_eff->SetMinimum(0);
   ks_dNdy_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[0][i]);
   }


//pt

   ks_dNdpT_rec->Sumw2();
   ks_dNdpT_gen->Sumw2();

   TH1F* ks_dNdpT_rec_scale = ks_dNdpT_rec->Clone("ks_dNdpT_rec_scale");
   ks_dNdpT_rec_scale->Divide(ks_ptscale);
   TH1F* ks_dNdpT_gen_scale = ks_dNdpT_gen->Clone("ks_dNdpT_gen_scale");
   ks_dNdpT_gen_scale->Divide(ks_ptscale);

   ks_dNdpT_eff->Divide(ks_dNdpT_rec,ks_dNdpT_gen,1,1,"B");

   ks_dNdpT_rec_scale->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
   ks_dNdpT_rec_scale->SetYTitle("Candidates / 0.1 GeV/c");
   ks_dNdpT_rec_scale->SetMinimum(0);
   ks_dNdpT_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[1][i]);
   }

   ks_dNdpT_gen_scale->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
   ks_dNdpT_gen_scale->SetYTitle("Candidates / 0.1 GeV/c");
   ks_dNdpT_gen_scale->SetMinimum(0);
   ks_dNdpT_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[1][i]);
   }

   ks_dNdpT_eff->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
   ks_dNdpT_eff->SetYTitle("Efficiency");
   ks_dNdpT_eff->SetMinimum(0);
   ks_dNdpT_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[1][i]);
   }


// lifetime

   ks_dNdlife_rec->Sumw2();
   ks_dNdlife_gen->Sumw2();

   TH1F* ks_dNdlife_rec_scale = ks_dNdlife_rec->Clone("ks_dNdlife_rec_scale");
   ks_dNdlife_rec_scale->Divide(ks_lifescale);
   TH1F* ks_dNdlife_gen_scale = ks_dNdlife_gen->Clone("ks_dNdlife_gen_scale");
   ks_dNdlife_gen_scale->Divide(ks_lifescale);

   ks_dNdlife_eff->Divide(ks_dNdlife_rec,ks_dNdlife_gen,1,1,"B");

   ks_dNdlife_rec_scale->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
   ks_dNdlife_rec_scale->SetYTitle("Candidates / 0.1 cm");
   ks_dNdlife_rec_scale->SetMinimum(0);
   ks_dNdlife_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[2][i]);
   }

   ks_dNdlife_gen_scale->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
   ks_dNdlife_gen_scale->SetYTitle("Candidates / 0.1 cm");
   ks_dNdlife_gen_scale->SetMinimum(0);
   ks_dNdlife_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[2][i]);
   }

   ks_dNdlife_eff->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
   ks_dNdlife_eff->SetYTitle("Efficiency");
   ks_dNdlife_eff->SetMinimum(0);
   ks_dNdlife_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[2][i]);
   }


 } else {


   TCanvas* c1 = new TCanvas("c1", "Lambda results");
   
   lam_dNdy_rec->Sumw2();
   lam_dNdy_gen->Sumw2();
   lam_dNdy_eff->Divide(lam_dNdy_rec,lam_dNdy_gen,1,1,"B");

   TH1F* lam_dNdy_rec_copy = lam_dNdy_rec->Clone("lam_dNdy_rec_copy");
   TH1F* lam_dNdy_gen_copy = lam_dNdy_gen->Clone("lam_dNdy_gen_copy");

   lam_dNdy_rec->SetXTitle("#Lambda^{0} rapidity");
   lam_dNdy_rec->SetYTitle("Candidates / 0.1");
   lam_dNdy_rec->SetNdivisions(506);
   lam_dNdy_rec->SetMinimum(0);
   lam_dNdy_rec->Draw();
   labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[0][i]);
   }

   lam_dNdy_gen->SetXTitle("#Lambda^{0} rapidity");
   lam_dNdy_gen->SetYTitle("Candidates / 0.1");
   lam_dNdy_gen->SetNdivisions(506);
   lam_dNdy_gen->SetMinimum(0);
   lam_dNdy_gen->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[0][i]);
   }

   lam_dNdy_eff->SetXTitle("#Lambda^{0} rapidity");
   lam_dNdy_eff->SetYTitle("Efficiency");
   lam_dNdy_eff->SetNdivisions(506);
   lam_dNdy_eff->SetMinimum(0);
   lam_dNdy_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[0][i]);
   }


//pt

   lam_dNdpT_rec->Sumw2();
   lam_dNdpT_gen->Sumw2();

   TH1F* lam_dNdpT_rec_scale = lam_dNdpT_rec->Clone("lam_dNdpT_rec_scale");
   lam_dNdpT_rec_scale->Divide(lam_ptscale);
   TH1F* lam_dNdpT_gen_scale = lam_dNdpT_gen->Clone("lam_dNdpT_gen_scale");
   lam_dNdpT_gen_scale->Divide(lam_ptscale);

   lam_dNdpT_eff->Divide(lam_dNdpT_rec,lam_dNdpT_gen,1,1,"B");

   lam_dNdpT_rec_scale->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
   lam_dNdpT_rec_scale->SetYTitle("Candidates / 0.1 GeV/c");
   lam_dNdpT_rec_scale->SetMinimum(0);
   lam_dNdpT_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[1][i]);
   }

   lam_dNdpT_gen_scale->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
   lam_dNdpT_gen_scale->SetYTitle("Candidates / 0.1 GeV/c");
   lam_dNdpT_gen_scale->SetMinimum(0);
   lam_dNdpT_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[1][i]);
   }

   lam_dNdpT_eff->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
   lam_dNdpT_eff->SetYTitle("Efficiency");
   lam_dNdpT_eff->SetMinimum(0);
   lam_dNdpT_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[1][i]);
   }


// lifetime

   lam_dNdlife_rec->Sumw2();
   lam_dNdlife_gen->Sumw2();

   TH1F* lam_dNdlife_rec_scale = lam_dNdlife_rec->Clone("lam_dNdlife_rec_scale");
   lam_dNdlife_rec_scale->Divide(lam_lifescale);
   TH1F* lam_dNdlife_gen_scale = lam_dNdlife_gen->Clone("lam_dNdlife_gen_scale");
   lam_dNdlife_gen_scale->Divide(lam_lifescale);

   lam_dNdlife_eff->Divide(lam_dNdlife_rec,lam_dNdlife_gen,1,1,"B");

   lam_dNdlife_rec_scale->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
   lam_dNdlife_rec_scale->SetYTitle("Candidates / 0.1 cm");
   lam_dNdlife_rec_scale->SetMinimum(0);
   lam_dNdlife_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[2][i]);
   }

   lam_dNdlife_gen_scale->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
   lam_dNdlife_gen_scale->SetYTitle("Candidates / 0.1 cm");
   lam_dNdlife_gen_scale->SetMinimum(0);
   lam_dNdlife_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[2][i]);
   }

   lam_dNdlife_eff->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
   lam_dNdlife_eff->SetYTitle("Efficiency");
   lam_dNdlife_eff->SetMinimum(0);
   lam_dNdlife_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[2][i]);
   }


 }

  fout->Write();


}
