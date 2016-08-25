{
#include <tdrstyle.C>
#include <sstream>
//#include <iostream>
#include <iomanip.h>
#include <string>
#include <vector>
//#include "ks_efficiency.C"
//#include "lam_efficiency.C"
#include "constants.h"
//#include "ntrk_wgt.C"
//#include "kspt_wgt.C"
//#include "lampt_wgt.C"
//#include "ksy_wgt.C"
//#include "lamy_wgt.C"

//  .L ks_efficiency.C+
//  .L lam_efficiency.C+

using namespace std;
 cout << "TEst" << endl;
Int_t status_ks = gSystem->Load("./ks_efficiency_C.so");
Int_t status_lam = gSystem->Load("./lam_efficiency_subgeant2_C.so");

cout << "stat = " << status_ks << " and " << status_lam << endl;
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

  TString part_output_rec_plots[3][4];
  TString part_output_gen_plots[3][4];
  TString part_output_eff_plots[3][4];

  TString anti_output_rec_plots[3][4];
  TString anti_output_gen_plots[3][4];
  TString anti_output_eff_plots[3][4];

  TString mctype = "def";
  TString energy = "7TeV";
  if ( mc_choice >= 10 ) energy = "900GeV";
  TString v0type = "ks";
  if (!kshort) v0type = "lam";

  if ( mc_choice == 1 ) { // D6T tune
    mctype = "tuneD6T";
    TFile* fin = TFile::Open("/nfs/data37/cms/v0ntuple/mc_7TeV_tuneD6T_50M_START3X_V26B_june4fixesGood.root");
    weight_choice = 101;
    TString weightchoice = "101";

  } else if ( mc_choice == 2 ) { // pythia8
    mctype = "pythia8";
    TFile* fin = TFile::Open("/nfs/data37/cms/v0ntuple/mc_7TeV_pythia8_START3X_V26B_june4fixesGood.root");
    weight_choice = 0;
    TString weightchoice = "0";

  } else if ( mc_choice == 3 ) { // p0
    mctype = "tuneP0";
    TFile* fin = TFile::Open("/nfs/data37/cms/v0ntuple/mc_7TeV_tuneP0_START3X_V26B_june4fixesGood.root");
    weight_choice = 103;
    TString weightchoice = "103";

  } else if ( mc_choice == 11 ) { // 900 GeV D6T tune
    mctype = "tuneD6T";
    TFile* fin = TFile::Open("/nfs/data37/cms/v0ntuple/mc_900GeV_START3X_V22A_june4fixesGood_PARTIAL.root");
    weight_choice = 111;
    TString weightchoice = "111";
  }

  for ( Int_t i = 0; i < variables.size(); ++i ) {
    for ( Int_t j = 0; j < formats.size(); ++j ) {
      output_rec_plots[i][j] = "plots/" + v0type + "_dNd" + variables[i] + "_rec_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
      output_gen_plots[i][j] = "plots/" + v0type + "_dNd" + variables[i] + "_gen_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
      output_eff_plots[i][j] = "plots/" + v0type + "_dNd" + variables[i] + "_eff_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
    }
  }

if ( v0type == "lam" ) {

  for ( Int_t i = 0; i < variables.size(); ++i ) {
    for ( Int_t j = 0; j < formats.size(); ++j ) {
      part_output_rec_plots[i][j] = "plots/" + v0type + "Reg_dNd" + variables[i] + "_rec_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
      part_output_gen_plots[i][j] = "plots/" + v0type + "Reg_dNd" + variables[i] + "_gen_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
      part_output_eff_plots[i][j] = "plots/" + v0type + "Reg_dNd" + variables[i] + "_eff_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
      anti_output_rec_plots[i][j] = "plots/" + v0type + "Bar_dNd" + variables[i] + "_rec_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
      anti_output_gen_plots[i][j] = "plots/" + v0type + "Bar_dNd" + variables[i] + "_gen_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
      anti_output_eff_plots[i][j] = "plots/" + v0type + "Bar_dNd" + variables[i] + "_eff_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2." + formats[j];
    }
  }
}

  TString outFileName = "eff_dNdX_" + v0type + "_" + mctype + "_" + energy + "_weight" + weightchoice + "_subgeant2.root";

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


 gStyle->SetOptTitle(0);

 // Find number of events.
 TH1F* veematch = new TH1F("veematch", "veematch", 5, 0., 5);
 if ( kshort ) {
   cout << "Doing Kshort..." << endl;

   TH1F* ks_dNdy_rec = new TH1F("ks_dNdy_rec", "dN/dy rec", ksyNbins, ksyXmin, ksyXmax);
   TH1F* ks_dNdy_gen = new TH1F("ks_dNdy_gen", "dN/dy gen", ksyNbins, ksyXmin, ksyXmax);
   TH1F* ks_dNdy_eff = new TH1F("ks_dNdy_eff", "dN/dy eff", ksyNbins, ksyXmin, ksyXmax);
   ks_dNdy_rec->Sumw2();
   ks_dNdy_gen->Sumw2();
   ks_dNdy_rec->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_rec->SetYTitle("K^{0}_{S} / 0.1");
   ks_dNdy_rec->SetNdivisions(506);
   ks_dNdy_rec->SetMinimum(0);
   ks_dNdy_gen->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_gen->SetYTitle("K^{0}_{S} / 0.1");
   ks_dNdy_gen->SetNdivisions(506);
   ks_dNdy_gen->SetMinimum(0);
   ks_dNdy_eff->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_eff->SetYTitle("Efficiency");
   ks_dNdy_eff->SetNdivisions(506);
   ks_dNdy_eff->SetMinimum(0);

   TH1F* ks_dNdy_neg_rec = new TH1F("ks_dNdy_neg_rec", "dN/dy neg rec", ksyNbins, ksyXmin, ksyXmax);
   TH1F* ks_dNdy_neg_gen = new TH1F("ks_dNdy_neg_gen", "dN/dy neg gen", ksyNbins, ksyXmin, ksyXmax);
   TH1F* ks_dNdy_neg_eff = new TH1F("ks_dNdy_neg_eff", "dN/dy neg eff", ksyNbins, ksyXmin, ksyXmax);
   ks_dNdy_neg_rec->Sumw2();
   ks_dNdy_neg_gen->Sumw2();
   ks_dNdy_neg_rec->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_neg_rec->SetYTitle("K^{0}_{S} / 0.1");
   ks_dNdy_neg_rec->SetNdivisions(506);
   ks_dNdy_neg_rec->SetMinimum(0);
   ks_dNdy_neg_gen->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_neg_gen->SetYTitle("K^{0}_{S} / 0.1");
   ks_dNdy_neg_gen->SetNdivisions(506);
   ks_dNdy_neg_gen->SetMinimum(0);
   ks_dNdy_neg_eff->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_neg_eff->SetYTitle("Efficiency");
   ks_dNdy_neg_eff->SetNdivisions(506);
   ks_dNdy_neg_eff->SetMinimum(0);

   TH1F* ks_dNdy_pos_rec = new TH1F("ks_dNdy_pos_rec", "dN/dy pos rec", ksyNbins, ksyXmin, ksyXmax);
   TH1F* ks_dNdy_pos_gen = new TH1F("ks_dNdy_pos_gen", "dN/dy pos gen", ksyNbins, ksyXmin, ksyXmax);
   TH1F* ks_dNdy_pos_eff = new TH1F("ks_dNdy_pos_eff", "dN/dy pos eff", ksyNbins, ksyXmin, ksyXmax);
   ks_dNdy_pos_rec->Sumw2();
   ks_dNdy_pos_gen->Sumw2();
   ks_dNdy_pos_rec->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_pos_rec->SetYTitle("K^{0}_{S} / 0.1");
   ks_dNdy_pos_rec->SetNdivisions(506);
   ks_dNdy_pos_rec->SetMinimum(0);
   ks_dNdy_pos_gen->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_pos_gen->SetYTitle("K^{0}_{S} / 0.1");
   ks_dNdy_pos_gen->SetNdivisions(506);
   ks_dNdy_pos_gen->SetMinimum(0);
   ks_dNdy_pos_eff->SetXTitle("K^{0}_{S} rapidity");
   ks_dNdy_pos_eff->SetYTitle("Efficiency");
   ks_dNdy_pos_eff->SetNdivisions(506);
   ks_dNdy_pos_eff->SetMinimum(0);

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
   ks_dNdpT_rec->Sumw2();
   ks_dNdpT_gen->Sumw2();
   ks_dNdpT_rec->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
   ks_dNdpT_rec->SetYTitle("K^{0}_{S} / 0.1 GeV/c");
   ks_dNdpT_rec->SetMinimum(0);
   ks_dNdpT_gen->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
   ks_dNdpT_gen->SetYTitle("K^{0}_{S} / 0.1 GeV/c");
   ks_dNdpT_gen->SetMinimum(0);
   ks_dNdpT_eff->SetXTitle("K^{0}_{S} transverse momentum [GeV/c]");
   ks_dNdpT_eff->SetYTitle("Efficiency");
   ks_dNdpT_eff->SetMinimum(0);

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
   ks_dNdlife_rec->Sumw2();
   ks_dNdlife_gen->Sumw2();
   ks_dNdlife_rec->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
   ks_dNdlife_rec->SetYTitle("K^{0}_{S} / 0.1 cm");
   ks_dNdlife_rec->SetMinimum(0);
   ks_dNdlife_gen->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
   ks_dNdlife_gen->SetYTitle("K^{0}_{S} / 0.1 cm");
   ks_dNdlife_gen->SetMinimum(0);
   ks_dNdlife_eff->SetXTitle("K^{0}_{S} lifetime (ct)  [cm]");
   ks_dNdlife_eff->SetYTitle("Efficiency");
   ks_dNdlife_eff->SetMinimum(0);

 } 
 else {
   cout << "Doing Lambda..." << endl;
 
   TH1F* lam_dNdy_rec = new TH1F("lam_dNdy_rec", "dN/dy rec", lamyNbins, lamyXmin, lamyXmax);
   TH1F* lam_dNdy_gen = new TH1F("lam_dNdy_gen", "dN/dy gen", lamyNbins, lamyXmin, lamyXmax);
   TH1F* lam_dNdy_eff = new TH1F("lam_dNdy_eff", "dN/dy eff", lamyNbins, lamyXmin, lamyXmax);
   lam_dNdy_rec->Sumw2();
   lam_dNdy_gen->Sumw2();
   lam_dNdy_rec->SetXTitle("#Lambda^{0} rapidity");
   lam_dNdy_rec->SetYTitle("#Lambda^{0} / 0.1");
   lam_dNdy_rec->SetNdivisions(506);
   lam_dNdy_rec->SetMinimum(0);
   lam_dNdy_gen->SetXTitle("#Lambda^{0} rapidity");
   lam_dNdy_gen->SetYTitle("#Lambda^{0} / 0.1");
   lam_dNdy_gen->SetNdivisions(506);
   lam_dNdy_gen->SetMinimum(0);
   lam_dNdy_eff->SetXTitle("#Lambda^{0} rapidity");
   lam_dNdy_eff->SetYTitle("Efficiency");
   lam_dNdy_eff->SetNdivisions(506);
   lam_dNdy_eff->SetMinimum(0);

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
   lam_dNdpT_rec->Sumw2();
   lam_dNdpT_gen->Sumw2();
   lam_dNdpT_rec->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
   lam_dNdpT_rec->SetYTitle("#Lambda^{0} / 0.1 GeV/c");
   lam_dNdpT_rec->SetMinimum(0);
   lam_dNdpT_gen->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
   lam_dNdpT_gen->SetYTitle("#Lambda^{0} / 0.1 GeV/c");
   lam_dNdpT_gen->SetMinimum(0);
   lam_dNdpT_eff->SetXTitle("#Lambda^{0} transverse momentum [GeV/c]");
   lam_dNdpT_eff->SetYTitle("Efficiency");
   lam_dNdpT_eff->SetMinimum(0);


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
   lam_dNdlife_rec->Sumw2();
   lam_dNdlife_gen->Sumw2();
   lam_dNdlife_rec->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
   lam_dNdlife_rec->SetYTitle("#Lambda^{0} / 0.1 cm");
   lam_dNdlife_rec->SetMinimum(0);
   lam_dNdlife_gen->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
   lam_dNdlife_gen->SetYTitle("#Lambda^{0} / 0.1 cm");
   lam_dNdlife_gen->SetMinimum(0);
   lam_dNdlife_eff->SetXTitle("#Lambda^{0} lifetime (ct)  [cm]");
   lam_dNdlife_eff->SetYTitle("Efficiency");
   lam_dNdlife_eff->SetMinimum(0);

   // Lambda (particle)

   TH1F* lamReg_dNdy_rec = new TH1F("lamReg_dNdy_rec", "dN/dy part rec", lamyNbins, lamyXmin, lamyXmax);
   TH1F* lamReg_dNdy_gen = new TH1F("lamReg_dNdy_gen", "dN/dy part gen", lamyNbins, lamyXmin, lamyXmax);
   TH1F* lamReg_dNdy_eff = new TH1F("lamReg_dNdy_eff", "dN/dy part eff", lamyNbins, lamyXmin, lamyXmax);
   lamReg_dNdy_rec->Sumw2();
   lamReg_dNdy_gen->Sumw2();
   lamReg_dNdy_rec->SetXTitle("#Lambda^{0} (part) rapidity");
   lamReg_dNdy_rec->SetYTitle("#Lambda^{0} / 0.1");
   lamReg_dNdy_rec->SetNdivisions(506);
   lamReg_dNdy_rec->SetMinimum(0);
   lamReg_dNdy_gen->SetXTitle("#Lambda^{0} (part) rapidity");
   lamReg_dNdy_gen->SetYTitle("#Lambda^{0} / 0.1");
   lamReg_dNdy_gen->SetNdivisions(506);
   lamReg_dNdy_gen->SetMinimum(0);
   lamReg_dNdy_eff->SetXTitle("#Lambda^{0} (part) rapidity");
   lamReg_dNdy_eff->SetYTitle("Efficiency");
   lamReg_dNdy_eff->SetNdivisions(506);
   lamReg_dNdy_eff->SetMinimum(0);

   TH1F* lamReg_dNdpT_rec = new TH1F("lamReg_dNdpT_rec", "dN/dp_{T} part rec", lamptNbins, &lamptBinMinima[0]);
   TH1F* lamReg_dNdpT_gen = new TH1F("lamReg_dNdpT_gen", "dN/dp_{T} part gen", lamptNbins, &lamptBinMinima[0]);
   TH1F* lamReg_dNdpT_eff = new TH1F("lamReg_dNdpT_eff", "dN/dp_{T} part eff", lamptNbins, &lamptBinMinima[0]);
   lamReg_dNdpT_rec->Sumw2();
   lamReg_dNdpT_gen->Sumw2();
   lamReg_dNdpT_rec->SetXTitle("#Lambda^{0} (part) transverse momentum [GeV/c]");
   lamReg_dNdpT_rec->SetYTitle("#Lambda^{0} / 0.1 GeV/c");
   lamReg_dNdpT_rec->SetMinimum(0);
   lamReg_dNdpT_gen->SetXTitle("#Lambda^{0} (part) transverse momentum [GeV/c]");
   lamReg_dNdpT_gen->SetYTitle("#Lambda^{0} / 0.1 GeV/c");
   lamReg_dNdpT_gen->SetMinimum(0);
   lamReg_dNdpT_eff->SetXTitle("#Lambda^{0} (part) transverse momentum [GeV/c]");
   lamReg_dNdpT_eff->SetYTitle("Efficiency");
   lamReg_dNdpT_eff->SetMinimum(0);


   TH1F* lamReg_dNdlife_rec = new TH1F("lamReg_dNdlife_rec", "dN/dct part rec", lamlifeNbins, &lamlifeBinMinima[0]);
   TH1F* lamReg_dNdlife_gen = new TH1F("lamReg_dNdlife_gen", "dN/dct part gen", lamlifeNbins, &lamlifeBinMinima[0]);
   TH1F* lamReg_dNdlife_eff = new TH1F("lamReg_dNdlife_eff", "dN/dct part eff", lamlifeNbins, &lamlifeBinMinima[0]);
   lamReg_dNdlife_rec->Sumw2();
   lamReg_dNdlife_gen->Sumw2();
   lamReg_dNdlife_rec->SetXTitle("#Lambda^{0} (part) lifetime (ct)  [cm]");
   lamReg_dNdlife_rec->SetYTitle("#Lambda^{0} / 0.1 cm");
   lamReg_dNdlife_rec->SetMinimum(0);
   lamReg_dNdlife_gen->SetXTitle("#Lambda^{0} (part) lifetime (ct)  [cm]");
   lamReg_dNdlife_gen->SetYTitle("#Lambda^{0} / 0.1 cm");
   lamReg_dNdlife_gen->SetMinimum(0);
   lamReg_dNdlife_eff->SetXTitle("#Lambda^{0} (part) lifetime (ct)  [cm]");
   lamReg_dNdlife_eff->SetYTitle("Efficiency");
   lamReg_dNdlife_eff->SetMinimum(0);


   TH1F* lamBar_dNdy_rec = new TH1F("lamBar_dNdy_rec", "dN/dy anti rec", lamyNbins, lamyXmin, lamyXmax);
   TH1F* lamBar_dNdy_gen = new TH1F("lamBar_dNdy_gen", "dN/dy anti gen", lamyNbins, lamyXmin, lamyXmax);
   TH1F* lamBar_dNdy_eff = new TH1F("lamBar_dNdy_eff", "dN/dy anti eff", lamyNbins, lamyXmin, lamyXmax);
   lamBar_dNdy_rec->Sumw2();
   lamBar_dNdy_gen->Sumw2();
   lamBar_dNdy_rec->SetXTitle("#Lambda^{0} (anti) rapidity");
   lamBar_dNdy_rec->SetYTitle("#Lambda^{0} / 0.1");
   lamBar_dNdy_rec->SetNdivisions(506);
   lamBar_dNdy_rec->SetMinimum(0);
   lamBar_dNdy_gen->SetXTitle("#Lambda^{0} (anti) rapidity");
   lamBar_dNdy_gen->SetYTitle("#Lambda^{0} / 0.1");
   lamBar_dNdy_gen->SetNdivisions(506);
   lamBar_dNdy_gen->SetMinimum(0);
   lamBar_dNdy_eff->SetXTitle("#Lambda^{0} (anti) rapidity");
   lamBar_dNdy_eff->SetYTitle("Efficiency");
   lamBar_dNdy_eff->SetNdivisions(506);
   lamBar_dNdy_eff->SetMinimum(0);

   TH1F* lamBar_dNdpT_rec = new TH1F("lamBar_dNdpT_rec", "dN/dp_{T} anti rec", lamptNbins, &lamptBinMinima[0]);
   TH1F* lamBar_dNdpT_gen = new TH1F("lamBar_dNdpT_gen", "dN/dp_{T} anti gen", lamptNbins, &lamptBinMinima[0]);
   TH1F* lamBar_dNdpT_eff = new TH1F("lamBar_dNdpT_eff", "dN/dp_{T} anti eff", lamptNbins, &lamptBinMinima[0]);
   lamBar_dNdpT_rec->Sumw2();
   lamBar_dNdpT_gen->Sumw2();
   lamBar_dNdpT_rec->SetXTitle("#Lambda^{0} (anti) transverse momentum [GeV/c]");
   lamBar_dNdpT_rec->SetYTitle("#Lambda^{0} / 0.1 GeV/c");
   lamBar_dNdpT_rec->SetMinimum(0);
   lamBar_dNdpT_gen->SetXTitle("#Lambda^{0} (anti) transverse momentum [GeV/c]");
   lamBar_dNdpT_gen->SetYTitle("#Lambda^{0} / 0.1 GeV/c");
   lamBar_dNdpT_gen->SetMinimum(0);
   lamBar_dNdpT_eff->SetXTitle("#Lambda^{0} (anti) transverse momentum [GeV/c]");
   lamBar_dNdpT_eff->SetYTitle("Efficiency");
   lamBar_dNdpT_eff->SetMinimum(0);


   TH1F* lamBar_dNdlife_rec = new TH1F("lamBar_dNdlife_rec", "dN/dct anti rec", lamlifeNbins, &lamlifeBinMinima[0]);
   TH1F* lamBar_dNdlife_gen = new TH1F("lamBar_dNdlife_gen", "dN/dct anti gen", lamlifeNbins, &lamlifeBinMinima[0]);
   TH1F* lamBar_dNdlife_eff = new TH1F("lamBar_dNdlife_eff", "dN/dct anti eff", lamlifeNbins, &lamlifeBinMinima[0]);
   lamBar_dNdlife_rec->Sumw2();
   lamBar_dNdlife_gen->Sumw2();
   lamBar_dNdlife_rec->SetXTitle("#Lambda^{0} (anti) lifetime (ct)  [cm]");
   lamBar_dNdlife_rec->SetYTitle("#Lambda^{0} / 0.1 cm");
   lamBar_dNdlife_rec->SetMinimum(0);
   lamBar_dNdlife_gen->SetXTitle("#Lambda^{0} (anti) lifetime (ct)  [cm]");
   lamBar_dNdlife_gen->SetYTitle("#Lambda^{0} / 0.1 cm");
   lamBar_dNdlife_gen->SetMinimum(0);
   lamBar_dNdlife_eff->SetXTitle("#Lambda^{0} (anti) lifetime (ct)  [cm]");
   lamBar_dNdlife_eff->SetYTitle("Efficiency");
   lamBar_dNdlife_eff->SetMinimum(0);

 }



 fin->cd("matchVeesByDeltaR");
 TTree* ksTree;
 ks_efficiency *keff;
 TTree* lamTree;
 lam_efficiency_subgeant2 *leff;

 if( kshort ) {
   ksTree = (TTree*)gDirectory->Get("kseff");
   keff = new ks_efficiency(ksTree);
   keff->Loop(weight_choice, veematch, ks_dNdy_gen, ks_dNdy_rec, ks_dNdpT_gen, ks_dNdpT_rec, ks_dNdlife_gen, ks_dNdlife_rec, ks_dNdy_neg_gen, ks_dNdy_neg_rec, ks_dNdy_pos_gen, ks_dNdy_pos_rec);
   delete keff;
 }
 else {
   lamTree = (TTree*)gDirectory->Get("lameff");
   leff = new lam_efficiency_subgeant2(lamTree);
   leff->Loop(weight_choice, veematch, lam_dNdy_gen, lam_dNdy_rec, lam_dNdpT_gen, lam_dNdpT_rec, lam_dNdlife_gen, lam_dNdlife_rec, lamReg_dNdy_gen, lamReg_dNdy_rec, lamReg_dNdpT_gen, lamReg_dNdpT_rec, lamReg_dNdlife_gen, lamReg_dNdlife_rec, lamBar_dNdy_gen, lamBar_dNdy_rec, lamBar_dNdpT_gen, lamBar_dNdpT_rec, lamBar_dNdlife_gen, lamBar_dNdlife_rec);
   delete leff;
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
   
   ks_dNdy_eff->Divide(ks_dNdy_rec,ks_dNdy_gen);
   TH1F* ks_dNdy_rec_copy = ks_dNdy_rec->Clone("ks_dNdy_rec_copy");
   TH1F* ks_dNdy_gen_copy = ks_dNdy_gen->Clone("ks_dNdy_gen_copy");

   ks_dNdy_rec->Draw();
   labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[0][i]);
   }
   ks_dNdy_gen->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[0][i]);
   }
   ks_dNdy_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[0][i]);
   }

   ks_dNdy_neg_eff->Divide(ks_dNdy_neg_rec,ks_dNdy_neg_gen);
   TH1F* ks_dNdy_neg_rec_copy = ks_dNdy_neg_rec->Clone("ks_dNdy_neg_rec_copy");
   TH1F* ks_dNdy_neg_gen_copy = ks_dNdy_neg_gen->Clone("ks_dNdy_neg_gen_copy");

   ks_dNdy_neg_rec->Draw();
   labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[0][i]);
   }
   ks_dNdy_neg_gen->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[0][i]);
   }
   ks_dNdy_neg_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[0][i]);
   }

   ks_dNdy_pos_eff->Divide(ks_dNdy_pos_rec,ks_dNdy_pos_gen);
   TH1F* ks_dNdy_pos_rec_copy = ks_dNdy_pos_rec->Clone("ks_dNdy_pos_rec_copy");
   TH1F* ks_dNdy_pos_gen_copy = ks_dNdy_pos_gen->Clone("ks_dNdy_pos_gen_copy");

   ks_dNdy_pos_rec->Draw();
   labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[0][i]);
   }
   ks_dNdy_pos_gen->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[0][i]);
   }
   ks_dNdy_pos_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[0][i]);
   }

//pt

   TH1F* ks_dNdpT_rec_scale = ks_dNdpT_rec->Clone("ks_dNdpT_rec_scale");
   ks_dNdpT_rec_scale->Divide(ks_ptscale);
   TH1F* ks_dNdpT_gen_scale = ks_dNdpT_gen->Clone("ks_dNdpT_gen_scale");
   ks_dNdpT_gen_scale->Divide(ks_ptscale);
   ks_dNdpT_eff->Divide(ks_dNdpT_rec,ks_dNdpT_gen);

   ks_dNdpT_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[1][i]);
   }
   ks_dNdpT_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[1][i]);
   }
   ks_dNdpT_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[1][i]);
   }


// lifetime

   TH1F* ks_dNdlife_rec_scale = ks_dNdlife_rec->Clone("ks_dNdlife_rec_scale");
   ks_dNdlife_rec_scale->Divide(ks_lifescale);
   TH1F* ks_dNdlife_gen_scale = ks_dNdlife_gen->Clone("ks_dNdlife_gen_scale");
   ks_dNdlife_gen_scale->Divide(ks_lifescale);
   ks_dNdlife_eff->Divide(ks_dNdlife_rec,ks_dNdlife_gen);

   ks_dNdlife_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[2][i]);
   }
   ks_dNdlife_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[2][i]);
   }
   ks_dNdlife_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[2][i]);
   }
  delete c1;


 } else {


   TCanvas* c1 = new TCanvas("c1", "Lambda results");
   
   lam_dNdy_eff->Divide(lam_dNdy_rec,lam_dNdy_gen);
   TH1F* lam_dNdy_rec_copy = lam_dNdy_rec->Clone("lam_dNdy_rec_copy");
   TH1F* lam_dNdy_gen_copy = lam_dNdy_gen->Clone("lam_dNdy_gen_copy");

   lam_dNdy_rec->Draw();
   labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[0][i]);
   }
   lam_dNdy_gen->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[0][i]);
   }
   lam_dNdy_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[0][i]);
   }

//pt

   TH1F* lam_dNdpT_rec_scale = lam_dNdpT_rec->Clone("lam_dNdpT_rec_scale");
   lam_dNdpT_rec_scale->Divide(lam_ptscale);
   TH1F* lam_dNdpT_gen_scale = lam_dNdpT_gen->Clone("lam_dNdpT_gen_scale");
   lam_dNdpT_gen_scale->Divide(lam_ptscale);
   lam_dNdpT_eff->Divide(lam_dNdpT_rec,lam_dNdpT_gen);

   lam_dNdpT_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[1][i]);
   }
   lam_dNdpT_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[1][i]);
   }
   lam_dNdpT_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[1][i]);
   }


// lifetime

   TH1F* lam_dNdlife_rec_scale = lam_dNdlife_rec->Clone("lam_dNdlife_rec_scale");
   lam_dNdlife_rec_scale->Divide(lam_lifescale);
   TH1F* lam_dNdlife_gen_scale = lam_dNdlife_gen->Clone("lam_dNdlife_gen_scale");
   lam_dNdlife_gen_scale->Divide(lam_lifescale);
   lam_dNdlife_eff->Divide(lam_dNdlife_rec,lam_dNdlife_gen);

   lam_dNdlife_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_rec_plots[2][i]);
   }
   lam_dNdlife_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_gen_plots[2][i]);
   }
   lam_dNdlife_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(output_eff_plots[2][i]);
   }


   //Particle

   lamReg_dNdy_eff->Divide(lamReg_dNdy_rec,lamReg_dNdy_gen);
   TH1F* lamReg_dNdy_rec_copy = lamReg_dNdy_rec->Clone("lamReg_dNdy_rec_copy");
   TH1F* lamReg_dNdy_gen_copy = lamReg_dNdy_gen->Clone("lamReg_dNdy_gen_copy");

   lamReg_dNdy_rec->Draw();
   labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_rec_plots[0][i]);
   }
   lamReg_dNdy_gen->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_gen_plots[0][i]);
   }
   lamReg_dNdy_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_eff_plots[0][i]);
   }

//pt

   TH1F* lamReg_dNdpT_rec_scale = lamReg_dNdpT_rec->Clone("lamReg_dNdpT_rec_scale");
   lamReg_dNdpT_rec_scale->Divide(lam_ptscale);
   TH1F* lamReg_dNdpT_gen_scale = lamReg_dNdpT_gen->Clone("lamReg_dNdpT_gen_scale");
   lamReg_dNdpT_gen_scale->Divide(lam_ptscale);
   lamReg_dNdpT_eff->Divide(lamReg_dNdpT_rec,lamReg_dNdpT_gen);

   lamReg_dNdpT_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_rec_plots[1][i]);
   }
   lamReg_dNdpT_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_gen_plots[1][i]);
   }
   lamReg_dNdpT_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_eff_plots[1][i]);
   }


// lifetime

   TH1F* lamReg_dNdlife_rec_scale = lamReg_dNdlife_rec->Clone("lamReg_dNdlife_rec_scale");
   lamReg_dNdlife_rec_scale->Divide(lam_lifescale);
   TH1F* lamReg_dNdlife_gen_scale = lamReg_dNdlife_gen->Clone("lamReg_dNdlife_gen_scale");
   lamReg_dNdlife_gen_scale->Divide(lam_lifescale);
   lamReg_dNdlife_eff->Divide(lamReg_dNdlife_rec,lamReg_dNdlife_gen);

   lamReg_dNdlife_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_rec_plots[2][i]);
   }
   lamReg_dNdlife_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_gen_plots[2][i]);
   }
   lamReg_dNdlife_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(part_output_eff_plots[2][i]);
   }


   //Anti Particle

   lamBar_dNdy_eff->Divide(lamBar_dNdy_rec,lamBar_dNdy_gen);
   TH1F* lamBar_dNdy_rec_copy = lamBar_dNdy_rec->Clone("lamBar_dNdy_rec_copy");
   TH1F* lamBar_dNdy_gen_copy = lamBar_dNdy_gen->Clone("lamBar_dNdy_gen_copy");

   lamBar_dNdy_rec->Draw();
   labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_rec_plots[0][i]);
   }
   lamBar_dNdy_gen->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_gen_plots[0][i]);
   }
   lamBar_dNdy_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_eff_plots[0][i]);
   }

//pt

   TH1F* lamBar_dNdpT_rec_scale = lamBar_dNdpT_rec->Clone("lamBar_dNdpT_rec_scale");
   lamBar_dNdpT_rec_scale->Divide(lam_ptscale);
   TH1F* lamBar_dNdpT_gen_scale = lamBar_dNdpT_gen->Clone("lamBar_dNdpT_gen_scale");
   lamBar_dNdpT_gen_scale->Divide(lam_ptscale);
   lamBar_dNdpT_eff->Divide(lamBar_dNdpT_rec,lamBar_dNdpT_gen);

   lamBar_dNdpT_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_rec_plots[1][i]);
   }
   lamBar_dNdpT_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_gen_plots[1][i]);
   }
   lamBar_dNdpT_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_eff_plots[1][i]);
   }


// lifetime

   TH1F* lamBar_dNdlife_rec_scale = lamBar_dNdlife_rec->Clone("lamBar_dNdlife_rec_scale");
   lamBar_dNdlife_rec_scale->Divide(lam_lifescale);
   TH1F* lamBar_dNdlife_gen_scale = lamBar_dNdlife_gen->Clone("lamBar_dNdlife_gen_scale");
   lamBar_dNdlife_gen_scale->Divide(lam_lifescale);
   lamBar_dNdlife_eff->Divide(lamBar_dNdlife_rec,lamBar_dNdlife_gen);

   lamBar_dNdlife_rec_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_rec_plots[2][i]);
   }
   lamBar_dNdlife_gen_scale->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_gen_plots[2][i]);
   }
   lamBar_dNdlife_eff->Draw();
   //labl->Draw("same");
   for ( Int_t i = 0; i < formats.size(); ++i ) {
     c1->SaveAs(anti_output_eff_plots[2][i]);
   }




  delete c1;
 }

  fout->Write();


}
