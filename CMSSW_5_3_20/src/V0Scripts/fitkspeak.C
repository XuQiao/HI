{
#include <tdrstyle.C>
#include <sstream>
//#include <iostream>
#include <iomanip.h>
#include <string>
#include <vector>

#include "ksfit.C"
#include "lamfit.C"
#include "kshortrec_masspeak.C"

  // Switches
  bool isMC = false;
  bool debug = false;
  bool kshort = true;

  ksfit* ksfitobj = new ksfit();
  lamfit* lamfitobj = new lamfit();

  // Style and other ROOT commands
  setTDRStyle();
  gStyle->SetOptStat(1);
  TGaxis::SetMaxDigits(3);
  gROOT->SetBatch(kTRUE);

  TString fin;
  TString outFileName;

  if( !isMC ) {
    //fin = "/nfs/data35/cms/v0ntuple/data_900GeV_may22fixes.root";
    fin = "/nfs/data37/cms/v0ntuple/data_7TeV_june3fixes_PARTIAL.root";
    //    fin = "../data_7TeV_allNtuples_400Jobs.root";
    if ( kshort ) {
      outFileName = "data_kshort_7000.root";
    } 
    else {
      outFileName = "data_lambda_7000.root";
    }
  }
  else {
    fin = "/nfs/data35/cms/v0ntuple/mc_7TeV_pythia8_may22fixes.root";
    //    fin = "/nfs/data35/cms/v0ntuple/mc_500kEvents.root";
    //    fin = "../mc_7TeV_pythia8.root";
	//fin = "../mc_500kEvents.root";
    if (kshort) {
      outFileName = "sim_kshort_pythia8_7000.root";
    } 
    else {
      outFileName = "sim_lambda_pythia8_7000.root";
    }
  }
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

  TPaveText* labl = new TPaveText(0.25, 0.88, 0.5, 0.95, "brNDC");
  labl->SetBorderSize(0);
  labl->SetFillColor(0);
  labl->SetFillStyle(0);
  labl->SetTextSize(0.035);
  TPaveText* labl2 = new TPaveText(0.25, 0.79, 0.55, 0.92, "brNDC");
  labl2->SetBorderSize(0);
  labl2->SetFillColor(0);
  labl2->SetFillStyle(0);
  labl2->SetTextSize(0.035);

  if( isMC ) {
    labltxt1_1 = labl->AddText("CMS Simulation pythia8");
    labltxt1_2 = labl->AddText("#sqrt{s} = 7 TeV");
    labltxt2_1 = labl2->AddText("CMS Simulation pythia8");
    labltxt2_2 = labl2->AddText("#sqrt{s} = 7 TeV");
  }
  else {
    labltxt1_1 = labl->AddText("CMS Preliminary");
    labltxt1_2 = labl->AddText("#sqrt{s} = 7 TeV");
    labltxt2_1 = labl2->AddText("CMS Preliminary");
    labltxt2_2 = labl2->AddText("#sqrt{s} = 7 TeV");
  }

  // Histogram quantities
  Double_t ksMassXmin = 0.340;
  Double_t ksMassXmax = 0.740;
  int ksMassNbins = 400;
  Double_t ksMassBinWidth = (ksMassXmax - ksMassXmin) / ksMassNbins;
  Double_t ksMassBinWidth_DG = ksMassBinWidth / sqrt(2*TMath::Pi());

  Double_t lamMassXmin = 1.09;
  Double_t lamMassXmax = 1.19;
  int lamMassNbins = 200;
  Double_t lamMassBinWidth = (lamMassXmax - lamMassXmin) / lamMassNbins;
  Double_t lamMassBinWidth_DG = lamMassBinWidth / sqrt(2*TMath::Pi());

  double ntrkarray[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,84,88,92,96,100,110,120,140,160,200};
  int ntrkNbins = sizeof(ntrkarray)/sizeof(double) -1;
  vector<double> ntrkBinMinima;
  for ( int ndx = 0; ndx <= ntrkNbins; ndx++ ) {
    ntrkBinMinima.push_back( ntrkarray[ndx] );
  }
  Double_t ntrk_default_binWidth = 0.1;

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
    TH1F* vee_cut1 = new TH1F("vee_cut1", "vee_cut1", ksMassNbins, ksMassXmin, ksMassXmax);
    TH1F* vee_cut2 = new TH1F("vee_cut2", "vee_cut2", ksMassNbins, ksMassXmin, ksMassXmax);
    TH1F* vee_cut3 = new TH1F("vee_cut3", "vee_cut3", ksMassNbins, ksMassXmin, ksMassXmax);
    TH1F* vee_cut4 = new TH1F("vee_cut4", "vee_cut4", ksMassNbins, ksMassXmin, ksMassXmax);
    TH1F* vee_cut5 = new TH1F("vee_cut5", "vee_cut5", ksMassNbins, ksMassXmin, ksMassXmax);
    vee_nocuts->SetXTitle("#pi^{+}#pi^{-} invariant mass  [GeV/c^{2}]");
    vee_nocuts->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut1->SetXTitle("#pi^{+}#pi^{-} invariant mass  [GeV/c^{2}]");
    vee_cut1->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut2->SetXTitle("#pi^{+}#pi^{-} invariant mass  [GeV/c^{2}]");
    vee_cut2->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut3->SetXTitle("#pi^{+}#pi^{-} invariant mass  [GeV/c^{2}]");
    vee_cut3->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut4->SetXTitle("#pi^{+}#pi^{-} invariant mass  [GeV/c^{2}]");
    vee_cut4->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut5->SetXTitle("#pi^{+}#pi^{-} invariant mass  [GeV/c^{2}]");
    vee_cut5->SetYTitle("Candidates / 0.1 GeV/c^{2}");
  } 
  else {
    TH1F* vee_nocuts = new TH1F("vee_nocuts", "vee_nocuts", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* vee_cut1 = new TH1F("vee_cut1", "vee_cut1", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* vee_cut2 = new TH1F("vee_cut2", "vee_cut2", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* vee_cut3 = new TH1F("vee_cut3", "vee_cut3", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* vee_cut4 = new TH1F("vee_cut4", "vee_cut4", lamMassNbins, lamMassXmin, lamMassXmax);
    TH1F* vee_cut5 = new TH1F("vee_cut5", "vee_cut5", lamMassNbins, lamMassXmin, lamMassXmax);
    vee_nocuts->SetXTitle("p^{+}#pi^{-} (+ c.c.) invariant mass  [GeV/c^{2}]");
    vee_nocuts->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut1->SetXTitle("p^{+}#pi^{-} (+ c.c.) invariant mass  [GeV/c^{2}]");
    vee_cut1->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut2->SetXTitle("p^{+}#pi^{-} (+ c.c.) invariant mass  [GeV/c^{2}]");
    vee_cut2->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut3->SetXTitle("p^{+}#pi^{-} (+ c.c.) invariant mass  [GeV/c^{2}]");
    vee_cut3->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut4->SetXTitle("p^{+}#pi^{-} (+ c.c.) invariant mass  [GeV/c^{2}]");
    vee_cut4->SetYTitle("Candidates / 0.1 GeV/c^{2}");
    vee_cut5->SetXTitle("p^{+}#pi^{-} (+ c.c.) invariant mass  [GeV/c^{2}]");
    vee_cut5->SetYTitle("Candidates / 0.1 GeV/c^{2}");
  }
  vee_nocuts->SetMinimum(0);
  vee_cut1->SetMinimum(0);
  vee_cut2->SetMinimum(0);
  vee_cut3->SetMinimum(0);
  vee_cut4->SetMinimum(0);
  vee_cut5->SetMinimum(0);


  // Find number of events.
  TH1F* nV0sH = new TH1F("nV0sH", "nV0sH", 100, 0., 100);
  TH1F* numvee = new TH1F("numvee", "numvee", 100, 0., 100);

  kshortrec_masspeak t;
  t.Loop();

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
    c1->SaveAs("ntrk_MC_pythia8_7000.png");
    c1->SaveAs("ntrk_MC_pythia8_7000.eps");
    c1->SaveAs("ntrk_MC_pythia8_7000.root");
    c1->SaveAs("ntrk_MC_pythia8_7000.C");
  }
  else {
    c1->SaveAs("ntrk_7000.png");
    c1->SaveAs("ntrk_7000.eps");
    c1->SaveAs("ntrk_7000.root");
    c1->SaveAs("ntrk_7000.C");
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
  // double height = 0; double heighte = 0;
  // double avgsigma = 0;
  // double yield = 0; double yielde = 0;
  // int fit0; int fitstat;
  // double frac; double sigma1; double sigma2; bool sigma1gsigma2;
  // double ston; double stone;

  gStyle->SetOptTitle(1);

  // Now we only need to create this object once, and we can refit with it
  //   as many times as we like.  Just pass the TH1F* to the Fit method.
  gStyle->SetOptTitle(0);
  fout->cd();
  if (kshort) {
    cout << "fitting full Kshort histogram\n";
    int gfitresult = ksfitobj->Fit(vee_nocuts);
    cout << "Total yield nocuts = " << ksfitobj->sig() << " ± " << ksfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "K0S mass plot no cut");
    int gplotresult = ksfitobj->Plot(vee_nocuts);
    if( isMC ) {
      c1->SaveAs("plots/ks_yldfit_nocut_MC_pythia8_7000.png");
      c1->SaveAs("plots/ks_yldfit_nocut_MC_pythia8_7000.eps");
      c1->SaveAs("plots/ks_yldfit_nocut_MC_pythia8_7000.root");
      c1->SaveAs("plots/ks_yldfit_nocut_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/ks_yldfit_nocut_7000.png");
      c1->SaveAs("plots/ks_yldfit_nocut_7000.eps");
      c1->SaveAs("plots/ks_yldfit_nocut_7000.root");
      c1->SaveAs("plots/ks_yldfit_nocut_7000.C");
    }
    cout << "fitting full Kshort histogram\n";
    int gfitresult = ksfitobj->Fit(vee_cut1);
    cout << "Total yield cut1 = " << ksfitobj->sig() << " ± " << ksfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "K0S mass plot cut1");
    int gplotresult = ksfitobj->Plot(vee_cut1);
    if( isMC ) {
      c1->SaveAs("plots/ks_yldfit_cut1_MC_pythia8_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut1_MC_pythia8_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut1_MC_pythia8_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut1_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/ks_yldfit_cut1_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut1_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut1_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut1_7000.C");
    }
    cout << "fitting full Kshort histogram\n";
    int gfitresult = ksfitobj->Fit(vee_cut2);
    cout << "Total yield cut2 = " << ksfitobj->sig() << " ± " << ksfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "K0S mass plot cut2");
    int gplotresult = ksfitobj->Plot(vee_cut2);
    if( isMC ) {
      c1->SaveAs("plots/ks_yldfit_cut2_MC_pythia8_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut2_MC_pythia8_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut2_MC_pythia8_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut2_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/ks_yldfit_cut2_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut2_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut2_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut2_7000.C");
    }
    cout << "fitting full Kshort histogram\n";
    int gfitresult = ksfitobj->Fit(vee_cut3);
    cout << "Total yield cut3 = " << ksfitobj->sig() << " ± " << ksfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "K0S mass plot cut3");
    int gplotresult = ksfitobj->Plot(vee_cut3);
    if( isMC ) {
      c1->SaveAs("plots/ks_yldfit_cut3_MC_pythia8_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut3_MC_pythia8_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut3_MC_pythia8_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut3_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/ks_yldfit_cut3_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut3_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut3_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut3_7000.C");
    }
    cout << "fitting full Kshort histogram\n";
    int gfitresult = ksfitobj->Fit(vee_cut4);
    cout << "Total yield cut4 = " << ksfitobj->sig() << " ± " << ksfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "K0S mass plot cut4");
    int gplotresult = ksfitobj->Plot(vee_cut4);
    if( isMC ) {
      c1->SaveAs("plots/ks_yldfit_cut4_MC_pythia8_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut4_MC_pythia8_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut4_MC_pythia8_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut4_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/ks_yldfit_cut4_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut4_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut4_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut4_7000.C");
    }
    cout << "fitting full Kshort histogram\n";
    int gfitresult = ksfitobj->Fit(vee_cut5);
    cout << "Total yield cut5 = " << ksfitobj->sig() << " ± " << ksfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "K0S mass plot cut5");
    int gplotresult = ksfitobj->Plot(vee_cut5);
    if( isMC ) {
      c1->SaveAs("plots/ks_yldfit_cut5_MC_pythia8_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut5_MC_pythia8_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut5_MC_pythia8_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut5_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/ks_yldfit_cut5_7000.png");
      c1->SaveAs("plots/ks_yldfit_cut5_7000.eps");
      c1->SaveAs("plots/ks_yldfit_cut5_7000.root");
      c1->SaveAs("plots/ks_yldfit_cut5_7000.C");
    }
  } else {
    cout << "fitting full Lambda histogram\n";
    int gfitresult = lamfitobj->Fit(vee_nocuts);
    cout << "Total yield nocuts = " << lamfitobj->sig() << " ± " << lamfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "Lambda mass plot no cut");
    int gplotresult = lamfitobj->Plot(vee_nocuts);
    if( isMC ) {
      c1->SaveAs("plots/lam_yldfit_nocut_MC_pythia8_7000.png");
      c1->SaveAs("plots/lam_yldfit_nocut_MC_pythia8_7000.eps");
      c1->SaveAs("plots/lam_yldfit_nocut_MC_pythia8_7000.root");
      c1->SaveAs("plots/lam_yldfit_nocut_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/lam_yldfit_nocut_7000.png");
      c1->SaveAs("plots/lam_yldfit_nocut_7000.eps");
      c1->SaveAs("plots/lam_yldfit_nocut_7000.root");
      c1->SaveAs("plots/lam_yldfit_nocut_7000.C");
    }
    cout << "fitting full Lambda histogram\n";
    int gfitresult = lamfitobj->Fit(vee_cut1);
    cout << "Total yield cut1 = " << lamfitobj->sig() << " ± " << lamfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "Lambda mass plot cut1");
    int gplotresult = lamfitobj->Plot(vee_cut1);
    if( isMC ) {
      c1->SaveAs("plots/lam_yldfit_cut1_MC_pythia8_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut1_MC_pythia8_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut1_MC_pythia8_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut1_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/lam_yldfit_cut1_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut1_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut1_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut1_7000.C");
    }
    cout << "fitting full Lambda histogram\n";
    int gfitresult = lamfitobj->Fit(vee_cut2);
    cout << "Total yield cut2 = " << lamfitobj->sig() << " ± " << lamfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "Lambda mass plot cut2");
    int gplotresult = lamfitobj->Plot(vee_cut2);
    if( isMC ) {
      c1->SaveAs("plots/lam_yldfit_cut2_MC_pythia8_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut2_MC_pythia8_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut2_MC_pythia8_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut2_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/lam_yldfit_cut2_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut2_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut2_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut2_7000.C");
    }
    cout << "fitting full Lambda histogram\n";
    int gfitresult = lamfitobj->Fit(vee_cut3);
    cout << "Total yield cut3 = " << lamfitobj->sig() << " ± " << lamfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "Lambda mass plot cut3");
    int gplotresult = lamfitobj->Plot(vee_cut3);
    if( isMC ) {
      c1->SaveAs("plots/lam_yldfit_cut3_MC_pythia8_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut3_MC_pythia8_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut3_MC_pythia8_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut3_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/lam_yldfit_cut3_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut3_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut3_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut3_7000.C");
    }
    cout << "fitting full Lambda histogram\n";
    int gfitresult = lamfitobj->Fit(vee_cut4);
    cout << "Total yield cut4 = " << lamfitobj->sig() << " ± " << lamfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "Lambda mass plot cut4");
    int gplotresult = lamfitobj->Plot(vee_cut4);
    if( isMC ) {
      c1->SaveAs("plots/lam_yldfit_cut4_MC_pythia8_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut4_MC_pythia8_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut4_MC_pythia8_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut4_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/lam_yldfit_cut4_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut4_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut4_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut4_7000.C");
    }
    cout << "fitting full Lambda histogram\n";
    int gfitresult = lamfitobj->Fit(vee_cut5);
    cout << "Total yield cut5 = " << lamfitobj->sig() << " ± " << lamfitobj->sig_err() << endl;
    TCanvas* c1 = new TCanvas("c1", "Lambda mass plot cut5");
    int gplotresult = lamfitobj->Plot(vee_cut5);
    if( isMC ) {
      c1->SaveAs("plots/lam_yldfit_cut5_MC_pythia8_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut5_MC_pythia8_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut5_MC_pythia8_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut5_MC_pythia8_7000.C");
    }
    else {
      c1->SaveAs("plots/lam_yldfit_cut5_7000.png");
      c1->SaveAs("plots/lam_yldfit_cut5_7000.eps");
      c1->SaveAs("plots/lam_yldfit_cut5_7000.root");
      c1->SaveAs("plots/lam_yldfit_cut5_7000.C");
    }


  }

  delete c1;

    //    for( int ndx = 0; ndx < lamlifeNbins; ndx++ ) {
    //      delete lamlifePtrs[ndx];
    //    }
    //    delete lamlifePtrs[];

  fout->Write();

}
