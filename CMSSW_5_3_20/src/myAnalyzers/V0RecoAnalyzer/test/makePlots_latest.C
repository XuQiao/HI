{
#include <tdrstyle.C>
#include <sstream>
//#include <iostream>
#include <string>
#include <vector>



gROOT->Reset();
setTDRStyle();
gStyle->SetLabelOffset(0.017, "X");
gStyle->SetPadTopMargin(0.06);
gStyle->SetPadBottomMargin(0.12);
gStyle->SetStatFormat("5.5g");
gStyle->SetStatY(0.94);
gStyle->SetFitFormat("3.3g");
gStyle->SetStatFontSize(0.042);
gStyle->SetStatFont(62);
gStyle->SetStatH(0.3);
gStyle->SetStatW(0.2);
//gStyle->SetStatFormat("NELU",AutoPrecision(3));

TFile* fin = TFile::Open("v0analysis_test.root", "READ");
TFile* fullmcfin = TFile::Open("v0analysis_MC.root", "READ");
TFile* fmassin = TFile::Open("ksMass_histo.root", "READ");
TFile* lifein = TFile::Open("ksLifetime_MC.root", "READ");
TFile* fout = TFile::Open("ksLifetime.root", "RECREATE");
fin->cd("analyzeKshort");

int choice;
cout << "For Monte Carlo, press 1.  For data, press 2.  Then press ENTER." 
     << endl;
cin >> choice;

cout << "Refit the K0S mass plot with a double gaussian?  1 for yes." << endl;
int refit = 0;
cin >> refit;


string mcordata;

if( choice == 1 ) {
  mcordata = "Monte Carlo";
}
else mcordata = "CMS Preliminary";

cout << mcordata << endl;

// Some constants
double piMass = 0.13957018;
double protonMass = 0.93827203;

TGaxis::SetMaxDigits(3);

labl = new TPaveText(0.5, 0.79, 0.8, 0.92, "brNDC");
labl->SetBorderSize(0);
labl->SetFillColor(0);
labl->SetTextSizePixels(20);
if( choice == 1 ) {
  labl->SetTextColor(13);
}
TText* labtxt = labl->AddText(mcordata.c_str());
TText* labtxt2;
if(mcordata == "CMS Preliminary") labtxt2 = labl->AddText("#sqrt{s} = 900 GeV");

lablMC = new TPaveText(0.5, 0.79, 0.8, 0.92, "brNDC");
lablMC->SetBorderSize(0);
lablMC->SetFillColor(0);
lablMC->SetTextSizePixels(20);
if( choice == 1 ) {
  lablMC->SetTextColor(13);
}
TText* labtxtMC = lablMC->AddText(//mcordata.c_str());
			      "CMS Simulation");
TText* labtxt2MC;
if(mcordata == "CMS Preliminary") labtxt2MC = lablMC->AddText("#sqrt{s} = 900 GeV");

labl2_0 = new TPaveText(0.25, 0.79, 0.55, 0.92, "brNDC");
labl2_0->SetBorderSize(0);
labl2_0->SetFillColor(0);
labl2_0->SetTextSizePixels(20);
if( choice == 1 ) {
  labl2_0->SetTextColor(13);
}
TText* labtxt2_0 = labl2_0->AddText(mcordata.c_str());
TText* labtxt2_2;
if(mcordata == "CMS Preliminary") labtxt2_2 = labl2_0->AddText("#sqrt{s} = 900 GeV");

ksLifetimeNbins = 25;
ksLifetimeXmin = 0.;
ksLifetimeXmax = 15.;
ksLifetimeBinWidth = (ksLifetimeXmax - ksLifetimeXmin) / ksLifetimeNbins; 

// This crap is a cross check using GenParticle decay distances.
/*
TFile *genin = TFile::Open("v0GenAnalysis.root", "READ");

genin->cd("ksGenAnalyzer");

ksGenBGCT = new TH1F("ksGenBGCT", "K^{0}_{S} gen #beta#gamma c#tau",
		     ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);

fullmcfin->cd("analyzeKshort");
ksMCBGCT = new TH1F("ksMCBGCT", "K^{0}_{S} MC #beta#gamma c#tau",
		    ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);

ksMCBGCT_BG = new TH1F("ksMCBGCT_BG", "Background",
		       ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);

genin->cd("ksGenAnalyzer");

gentuple->Project("ksGenBGCT", "fltDist", "");

ksGenBGCT->Draw();
c1->SaveAs("GenTest1.png");

fullmcfin->cd("analyzeKshort");

ntuple->Project("ksMCBGCT", "v0VtxPriDist", 
		"v0VtxSig > 15. && abs(v0CandMass - 0.498) < 0.02");
ntuple->Project("ksMCBGCT_BG", "v0VtxPriDist", 
		"v0VtxSig > 15 && abs(v0CandMass - (0.498 - 0.04)) < 0.02 && abs(v0CandMass - (0.498 + 0.04)) < 0.02");

ksMCBGCT->Draw();
c1->SaveAs("MCtest.png");

ksMCBGCT->Add( ksMCBGCT_BG, -0.5 );

ksMCBGCT->Divide( ksGenBGCT );

ksMCBGCT->Draw();
c1->SaveAs("GenDivided.png");

fin->cd("analyzeKshort");

//std::vector<TH1F*> ksDataCtauHistos;
//std::vector<TH1F*> ksDataCtauHistos_BG;
//std::vector<double> bgctBinContents;

TH1F* correctedKsCtau = new TH1F( "ksCtauGen", "K^{0}_{S} ctau",
				  ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax );
TH1F* correctedKsCtauBG = new TH1F( "ksCtauGenBG", "K^{0}_{S} ctau",
				    ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax );
for( int ndx = 1; ndx <= ksLifetimeNbins; ndx++ ) {
  double low = ndx*ksLifetimeBinWidth - ksLifetimeBinWidth;
  double high = ndx*ksLifetimeBinWidth;
  ostringstream signalss;
  ostringstream bgss;
  ostringstream titless;
  ostringstream titlessbg;
  signalss << "v0VtxSig > 15."
	   << " && v0VtxPriDist > " << low << " && v0VtxPriDist < " << high
	   << " && abs(v0CandMass - 0.498) < 0.02";
  bgss << "v0VtxSig > 15."
       << " && v0VtxPriDist > " << low << " && v0VtxPriDist < " << high
       << " && abs(v0CandMass - (0.498 - 0.04)) < 0.02"
       << " && abs(v0CandMass - (0.498 + 0.01)) < 0.02";
  titless << "ksCtau_" << ndx;
  titlessbg << "ksCtau_BG_" << ndx;


  TH1F* temp = new TH1F(titless.str().c_str(), "K^{0}_{S} c#tau",
			ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);
  TH1F* tempBG = new TH1F(titlessbg.str().c_str(), "K^{0}_{S} c#tau BG",
			  ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);

  ntuple->Project( titless.str().c_str(), "v03DLifetime",
		   signalss.str().c_str() );
  ntuple->Project( titlessbg.str().c_str(), "v03DLifetime",
		   bgss.str().c_str() );

  //ksDataCtauHistos.push_back( temp );
  //ksDataCtauHistos_BG.push_back( tempBG );
  //bgctBinContents.push_back( ksMCBGCT->GetBinContent(ndx) );
  cout << "ndx = " << ndx << endl;
  if( ndx == 1 ) {
    ntuple->Project( "ksCtauGen", "v03DLifetime",
		     signalss.str().c_str() );
    ntuple->Project( "ksCtauGenBG", "v03DLifetime",
		     bgss.str().c_str() );

    correctedKsCtau->Draw();
    correctedKsCtauBG->Draw();

    correctedKsCtau->Scale( ksMCBGCT->GetBinContent(ndx) );
    correctedKsCtauBG->Scale( ksMCBGCT->GetBinContent(ndx) );
  }
  else {
    correctedKsCtau->Add( temp, ksMCBGCT->GetBinContent(ndx) );
    correctedKsCtauBG->Add( tempBG, ksMCBGCT->GetBinContent(ndx) );
  }
  delete temp;
  delete tempBG;
}
correctedKsCtau->Draw();
c1->SaveAs("correctedKsCtau_Signal.png");

correctedKsCtauBG->Draw();
c1->SaveAs("correctedKsCtau_BG.png");

*/

gStyle->SetOptStat(0);
if(1) {
  ksLifetimeNbins = 20;
  ksLifetimeXmin = 0.;
  ksLifetimeXmax = 7.;
  ksCtauNbins = 10;

  //gStyle->SetErrorX((ksLifetimeXmax - ksLifetimeXmin) / ksCtauNbins);
  gStyle->SetErrorX(0.5);

  fin->cd("analyzeKshort");

  // Let's try doing mass plots in bins of ctau, fit, and use the yields to calculate
  //  the lifetime
  std::vector<double> fitYields;
  std::vector<double> fitYields_DG;
  std::vector<double> fitYieldErrors;
  std::vector<double> fitYieldErrors_DG;
  int ksCanvasX = ksCtauNbins / 2;
  int ksCanvasY = 2;
  TCanvas* ksCanvas = new TCanvas( "K0S2Dmasses", 
				   "K^{0}_{S} mass fits in c#tau bins",
				   2000, 800 );
  gStyle->SetOptTitle(1);
  ksCanvas->Divide( ksCanvasX, ksCanvasY );
  TPaveText* paves[10];
  TText* texts[10];
  for( int ndx = 0; ndx < ksCtauNbins; ndx++ ) {
    ksCanvas->cd( ndx + 1 );
    ostringstream ksHisto;
    ksHisto << "ksMTemp" << ndx;

    cout << ksHisto.str() << endl;
    ksCtauBinWidth = (ksLifetimeXmax - ksLifetimeXmin) / ksCtauNbins;
    cout << "BIN WIDTH: " << ksCtauBinWidth << endl << endl;
    //ksCtauBinWidth_DG = ksCtauBinWidth/sqrt(2*3.141592654);
    ksMNbins = 80;
    ksMXmin = 0.351;
    ksMXmax = 0.751;
    ksMBinWidth = (ksMXmax - ksMXmin) / ksMNbins;
    ksMBinWidth_DG = ksMBinWidth / sqrt(2*3.141592654);
    
    ostringstream ksselection;
    ksselection << "priVtxChi2 > 0. && v0VtxSig >= 15. && v0Lifetime > "
		<< ndx*ksCtauBinWidth << " && v0Lifetime < "
		<< (ndx+1)*ksCtauBinWidth
		<< " && abs(v0OtherCandMass - 1.1158) > 0. && v0DauNormChi2 < 5.";
    cout << ksselection.str() << endl;
    ostringstream labloss;
    labloss << ndx*ksCtauBinWidth << " < c#tau < " << (ndx + 1)*ksCtauBinWidth
	    << " cm";

    paves[ndx] = new TPaveText(0.6, 0.5, 0.8, 0.62, "brNDC");
    paves[ndx]->SetBorderSize(0);
    paves[ndx]->SetFillColor(0);
    paves[ndx]->SetTextColor(13);
    //TText* fitLabelText = fitLabel->AddText(labloss.str().c_str());
    texts[ndx] = paves[ndx]->AddText(labloss.str().c_str());
    ksMTemp = new TH1F(ksHisto.str().c_str(), //"K^{0}_{S} mass",
		       labloss.str().c_str(),
		       ksMNbins, ksMXmin, ksMXmax);
    ntuple->Project(ksHisto.str().c_str(), "v0CandMass", ksselection.str().c_str());
    ksMTemp->SetXTitle("m_{#pi^{+}#pi^{-}} (GeV/c^{2})");
    ksMTemp->SetYTitle("Entries/0.05 GeV/c^{2}");
    ksMTemp->SetNdivisions(506);

    ostringstream ksfitoss;
    ksfitoss << ksMBinWidth << "*gausn(0) + [3] + [4]*(x-[1])"
	     << " + [5]*(x-[1])^2";
    ostringstream ksfitoss_DG;
    ksfitoss_DG << ksMBinWidth_DG
		<< "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		<< ksMBinWidth_DG 
		<< "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		<< "[5] + [6]*(x-[1]) + [7]*(x-[1])^2";

    TF1* ksFit_SG = new TF1("ksFit_SG", ksfitoss.str().c_str(), ksMXmin, ksMXmax);
    ksFit_SG->SetParName( 0, "yield" );
    ksFit_SG->SetParName( 1, "mean" );
    ksFit_SG->SetParName( 2, "sigma" );
    ksFit_SG->SetParName( 3, "const" );
    ksFit_SG->SetParName( 4, "slope" );
    ksFit_SG->SetParName( 5, "quadconst" );
    
    ksFit_SG->SetParameter( 0, 4.4444*ksMTemp->GetMaximum() );
    ksFit_SG->SetParameter( 1, 0.49767 );
    ksFit_SG->SetParameter( 2, 0.005 );
    ksFit_SG->SetParLimits( 2, 0., 0.010 );
    ksFit_SG->SetParameter( 3, 0.1*ksMTemp->GetMaximum() );
    ksFit_SG->SetParameter( 4, 1. );
    ksFit_SG->SetParameter( 5, -1. );

    TF1* ksFit_DG = new TF1("ksFit_DG",
			    ksfitoss_DG.str().c_str(),
			    ksMXmin, ksMXmax);
    
    ksFit_DG->SetParName( 0, "yield" );
    ksFit_DG->SetParName( 1, "mean" );
    ksFit_DG->SetParName( 2, "sigma1" );
    ksFit_DG->SetParName( 3, "fraction" );
    ksFit_DG->SetParLimits( 3, 0., 1. );
    ksFit_DG->SetParName( 4, "sigma2" );
    ksFit_DG->SetParName( 5, "const" );
    ksFit_DG->SetParName( 6, "slope" );
    ksFit_DG->SetParName( 7, "quadconst" );

    ksFit_DG->SetParameter( 0, 4.4444*ksMTemp->GetMaximum() );
    ksFit_DG->SetParameter( 1, 0.49767 );
    ksFit_DG->SetParameter( 2, 0.0045 );
    ksFit_DG->SetParameter( 3, 0.75 );
    ksFit_DG->SetParameter( 4, 0.011 );
    ksFit_DG->SetParameter( 5, 0.1*ksMTemp->GetMaximum() );
    ksFit_DG->SetParameter( 6, 1. );
    ksFit_DG->SetParameter( 7, -1. );

    int fitstat = ksMTemp->Fit("ksFit_SG", "RLE");

    ostringstream screwy;
    screwy << "ksMTemp_" << ndx << ".png";
    ksMTemp->Draw();
    paves[ndx]->Draw("same");
    //gPad->Update();
    //string fname = string("ksMTemp_") + string(ndx) + string(".png");
    //c1->SaveAs(screwy.str().c_str());
    if( !fitstat ) {
      fitYields->push_back( ksFit_SG->GetParameter(0) );
      fitYieldErrors->push_back( ksFit_SG->GetParError(0) );
    }
    else {
      fitYields->push_back(0.);
      fitYieldErrors->push_back( 0. );
    }

    fitstat = ksMTemp->Fit("ksFit_DG", "RLE");
    if( !fitstat ) {
      fitYields_DG->push_back( ksFit_DG->GetParameter(0) );
      fitYieldErrors_DG->push_back( ksFit_DG->GetParError(0) );
    }
    else {
      fitYields_DG->push_back( 0. );
      fitYieldErrors_DG->push_back( 0. );
    }
    delete ksFit_DG;
    delete ksFit_SG;
    //delete fitLabel;
    //delete ksMTemp;
  }
  ksCanvas->Print("KsMassFitsInCtauBins2D.png");
  ksCanvas->Print("KsMassFitsInCtauBins2D.root");
  delete ksCanvas;
  gStyle->SetOptTitle(0);

  // Now divide by real lifetime, etc.
  lifein->cd();
  TF1* ksPDGCtau = new TF1(*ksCtauReal);
  TH1F* ksCtauFromYields_MC = new TH1F(*ksCtauFromYields);
  fin->cd("analyzeKshort");

  ksCtauFromYields = new TH1F("ksCtauFromYields", "K^{0}_{S} c#tau",
			      ksCtauNbins, ksLifetimeXmin, ksLifetimeXmax);
  ksCtauFromYields_DG = new TH1F("ksCtauFromYields_DG", "K^{0}_{S} c#tau",
				 ksCtauNbins, ksLifetimeXmin, ksLifetimeXmax);

  ksCtauFromYields->Sumw2();
  ksCtauFromYields_DG->Sumw2();

  for( int ndx = 0; ndx < fitYields.size(); ndx++ ) {
    ksCtauFromYields->SetBinContent( ndx + 1, fitYields[ndx] );
    ksCtauFromYields->SetBinError( ndx + 1, fitYieldErrors[ndx] );
  }

  ksCtauFromYields->GetXaxis()->SetTitle("K^{0}_{S} c#tau (cm)");
  ksCtauFromYields->GetYaxis()->SetTitle("K^{0}_{S} mass fit yield");
  ksCtauFromYields->GetXaxis()->SetTitleOffset(0.86);
  ksCtauFromYields->Draw();
  labl->Draw("same");
  c1->cd(1);
  gPad->SetLogy(1);
  c1->SaveAs("ksCtauFromYields_noScaling.png");
  c1->SaveAs("ksCtauFromYields_noScaling.root");
  c1->SaveAs("ksCtauFromYields_noScaling.pdf");

  for( int ndx = 0; ndx < fitYields_DG.size(); ndx++ ) {
    ksCtauFromYields_DG->SetBinContent( ndx + 1, fitYields_DG[ndx] );
    ksCtauFromYields_DG->SetBinError( ndx + 1, fitYieldErrors_DG[ndx] );
  }

  ksCtauFromYields_DG->Draw();
  c1->SaveAs("ksCtauFromYields_DG_noScaling.png");

  //Double_t c_const = 2.99792e10;// in cm/s
  Double_t c_const = 2.99792e-2;//in cm/ps, should put tau in ps
  ostringstream ksctauoss_;
  ksctauoss_ << "exp([0])*exp(-x/(" << c_const << "*[1]))";

  ksCtauFromYieldsFit = new TF1("ksCtauFromYieldsFit", ksctauoss_.str().c_str(),
				ksLifetimeXmin, ksLifetimeXmax);

  Double_t exp_p0_y = log( ksCtauFromYields->GetMaximum() );
  Double_t rise_y = log( ksCtauFromYields->GetMaximum() 
			 - ksCtauFromYields->GetBinContent(ksCtauNbins));
  Double_t exp_p1_y = (rise_y / ksLifetimeXmax);
  exp_p1_y = 1/c_const/exp_p1_y;

  ksCtauFromYieldsFit->SetParameter(0, exp_p0_y);
  ksCtauFromYieldsFit->SetParameter(1, exp_p1_y);
  ksCtauFromYieldsFit->SetParName(1, "#tau (ps)");

  ksCtauFromYieldsFit->Draw();
  c1->SaveAs("fitFunc.png");

  Double_t realCtau_y = 2.6786;
  Double_t realCtau_p1_y = realCtau_y / c_const;
  ksCtauReal_y = new TF1("ksCtauReal_y", ksctauoss_.str().c_str(),
			 ksLifetimeXmin, ksLifetimeXmax);
  Double_t realCtau_p0_y = //ksCtauFromYieldsFit->GetParameter(0);
    log(ksCtauFromYields_MC->GetMaximum());
  ksCtauReal_y->SetParameter(0, realCtau_p0_y);
  ksCtauReal_y->SetParameter(1, realCtau_p1_y);

  ksCtauFromYields_MC->Divide( ksCtauReal_y );
  gPad->SetLogy(0);
  gPad->Update();
  ksCtauFromYields_MC->GetYaxis()->SetTitle("Correction factor");
  ksCtauFromYields_MC->GetXaxis()->SetTitle("K^{0}_{S} c#tau (cm)");
  ksCtauFromYields_MC->SetMinimum(0.);
  ksCtauFromYields_MC->GetXaxis()->SetTitleOffset(0.86);
  ksCtauFromYields_MC->Draw();
  lablMC->Draw("same");
  c1->cd(1);
  //gPad->SetLogy(0);
  //gPad->Update();
  c1->SaveAs("KsYieldsRelEff.png");
  c1->SaveAs("KsYieldsRelEff.root");
  c1->SaveAs("KsYieldsRelEff.pdf");
  gPad->SetLogy(1);
  gPad->Update();

  ksCtauFromYields->Divide(ksCtauFromYields_MC);
  ksCtauFromYields->Fit("ksCtauFromYieldsFit", "REI");
  legen_ = new TLegend(0.2, 0.2, 0.45, 0.39);
  legen_->AddEntry(ksCtauFromYields, "Corrected data", "lp");
  legen_->AddEntry(ksCtauFromYieldsFit, "Exponential fit", "l");
  legen_->SetFillColor(0);
  legen_->SetBorderSize(0);
  legen_->SetTextSize(0.035);
  ksCtauFromYields->GetYaxis()->SetTitle("Corrected K^{0}_{S} mass fit yield");
  ksCtauFromYields->GetXaxis()->SetTitle("K^{0}_{S} c#tau (cm)");
  ksCtauFromYields->GetXaxis()->SetTitleOffset(0.86);
  //ksCtauFromYields->Draw();
  //c1->cd(1);
  //gPad->SetLogy(1);
  //gPad->Update();
  //c1->SaveAs("YieldsScaledCtau_noFit.png");

  ksCtauFromYields->Draw();
  legen_->Draw("same");
  labl2_0->Draw("same");
  //delete legen_;

  c1->SaveAs("ksCtauFromYields_withFit.png");
  c1->SaveAs("ksCtauFromYields_withFit.root");
  c1->SaveAs("ksCtauFromYields_withFit.pdf");
  delete legen_;

  delete ksCtauFromYields;

  gPad->SetLogy(0);

  lamLifetimeNbins = 20;
  lamLifetimeXmin = 0.;
  lamLifetimeXmax = 16.;
  lamCtauNbins = 8;

  //gStyle->SetErrorX((lamLifetimeXmax - lamLifetimeXmin) / lamCtauNbins);

  fin->cd("analyzeLambda");

  std::vector<double> lamFitYields;
  std::vector<double> lamFitYieldErrors;
  int lamCanvasX = lamCtauNbins / 2;
  int lamCanvasY = 2;
  TCanvas* lamCanvas = new TCanvas( "Lam2Dmasses",
				    "#Lambda^{0} mass fits in c#tau bins",
				    1600, 800 );
  gStyle->SetOptTitle(1);
  lamCanvas->Divide( lamCanvasX, lamCanvasY );
  /*lamMFull = new TH1F("lamMFull", "#Lambda",
    50, 1.08, 1.18);

    ntuple->Project("lamMFull", "v0CandMass", "priVtxChi2 > 0. && v0VtxSig > 15. && v0Lifetime > 0. && v0Lifetime < 20. && abs(v0OtherCandMass - 0.498) > 0.02");
    lamMFull->Draw();
    c1->SaveAs("LamMassTest_0_ctau_20.png");*/
  gPad->SetLogy(0);
  for( int ndx = 0; ndx < lamCtauNbins; ndx++ ) {
    lamCanvas->cd(ndx + 1);
    ostringstream lamHisto;
    lamHisto << "lamMTemp" << ndx;
    lamCtauBinWidth = (lamLifetimeXmax - lamLifetimeXmin) / lamCtauNbins;
    lamMNbins = 50;
    lamMXmin = 1.08;
    lamMXmax = 1.18;
    lamMBinWidth = (lamMXmax - lamMXmin) / lamMNbins;

    ostringstream lamselection;
    lamselection << "priVtxChi2 > 0. && v0VtxSig > 15. && v0Lifetime > "
		 << ndx*lamCtauBinWidth << " && v0Lifetime < "
		 << (ndx+1)*lamCtauBinWidth
		 << " && abs(v0OtherCandMass - 0.498) > 0. && v0DauNormChi2 < 5.";
    ostringstream labloss;
    labloss << ndx*lamCtauBinWidth << " < c#tau < " << (ndx+1)*lamCtauBinWidth
	    << " cm";
    cout << lamselection.str() << endl;
    lamMTemp = //new TH1F("lamMTemp", "#Lambda^{0} mass",
      new TH1F(lamHisto.str().c_str(), labloss.str().c_str(),
	       lamMNbins, lamMXmin, lamMXmax);

    ntuple->Project(lamHisto.str().c_str(), "v0CandMass", lamselection.str().c_str());
    lamMTemp->SetXTitle("m_{p^{+}#pi^{-}} ( + c.c.) (GeV/c^{2})");
    lamMTemp->SetYTitle("Entries/0.02 GeV/c^{2}");
    lamMTemp->SetNdivisions(506);

    ostringstream lamfitoss;
    lamfitoss << lamMBinWidth << "*gausn(0) + [3]*(x - "
	      << piMass + protonMass << ")^(1/2) + [4]*(x - "
	      << piMass + protonMass << ")^(3/2)";
    TF1* lamFit_SG = new TF1("lamFit_SG", lamfitoss.str().c_str(), lamMXmin, lamMXmax);
    lamFit_SG->SetParName( 0, "yield" );
    lamFit_SG->SetParName( 1, "GausMean" );
    lamFit_SG->SetParName( 2, "GausSigma" );
    lamFit_SG->SetParName( 3, "c1" );
    lamFit_SG->SetParName( 4, "c2" );

    lamFit_SG->SetParameter( 0, 2.6*lamMTemp->GetMaximum() );
    lamFit_SG->SetParameter( 1, 1.1156 );
    lamFit_SG->SetParameter( 2, 0.002 );
    lamFit_SG->SetParameter( 3, 0.75*lamMTemp->GetMaximum() );
    lamFit_SG->SetParameter( 4, 1.85*lamMTemp->GetMaximum() );

    int fitstat = lamMTemp->Fit("lamFit_SG", "RLE");

    ostringstream hmmm;
    hmmm << "lamMTemp_" << ndx << ".png";
    lamMTemp->Draw();
    //gPad->SetLogy(0);
    //c1->SaveAs(hmmm.str().c_str());
    cout << "Fit status: " << fitstat << endl;
    if( !fitstat ) {
      lamFitYields->push_back( lamFit_SG->GetParameter(0) );
      lamFitYieldErrors->push_back( lamFit_SG->GetParError(0) );
    }
    else {
      lamFitYields->push_back( 0. );
      lamFitYieldErrors->push_back( 0. );
    }
    delete lamFit_SG;
    //delete lamMTemp;
  }

  lamCanvas->Print("LamMassFitsInCtauBins2D.png");
  lamCanvas->Print("LamMassFitsInCtauBins2D.root");
  delete lamCanvas;
  gStyle->SetOptTitle(0);

  lifein->cd();
  TF1* lamPDGCtau = new TF1(*lamCtauReal);
  TH1F* lamCtauFromYields_MC = new TH1F(*lamCtauFromYields);
  fin->cd("analyzeLambda");

  lamCtauFromYields = new TH1F("lamCtauFromYields", "#Lambda c#tau",
			       lamCtauNbins, lamLifetimeXmin, lamLifetimeXmax);

  lamCtauFromYields->Sumw2();

  for( int ndx = 0; ndx < fitYields.size(); ndx++ ) {
    lamCtauFromYields->SetBinContent( ndx + 1, lamFitYields[ndx] );
    lamCtauFromYields->SetBinError( ndx + 1, lamFitYieldErrors[ndx] );
  }

  lamCtauFromYields->Draw();
  labl->Draw("same");
  lamCtauFromYields->SetXTitle("#Lambda^{0} c#tau (cm)");
  lamCtauFromYields->SetYTitle("#Lambda^{0} mass fit yield");
  lamCtauFromYields->GetXaxis()->SetRangeUser(lamLifetimeXmin, lamLifetimeXmax + 0.1);
  lamCtauFromYields->GetXaxis()->SetTitleOffset(1);
  gPad->SetLogy(1);
  gPad->Update();
  c1->SaveAs("lamCtauFromYields_noScaling.png");
  c1->SaveAs("lamCtauFromYields_noScaling.root");
  c1->SaveAs("lamCtauFromYields_noScaling.pdf");
  gPad->SetLogy(0);

  ostringstream lamctauoss_;
  lamctauoss_ << "exp([0])*exp(-x/(" << c_const << "*[1]))";

  lamCtauFromYieldsFit = new TF1("lamCtauFromYieldsFit", lamctauoss_.str().c_str(),
				 lamLifetimeXmin, lamLifetimeXmax);

  Double_t exp_p0_ly = log( lamCtauFromYields->GetMaximum() );
  Double_t rise_ly = log( lamCtauFromYields->GetMaximum()
			  - lamCtauFromYields->GetBinContent(lamCtauNbins) );
  Double_t exp_p1_ly = (rise_ly / lamLifetimeXmax);
  exp_p1_ly = 1/c_const/exp_p1_ly;

  lamCtauFromYieldsFit->SetParameter(0, exp_p0_ly);
  lamCtauFromYieldsFit->SetParameter(1, exp_p1_ly);
  lamCtauFromYieldsFit->SetParName(1, "#tau (ps)");

  lamCtauFromYieldsFit->Draw();
  c1->SaveAs("lamFitFunc.png");

  Double_t realCtau_ly = 7.89;
  Double_t realCtau_p1_ly = realCtau_ly/c_const;
  lamCtauReal_ly = new TF1("lamCtauReal_ly", lamctauoss_.str().c_str(),
			   lamLifetimeXmin, lamLifetimeXmax);
  Double_t realCtau_p0_ly = //lamCtauFromYieldsFit->GetParameter(0);
    log(lamCtauFromYields_MC->GetMaximum());
  lamCtauReal_ly->SetParameter(0, realCtau_p0_ly);
  lamCtauReal_ly->SetParameter(1, realCtau_p1_ly);

  lamCtauFromYields_MC->Divide( lamCtauReal_ly );
  gPad->SetLogy(0);
  gPad->Update();
  lamCtauFromYields_MC->SetXTitle("#Lambda^{0} c#tau (cm)");
  lamCtauFromYields_MC->GetYaxis()->SetTitle("Correction factor");
  lamCtauFromYields_MC->SetMinimum(0.);
  lamCtauFromYields_MC->GetXaxis()->SetLimits(lamLifetimeXmin, lamLifetimeXmax + 0.1);
  lamCtauFromYields_MC->GetXaxis()->SetTitleOffset(1);
  lamCtauFromYields_MC->Draw();
  labl->Draw("same");
  //c1->cd(1);
  //gPad->SetLogy(0);
  //gPad->Update();
  c1->SaveAs("LamYieldsRelEff.png");
  c1->SaveAs("LamYieldsRelEff.root");
  c1->SaveAs("LamYieldsRelEff.pdf");
  gPad->SetLogy(1);
  gPad->Update();
  gStyle->SetFitFormat("3.0f");

  lamCtauFromYields->Divide(lamCtauFromYields_MC);
  lamCtauFromYields->Fit("lamCtauFromYieldsFit", "REI");
  legen_ = new TLegend(0.2, 0.2, 0.45, 0.39);
  legen_->AddEntry(lamCtauFromYields, "Corrected data", "lp");
  legen_->AddEntry(lamCtauFromYieldsFit, "Exponential fit", "l");
  legen_->SetFillColor(0);
  legen_->SetBorderSize(0);
  legen_->SetTextSize(0.035);
  lamCtauFromYields->GetYaxis()->SetTitle("Corrected #Lambda^{0} mass fit yield");
  lamCtauFromYields->GetXaxis()->SetTitle("c#tau (cm)");
  lamCtauFromYields->GetXaxis()->SetTitleOffset(1);

  lamCtauFromYields->Draw();
  legen_->Draw("same");
  labl2_0->Draw("same");

  //delete legen_;

  c1->SaveAs("lamCtauFromYields_withFit.png");
  c1->SaveAs("lamCtauFromYields_withFit.root");
  c1->SaveAs("lamCtauFromYields_withFit.pdf");

  delete legen_;
  gStyle->SetFitFormat("3.3g");

}

if(1) {
  ksLifetimeNbins = 20;
  ksLifetimeXmin = 0.;
  ksLifetimeXmax = 7.;
  ksCtauNbins = 10;

  //gStyle->SetErrorX((ksLifetimeXmax - ksLifetimeXmin) / ksCtauNbins);

  fin->cd("analyzeKshort");
  gPad->SetLogy(0);

  // Let's try doing mass plots in bins of ctau, fit, and use the yields to calculate
  //  the lifetime
  std::vector<double> fitYields;
  std::vector<double> fitYields_DG;
  std::vector<double> fitYieldErrors;
  std::vector<double> fitYieldErrors_DG;
  int ksCanvasX = ksCtauNbins / 2;
  int ksCanvasY = 2;
  TCanvas* ksCanvas = new TCanvas( "K0S3Dmasses",
				   "K^{0}_{S} mass fits in c#tau bins",
				   2000, 800 );
  gStyle->SetOptTitle(1);
  ksCanvas->Divide( ksCanvasX, ksCanvasY );
  for( int ndx = 0; ndx < ksCtauNbins; ndx++ ) {
    ksCanvas->cd( ndx + 1 );
    ostringstream ksHisto;
    ksHisto << "ksMTemp" << ndx;
    ksCtauBinWidth = (ksLifetimeXmax - ksLifetimeXmin) / ksCtauNbins;
    cout << "BIN WIDTH: " << ksCtauBinWidth << endl << endl;
    //ksCtauBinWidth_DG = ksCtauBinWidth/sqrt(2*3.141592654);
    ksMNbins = 80;
    ksMXmin = 0.351;
    ksMXmax = 0.751;
    ksMBinWidth = (ksMXmax - ksMXmin) / ksMNbins;
    ksMBinWidth_DG = ksMBinWidth / sqrt(2*3.141592654);
    
    ostringstream ksselection;
    ksselection << "priVtxChi2 > 0. && v0VtxSig >= 15. && v03DLifetime > "
		<< ndx*ksCtauBinWidth << " && v03DLifetime < "
		<< (ndx+1)*ksCtauBinWidth
		<< " && abs(v0OtherCandMass - 1.1158) > 0. && v0DauNormChi2 < 5.";
    cout << ksselection.str() << endl;
    ostringstream labloss;
    labloss << ndx*ksCtauBinWidth << " < c#tau < " << (ndx + 1)*ksCtauBinWidth
	    << " cm";

    ksMTemp = new TH1F(ksHisto.str().c_str(),// "K^{0}_{S} mass",
		       labloss.str().c_str(),
		       ksMNbins, ksMXmin, ksMXmax);
    ntuple->Project(ksHisto.str().c_str(), "v0CandMass", ksselection.str().c_str());
    ksMTemp->SetXTitle("m_{#pi^{+}#pi^{-}} (GeV/c^{2})");
    ksMTemp->SetYTitle("Entries/0.05 GeV/c^{2}");
    ksMTemp->SetNdivisions(506);

    ostringstream ksfitoss;
    ksfitoss << ksMBinWidth << "*gausn(0) + [3] + [4]*(x-[1])"
	     << " + [5]*(x-[1])^2";
    ostringstream ksfitoss_DG;
    ksfitoss_DG << ksMBinWidth_DG
		<< "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		<< ksMBinWidth_DG 
		<< "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		<< "[5] + [6]*(x-[1]) + [7]*(x-[1])^2";

    TF1* ksFit_SG = new TF1("ksFit_SG", ksfitoss.str().c_str(), ksMXmin, ksMXmax);
    ksFit_SG->SetParName( 0, "yield" );
    ksFit_SG->SetParName( 1, "mean" );
    ksFit_SG->SetParName( 2, "sigma" );
    ksFit_SG->SetParName( 3, "const" );
    ksFit_SG->SetParName( 4, "slope" );
    ksFit_SG->SetParName( 5, "quadconst" );
    
    ksFit_SG->SetParameter( 0, 4.4444*ksMTemp->GetMaximum() );
    ksFit_SG->SetParameter( 1, 0.49767 );
    ksFit_SG->SetParameter( 2, 0.005 );
    ksFit_SG->SetParLimits( 2, 0., 0.010 );
    ksFit_SG->SetParameter( 3, 0.1*ksMTemp->GetMaximum() );
    ksFit_SG->SetParameter( 4, 1. );
    ksFit_SG->SetParameter( 5, -1. );

    TF1* ksFit_DG = new TF1("ksFit_DG",
			    ksfitoss_DG.str().c_str(),
			    ksMXmin, ksMXmax);
    
    ksFit_DG->SetParName( 0, "yield" );
    ksFit_DG->SetParName( 1, "mean" );
    ksFit_DG->SetParName( 2, "sigma1" );
    ksFit_DG->SetParName( 3, "fraction" );
    ksFit_DG->SetParLimits( 3, 0., 1. );
    ksFit_DG->SetParName( 4, "sigma2" );
    ksFit_DG->SetParName( 5, "const" );
    ksFit_DG->SetParName( 6, "slope" );
    ksFit_DG->SetParName( 7, "quadconst" );

    ksFit_DG->SetParameter( 0, 4.4444*ksMTemp->GetMaximum() );
    ksFit_DG->SetParameter( 1, 0.49767 );
    ksFit_DG->SetParameter( 2, 0.0045 );
    ksFit_DG->SetParameter( 3, 0.75 );
    ksFit_DG->SetParameter( 4, 0.011 );
    ksFit_DG->SetParameter( 5, 0.1*ksMTemp->GetMaximum() );
    ksFit_DG->SetParameter( 6, 1. );
    ksFit_DG->SetParameter( 7, -1. );

    int fitstat = ksMTemp->Fit("ksFit_SG", "RLE");

    ostringstream screwy;
    screwy << "ksMTemp3D_" << ndx << ".png";
    ksMTemp->Draw();
    //string fname = string("ksMTemp_") + string(ndx) + string(".png");
    //c1->SaveAs(screwy.str().c_str());
    if( !fitstat ) {
      fitYields->push_back( ksFit_SG->GetParameter(0) );
      fitYieldErrors->push_back( ksFit_SG->GetParError(0) );
    }
    else {
      fitYields->push_back(0.);
      fitYieldErrors->push_back( 0. );
    }

    fitstat = ksMTemp->Fit("ksFit_DG", "RLE");
    if( !fitstat ) {
      fitYields_DG->push_back( ksFit_DG->GetParameter(0) );
      fitYieldErrors_DG->push_back( ksFit_DG->GetParError(0) );
    }
    else {
      fitYields_DG->push_back( 0. );
      fitYieldErrors_DG->push_back( 0. );
    }
    delete ksFit_DG;
    delete ksFit_SG;
    //delete ksMTemp;
  }
  ksCanvas->Print("KsMassFitsInCtauBins3D.png");
  ksCanvas->Print("KsMassFitsInCtauBins3D.root");
  delete ksCanvas;
  gStyle->SetOptTitle(0);

  // Now divide by real lifetime, etc.
  lifein->cd();
  TF1* ksPDGCtau3D = new TF1(*ksCtauReal);
  TH1F* ksCtauFromYields_MC3D = new TH1F(*ksCtauFromYields3D);
  ksCtauFromYields_MC3D->Draw();
  c1->cd(1);
  gPad->SetLogy(1);
  c1->SaveAs("ksCtauFromYields_MC_noScaling3D.png");
  fin->cd("analyzeKshort");

  ksCtauFromYields3D_ = new TH1F("ksCtauFromYields3D_", "K^{0}_{S} c#tau",
			      ksCtauNbins, ksLifetimeXmin, ksLifetimeXmax);
  ksCtauFromYields_DG3D = new TH1F("ksCtauFromYields_DG3D", "K^{0}_{S} c#tau",
				 ksCtauNbins, ksLifetimeXmin, ksLifetimeXmax);

  ksCtauFromYields3D_->Sumw2();
  ksCtauFromYields_DG3D->Sumw2();

  for( int ndx = 0; ndx < fitYields.size(); ndx++ ) {
    cout << "Filling histo from vector..." << endl;
    cout << fitYields[ndx] << endl;
    ksCtauFromYields3D_->SetBinContent( ndx + 1, fitYields[ndx] );
    ksCtauFromYields3D_->SetBinError( ndx + 1, fitYieldErrors[ndx] );
  }

  ksCtauFromYields3D_->GetXaxis()->SetTitle("K^{0}_{S} c#tau (cm)");
  ksCtauFromYields3D_->GetYaxis()->SetTitle("K^{0}_{S} mass fit yield");
  ksCtauFromYields3D_->GetXaxis()->SetTitleOffset(0.86);
  ksCtauFromYields3D_->Draw();
  labl->Draw("same");
  //c1->cd(1);
  //gPad->SetLogy(1);
  c1->SaveAs("ksCtauFromYields_noScaling3D.png");
  c1->SaveAs("ksCtauFromYields_noScaling3D.root");
  c1->SaveAs("ksCtauFromYields_noScaling3D.pdf");

  for( int ndx = 0; ndx < fitYields_DG.size(); ndx++ ) {
    ksCtauFromYields_DG3D->SetBinContent( ndx + 1, fitYields_DG[ndx] );
    ksCtauFromYields_DG3D->SetBinError( ndx + 1, fitYieldErrors_DG[ndx] );
  }

  ksCtauFromYields_DG3D->Draw();
  c1->SaveAs("ksCtauFromYields_DG_noScaling3D.png");

  //Double_t c_const = 2.99792e10;//in cm/s
  Double_t c_const = 2.99792e-2;//in cm/ps, should put tau in ps
  ostringstream ksctauoss_;
  ksctauoss_ << "exp([0])*exp(-x/(" << c_const << "*[1]))";

  ksCtauFromYieldsFit3D_ = new TF1("ksCtauFromYieldsFit3D_", ksctauoss_.str().c_str(),
				ksLifetimeXmin, ksLifetimeXmax);

  Double_t exp_p0_y = log( ksCtauFromYields3D_->GetMaximum() );
  Double_t rise_y = log( ksCtauFromYields3D_->GetMaximum() 
			 - ksCtauFromYields3D_->GetBinContent(ksCtauNbins));
  Double_t exp_p1_y = (rise_y / ksLifetimeXmax);
  exp_p1_y = 1/c_const/exp_p1_y;

  ksCtauFromYieldsFit3D_->SetParameter(0, exp_p0_y);
  ksCtauFromYieldsFit3D_->SetParameter(1, exp_p1_y);
  ksCtauFromYieldsFit3D_->SetParName(1, "#tau (ps)");

  ksCtauFromYieldsFit3D_->Draw();
  c1->SaveAs("fitFunc3D.png");

  Double_t realCtau_y = 2.6786;
  Double_t realCtau_p1_y = realCtau_y / c_const;
  cout << realCtau_p1_y << endl;
  ksCtauReal_y3D = new TF1("ksCtauReal_y3D", ksctauoss_.str().c_str(),
			 ksLifetimeXmin, ksLifetimeXmax);
  Double_t realCtau_p0_y = //ksCtauFromYieldsFit3D->GetParameter(0);
    log(ksCtauFromYields_MC3D->GetMaximum());
  ksCtauReal_y3D->SetParameter(0, realCtau_p0_y);
  ksCtauReal_y3D->SetParameter(1, realCtau_p1_y);

  //ksCtauFromYields_MC3D->Draw();
  //c1->SaveAs("before.png");
  //ksCtauReal_y3D->Draw();
  //c1->SaveAs("function_before.png");
  ksCtauFromYields_MC3D->Divide( ksCtauReal_y3D );
  //ksCtauFromYields_MC3D->Draw();
  //c1->SaveAs("after.png");
  c1->cd(1);
  gPad->SetLogy(0);
  gPad->Update();
  ksCtauFromYields_MC3D->GetYaxis()->SetTitle("Correction factor");
  ksCtauFromYields_MC3D->GetXaxis()->SetTitle("K^{0}_{S} c#tau (cm)");
  ksCtauFromYields_MC3D->GetXaxis()->SetTitleOffset(0.86);
  ksCtauFromYields_MC3D->SetMinimum(0.);
  ksCtauFromYields_MC3D->Draw();
  cout << "Here?" << endl;
  lablMC->Draw("same");

  c1->SaveAs("KsYieldsRelEff3D.png");
  c1->SaveAs("KsYieldsRelEff3D.root");
  c1->SaveAs("KsYieldsRelEff3D.pdf");
  gPad->SetLogy(1);
  gPad->Update();

  ksCtauFromYields3D_->Divide(ksCtauFromYields_MC3D);
  ksCtauFromYields3D_->Fit("ksCtauFromYieldsFit3D_", "REI");
  legen_ = new TLegend(0.2, 0.2, 0.5, 0.39);
  legen_->AddEntry(ksCtauFromYields3D_, "Corrected data", "lp");
  legen_->AddEntry(ksCtauFromYieldsFit3D_, "Exponential fit", "l");
  legen_->SetFillColor(0);
  legen_->SetBorderSize(0);
  legen_->SetTextSize(0.035);
  ksCtauFromYields3D_->GetYaxis()->SetTitle("Corrected K^{0}_{S} mass fit yield");
  ksCtauFromYields3D_->GetXaxis()->SetTitle("K^{0}_{S} c#tau (cm)");
  ksCtauFromYields3D_->GetXaxis()->SetTitleOffset(0.86);
  //ksCtauFromYields->Draw();
  //c1->cd(1);
  //gPad->SetLogy(1);
  //gPad->Update();
  //c1->SaveAs("YieldsScaledCtau_noFit.png");

  ksCtauFromYields3D_->Draw();
  legen_->Draw("same");
  labl2_0->Draw("same");

  c1->SaveAs("ksCtauFromYields_withFit3D.png");
  c1->SaveAs("ksCtauFromYields_withFit3D.root");
  c1->SaveAs("ksCtauFromYields_withFit3D.pdf");

  delete legen_;

  gPad->SetLogy(0);

  lamLifetimeNbins = 20;
  lamLifetimeXmin = 0.;
  lamLifetimeXmax = 16.;
  lamCtauNbins = 8;

  //gStyle->SetErrorX((lamLifetimeXmax - lamLifetimeXmin) / lamCtauNbins);

  fin->cd("analyzeLambda");

  std::vector<double> lamFitYields;
  std::vector<double> lamFitYieldErrors;
  int lamCanvasX = lamCtauNbins / 2;
  int lamCanvasY = 2;
  TCanvas* lamCanvas = new TCanvas( "Lam3Dmasses",
				    "#Lambda^{0} mass fits in c#tau bins",
				    1600, 800 );
  gStyle->SetOptTitle(1);
  lamCanvas->Divide( lamCanvasX, lamCanvasY );

  /*lamMFull = new TH1F("lamMFull", "#Lambda",
    50, 1.08, 1.18);

    ntuple->Project("lamMFull", "v0CandMass", "priVtxChi2 > 0. && v0VtxSig > 15. && v0Lifetime > 0. && v0Lifetime < 20. && abs(v0OtherCandMass - 0.498) > 0.02");
    lamMFull->Draw();
    c1->SaveAs("LamMassTest_0_ctau_20.png");*/

  for( int ndx = 0; ndx < lamCtauNbins; ndx++ ) {
    lamCanvas->cd(ndx + 1);
    ostringstream lamHisto;
    lamHisto << "lamMTemp" << ndx;
    lamCtauBinWidth = (lamLifetimeXmax - lamLifetimeXmin) / lamCtauNbins;
    lamMNbins = 50;
    lamMXmin = 1.08;
    lamMXmax = 1.18;
    lamMBinWidth = (lamMXmax - lamMXmin) / lamMNbins;

    ostringstream lamselection;
    lamselection << "priVtxChi2 > 0. && v0VtxSig > 15. && v03DLifetime > "
		 << ndx*lamCtauBinWidth << " && v03DLifetime < "
		 << (ndx+1)*lamCtauBinWidth
		 << " && abs(v0OtherCandMass - 0.498) > 0. && v0DauNormChi2 < 5.";
    cout << lamselection.str() << endl;
    lamMTemp = new TH1F(//"lamMTemp", "#Lambda^{0} mass",
			lamHisto.str().c_str(), labloss.str().c_str(),
			lamMNbins, lamMXmin, lamMXmax);

    ntuple->Project(lamHisto.str().c_str(), "v0CandMass", lamselection.str().c_str());
    lamMTemp->SetXTitle("m_{p^{+}#pi^{-}} ( + c.c.) (GeV/c^{2})");
    lamMTemp->SetYTitle("Entries/0.02 GeV/c^{2}");
    lamMTemp->SetNdivisions(506);

    ostringstream lamfitoss;
    lamfitoss << lamMBinWidth << "*gausn(0) + [3]*(x - "
	      << piMass + protonMass << ")^(1/2) + [4]*(x - "
	      << piMass + protonMass << ")^(3/2)";
    TF1* lamFit_SG = new TF1("lamFit_SG", lamfitoss.str().c_str(), lamMXmin, lamMXmax);
    lamFit_SG->SetParName( 0, "yield" );
    lamFit_SG->SetParName( 1, "GausMean" );
    lamFit_SG->SetParName( 2, "GausSigma" );
    lamFit_SG->SetParName( 3, "c1" );
    lamFit_SG->SetParName( 4, "c2" );

    lamFit_SG->SetParameter( 0, 2.6*lamMTemp->GetMaximum() );
    lamFit_SG->SetParameter( 1, 1.1156 );
    lamFit_SG->SetParameter( 2, 0.002 );
    lamFit_SG->SetParameter( 3, 0.75*lamMTemp->GetMaximum() );
    lamFit_SG->SetParameter( 4, 1.85*lamMTemp->GetMaximum() );

    int fitstat = lamMTemp->Fit("lamFit_SG", "RLE");

    ostringstream hmmm;
    hmmm << "lamMTemp3D_" << ndx << ".png";
    lamMTemp->Draw();
    //c1->SaveAs(hmmm.str().c_str());
    cout << "Fit status: " << fitstat << endl;
    if( !fitstat ) {
      lamFitYields->push_back( lamFit_SG->GetParameter(0) );
      lamFitYieldErrors->push_back( lamFit_SG->GetParError(0) );
    }
    else {
      lamFitYields->push_back( 0. );
      lamFitYieldErrors->push_back( 0. );
    }
    delete lamFit_SG;
    //delete lamMTemp;
  }

  lamCanvas->Print("LamMassFitsInCtauBins3D.png");
  lamCanvas->Print("LamMassFitsInCtauBins3D.root");
  delete lamCanvas;
  gStyle->SetOptTitle(0);

  lifein->cd();
  TF1* lamPDGCtau3D = new TF1(*lamCtauReal);
  TH1F* lamCtauFromYields_MC3D = new TH1F(*lamCtauFromYields3D);
  fin->cd("analyzeLambda");

  lamCtauFromYields3D = new TH1F("lamCtauFromYields3D", "#Lambda c#tau",
			       lamCtauNbins, lamLifetimeXmin, lamLifetimeXmax);

  lamCtauFromYields3D->Sumw2();

  for( int ndx = 0; ndx < fitYields.size(); ndx++ ) {
    lamCtauFromYields3D->SetBinContent( ndx + 1, lamFitYields[ndx] );
    lamCtauFromYields3D->SetBinError( ndx + 1, lamFitYieldErrors[ndx] );
  }

  lamCtauFromYields3D->GetXaxis()->SetLimits(lamLifetimeXmin, lamLifetimeXmax + 0.1);

  lamCtauFromYields3D->Draw();
  labl->Draw("same");
  lamCtauFromYields3D->SetXTitle("#Lambda^{0} c#tau (cm)");
  lamCtauFromYields3D->SetYTitle("#Lambda^{0} mass fit yield");
  lamCtauFromYields3D->GetXaxis()->SetTitleOffset(1);
  gPad->SetLogy(1);
  gPad->Update();
  c1->SaveAs("lamCtauFromYields_noScaling3D.png");
  c1->SaveAs("lamCtauFromYields_noScaling3D.root");
  c1->SaveAs("lamCtauFromYields_noScaling3D.pdf");
  gPad->SetLogy(0);

  ostringstream lamctauoss_;
  lamctauoss_ << "exp([0])*exp(-x/(" << c_const << "*[1]))";

  lamCtauFromYieldsFit3D = new TF1("lamCtauFromYieldsFit3D", lamctauoss_.str().c_str(),
				 lamLifetimeXmin, lamLifetimeXmax);

  Double_t exp_p0_ly = log( lamCtauFromYields3D->GetMaximum() );
  Double_t rise_ly = log( lamCtauFromYields3D->GetMaximum()
			  - lamCtauFromYields3D->GetBinContent(lamCtauNbins) );
  Double_t exp_p1_ly = (rise_ly / lamLifetimeXmax);
  exp_p1_ly = 1/c_const/exp_p1_ly;

  lamCtauFromYieldsFit3D->SetParameter(0, exp_p0_ly);
  lamCtauFromYieldsFit3D->SetParameter(1, exp_p1_ly);
  lamCtauFromYieldsFit3D->SetParName(1, "#tau (ps)");

  lamCtauFromYieldsFit3D->Draw();
  c1->SaveAs("lamFitFunc3D.png");

  Double_t realCtau_ly = 7.89;
  Double_t realCtau_p1_ly = realCtau_ly/c_const;
  lamCtauReal_ly3D = new TF1("lamCtauReal_ly3D", lamctauoss_.str().c_str(),
			   lamLifetimeXmin, lamLifetimeXmax);
  Double_t realCtau_p0_ly = //lamCtauFromYieldsFit3D->GetParameter(0);
    log(lamCtauFromYields_MC3D->GetMaximum());
  lamCtauReal_ly3D->SetParameter(0, realCtau_p0_ly);
  lamCtauReal_ly3D->SetParameter(1, realCtau_p1_ly);

  lamCtauFromYields_MC3D->Divide( lamCtauReal_ly3D );
  c1->cd(1);
  gPad->SetLogy(0);
  gPad->Update();
  lamCtauFromYields_MC3D->GetYaxis()->SetTitle("Correction factor");
  lamCtauFromYields_MC3D->SetXTitle("#Lambda^{0} c#tau (cm)");
  lamCtauFromYields_MC3D->GetXaxis()->SetTitleOffset(1);
  lamCtauFromYields_MC3D->SetMinimum(0.);
  lamCtauFromYields_MC3D->Draw();
  lablMC->Draw("same");

  c1->SaveAs("LamYieldsRelEff3D.png");
  c1->SaveAs("LamYieldsRelEff3D.root");
  c1->SaveAs("LamYieldsRelEff3D.pdf");
  gPad->SetLogy(1);
  gPad->Update();
  gStyle->SetFitFormat("3.0f");

  lamCtauFromYields3D->Divide(lamCtauFromYields_MC3D);
  lamCtauFromYields3D->Fit("lamCtauFromYieldsFit3D", "REI");
  legen_ = new TLegend(0.2, 0.2, 0.5, 0.39);
  legen_->AddEntry(lamCtauFromYields3D, "Corrected data", "lp");
  legen_->AddEntry(lamCtauFromYieldsFit3D, "Exponential fit", "l");
  legen_->SetFillColor(0);
  legen_->SetBorderSize(0);
  legen_->SetTextSize(0.035);
  lamCtauFromYields3D->GetYaxis()->SetTitle("Corrected #Lambda^{0} mass fit yield");
  lamCtauFromYields3D->GetXaxis()->SetTitle("#Lambda^{0} c#tau (cm)");
  lamCtauFromYields3D->GetXaxis()->SetTitleOffset(1);

  lamCtauFromYields3D->Draw();
  legen_->Draw("same");
  labl2_0->Draw("same");
  //delete legen_;

  c1->SaveAs("lamCtauFromYields_withFit3D.png");
  c1->SaveAs("lamCtauFromYields_withFit3D.root");
  c1->SaveAs("lamCtauFromYields_withFit3D.pdf");
  delete legen_;
  gStyle->SetFitFormat("5.5g");

}

if(0) {
.q;

gStyle->SetErrorX(0);

fin->cd("analyzeKshort");
// ks lifetime, using MC Scaling with correct lifetime

ksLifetimeNbins = 18;
ksLifetimeXmin = 0.;
ksLifetimeXmax = 9.;

//gStyle->SetOptStat(0);
ksLifetimeLoc = new TH1F("ksCtauLoc", "K^{0}_{S} c#tau",
		      ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);
ksLifetimeLoc->Sumw2();
ntuple->Project("ksCtauLoc", "v0Lifetime", 
		"v0VtxSig > 15. && v0Lifetime > 0.  && abs(v0CandMass - 0.498) < 0.02 && priVtxChi2 > 0. && abs(v0OtherCandMass - 1.1158) > 0. && v0DauNormChi2 < 5.");
ksLifetimeBG = new TH1F("ksCtauBG", "K^{0}_{S} c#tau background",
			ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);
ksLifetimeBG->Sumw2();
ntuple->Project("ksCtauBG", "v0Lifetime",
		"v0VtxSig > 15. && v0Lifetime > 0.  && (abs(v0CandMass - (0.498 - 0.06)) < 0.02 || abs(v0CandMass - (0.498 + 0.06)) < 0.02) && priVtxChi2 > 0. && abs(v0OtherCandMass - 1.1158) > 0. && v0DauNormChi2 < 5.");
ksLifetimeLoc->Add( ksLifetimeBG, -0.5 );

ksLifetimeLoc->GetXaxis()->SetTitle("K^{0}_{S} c#tau (cm)");
ksLifetimeLoc->GetYaxis()->SetTitle("Entries/0.5 cm");
ksLmax = ksLifetimeLoc->GetMaximum();
ksLmin = ksLifetimeLoc->GetBinContent( ksLifetimeNbins );

Double_t spdolight = 2.99792e10;
ostringstream ksctauoss;
//ksctauoss << "expo(0)";
ksctauoss << "exp([0])*exp(-x/(" << spdolight << "*[1]))";
cout << ksctauoss.str() << endl;
ksCtauFitLoc = new TF1("ksCtauFitLoc", ksctauoss.str().c_str(), 
		    ksLifetimeXmin, ksLifetimeXmax );

Double_t exp_p0 = log( ksLmax );
Double_t rise = log( ksLmax - ksLmin );
Double_t exp_p1 = ( rise / ksLifetimeXmax );
exp_p1 = 1/ spdolight / exp_p1;
cout << exp_p1 << endl;;
//ksCtauFit->SetParameter( 0, 11.56252 );
//ksCtauFit->SetParameter( 1, -1 );
ksCtauFitLoc->SetParameter( 0, exp_p0 );
ksCtauFitLoc->SetParameter( 1, exp_p1 );
ksCtauFitLoc->SetParName(1, "lifetime");

ksLifetimeLoc->Fit("ksCtauFitLoc", "RE");

Double_t realCtau = 2.6786;
//Double_t realCtau_p1 = -1./realCtau;
Double_t realCtau_p1 = realCtau/spdolight;
ksCtauRealLoc = new TF1("ksCtauRealLoc", ksctauoss.str().c_str(),
		     ksLifetimeXmin, ksLifetimeXmax);
Double_t realCtau_p0 = ksCtauFitLoc->GetParameter(0);
ksCtauRealLoc->SetParameter( 0, realCtau_p0 );
ksCtauRealLoc->SetParameter( 1, realCtau_p1 );

//bin20 = ksLifetime->GetBinContent(ksLifetimeNbins);
//cout << "bin20 = " << bin20 << endl;

fout->cd();
ksLifetimeLoc->Write();
ksCtauRealLoc->Write();
ksCtauFitLoc->Write();
fout->Write();
fin->cd("analyzeKshort");

ksLifetimeLoc->Draw();
labl2_0->Draw("same");
//ksCtauReal->Draw("same");

c1->cd(1);
gPad->SetLogy(1);

c1->SaveAs("ksCtau_noScaling.png");
c1->SaveAs("ksCtau_noScaling.root");
//c1->SaveAs("ksCtau.ps");

lifein->cd();
//ksCtauFit->Draw();
TF1* ksLTFit = new TF1(*ksCtauFit);
ksLTFit->Draw();
c1->SaveAs("ksLTFit.png");
TF1* ksLTReal = new TF1(*ksCtauReal);
//fin->cd("analyzeKshort");
TH1F* ksLTHist = new TH1F(*ksCtau);
ksLTHist->Draw();
c1->SaveAs("ltHistTest.png");
ksLTReal->Draw();
c1->SaveAs("ltFunc.png");
ksLTHist->Divide(ksLTReal);
ksLTHist->GetYaxis()->SetTitle("Correction factor");
ksLTHist->SetXTitle("K^{0}_{S} c#tau (cm)");
ksLTHist->Draw();
labl2_0->Draw("same");
gPad->SetLogy(0);
c1->SaveAs("KsRelEff.png");
c1->SaveAs("KsRelEff.root");
gPad->SetLogy(0);
fin->cd("analyzeKshort");
ksLifetimeLoc->Divide(ksLTHist);
//ksLifetime->Multiply(ksLTReal);
ksLifetimeLoc->Fit("ksCtauFitLoc", "REI");
ksLifetimeLoc->SetYTitle("Corrected entries/0.5 cm");
ksLifetimeLoc->Draw();
labl2_0->Draw("same");
//ksLTReal->Draw("same");
gPad->SetLogy(1);

c1->SaveAs("ksCtau_scaled.png");
c1->SaveAs("ksCtau_scaled.root");

gPad->SetLogy(0);
c1->SaveAs("ksCtau_scaled_noLogScale.png");

// 3D Lifetime:
//***************************************************

//gStyle->SetOptStat(0);
ksLifetime3D = new TH1F("ksCtau3D", "K^{0}_{S} c#tau 3D",
		      ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);
ksLifetime3D->Sumw2();
ntuple->Project("ksCtau3D", "v03DLifetime", 
		"v0VtxSig > 15. && v03DLifetime > 0.  && abs(v0CandMass - 0.498) < 0.02  && priVtxChi2 > 0. && abs(v0OtherCandMass - 1.1158) > 0. && v0DauNormChi2 < 5.");
ksLifetimeBG3D = new TH1F("ksCtauBG3D", "K^{0}_{S} c#tau background",
			ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);
ksLifetimeBG3D->Sumw2();
ntuple->Project("ksCtauBG3D", "v03DLifetime",
		"v0VtxSig > 15. && v03DLifetime > 0.  && (abs(v0CandMass - (0.498 - 0.06)) < 0.02 || abs(v0CandMass - (0.498 + 0.06)) < 0.02)  && priVtxChi2 > 0. && abs(v0OtherCandMass - 1.1158) > 0. && v0DauNormChi2 < 5.");
ksLifetime3D->Add( ksLifetimeBG3D, -0.5 );

ksLifetime3D->GetXaxis()->SetTitle("K^{0}_{S} c#tau (cm)");
ksLifetime3D->GetYaxis()->SetTitle("Entries/0.5 cm");
ksLmax3d = ksLifetime3D->GetMaximum();
ksLmin3d = ksLifetime3D->GetBinContent( ksLifetimeNbins );

//Double_t spdolight = 2.99792e10;
ostringstream ksctauoss3d;
//ksctauoss << "expo(0)";
ksctauoss3d << "exp([0])*exp(-x/(" << spdolight << "*[1]))";
cout << ksctauoss3d.str() << endl;
ksCtauFit3D = new TF1("ksCtauFit3D", ksctauoss3d.str().c_str(), 
		    ksLifetimeXmin, ksLifetimeXmax );

Double_t exp_p03d = log( ksLmax3d );
Double_t rise3d = log( ksLmax3d - ksLmin3d );
Double_t exp_p13d = ( rise3d / ksLifetimeXmax );
exp_p13d = 1/ spdolight / exp_p13d;
cout << exp_p13d << endl;;
//ksCtauFit->SetParameter( 0, 11.56252 );
//ksCtauFit->SetParameter( 1, -1 );
ksCtauFit3D->SetParameter( 0, exp_p03d );
ksCtauFit3D->SetParameter( 1, exp_p13d );
ksCtauFit3D->SetParName(1, "lifetime3d");

ksLifetime3D->Fit("ksCtauFit3D", "RE");

//Double_t realCtau = 2.6786;
//Double_t realCtau_p1 = -1./realCtau;
Double_t realCtau_p13d = realCtau/spdolight;
ksCtauReal3D = new TF1("ksCtauReal3D", ksctauoss3d.str().c_str(),
		     ksLifetimeXmin, ksLifetimeXmax);
Double_t realCtau_p03d = ksCtauFit3D->GetParameter(0);
ksCtauReal3D->SetParameter( 0, realCtau_p03d );
ksCtauReal3D->SetParameter( 1, realCtau_p13d );

//bin20 = ksLifetime->GetBinContent(ksLifetimeNbins);
//cout << "bin20 = " << bin20 << endl;

fout->cd();
ksLifetime3D->Write();
ksCtauReal3D->Write();
ksCtauFit3D->Write();
fout->Write();
fin->cd("analyzeKshort");

ksLifetime3D->Draw();
labl2_0->Draw("same");
//ksCtauReal->Draw("same");

c1->cd(1);
gPad->SetLogy(1);

c1->SaveAs("ksCtau_noScaling3D.png");
c1->SaveAs("ksCtau_noScaling3D.root");
//c1->SaveAs("ksCtau.ps");

lifein->cd();
//ksCtauFit->Draw();
TF1* ksLTFit3D = new TF1(*ksCtauFit3D);
ksLTFit3D->Draw();
c1->SaveAs("ksLTFit3D.png");
TF1* ksLTReal3D = new TF1(*ksCtauReal3D);
//fin->cd("analyzeKshort");
TH1F* ksLTHist3D = new TH1F(*ksCtau3D);
ksLTHist3D->Draw();
c1->SaveAs("ltHistTest3D.png");
ksLTReal3D->Draw();
c1->SaveAs("ltFunc3D.png");
ksLTHist3D->Divide(ksLTReal3D);
ksLTHist3D->GetYaxis()->SetTitle("Correction factor");
ksLTHist3D->SetXTitle("K^{0}_{S} c#tau (cm)");
ksLTHist3D->Draw();
labl2_0->Draw("same");
gPad->SetLogy(0);
c1->SaveAs("KsRelEff3D.png");
c1->SaveAs("KsRelEff3D.root");
gPad->SetLogy(0);
fin->cd("analyzeKshort");
ksLifetime3D->Divide(ksLTHist3D);
//ksLifetime->Multiply(ksLTReal);
ksLifetime3D->Fit("ksCtauFit3D", "REI");
ksLifetime3D->SetYTitle("Corrected entries/0.5 cm");
ksLifetime3D->Draw();
labl2_0->Draw("same");
//ksLTReal->Draw("same");
gPad->SetLogy(1);

c1->SaveAs("ksCtau_scaled3D.png");
c1->SaveAs("ksCtau_scaled3D.root");

gPad->SetLogy(0);
c1->SaveAs("ksCtau_scaled_noLogScale3D.png");

// END 3D LIFETIME

// Lambda lifetime
// ******************************************

fin->cd("analyzeLambda");

lamLifetimeNbins = 20;//line number 280
lamLifetimeXmin = 0.;
lamLifetimeXmax = 10.;

//gStyle->SetOptStat(0);
lamLifetimeLoc = new TH1F("lamCtauLoc", "#Lambda^{0} c#tau",
		      lamLifetimeNbins, lamLifetimeXmin, lamLifetimeXmax);
lamLifetimeLoc->Sumw2();
ntuple->Project("lamCtauLoc", "v0Lifetime", 
		"v0VtxSig > 15. && v0Lifetime > 0.  && abs(v0CandMass - 1.1158) < 0.0075 && abs(v0OtherCandMass - 0.498) > 0. && v0DauNormChi2 < 5.");

lamLifetimeBG = new TH1F("lamCtauBG", "#Lambda^{0} c#tau background",
			lamLifetimeNbins, lamLifetimeXmin, lamLifetimeXmax);
lamLifetimeBG->Sumw2();
ntuple->Project("lamCtauBG", "v0Lifetime",
		"v0VtxSig > 15. && v0Lifetime > 0.  && (abs(v0CandMass - 1.14) < 0.0075 || abs(v0CandMass - 1.095) < 0.0075) && abs(v0OtherCandMass - 0.498) > 0. && v0DauNormChi2 < 5.");
lamLifetimeLoc->Add( lamLifetimeBG, -0.5 );

lamLifetimeLoc->GetXaxis()->SetTitle("#Lambda^{0} c#tau (cm)");
lamLifetimeLoc->GetYaxis()->SetTitle("Entries/0.5 cm");
lamLifetimeLoc->GetXaxis()->SetLimits(lamLifetimeXmin, lamLifetimeXmax + 0.1);
lamLmax = lamLifetimeLoc->GetMaximum();
lamLmin = lamLifetimeLoc->GetBinContent( lamLifetimeNbins );

spdolight = 2.99792e10;
ostringstream lamctauoss;
//lamctauoss << "expo(0)";
lamctauoss << "exp([0])*exp(-x/(" << spdolight << "*[1]))";
cout << lamctauoss.str() << endl;
lamCtauFitLoc = new TF1("lamCtauFitLoc", lamctauoss.str().c_str(), 
		    lamLifetimeXmin, lamLifetimeXmax );

Double_t lamexp_p0 = log( lamLmax );
Double_t lamrise = log( lamLmax - lamLmin );
Double_t lamexp_p1 = ( lamrise / lamLifetimeXmax );
lamexp_p1 = 1/ spdolight / lamexp_p1;
cout << lamexp_p1 << endl;;
//lamCtauFit->SetParameter( 0, 11.56252 );
//lamCtauFit->SetParameter( 1, -1 );
lamCtauFitLoc->SetParameter( 0, lamexp_p0 );
lamCtauFitLoc->SetParameter( 1, lamexp_p1 );
lamCtauFitLoc->SetParName(1, "lifetime");

lamLifetimeLoc->Fit("lamCtauFitLoc", "RE");

Double_t lamrealCtau = 7.89;
//Double_t realCtau_p1 = -1./realCtau;
Double_t lamrealCtau_p1 = lamrealCtau/spdolight;
lamCtauRealLoc = new TF1("lamCtauRealLoc", lamctauoss.str().c_str(),
		     lamLifetimeXmin, lamLifetimeXmax);
Double_t lamrealCtau_p0 = lamCtauFitLoc->GetParameter(0);
lamCtauRealLoc->SetParameter( 0, lamrealCtau_p0 );
lamCtauRealLoc->SetParameter( 1, lamrealCtau_p1 );

//bin20 = lamLifetime->GetBinContent(lamLifetimeNbins);
//cout << "bin20 = " << bin20 << endl;

fout->cd();
lamLifetimeLoc->Write();
lamCtauRealLoc->Write();
lamCtauFitLoc->Write();
fout->Write();
fin->cd("analyzeLambda");

lamLifetimeLoc->Draw();
labl2_0->Draw("same");
//lamCtauReal->Draw("same");

c1->cd(1);
gPad->SetLogy(1);

c1->SaveAs("lamCtau_noScaling.png");
c1->SaveAs("lamCtau_noScaling.root");
//c1->SaveAs("lamCtau.ps");

lifein->cd();
//lamCtauFit->Draw();
TF1* lamLTFit = new TF1(*lamCtauFit);
lamLTFit->Draw();
c1->SaveAs("lamLTFit.png");
TF1* lamLTReal = new TF1(*lamCtauReal);
//fin->cd("analyzeLambda");
TH1F* lamLTHist = new TH1F(*lamCtau);
lamLTHist->Draw();
c1->SaveAs("lamltHistTest.png");
lamLTReal->Draw();
c1->SaveAs("lamltFunc.png");
lamLTHist->Divide(lamLTReal);
lamLTHist->GetYaxis()->SetTitle("Correction factor");
lamLTHist->SetXTitle("#Lambda^{0} c#tau (cm)");
lamLTHist->GetXaxis()->SetLimits(lamLifetimeXmin, lamLifetimeXmax + 0.1);
lamLTHist->Draw();
labl2_0->Draw("same");
gPad->SetLogy(0);
c1->SaveAs("LamRelEff.png");
c1->SaveAs("LamRelEff.root");
gPad->SetLogy(0);
fin->cd("analyzeLambda");
lamLifetimeLoc->Divide(lamLTHist);
//lamLifetime->Multiply(lamLTReal);

lamLmax = lamLifetimeLoc->GetMaximum();
lamLmin = lamLifetimeLoc->GetBinContent( lamLifetimeNbins );

Double_t lamexp_p0x = log( lamLmax );
Double_t lamrisex = log( lamLmax - lamLmin );
Double_t lamexp_p1x = ( lamrisex / lamLifetimeXmax );
lamexp_p1x = 1/ spdolight / lamexp_p1x;
cout << lamexp_p1x << endl;;
//lamCtauFit->SetParameter( 0, 11.56252 );
//lamCtauFit->SetParameter( 1, -1 );
lamCtauFitLoc->SetParameter( 0, lamexp_p0x );
lamCtauFitLoc->SetParameter( 1, lamexp_p1x );
//lamCtauFitLoc->SetParName(1, "lifetime");

lamLifetimeLoc->Fit("lamCtauFitLoc", "REI");
lamLifetimeLoc->SetYTitle("Corrected entries/0.5 cm");
lamLifetimeLoc->Draw();
labl2_0->Draw("same");
//lamLTReal->Draw("same");
gPad->SetLogy(1);

c1->SaveAs("lamCtau_scaled.png");
c1->SaveAs("lamCtau_scaled.root");

gPad->SetLogy(0);
c1->SaveAs("lamCtau_scaled_noLogScale.png");

// 3D Lifetime:
//***************************************************

lamLifetimeXmin = 1.;

//gStyle->SetOptStat(0);
lamLifetime3D = new TH1F("lamCtau3D", "#Lambda^{0} c#tau 3D",
		      lamLifetimeNbins, lamLifetimeXmin, lamLifetimeXmax);
lamLifetime3D->Sumw2();
ntuple->Project("lamCtau3D", "v03DLifetime", 
		"v0VtxSig > 15. && v03DLifetime > 0.  && abs(v0CandMass - 1.1158) < 0.0075 && abs(v0OtherCandMass - 0.498) > 0. && v0DauNormChi2 < 5.");
lamLifetimeBG3D = new TH1F("lamCtauBG3D", "#Lambda^{0} c#tau background",
			lamLifetimeNbins, lamLifetimeXmin, lamLifetimeXmax);
lamLifetimeBG3D->Sumw2();
ntuple->Project("lamCtauBG3D", "v03DLifetime",
		"v0VtxSig > 15. && v03DLifetime > 0.  && (abs(v0CandMass - 1.14) < 0.0075 || abs(v0CandMass - 1.095) < 0.0075) && abs(v0OtherCandMass - 0.498) > 0. && v0DauNormChi2 < 5.");
lamLifetime3D->Add( lamLifetimeBG3D, -0.5 );

lamLifetime3D->GetXaxis()->SetTitle("#Lambda^{0} c#tau (cm)");
lamLifetime3D->GetYaxis()->SetTitle("Entries/0.5 cm");
lamLifetime3D->GetXaxis()->SetLimits(lamLifetimeXmin, lamLifetimeXmax + 0.1);
lamLmax3d = lamLifetime3D->GetMaximum();
lamLmin3d = lamLifetime3D->GetBinContent( lamLifetimeNbins );

//Double_t spdolight = 2.99792e10;
ostringstream lamctauoss3d;
//lamctauoss << "expo(0)";
lamctauoss3d << "exp([0])*exp(-x/(" << spdolight << "*[1]))";
cout << lamctauoss3d.str() << endl;
lamCtauFit3D = new TF1("lamCtauFit3D", lamctauoss3d.str().c_str(), 
		    lamLifetimeXmin, lamLifetimeXmax );

Double_t lamexp_p03d = log( lamLmax3d );
Double_t lamrise3d = log( lamLmax3d - lamLmin3d );
Double_t lamexp_p13d = ( lamrise3d / lamLifetimeXmax );
lamexp_p13d = 1/ spdolight / lamexp_p13d;
cout << lamexp_p13d << endl;;
//lamCtauFit->SetParameter( 0, 11.56252 );
//lamCtauFit->SetParameter( 1, -1 );
lamCtauFit3D->SetParameter( 0, lamexp_p03d );
lamCtauFit3D->SetParameter( 1, lamexp_p13d );
lamCtauFit3D->SetParName(1, "lifetime3d");

lamLifetime3D->Fit("lamCtauFit3D", "RE");

//Double_t realCtau = 2.6786;
//Double_t realCtau_p1 = -1./realCtau;
Double_t lamrealCtau_p13d = lamrealCtau/spdolight;
lamCtauReal3D = new TF1("lamCtauReal3D", lamctauoss3d.str().c_str(),
		     lamLifetimeXmin, lamLifetimeXmax);
Double_t lamrealCtau_p03d = lamCtauFit3D->GetParameter(0);
lamCtauReal3D->SetParameter( 0, lamrealCtau_p03d );
lamCtauReal3D->SetParameter( 1, lamrealCtau_p13d );

//bin20 = lamLifetime->GetBinContent(lamLifetimeNbins);
//cout << "bin20 = " << bin20 << endl;

fout->cd();
lamLifetime3D->Write();
lamCtauReal3D->Write();
lamCtauFit3D->Write();
fout->Write();
fin->cd("analyzeLambda");

//lamLifetime3D->SetYTitle("Entries/0.5 cm");
lamLifetime3D->Draw();
labl2_0->Draw("same");
//lamCtauReal->Draw("same");

c1->cd(1);
gPad->SetLogy(1);

c1->SaveAs("lamCtau_noScaling3D.png");
c1->SaveAs("lamCtau_noScaling3D.root");
//c1->SaveAs("lamCtau.ps");

lifein->cd();
//lamCtauFit->Draw();
TF1* lamLTFit3D = new TF1(*lamCtauFit3D);
lamLTFit3D->Draw();
c1->SaveAs("lamLTFit3D.png");
TF1* lamLTReal3D = new TF1(*lamCtauReal3D);
//fin->cd("analyzeLambda");
TH1F* lamLTHist3D = new TH1F(*lamCtau3D);
//gPad->SetLogy(0);
//gPad->Update();
lamLTHist3D->Draw();
c1->SaveAs("lamltHistTest3D.png");
lamLTReal3D->Draw();
c1->SaveAs("lamltFunc3D.png");
lamLTHist3D->Divide(lamLTReal3D);
lamLTHist3D->GetYaxis()->SetTitle("Correction factor");
lamLTHist3D->SetXTitle("#Lambda^{0} c#tau (cm)");
lamLTHist3D->GetXaxis()->SetLimits(lamLifetimeXmin, lamLifetimeXmax + 0.1);
lamLTHist3D->Draw();
labl2_0->Draw("same");
gPad->SetLogy(0);
c1->SaveAs("LamRelEff3D.png");
c1->SaveAs("LamRelEff3D.root");
gPad->SetLogy(0);
fin->cd("analyzeLambda");
lamLifetime3D->Divide(lamLTHist3D);
//lamLifetime->Multiply(lamLTReal);
lamLifetime3D->Fit("lamCtauFit3D", "REI");
lamLifetime3D->SetYTitle("Corrected entries/0.5 cm");
lamLifetime3D->Draw();
labl2_0->Draw("same");
//lamLTReal->Draw("same");
gPad->SetLogy(1);

c1->SaveAs("lamCtau_scaled3D.png");
c1->SaveAs("lamCtau_scaled3D.root");

gPad->SetLogy(0);
c1->SaveAs("lamCtau_scaled_noLogScale3D.png");

// END 3D LIFETIME
// END LAMBDA LIFETIME

//gStyle->SetOptStat("i");

ksMassNbins = 200;
ksMassXmin = 0.351;
ksMassXmax = 0.751;
ksMassBinWidth = (ksMassXmax-ksMassXmin)/ksMassNbins;
// Fit function:


// MASS:

TH1F *mcmass;
TH1F *lammcmass;
fmassin->cd();
//TH1F* mcmass = new TH1F(*ksmouthist);
mcmass = new TH1F(*ksMassTunable);
lammcmass = new TH1F(*lamMassTunable);
fin->cd("analyzeKshort");

//ksCandMass->Fit("ksFit_doubGaus0", "RLE");

ostringstream ksoss;
ksoss << ksMassBinWidth << "*gausn(0) + [3] + [4]*(x-[1])"
      << " + [5]*(x-[1])^2";
TF1 *ksFitSG = new TF1("ksFitSG", ksoss.str().c_str(), ksMassXmin, ksMassXmax);
ksFitSG->SetParName( 0, "yield" );
ksFitSG->SetParName( 1, "mean" );
ksFitSG->SetParName( 2, "sigma" );
ksFitSG->SetParName( 3, "const" );
ksFitSG->SetParName( 4, "slope" );
ksFitSG->SetParName( 5, "quadconst" );
//ksFitSG->SetParName( 6, "center" );

ksFitSG->SetParameter( 0, 4.4444*ksCandMass->GetMaximum() );
ksFitSG->SetParameter( 1, 0.49767 );
ksFitSG->SetParameter( 2, 0.005 );
ksFitSG->SetParLimits( 2, 0., 0.010 );
ksFitSG->SetParameter( 3, 0.1*ksCandMass->GetMaximum() );
ksFitSG->SetParameter( 4, 1. );
ksFitSG->SetParameter( 5, -1. );
//ksFitSG->SetParameter( 6, 0.49767 );

ksMassTunable = new TH1F("ksMassTunable", "K^{0}_{S} mass",
			 ksMassNbins, ksMassXmin, ksMassXmax);
ksMassTunable2 = new TH1F("ksMassTunable2", "K^{0}_{S} mass",
			 ksMassNbins, ksMassXmin, ksMassXmax);
ntuple->Project("ksMassTunable", "v0CandMass", "v0VtxSig > 15. && v0DauNormChi2 < 5. && priVtxChi2 > 0. && abs(v0Eta) < 1.");
//ntuple->Project("ksMassTunable", "v0CandMass", "");
//ntuple->Project("ksMassTunable2", "v0CandMass", "nTracks < 150 && priVtxChi2 > 0.");



ksMassBinWidth_DG = ksMassBinWidth/sqrt(2*3.141592654);
ostringstream ksoss_doubgaus0;

ksoss_doubgaus0 << ksMassBinWidth_DG 	
                << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		<< ksMassBinWidth_DG 
		<< "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		<< "[5] + [6]*(x-[1]) + [7]*(x-[1])^2";
cout << ksoss_doubgaus0.str() << endl;
TF1 *ksFit_doubGaus0 = new TF1("ksFit_doubGaus0", 
			       ksoss_doubgaus0.str().c_str(),
			       ksMassXmin, ksMassXmax);
ksFit_doubGaus0->SetParName( 0, "yield" );
ksFit_doubGaus0->SetParName( 1, "mean" );
ksFit_doubGaus0->SetParName( 2, "sigma1" );
ksFit_doubGaus0->SetParName( 3, "fraction" );
ksFit_doubGaus0->SetParLimits( 3, 0., 1. );
ksFit_doubGaus0->SetParName( 4, "sigma2" );
ksFit_doubGaus0->SetParName( 5, "const" );
ksFit_doubGaus0->SetParName( 6, "slope" );
ksFit_doubGaus0->SetParName( 7, "quadconst" );

ksFit_doubGaus0->SetParameter( 0, 4.4444*ksMassTunable->GetMaximum() );
ksFit_doubGaus0->SetParameter( 1, 0.49767 );
ksFit_doubGaus0->SetParameter( 2, 0.0045 );
ksFit_doubGaus0->SetParameter( 3, 0.75 );
ksFit_doubGaus0->SetParameter( 4, 0.011 );
ksFit_doubGaus0->SetParameter( 5, 0.1*ksMassTunable->GetMaximum() );
ksFit_doubGaus0->SetParameter( 6, 1. );
ksFit_doubGaus0->SetParameter( 7, -1. );
//ntuple->Project("ksMassTunable", "v0CandMass");

ksMassTunable->GetXaxis()->SetTitle("#pi^{+}#pi^{-} invariant mass (GeV/c^{2})");
ksMassTunable->GetYaxis()->SetTitle("Entries/.002 GeV/c^{2}");
//ksMassTunable2->GetXaxis()->SetTitle("#pi^{+}#pi^{-} invariant mass (GeV/c^{2})");
//ksMassTunable2->GetYaxis()->SetTitle("Entries/.002 GeV/c^{2}");
ksMassTunable->Fit("ksFitSG", "RLE");
ksxmint = ksMassTunable->GetXaxis()->GetXmin();
ksxmaxt = ksMassTunable->GetXaxis()->GetXmax();
ksxmaxt -= 0.0001;
ksymaxt = ksMassTunable->GetMaximum();
ksymaxt += 0.1*ksymaxt;
ksMassTunable->SetAxisRange(0., ksymaxt, "Y");
ksMassTunable->GetXaxis()->SetLimits(ksxmint, ksxmaxt);
ksMassTunable->GetXaxis()->SetNdivisions(506);

ksMassTunable->Draw();

// PDG mass
kspdgmasst = new TLine(0.49767, 0., 0.49767, ksymaxt);
kspdgmasst->SetLineColor(9);
kspdgmasst->SetLineWidth(2);
kspdgmasst->SetLineStyle(2);
kspdgmasst->Draw("same");

//ksconstt = ksMassTunable->GetFunction("ksFitSG")->GetParameter(3);
ksconstt = ksFitSG->GetParameter(3);

//ksslopet = ksMassTunable->GetFunction("ksFitSG")->GetParameter(4);
ksslopet = ksFitSG->GetParameter(4);
//ksquadt = ksMassTunable->GetFunction("ksFitSG")->GetParameter(5);
ksquadt = ksFitSG->GetParameter(5);
cout << "ksconstt: " << ksconstt << ", ksslopet: "<< ksslopet << ", ksquadt: "
<< ksquadt << endl;
Double_t range_mint = 0.;
Double_t range_maxt = 0.;
ksMassTunable->GetFunction("ksFitSG")->GetRange(range_mint, range_maxt);
//TF1* ksbgt = new TF1("ksBgSG", "[0] + [1]*(x-0.49767) + [2]*(x-0.49767)^2", range_mint, range_maxt);
TF1* ksbgt = new TF1("ksBgSG", "[0] + [1]*(x-0.49767) + [2]*(x-0.49767)^2", ksMassXmin, ksMassXmax);
ksbgt->SetParName(0, "const");
ksbgt->SetParName(1, "slope");
ksbgt->SetParameter(0, ksconstt);
ksbgt->SetParameter(1, ksslopet);
ksbgt->SetParName(2, "quadconst");
ksbgt->SetParameter(2, ksquadt);
ksbgt->SetLineColor(2);
ksbgt->SetLineWidth(2);
ksbgt->SetLineStyle(3);
ksbgt->Draw("same");

labl->Draw("same");

c1->SaveAs("ksMass_tuned.png");


/*ksxmint2 = ksMassTunable->GetXaxis()->GetXmin();
ksxmaxt2 = ksMassTunable->GetXaxis()->GetXmax();
ksxmaxt2 -= 0.0001;
ksymaxt2 = ksMassTunable->GetMaximum();
ksymaxt2 += 0.1*ksymaxt;
ksMassTunable2->SetAxisRange(0., ksymaxt2, "Y");
ksMassTunable2->GetXaxis()->SetLimits(ksxmint2, ksxmaxt2);
ksMassTunable2->GetXaxis()->SetNdivisions(506);*/
ksMassTunable->Fit("ksFit_doubGaus0", "RLE");
ksMassTunable->Draw();
ksconstt2 = ksMassTunable->GetFunction("ksFit_doubGaus0")->GetParameter(5);
ksslopet2 = ksMassTunable->GetFunction("ksFit_doubGaus0")->GetParameter(6);
ksqc2 = ksMassTunable->GetFunction("ksFit_doubGaus0")->GetParameter(7);
TF1* ksbgt2 = new TF1("ksBgSG2", "[0] + [1]*(x-0.49767) + [2]*(x-0.49767)^2", ksMassXmin, ksMassXmax);
ksbgt2->SetParName(0, "const");
ksbgt2->SetParName(1, "slope");
ksbgt2->SetParameter(0, ksconstt2);
ksbgt2->SetParameter(1, ksslopet2);
ksbgt2->SetParName(2, "quadconst");
ksbgt2->SetParameter(2, ksqc2);
ksbgt2->SetLineColor(2);
ksbgt2->SetLineWidth(2);
ksbgt2->SetLineStyle(3);
ksbgt2->Draw("same");
kspdgmasst->Draw("same");
labl->Draw("same");
//xksbgt->Draw("same");

c1->SaveAs("ksMass_tuned_doubGaus.png");

// MAKE A PLOT WITH SUPERIMPOSED MASSES
//TCanvas *c2 = new TCanvas("c1","different scales hists",600,400); //create, fill and draw h1

leg = new TLegend(0.7,0.5,0.92,0.69); 
leg->AddEntry(mcmass,"Monte Carlo","l"); 
leg->AddEntry(ksMassTunable,"Data","l"); 
leg->AddEntry(ksFit_doubGaus0, "Data Fit", "l");
leg->SetFillColor(0);
Float_t rightmax = mcmass->Integral();
Float_t scale = ksMassTunable->Integral()/rightmax;

cout << "Integral(): " << rightmax << ", Integral(bin1, binN): "
<< mcmass->Integral(1, ksMassNbins) << endl << endl;

mcmass->SetLineColor(kGreen);
//mcmass->Scale(scale*1.5); 
scale = 293446./1000000.;
mcmass->Scale(scale);
mcmass->GetXaxis()->SetNdivisions(506);
mcmass->Draw(); 
ksMassTunable->Fit("ksFit_doubGaus0", "RLE");
ksMassTunable->Draw("same");
kspdgmasst->Draw("same");
ksbgt2->Draw("same");
labl->Draw("same");

c1->Update();
//draw an axis on the right side
TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,508,"+L");
//axis->SetLineColor(30;)
axis->SetLabelColor(50);
axis->SetLabelOffset(-0.5);
//axis->Draw();

//c1->SaveAs("ksMass_superimposed.png");

ksEtaSignal->GetXaxis()->SetTitle("K^{0}_{S} #eta");
ksEtaSignal->GetYaxis()->SetTitle("Entries");
ksEtaSignal->Draw();

c1->SaveAs("ksEtaSignal.png");

//gStyle->SetOptFit(0);
//c1->SaveAs("ksMass_superimposed_noFitBox.png");
ksMassTunable->Draw(); 
ksymax1 = mcmass->GetMaximum();
ksymax1 += 0.1*ksymax1;
ksMassTunable->SetAxisRange(0., ksymax1, "Y");
//ksMassTunable->GetYaxis()->SetLimits(0., ksymax1);

kspdgmass1 = new TLine(0.49767, 0., 0.49767, ksymax1);
kspdgmass1->SetLineColor(9);
kspdgmass1->SetLineWidth(1);
kspdgmass1->SetLineStyle(2);
//kspdgmass1->Draw("same");

//ksMassTunable->Draw("same");
mcmass->Draw("same");
kspdgmass1->Draw("same");
ksbgt2->Draw("same");
labl->Draw("same");
leg->Draw("same");
c1->Update();

c1->SaveAs("ksMass_superimposed.png");

//gStyle->SetOptStat("i");


ksCandMass->GetXaxis()->SetTitle("#pi#pi invariant mass (GeV)");
ksCandMass->GetYaxis()->SetTitle("Entries/7 MeV");

// Rerun fit with double gaussian...
if( refit == 1 ) {
  ksMassBinWidth = ksMassBinWidth/sqrt(2*3.141592654);
  ostringstream ksoss_doubgaus;

  ksoss_doubgaus << ksMassBinWidth 	
		 << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		 << ksMassBinWidth 
		 << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		 << "[5] + [6]*(x-[1]) + [7]*(x-[1])^2";
  cout << ksoss_doubgaus.str() << endl;
  TF1 *ksFit_doubGaus = new TF1("ksFit_doubGaus", ksoss_doubgaus.str().c_str(),
				0.3, 0.8);
  ksFit_doubGaus->SetParName( 0, "yield" );
  ksFit_doubGaus->SetParName( 1, "mean" );
  ksFit_doubGaus->SetParName( 2, "sigma1" );
  ksFit_doubGaus->SetParName( 3, "fraction" );
  ksFit_doubGaus->SetParLimits( 3, 0., 1. );
  ksFit_doubGaus->SetParName( 4, "sigma2" );
  ksFit_doubGaus->SetParName( 5, "const" );
  ksFit_doubGaus->SetParName( 6, "slope" );
  ksFit_doubGaus->SetParName( 7, "quadconst" );

  ksFit_doubGaus->SetParameter( 0, 4.4444*ksCandMass->GetMaximum() );
  ksFit_doubGaus->SetParameter( 1, 0.49767 );
  ksFit_doubGaus->SetParameter( 2, 0.0045 );
  ksFit_doubGaus->SetParameter( 3, 0.75 );
  ksFit_doubGaus->SetParameter( 4, 0.011 );
  ksFit_doubGaus->SetParameter( 5, 0.1*ksCandMass->GetMaximum() );
  ksFit_doubGaus->SetParameter( 6, 1. );
  ksFit_doubGaus->SetParameter( 7, -1. );

  ksCandMass->Fit("ksFit_doubGaus", "RLE");
}

ksxmin = ksCandMass->GetXaxis()->GetXmin();
ksxmax = ksCandMass->GetXaxis()->GetXmax();
ksxmax -= 0.0001;
ksymax = ksCandMass->GetMaximum();
ksymax += 0.1*ksymax;
ksCandMass->SetAxisRange(0., ksymax, "Y");
ksCandMass->GetXaxis()->SetLimits(ksxmin, ksxmax);
ksCandMass->GetXaxis()->SetNdivisions(508);
if( refit != 1 ) {
  ksCandMass->GetFunction("ksFit")->SetLineWidth(2);
  ksCandMass->GetFunction("ksFit")->SetLineColor(2);
}
if( refit == 1 ) {
  ksCandMass->GetFunction("ksFit_doubGaus")->SetLineWidth(2);
  ksCandMass->GetFunction("ksFit_doubGaus")->SetLineColor(2);
}
ksCandMass->Draw();

 
// TLine for PDG mass
kspdgmass = new TLine(0.49767, 0., 0.49767, ksymax);
kspdgmass->SetLineColor(9);
kspdgmass->SetLineWidth(2);
kspdgmass->SetLineStyle(2);
kspdgmass->Draw("same");


// TF1 for the background function
if( refit != 1 ) {
  ksconst = ksCandMass->GetFunction("ksFit")->GetParameter(3);
  ksslope = ksCandMass->GetFunction("ksFit")->GetParameter(4);
  ksquad = ksCandMass->GetFunction("ksFit")->GetParameter(5);
}
// For double gaussian:
if( refit == 1 ) {
  ksconst = ksCandMass->GetFunction("ksFit_doubGaus")->GetParameter(5);
  ksslope = ksCandMass->GetFunction("ksFit_doubGaus")->GetParameter(6);
  ksquad = ksCandMass->GetFunction("ksFit_doubGaus")->GetParameter(7);
}
Double_t range_min = 0.;
Double_t range_max = 0.;
if( refit != 1 ) {
  ksCandMass->GetFunction("ksFit")->GetRange(range_min, range_max);
}
if( refit == 1 ) {
  ksCandMass->GetFunction("ksFit_doubGaus")->GetRange(range_min, range_max);
}
TF1* ksbg;
//if( refit != 1 ) {
//  ksbg = new TF1("ksBg", "[0] + [1]*(x-0.49767)", range_min, range_max);
//}
//if( refit == 1 ) {
  ksbg = new TF1("ksBg", "[0] + [1]*(x-0.49767) + [2]*(x-0.49767)^2", range_min, range_max);
//}
ksbg->SetParName(0, "const");
ksbg->SetParName(1, "slope");
ksbg->SetParameter(0, ksconst);
ksbg->SetParameter(1, ksslope);
//if( refit == 1 ) {
  ksbg->SetParName(2, "quadconst");
  ksbg->SetParameter(2, ksquad);
//}
ksbg->SetLineColor(2);
ksbg->SetLineWidth(2);
ksbg->SetLineStyle(3);
ksbg->Draw("same");

//labl = new TPaveLabel(0.2, 0.8, 0.4, 0.9, mcordata.c_str(), "brNDC");

labl->Draw("same");

if( refit != 1 ) {
  c1->SaveAs("ksMass.png");
}
if( refit == 1 ) {
  c1->SaveAs("ksMass_refit.png");
}

//cin >> choice;

// Ks eta, plotting full and signal for now until better statistics
ksEtaFull->GetXaxis()->SetTitle("K^{0}_{S} #eta");
ksEtaFull->GetYaxis()->SetTitle("Entries");
ksEtaFull->Draw();

c1->SaveAs("ksEtaFull.png");

ksEtaSignal->GetXaxis()->SetTitle("K^{0}_{S} #eta");
ksEtaSignal->GetYaxis()->SetTitle("Entries");
ksEtaSignal->Draw();

c1->SaveAs("ksEtaSignal.png");

ksPhiFull->GetXaxis()->SetTitle("K^{0}_{S} #phi");
ksPhiFull->GetYaxis()->SetTitle("Entries");
ksPhiFull->Draw();

c1->SaveAs("ksPhiFull.png");

ksPhiSignal->GetXaxis()->SetTitle("K^{0}_{S} #phi");
ksPhiSignal->GetYaxis()->SetTitle("Entries");
ksPhiSignal->Draw();

c1->SaveAs("ksPhiSignal.png");

ksRFull->GetXaxis()->SetTitle("K^{0}_{S} decay radius R (cm)");
ksRFull->GetYaxis()->SetTitle("Entries");
ksRFull->Draw();

c1->SaveAs("ksRFull.png");

ksRSignal->GetXaxis()->SetTitle("K^{0}_{S} decay radius R (cm)");
ksRSignal->GetXaxis()->SetTitle("Entries");
ksRSignal->Draw();

c1->SaveAs("ksRSignal.png");

ksZvsRFull->GetXaxis()->SetTitle("K^{0}_{S} decay z position (cm)");
ksZvsRFull->GetYaxis()->SetTitle("K^{0}_{S} decay radius R (cm)");
ksZvsRFull->SetMarkerStyle(1);
ksZvsRFull->Draw();

c1->SaveAs("ksRvsZFull.png");

ksZvsRSignal->GetXaxis()->SetTitle("K^{0}_{S} decay z position (cm)");
ksZvsRSignal->GetYaxis()->SetTitle("K^{0}_{S} decay radius R (cm)");
ksZvsRSignal->SetMarkerStyle(1);
ksZvsRSignal->Draw();

c1->SaveAs("ksRvsZSignal.png");

//tdrStyle->SetPadLeftMargin(0.3);

//gStyle->SetOptStat(0);
KsMassBiasVsPhi->GetXaxis()->SetNdivisions(506, kFALSE);
KsMassBiasVsPhi->GetXaxis()->SetTitle("K^{0}_{S} #phi");
KsMassBiasVsPhi->GetYaxis()->SetTitle("Mass bias");
//KsMassBiasVsPhi->GetXaxis()->SetMaxDigits(2);
KsMassBiasVsPhi->GetYaxis()->SetTitleOffset(1.67);
pi_ = 3.141592654;
ksMBPhi = new TLine(-pi_, 0.49767, pi_, 0.49767);
ksMBPhi->SetLineStyle(2);
ksMBPhi->SetLineColor(9);
ksMBPhi->SetLineWidth(2);
KsMassBiasVsPhi->Draw();
ksMBPhi->Draw("same");
labl2_0->Draw("same");

c1->SetLeftMargin(0.2);

c1->SaveAs("ksMassBiasVsPhi.png");
//c1->SaveAs("ksMassBiasVsPhi.root");

KsMassBiasVsEta->GetXaxis()->SetNdivisions(506, kFALSE);
KsMassBiasVsEta->GetYaxis()->SetTitle("K^{0}_{S} #eta");
KsMassBiasVsEta->GetYaxis()->SetTitle("Mass bias");
KsMassBiasVsEta->GetYaxis()->SetTitleOffset(1.67);
//pi_ = 3.141592654;
ksMBEta = new TLine(0., 0.49767, 2.5, 0.49767);
ksMBEta->SetLineStyle(2);
ksMBEta->SetLineColor(9);
ksMBEta->SetLineWidth(2);
KsMassBiasVsEta->Draw();
ksMBEta->Draw("same");
labl2_0->Draw("same");

c1->SaveAs("ksMassBiasVsEta.png");

KsMassBiasVsPt->GetXaxis()->SetNdivisions(506, kFALSE);
KsMassBiasVsPt->GetXaxis()->SetTitle("K^{0}_{S} p_{T} (GeV/c)");
KsMassBiasVsPt->GetYaxis()->SetTitle("Mass bias");
KsMassBiasVsPt->GetYaxis()->SetTitleOffset(1.67);
//pi_ = 3.141592654;
ksMBPt = new TLine(0., 0.49767, 6., 0.49767);
ksMBPt->SetLineStyle(2);
ksMBPt->SetLineColor(9);
ksMBPt->SetLineWidth(2);
KsMassBiasVsPt->SetAxisRange(0.45, 0.52, "Y");
KsMassBiasVsPt->Draw();
ksMBPt->Draw("same");
labl2_0->Draw("same");

c1->SaveAs("ksMassBiasVsPt.png");
//c1->SaveAs("ksMassBiasVsPt.root");

c1->SetLeftMargin(0.16);

// ***************************************************************
// Lambda plots
fin->cd("analyzeLambda");

TGaxis::SetMaxDigits(3);

//gStyle->SetOptStat("i");

// Making plot from the ntuple


lamMassNbins = 100;
lamMassXmin = 1.08;
lamMassXmax = 1.18;
lamMassBinWidth = (lamMassXmax - lamMassXmin)/lamMassNbins;

ostringstream lambgoss;
lambgoss << "v0VtxSig > 15. && (abs(v0CandMass - 1.14) < 0.0075 || abs(v0CandMass - 1.095) < 0.0075) && abs(v0OtherCandMass - 0.498) > 0.015 && v0DauNormChi2 < 5.";
ostringstream lamsigoss;
lamsigoss << "v0VtxSig > 15. && abs(v0CandMass - 1.1158) < 0.0075 && abs(v0OtherCandMass - 0.498) > 0.015 && v0DauNormChi2 < 5.";
//ostringstream ksbgoss;
//ksbgoss << "v0VtxSig > 15. && (abs(v0CandMass - (0.498 - 0.06)) < 0.02 || abs(v0CandMass - (0.498 + 0.06)) < 0.02) && priVtxChi2 > 0. && abs(v0OtherCandMass - 1.1158) < 0.008";

lamMassWindow = new TH1F("lamMassWindow", "LAM", lamMassNbins, 
			 lamMassXmin, lamMassXmax);
lamMassWindow_ = new TH1F("lamMassWindow_", "LAM", lamMassNbins,
			  lamMassXmin, lamMassXmax);
lamMassWindowS = new TH1F("lamMassWindowS", "LAM", lamMassNbins,
			  lamMassXmin, lamMassXmax);
lamMassWindowBG = new TH1F("lamMassWindowBG", "LAM", lamMassNbins,
			   lamMassXmin, lamMassXmax);

ntuple->Project("lamMassWindow", "v0CandMass", "v0VtxSig > 15. && abs(v0OtherCandMass - 0.498) > 0.005 && v0DauNormChi2 < 5.");
ntuple->Project("lamMassWindow_", "v0CandMass", "v0VtxSig > 15.");
ntuple->Project("lamMassWindowS", "v0CandMass", lamsigoss.str().c_str());
ntuple->Project("lamMassWindowBG", "v0CandMass", lambgoss.str().c_str());

lamMassWindowS->SetFillColor(kBlue);
lamMassWindowBG->SetFillColor(kRed);

lamMassWindow->Draw();
lamMassWindowS->Draw("same");
lamMassWindowBG->Draw("same");
c1->SaveAs("LamMassWindows.png");

lamMassWindow_->SetLineColor(kBlue);
lamMassWindow_->Draw();
lamMassWindow->Draw("same");
c1->SaveAs("KsLamContaminationComparison.png");


// fit:
    ostringstream lamoss;
    lamoss << lamMassBinWidth << "*gausn(0) + [3]*(x - "
	   << piMass + protonMass << ")^(1/2) + [4]*(x - "
	   << piMass + protonMass << ")^(3/2)";
    TF1 *lamFitt = new TF1("lamFitt", lamoss.str().c_str(), piMass+protonMass, lamMassXmax);
    lamFitt->SetParName( 0, "yield" );
    lamFitt->SetParName( 1, "GausMean" );
    lamFitt->SetParName( 2, "GausSigma" );
    lamFitt->SetParName( 3, "c1" );
    lamFitt->SetParName( 4, "c2" );

    lamFitt->SetParameter( 0, 2.6*lamCandMass->GetMaximum() );
    lamFitt->SetParameter( 1, 1.1156 );
    lamFitt->SetParameter( 2, 0.002 );
    lamFitt->SetParameter( 3, 0.75*lamCandMass->GetMaximum() );
    lamFitt->SetParameter( 4, 1.85*lamCandMass->GetMaximum() );

TF1 *lamFitKsContam = new TF1("lamFitKsContam", lamoss.str().c_str(), piMass+protonMass, lamMassXmax);

lamFitKsContam->SetParName( 0, "yield");
lamFitKsContam->SetParName( 1, "GausMean");
lamFitKsContam->SetParName( 3, "c1");
lamFitKsContam->SetParName( 4, "c2");

lamFitKsContam->SetParameter( 0, 2.6*lamCandMass->GetMaximum() );
lamFitKsContam->SetParameter( 1, 1.1156 );
lamFitKsContam->SetParameter( 2, 0.002 );
lamFitKsContam->SetParameter( 3, 0.75*lamCandMass->GetMaximum() );
lamFitKsContam->SetParameter( 4, 1.85*lamCandMass->GetMaximum() );

lamMassTunable = new TH1F("lamMassTunable", "#Lambda^{0} mass",
			  lamMassNbins, lamMassXmin, lamMassXmax);
lamMassKsContam = new TH1F("lamMassKsContam", "#Lambda^{0} mass",
			   lamMassNbins, lamMassXmin, lamMassXmax);
ntuple->Project("lamMassTunable", "v0CandMass", 
		"v0VtxSig > 15. && priVtxChi2 > 0. && v0DauNormChi2 < 5.");
ntuple->Project("lamMassKsContam", "v0CandMass",
		"v0VtxSig > 15. && priVtxChi2 > 0. && abs(v0OtherCandMass - 0.498) < 0.015 && v0DauNormChi2 < 5.");

lamMassTunable->GetXaxis()->SetTitle("p#pi invariant mass (GeV/c^{2})");
lamMassTunable->GetYaxis()->SetTitle("Entries/.002 GeV/c^{2}");
lamMassTunable->Fit("lamFitt", "RLE");
lamMassKsContam->GetXaxis()->SetTitle("p#pi invariant mass (GeV/c^{2})");
lamMassKsContam->GetYaxis()->SetTitle("Entries/.002 GeV/c^{2}");
//lamMassKsContam->Fit("lamFitKsContam", "RLE");
lamxmint = lamMassTunable->GetXaxis()->GetXmin();
lamxmaxt = lamMassTunable->GetXaxis()->GetXmax();
lamxmaxt -= 0.00001;
lamymaxt = lamMassTunable->GetMaximum();
lamymaxt += 0.2*lamymaxt;
lamMassTunable->SetAxisRange(0., lamymaxt, "Y");
lamMassTunable->GetXaxis()->SetLimits(lamxmint, lamxmaxt);
lamMassTunable->GetXaxis()->SetNdivisions(506);
lamMassTunable->GetFunction("lamFitt")->SetLineWidth(2);
lamMassTunable->GetFunction("lamFitt")->SetLineColor(2);
lamMassTunable->Draw();

lampdgmasst = new TLine(1.1156, 0., 1.1156, lamymaxt);
lampdgmasst->SetLineColor(9);
lampdgmasst->SetLineWidth(2);
lampdgmasst->SetLineStyle(2);
lampdgmasst->Draw("same");

lamc1t = lamMassTunable->GetFunction("lamFitt")->GetParameter(3);
lamc2t = lamMassTunable->GetFunction("lamFitt")->GetParameter(4);
Double_t lamrange_mint = 0.;
Double_t lamrange_maxt = 0.;
lamMassTunable->GetFunction("lamFitt")->GetRange(lamrange_mint, lamrange_maxt);
ostringstream lamosst;
lamosst << "[0]*(x - " << 1.07784221 << ")^(1/2) + [1]*(x - "
        << 1.07784221 << ")^(3/2)";
TF1* lambgt = new TF1("lamBgt", lamosst.str().c_str(), 
		     lamrange_mint, lamrange_maxt);
lambgt->SetParName(0, "c1");
lambgt->SetParName(1, "c2");
lambgt->SetParameter(0, lamc1t);
lambgt->SetParameter(1, lamc2t);
lambgt->SetLineColor(2);
lambgt->SetLineWidth(2);
lambgt->SetLineStyle(3);
lambgt->Draw("same");

labl->Draw("same");

c1->SaveAs("lamMass_tuned.png");

lamMassKsContam->SetFillColor(kBlue);
lamMassKsContam->Draw("same");
//lampdgmasst->Draw("same");
//labl->Draw("same");
c1->SaveAs("lamMass_ksContamRemoved.png");

lamSignalRegion = new TH1F("lamSignalRegion", "lambda signal region",
			   lamMassNbins, lamMassXmin, lamMassXmax);
lamBackRegion = new TH1F("lamBackRegion", "lambda background region",
			 lamMassNbins, lamMassXmin, lamMassXmax);

ntuple->Project("lamSignalRegion", "v0CandMass",
		"v0VtxSig > 15. && priVtxChi2 > 0. && abs(v0CandMass - 1.1158) < 0.0075 && v0DauNormChi2 < 5.");
ntuple->Project("lamBackRegion", "v0CandMass",
		"v0VtxSig > 15. && priVtxChi2 > 0. && (abs(v0CandMass - 1.14) < 0.0075 || abs(v0CandMass - 1.095) < 0.0075) && v0DauNormChi2 < 5.");

lamMassTunable->Draw();
lampdgmasst->Draw("same");
labl->Draw("same");
lamSignalRegion->SetFillColor(kBlue - 1);
lamBackRegion->SetFillColor(kRed - 1);
lamSignalRegion->Draw("same");
lamBackRegion->Draw("same");
c1->SaveAs("LamCutRegions.png");

leg1 = new TLegend(0.7,0.5,0.92,0.69); 
leg1->AddEntry(lammcmass,"Monte Carlo","l"); 
leg1->AddEntry(lamMassTunable,"Data","l"); 
leg1->AddEntry(lamFitt, "Data Fit", "l");
leg1->SetFillColor(0);
Float_t rightmax1 = lammcmass->Integral();
Float_t scale1 = lamMassTunable->Integral()/rightmax1;

cout << "Integral(): " << rightmax1 << ", Integral(bin1, binN): "
<< lammcmass->Integral(1, ksMassNbins) << endl << endl;

lammcmass->SetLineColor(kCyan - 3);
//mcmass->Scale(scale*1.5); 
lammcmass->Scale(scale);
lammcmass->GetXaxis()->SetNdivisions(506);
lammcmass->Draw(); 
//lamMassTunable->Fit("ksFit_doubGaus0", "RLE");
lamMassTunable->Draw("same");
lampdgmasst->Draw("same");
lambgt->Draw("same");
labl->Draw("same");

c1->Update();

c1->SaveAs("lamMass_tempoverlay1.png");

//lammcmass->Draw(); 
lamMassTunable->Draw();
//lamymax1 = lammcmass->GetMaximum();
lamymax1 = lamMassTunable->GetMaximum();
lamymax1 += 0.1*lamymax1;
//lammcmass->SetAxisRange(0., lamymax1, "Y");
lamMassTunable->SetAxisRange(0., lamymax1, "Y");

lampdgmass1 = new TLine(0.49767, 0., 0.49767, lamymax1);
lampdgmass1->SetLineColor(9);
lampdgmass1->SetLineWidth(1);
lampdgmass1->SetLineStyle(2);
//kspdgmass1->Draw("same");

//lamMassTunable->Draw("same");
lammcmass->Draw("same");
lampdgmass1->Draw("same");
lambgt->Draw("same");
labl->Draw("same");
leg1->Draw("same");
c1->Update();

c1->SaveAs("lamMass_superimposed.png");

// From histo:
lamCandMass->GetXaxis()->SetTitle("p#pi invariant mass (GeV)");
lamCandMass->GetYaxis()->SetTitle("Entries/2 MeV");
lamxmin = lamCandMass->GetXaxis()->GetXmin();
lamxmax = lamCandMass->GetXaxis()->GetXmax();
lamxmax -= 0.00001;
lamymax = lamCandMass->GetMaximum();
lamymax += 0.1*lamymax;
lamCandMass->SetAxisRange(0., lamymax, "Y");
lamCandMass->GetXaxis()->SetLimits(lamxmin, lamxmax);
lamCandMass->GetXaxis()->SetNdivisions(506);
//ksCandMass->SetAxisRange(0., lamymaxt, "Y");
lamCandMass->GetFunction("lamFit")->SetLineWidth(2);
lamCandMass->GetFunction("lamFit")->SetLineColor(2);
lamCandMass->Draw();

//ksiourmen

lampdgmass = new TLine(1.1156, 0., 1.1156, lamymax);
lampdgmass->SetLineColor(9);
lampdgmass->SetLineWidth(2);
lampdgmass->SetLineStyle(2);
lampdgmass->Draw("same");

lamc1 = lamCandMass->GetFunction("lamFit")->GetParameter(3);
lamc2 = lamCandMass->GetFunction("lamFit")->GetParameter(4);
Double_t lamrange_min = 0.;
Double_t lamrange_max = 0.;
lamCandMass->GetFunction("lamFit")->GetRange(lamrange_min, lamrange_max);
ostringstream lamoss;
lamoss << "[0]*(x - " << 1.07784221 << ")^(1/2) + [1]*(x - "
       << 1.07784221 << ")^(3/2)";
TF1* lambg = new TF1("lamBg", lamoss.str().c_str(), 
		     lamrange_min, lamrange_max);
lambg->SetParName(0, "c1");
lambg->SetParName(1, "c2");
lambg->SetParameter(0, lamc1);
lambg->SetParameter(1, lamc2);
lambg->SetLineColor(2);
lambg->SetLineWidth(2);
lambg->SetLineStyle(3);
lambg->Draw("same");

labl->Draw("same");

c1->SaveAs("lamMass.png");

lamEtaFull->GetXaxis()->SetTitle("#Lambda^{0} #eta");
lamEtaFull->GetYaxis()->SetTitle("Entries");
lamEtaFull->Draw();

c1->SaveAs("lamEtaFull.png");

lamEtaSignal->GetXaxis()->SetTitle("#Lambda^{0} #eta");
lamEtaSignal->GetYaxis()->SetTitle("Entries");
lamEtaSignal->Draw();

c1->SaveAs("lamEtaSignal.png");

lamPhiFull->GetXaxis()->SetTitle("#Lambda^{0} #phi");
lamPhiFull->GetYaxis()->SetTitle("Entries");
lamPhiFull->Draw();

c1->SaveAs("lamPhiFull.png");

lamPhiSignal->GetXaxis()->SetTitle("#Lambda^{0} #phi");
lamPhiSignal->GetYaxis()->SetTitle("Entries");
lamPhiSignal->Draw();

c1->SaveAs("lamPhiSignal.png");

lamRFull->GetXaxis()->SetTitle("#Lambda^{0} decay radius R (cm)");
lamRFull->GetYaxis()->SetTitle("Entries");
lamRFull->Draw();

c1->SaveAs("lamRFull.png");

lamRSignal->GetXaxis()->SetTitle("#Lambda^{0} decay radius R (cm)");
lamRSignal->GetXaxis()->SetTitle("Entries");
lamRSignal->Draw();

c1->SaveAs("lamRSignal.png");

lamZvsRFull->GetXaxis()->SetTitle("#Lambda^{0} decay z position (cm)");
lamZvsRFull->GetYaxis()->SetTitle("#Lambda^{0} decay radius R (cm)");
lamZvsRFull->SetMarkerStyle(1);
lamZvsRFull->Draw();

c1->SaveAs("lamRvsZFull.png");

lamZvsRSignal->GetXaxis()->SetTitle("#Lambda^{0} decay z position (cm)");
lamZvsRSignal->GetYaxis()->SetTitle("#Lambda^{0} decay radius R (cm)");
lamZvsRSignal->SetMarkerStyle(1);
lamZvsRSignal->Draw();

c1->SaveAs("lamRvsZSignal.png");

//gStyle->SetOptStat(0);

LamMassBiasVsPhi->GetXaxis()->SetNdivisions(506, kFALSE);
LamMassBiasVsPhi->GetXaxis()->SetTitle("#Lambda^{0} #phi");
LamMassBiasVsPhi->GetYaxis()->SetTitle("Mass bias");
//LamMassBiasVsPhi->GetXaxis()->SetMaxDigits(2);
LamMassBiasVsPhi->GetYaxis()->SetTitleOffset(1.67);
pi_ = 3.141592654;
lamMBPhi = new TLine(-pi_, 0.49767, pi_, 0.49767);
lamMBPhi->SetLineStyle(2);
lamMBPhi->SetLineColor(9);
lamMBPhi->SetLineWidth(2);
LamMassBiasVsPhi->Draw();
lamMBPhi->Draw("same");
labl2_0->Draw("same");

c1->SetLeftMargin(0.2);

c1->SaveAs("lamMassBiasVsPhi.png");

LamMassBiasVsEta->GetXaxis()->SetNdivisions(506, kFALSE);
LamMassBiasVsEta->GetYaxis()->SetTitle("#Lambda^{0} #eta");
LamMassBiasVsEta->GetYaxis()->SetTitle("Mass bias");
LamMassBiasVsEta->GetYaxis()->SetTitleOffset(1.67);
//pi_ = 3.141592654;
lamMBEta = new TLine(0., 0.49767, 2.5, 0.49767);
lamMBEta->SetLineStyle(2);
lamMBEta->SetLineColor(9);
lamMBEta->SetLineWidth(2);
LamMassBiasVsEta->Draw();
lamMBEta->Draw("same");
labl2_0->Draw("same");

c1->SaveAs("lamMassBiasVsEta.png");

LamMassBiasVsPt->GetXaxis()->SetNdivisions(506, kFALSE);
LamMassBiasVsPt->GetXaxis()->SetTitle("#Lambda^{0} p_{T} (GeV/c)");
LamMassBiasVsPt->GetYaxis()->SetTitle("Mass bias");
LamMassBiasVsPt->GetYaxis()->SetTitleOffset(1.67);
//pi_ = 3.141592654;
lamMBPt = new TLine(0., 0.49767, 6., 0.49767);
lamMBPt->SetLineStyle(2);
lamMBPt->SetLineColor(9);
lamMBPt->SetLineWidth(2);
LamMassBiasVsPt->Draw();
lamMBPt->Draw("same");
labl2_0->Draw("same");

c1->SaveAs("lamMassBiasVsPt.png");

c1->SetLeftMargin(0.16);

fin->cd("analyzeKshort");

ksMassTunable2 = new TH1F("ksMassTunable2", "K^{0} mass",
			  ksMassNbins, ksMassXmin, ksMassXmax);
ksMassLamContam = new TH1F("ksMassLamContam", "K^{0} mass",
			   ksMassNbins, ksMassXmin, ksMassXmax);
ntuple->Project("ksMassTunable2", "v0CandMass", 
		"v0VtxSig > 15. && priVtxChi2 > 0. && v0DauNormChi2 < 5.");
ntuple->Project("ksMassLamContam", "v0CandMass",
		"v0VtxSig > 15. && priVtxChi2 > 0. && abs(v0OtherCandMass - 1.1158) < 0.008 && v0DauNormChi2 < 5.");

ksMassTunable2->GetXaxis()->SetTitle("#pi#pi invariant mass (GeV/c^{2})");
ksMassTunable2->GetYaxis()->SetTitle("Entries/? GeV/c^{2}");
ksMassLamContam->SetFillColor(kBlue);

ksMassTunable2->Fit("ksFitSG", "RLE");
ksMassTunable2->Draw();
ksMassLamContam->Draw("same");
kspdgmasst->Draw("same");
labl->Draw("same");

c1->SaveAs("ksMass_withContam.png");

ksMassLamContam->Draw();
c1->SaveAs("ksPureContamination.png");

ostringstream cutstring2;
cutstring2 << "(v0VtxSig > 15. && priVtxChi2 > 0. && abs(v0CandMass - 0.498) < 0.02)*1. && v0DauNormChi2 < 5.";
ostringstream bgcutstring2;
bgcutstring2 << "(v0VtxSig > 15. && priVtxChi2 > 0. && (abs(v0CandMass - (0.498-0.06)) < 0.02)|| abs(v0CandMass - (0.498+0.06)) < 0.02)*1. && v0DauNormChi2 < 5.";

invMass = new TH1F("ksMassFull", "Ks Mass", ksMassNbins, ksMassXmin, ksMassXmax);
invMassSignal = new TH1F("ksMassSignal", "Ks Mass", ksMassNbins, ksMassXmin, ksMassXmax);
invMassBG = new TH1F("ksMassBG", "Ks Mass", ksMassNbins, ksMassXmin, ksMassXmax);

ntuple->Project("ksMassFull", "v0CandMass", "v0VtxSig > 15. && priVtxChi2 > 0. && v0DauNormChi2 < 5.");
ntuple->Project("ksMassSignal", "v0CandMass", cutstring2.str().c_str());
ntuple->Project("ksMassBG", "v0CandMass", bgcutstring2.str().c_str());

invMassSignal->SetFillColor(kBlue - 1);
invMassBG->SetFillColor(kRed - 1);

invMass->GetXaxis()->SetTitle("#pi#pi invariant mass (GeV/c^{2})");
invMass->GetXaxis()->SetNdivisions(506);
invMass->GetYaxis()->SetTitle("Entries/.002 GeV/c^{2}");
invMass->Draw();
invMassSignal->Draw("same");
invMassBG->Draw("same");
c1->SaveAs("KsCutRegions.png");



/*
LamMassBiasVsPhi->GetXaxis()->SetNdivisions(506, kFALSE);
LamMassBiasVsPhi->GetXaxis()->SetTitle("#Lambda^{0} #phi");
LamMassBiasVsPhi->GetYaxis()->SetTitle("Mass bias");
lamphimin = LamMassBiasVsPhi->GetYaxis()->GetXmax();
cout << "lamphimin: " << lamphimin << endl;
pi_ = 3.141592654;
lamMBPhi = new TLine(-pi_, 1.1157, pi_, 1.1157);
lamMBPhi->SetLineStyle(2);
lamMBPhi->SetLineColor(9);
lamMBPhi->SetLineWidth(2);
LamMassBiasVsPhi->Draw();
lamMBPhi->Draw("same");

c1->SaveAs("lamMassBiasVsPhi.png");
*/

}
}
