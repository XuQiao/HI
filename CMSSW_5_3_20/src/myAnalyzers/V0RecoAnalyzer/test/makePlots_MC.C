{
#include <tdrstyle.C>
#include <sstream>
//#include <iostream>
#include <string>



gROOT->Reset();
setTDRStyle();
gStyle->SetLabelOffset(0.017, "X");
gStyle->SetPadTopMargin(0.06);
gStyle->SetPadBottomMargin(0.12);
gStyle->SetStatFormat("5.5g");
gStyle->SetStatY(0.94);
gStyle->SetFitFormat("5.5g");
//gStyle->SetStatFormat("NELU",AutoPrecision(3));

TFile* fin = TFile::Open("v0analysis_MC.root", "READ");
//TFile* lifein = TFile::Open("ksLifetime_MC.root", "READ");
TFile* fout = TFile::Open("ksLifetime_MC2.root", "RECREATE");
TFile* fmassout = TFile::Open("ksMass_histo.root", "RECREATE");
fin->cd("ksAnalyzer");

int choice;
/*
cout << "For Monte Carlo, press 1.  For data, press 2.  Then press ENTER." 
     << endl;
cin >> choice;

cout << "Refit the K0S mass plot with a double gaussian?  1 for yes." << endl;
int refit = 0;
cin >> refit;
*/
choice = 1;

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

labl = new TPaveText(0.19, 0.87, 0.4, 0.92, "brNDC");
labl->SetBorderSize(0);
labl->SetFillColor(0);
if( choice == 1 ) {
  labl->SetTextColor(13);
}
TText* labtxt = labl->AddText(mcordata.c_str());

labl2_0 = new TPaveText(0.23, 0.87, 0.8, 0.92, "brNDC");
labl2_0->SetBorderSize(0);
labl2_0->SetFillColor(0);
if( choice == 1 ) {
  labl2_0->SetTextColor(13);
}
TText* labtxt2_0 = labl2_0->AddText(mcordata.c_str());

ksLifetimeNbins = 20;
ksLifetimeXmin = 0.;
ksLifetimeXmax = 10.;

gStyle->SetOptStat(0);
ksLifetimeLoc = new TH1F("ksCtau", "K^{0}_{S} c#tau",
		      ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);
ksLifetimeLoc->Sumw2();
ntuple->Project("ksCtau", "v0Lifetime", 
		"v0Lifetime > 0. && v0pT < 1.5 && abs(v0CandMass - 0.498) < 0.01");
ksLifetimeBG = new TH1F("ksCtauBG", "K^{0}_{S} c#tau background",
			ksLifetimeNbins, ksLifetimeXmin, ksLifetimeXmax);
ksLifetimeBG->Sumw2();
ntuple->Project("ksCtauBG", "v0Lifetime",
		"v0Lifetime > 0. && v0pT < 1.5 && abs(v0CandMass - (0.498 - 0.04)) < 0.01 || abs(v0CandMass - (0.498 + 0.04)) < 0.01");
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
ksCtauFitLoc = new TF1("ksCtauFit", ksctauoss.str().c_str(), 
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

ksLifetimeLoc->Fit("ksCtauFit", "RE");

Double_t realCtau = 2.6786;
//Double_t realCtau_p1 = -1./realCtau;
Double_t realCtau_p1 = realCtau/spdolight;
ksCtauRealLoc = new TF1("ksCtauReal", ksctauoss.str().c_str(),
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
fin->cd("ksAnalyzer");
/*
ksLifetimeLoc->Draw();
labl->Draw("same");
//ksCtauReal->Draw("same");

c1->cd(1);
gPad->SetLogy(1);

c1->SaveAs("ksCtau.png");
c1->SaveAs("ksCtau.ps");

lifein->cd();
//ksCtauFit->Draw();
TF1* ksLTFit = new TF1(*ksCtauFit);
ksLTFit->Draw();
c1->SaveAs("ksLTFit.png");
TF1* ksLTReal = new TF1(*ksCtauReal);
//fin->cd("ksAnalyzer");
TH1F* ksLTHist = new TH1F(*ksCtau);
ksLTHist->Draw();
c1->SaveAs("ltHistTest.png");
ksLTReal->Draw();
c1->SaveAs("ltFunc.png");
ksLTHist->Divide(ksLTReal);
ksLTHist->GetYaxis()->SetTitle("Efficiency");
ksLTHist->Draw();
labl->Draw("same");
gPad->SetLogy(0);
c1->SaveAs("ltDivided.png");
gPad->SetLogy(0);
fin->cd("ksAnalyzer");
ksLifetimeLoc->Divide(ksLTHist);
//ksLifetime->Multiply(ksLTReal);
ksLifetimeLoc->Fit("ksCtauFitLoc", "RE");
ksLifetimeLoc->Draw();
labl->Draw("same");
//ksLTReal->Draw("same");
gPad->SetLogy(1);

c1->SaveAs("ksCtau_scaled.png");

gPad->SetLogy(0);
c1->SaveAs("ksCtau_scaled_noLogScale.png");
*/
gStyle->SetOptStat("i");

ksMassNbins = 200;
ksMassXmin = 0.351;
ksMassXmax = 0.751;
ksMassBinWidth = (ksMassXmax-ksMassXmin)/ksMassNbins;
// Fit function:


// MASS:

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
ntuple->Project("ksMassTunable", "v0CandMass", "");
//ntuple->Project("ksMassTunable2", "v0CandMass", "");

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

fmassout->cd();
TH1F* ksmouthist = new TH1F(*ksMassTunable);
ksmouthist->Write();
fmassout->Write();
fin->cd("ksAnalyzer");
/*

//ksMassTunable->

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
ksMassTunable2->GetXaxis()->SetNdivisions(506);
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

gStyle->SetOptStat(0);
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
c1->SaveAs("ksMassBiasVsPhi.root");

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
c1->SaveAs("ksMassBiasVsPt.root");

c1->SetLeftMargin(0.16);
*/

// ***************************************************************
// Lambda plots

fin->cd("lamAnalyzer");

TGaxis::SetMaxDigits(3);

gStyle->SetOptStat("i");

// Making plot from the ntuple


lamMassNbins = 100;
lamMassXmin = 1.08;
lamMassXmax = 1.18;
lamMassBinWidth = (lamMassXmax - lamMassXmin)/lamMassNbins;

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

lamMassTunable = new TH1F("lamMassTunable", "#Lambda^{0} mass",
			  lamMassNbins, lamMassXmin, lamMassXmax);
ntuple->Project("lamMassTunable", "v0CandMass", "");

lamMassTunable->GetXaxis()->SetTitle("p#pi invariant mass (GeV/c^{2})");
lamMassTunable->GetYaxis()->SetTitle("Entries/.002 GeV/c^{2}");
lamMassTunable->Fit("lamFitt", "RLE");
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

fmassout->cd();
TH1F* lammouthist = new TH1F(*lamMassTunable);
lammouthist->Write();
fmassout->Write();
fin->cd("lamAnalyzer");
/*

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

gStyle->SetOptStat(0);

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
*/

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
