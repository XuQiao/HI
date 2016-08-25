#include <sstream>
#include <vector>

void plot(){

  gROOT->SetBatch(kTRUE);

  bool do2D = true;
  bool do1DmassBias = true;
  bool drawHistos = true;
  double ksMass_eta_nBins = 20;
  double ksMass_phi_nBins = 12;
  double ksMass_eta_nBins_2D = 20;
  double ksMass_phi_nBins_2D = 6;

  double lamMass_eta_nBins = 10;
  double lamMass_phi_nBins = 6;

  Double_t piMass = 0.13957018;
  Double_t protonMass = 0.93827203;

  Double_t ksMassXmin = 0.351;
  Double_t ksMassXmax = 0.751;
  Double_t ksMassNbins = 200;
  Double_t ksMassBinWidth = (ksMassXmax-ksMassXmin)/ksMassNbins;
  Double_t ksMassBinWidth_DG = ksMassBinWidth/sqrt(2*TMath::Pi());

  Double_t lamMassXmin = 1.08;
  Double_t lamMassXmax = 1.18;
  Double_t lamMassNbins = 100;
  Double_t lamMassBinWidth = (lamMassXmax-lamMassXmin)/lamMassNbins;
  Double_t lamMassBinWidth_DG = lamMassBinWidth/sqrt(2*TMath::Pi());

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
//
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  gStyle->SetPadRightMargin(0.15);

  TFile* fin = TFile::Open("/nfs/data35/cms/v0ntuple/newmcout_wideCuts.root", "UPDATE");

  fin->cd("analyzeKshort");

  ostringstream ksoss_doubgaus;
  ksoss_doubgaus << ksMassBinWidth_DG 	
		 << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		 << ksMassBinWidth_DG 
		 << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		 << "[5] + [6]*(x-[1]) + [7]*(x-[1])^2";

  TF1 *ksFit_doubGaus = new TF1("ksFit_doubGaus", 
				ksoss_doubgaus.str().c_str(),
				ksMassXmin, ksMassXmax);
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
  ksFit_doubGaus->SetParameter( 2, 0.005 );
  ksFit_doubGaus->SetParameter( 3, 0.6 );
  ksFit_doubGaus->SetParameter( 4, 0.011 );
  ksFit_doubGaus->SetParameter( 5, 0.1*ksCandMass->GetMaximum() );
  ksFit_doubGaus->SetParameter( 6, 1. );
  ksFit_doubGaus->SetParameter( 7, -1. );

  ksCandMass->Fit("ksFit_doubGaus", "RLE");
  //ksCandMass->Draw();
  //ksCandMass->Write();

  fin->cd("analyzeLambda");

  ostringstream lamfitoss_doubGaus;
  lamfitoss_doubGaus << lamMassBinWidth_DG
		     << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		     << lamMassBinWidth_DG
		     << "*[0]*(1-[3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		     << "[5]*(x - " << piMass + protonMass << ")^(1/2) + [6]*(x-"
		     << piMass + protonMass << ")^(3/2)";
  TF1 *lamFit_doubGaus = new TF1("lamFit_doubGaus", lamfitoss_doubGaus.str().c_str(),
				 piMass+protonMass, lamMassXmax);
  lamFit_doubGaus->SetParName( 0, "yield" );
  lamFit_doubGaus->SetParName( 1, "mean" );
  lamFit_doubGaus->SetParName( 2, "sigma1" );
  lamFit_doubGaus->SetParName( 3, "fraction" );
  lamFit_doubGaus->SetParLimits( 3, 0., 1. );
  lamFit_doubGaus->SetParName( 4, "sigma2" );
  lamFit_doubGaus->SetParName( 5, "sqrt" );
  lamFit_doubGaus->SetParName( 6, "sqrtcube" );
  lamFit_doubGaus->SetParameter( 0, 8.*lamCandMass->GetMaximum() );
  lamFit_doubGaus->SetParameter( 1, 1.1159);
  lamFit_doubGaus->SetParameter( 2, 0.0016 );
  lamFit_doubGaus->SetParameter( 3, 0.44 );
  lamFit_doubGaus->SetParameter( 4, 0.0037);
  //lamFit_doubGaus->SetParameter( 5, 0.75*lamCandMass->GetMaximum() );
  //lamFit_doubGaus->SetParameter( 6, 1.85*lamCandMass->GetMaximum() );
  lamFit_doubGaus->SetParameter( 5, 0.001*lamCandMass->GetMaximum()/lamMassBinWidth );
  lamFit_doubGaus->SetParameter( 6, 0.001*lamCandMass->GetMaximum()/lamMassBinWidth );

  lamCandMass->Fit("lamFit_doubGaus", "RLE");
  //lamCandMass->Write();

  fin->cd("analyzeKshort");
  double ksMassConst = 0.497614;

  // Do mass bias additions and make a new histo to contain the yields
  //fin->cd();

  //gDirectory->ls();

  //std::vector<TH1F*> thePtrs;
  //thePtrs.push_back( (TH1F*) gROOT->FindObject("KsMassBiasVsEta"));

  //cout << "CREATED" << endl;
  //std::vector::iterator it = thePtrs.begin();
  //(*it)->Draw();
  //thePtrs[0]->Draw();


  if( do1DmassBias ) {
    TH1F* ksMassEta = (TH1F*) gROOT->FindObject("KsMassBiasVsEta");
    TH1F* ksMasspT = (TH1F*) gROOT->FindObject("KsMassBiasVsPt");
    TH1F* ksMassPhi = (TH1F*) gROOT->FindObject("KsMassBiasVsPhi");
    TH1F* ksMassBiasEtaVsPhi = (TH1F*) gROOT->FindObject("KsMassBiasEtaVsPhi");

    //ksMassEta->Draw();

    ostringstream ksoss;
    ksoss << ksMassBinWidth << "*gausn(0) + [3] + [4]*(x-0.49767)"
	  << " + [5]*(x-0.49767)^2";

    TF1 *ksFit = new TF1("ksFit", ksoss.str().c_str(), ksMassXmin, ksMassXmax);
    ksFit->SetParName( 0, "yield" );
    ksFit->SetParName( 1, "mean" );
    ksFit->SetParName( 2, "sigma" );
    ksFit->SetParName( 3, "const" );
    ksFit->SetParName( 4, "slope" );
    ksFit->SetParName( 5, "quadconst" );

    /*TH1F* etaHists[20];

    etaHists[0] = (TH1F*) gROOT->FindObject("ksMass_eta_-2.5");
    etaHists[1] = (TH1F*) gROOT->FindObject("ksMass_eta_-2.25");
    etaHists[2] = (TH1F*) gROOT->FindObject("ksMass_eta_-2");
    etaHists[3] = (TH1F*) gROOT->FindObject("ksMass_eta_-1.75");
    etaHists[4] = (TH1F*) gROOT->FindObject("ksMass_eta_-1.5");
    etaHists[5] = (TH1F*) gROOT->FindObject("ksMass_eta_-1.25");
    etaHists[6] = (TH1F*) gROOT->FindObject("ksMass_eta_-1");
    etaHists[7] = (TH1F*) gROOT->FindObject("ksMass_eta_-0.75");
    etaHists[8] = (TH1F*) gROOT->FindObject("ksMass_eta_-0.5");
    etaHists[9] = (TH1F*) gROOT->FindObject("ksMass_eta_-0.25");

    etaHists[10] = (TH1F*) gROOT->FindObject("ksMass_eta_0");
    etaHists[11] = (TH1F*) gROOT->FindObject("ksMass_eta_0.25");
    etaHists[12] = (TH1F*) gROOT->FindObject("ksMass_eta_0.5");
    etaHists[13] = (TH1F*) gROOT->FindObject("ksMass_eta_0.75");
    etaHists[14] = (TH1F*) gROOT->FindObject("ksMass_eta_1");
    etaHists[15] = (TH1F*) gROOT->FindObject("ksMass_eta_1.25");
    etaHists[16] = (TH1F*) gROOT->FindObject("ksMass_eta_1.5");
    etaHists[17] = (TH1F*) gROOT->FindObject("ksMass_eta_1.75");
    etaHists[18] = (TH1F*) gROOT->FindObject("ksMass_eta_2");
    etaHists[19] = (TH1F*) gROOT->FindObject("ksMass_eta_2.25");*/

    double ksMass_eta_binWidth = 5./ksMass_eta_nBins;
    unsigned int ksEtaCount = 1;
    //for( unsigned int ksEtaCount = 0; ksEtaCount < 20; ksEtaCount++ ) {
    for( double ksMass_eta_min = -2.5;
	 ksMass_eta_min < 2.5;
	 ksMass_eta_min += ksMass_eta_binWidth ) {
      fin->cd("analyzeKshort");
      ostringstream etaNameShort;
      etaNameShort << "ksMass_eta_" << ksMass_eta_min;
      TH1F* tmpHist = (TH1F*) gROOT->FindObject(etaNameShort.str().c_str());
      //etaHists[ksEtaCount];
      ksFit->SetParameter( 0, 4.4444*tmpHist->GetMaximum() );
      ksFit->SetParameter( 1, ksMassConst );
      ksFit->SetParameter( 2, 0.005 );
      ksFit->SetParameter( 3, 0.1*tmpHist->GetMaximum() );
      ksFit->SetParameter( 4, 1. );
      ksFit->SetParameter( 5, -1. );

      int fitStatus = tmpHist->Fit("ksFit", "RLE");//etaHists[ksEtaCount]->Fit("ksFit");
      if( !fitStatus ) {
	ksMassEta->SetBinContent( ksEtaCount, ksFit->GetParameter(1) - ksMassConst);
	ksMassEta->SetBinError( ksEtaCount, ksFit->GetParError(1) );
      }
      else {
	ksMassEta->SetBinContent( ksEtaCount, 0.);//ksMassConst );
	ksMassEta->SetBinError( ksEtaCount, 0. );
      }
      ksEtaCount++;
    }

    /*TH1F* phiHists[12];

    phiHists[0] = (TH1F*) gROOT->FindObject("ksMass_phi_-3.14159");
    phiHists[1] = (TH1F*) gROOT->FindObject("ksMass_phi_-2.61799");
    phiHists[2] = (TH1F*) gROOT->FindObject("ksMass_phi_-2.0944");
    phiHists[3] = (TH1F*) gROOT->FindObject("ksMass_phi_-1.5708");
    phiHists[4] = (TH1F*) gROOT->FindObject("ksMass_phi_-1.0472");
    phiHists[5] = (TH1F*) gROOT->FindObject("ksMass_phi_-0.523599");
    phiHists[6] = (TH1F*) gROOT->FindObject("ksMass_phi_-6.66134e-16");
    phiHists[7] = (TH1F*) gROOT->FindObject("ksMass_phi_0.523599");
    phiHists[8] = (TH1F*) gROOT->FindObject("ksMass_phi_1.0472");
    phiHists[9] = (TH1F*) gROOT->FindObject("ksMass_phi_1.5708");
    phiHists[10] = (TH1F*) gROOT->FindObject("ksMass_phi_2.0944");
    phiHists[11] = (TH1F*) gROOT->FindObject("ksMass_phi_2.61799");*/

    //for( unsigned int ndx = 0; ndx < 12; ndx++ ) {
    unsigned int phi_ndx = 1;
    for( double ksMass_phi_min = -TMath::Pi();
	 ksMass_phi_min < TMath::Pi() - 0.0001;
	 ksMass_phi_min += 2*TMath::Pi()/ksMass_phi_nBins) {
      fin->cd("analyzeKshort");
      ostringstream oss;
      oss << "ksMass_phi_" << ksMass_phi_min;
      //TH1F* tmpHist = phiHists[ndx];
      TH1F* tmpHist = (TH1F*) gROOT->FindObject(oss.str().c_str());
      ksFit->SetParameter( 0, 4.4444*tmpHist->GetMaximum() );
      ksFit->SetParameter( 1, ksMassConst );
      ksFit->SetParameter( 2, 0.005 );
      ksFit->SetParameter( 3, 0.1*tmpHist->GetMaximum() );
      ksFit->SetParameter( 4, 1. );
      ksFit->SetParameter( 5, -1. );

      int fitStatus = tmpHist->Fit("ksFit", "RLE");
      //cout << "Phi " << ndx << " fit status: " << fitStatus << endl;
      if( !fitStatus ) {
	//cout << "Setting bin content..." << endl;
	ksMassPhi->SetBinContent( phi_ndx, ksFit->GetParameter(1) - ksMassConst);
	ksMassPhi->SetBinError( phi_ndx, ksFit->GetParError(1) );
      }
      else {
	//cout << "Setting a ZERO bin content..." << endl;
	ksMassPhi->SetBinContent( phi_ndx, 0.);//ksMassConst );
	ksMassPhi->SetBinError( phi_ndx, 0. );
      }
      phi_ndx++;
    }

    /*TH2F* ksEtaVsPhi = new TH2F("ksEtaVsPhi", "K^{0}_{S} mass bias, #eta vs #phi",
      12, -TMath::Pi(), TMath::Pi(), 
      10, 0., 2.5);
      for( unsigned int etandx = 0; etandx < 10; etandx++ ) {
      for( unsigned int phindx = 0; phindx < 12; phindx++ ) {*/


    TH1F* ptHists[13];

    ptHists[0] = (TH1F*) gROOT->FindObject("ksMass_pt_0");
    ptHists[1] = (TH1F*) gROOT->FindObject("ksMass_pt_0.125");
    ptHists[2] = (TH1F*) gROOT->FindObject("ksMass_pt_0.25");
    ptHists[3] = (TH1F*) gROOT->FindObject("ksMass_pt_0.375");
    ptHists[4] = (TH1F*) gROOT->FindObject("ksMass_pt_0.5");
    ptHists[5] = (TH1F*) gROOT->FindObject("ksMass_pt_0.75");
    ptHists[6] = (TH1F*) gROOT->FindObject("ksMass_pt_1");
    ptHists[7] = (TH1F*) gROOT->FindObject("ksMass_pt_1.25");
    ptHists[8] = (TH1F*) gROOT->FindObject("ksMass_pt_1.5");
    ptHists[9] = (TH1F*) gROOT->FindObject("ksMass_pt_1.75");
    ptHists[10] = (TH1F*) gROOT->FindObject("ksMass_pt_2");
    ptHists[11] = (TH1F*) gROOT->FindObject("ksMass_pt_2.5");
    ptHists[12] = (TH1F*) gROOT->FindObject("ksMass_pt_3");

    for( unsigned int ndx = 0; ndx < 13; ndx++ ) {
      TH1F* tmpHist = ptHists[ndx];
      ksFit->SetParameter( 0, 4.4444*tmpHist->GetMaximum() );
      ksFit->SetParameter( 1, ksMassConst );
      ksFit->SetParameter( 2, 0.005 );
      ksFit->SetParameter( 3, 0.1*tmpHist->GetMaximum() );
      ksFit->SetParameter( 4, 1. );
      ksFit->SetParameter( 5, -1. );

      int fitStatus = tmpHist->Fit("ksFit", "RLE");
      if( !fitStatus ) {
	ksMasspT->SetBinContent( ndx + 1, ksFit->GetParameter(1) - ksMassConst);
	ksMasspT->SetBinError( ndx + 1, ksFit->GetParError(1) );
      }
      else {
	ksMasspT->SetBinContent( ndx + 1, 0.);//ksMassConst );
	ksMasspT->SetBinError( ndx + 1, 0. );
      }
    }

    fin->cd("analyzeLambda");
    double lamMassConst = 1.115683;
    double protonMass = 0.93827203;
    double piMass = 0.13957018;

    TH1F* lamMassEta = (TH1F*) gROOT->FindObject("LamMassBiasVsEta");
    TH1F* lamMasspT = (TH1F*) gROOT->FindObject("LamMassBiasVsPt");
    TH1F* lamMassPhi = (TH1F*) gROOT->FindObject("LamMassBiasVsPhi");

    cout << "eta, bin1: " << lamMassEta->GetBinContent(1) << endl;

    ostringstream lamoss;
    lamoss << lamMassBinWidth << "*gausn(0) + [3]*(x - "
	   << piMass + protonMass << ")^(1/2) + [4]*(x - "
	   << piMass + protonMass << ")^(3/2)";

    TF1 *lamFit = new TF1("lamFit", lamoss.str().c_str(), piMass+protonMass, lamMassXmax);
    lamFit->SetParName( 0, "yield" );
    lamFit->SetParName( 1, "GausMean" );
    lamFit->SetParName( 2, "GausSigma" );
    lamFit->SetParName( 3, "c1" );
    lamFit->SetParName( 4, "c2" );

    /*TH1F* lamEtaHists[10];
      lamEtaHists[0] = (TH1F*) gROOT->FindObject("lamMass_eta_-2.5");
      lamEtaHists[1] = (TH1F*) gROOT->FindObject("lamMass_eta_-2");
      lamEtaHists[2] = (TH1F*) gROOT->FindObject("lamMass_eta_-2.5");
      lamEtaHists[3] = (TH1F*) gROOT->FindObject("lamMass_eta_-1");
      lamEtaHists[4] = (TH1F*) gROOT->FindObject("lamMass_eta_-0.5");

      lamEtaHists[5] = (TH1F*) gROOT->FindObject("lamMass_eta_0");
      lamEtaHists[6] = (TH1F*) gROOT->FindObject("lamMass_eta_0.5");
      lamEtaHists[7] = (TH1F*) gROOT->FindObject("lamMass_eta_1");
      lamEtaHists[8] = (TH1F*) gROOT->FindObject("lamMass_eta_1.5");
      lamEtaHists[9] = (TH1F*) gROOT->FindObject("lamMass_eta_2");*/

    double lamMass_eta_binWidth = 5./lamMass_eta_nBins;
    unsigned int lamEtaCount = 1;
    //for( unsigned int ndx = 0; ndx < 10; ndx++ ) {
    for( double lamMass_eta_min = -2.5;
	 lamMass_eta_min < 2.5;
	 lamMass_eta_min += lamMass_eta_binWidth ) {
      fin->cd("analyzeLambda");
      ostringstream lamEtaNameShort;
      lamEtaNameShort << "lamMass_eta_" << lamMass_eta_min;
      //TH1F* tmpHist = lamEtaHists[ndx];
      TH1F* tmpHist = (TH1F*) gROOT->FindObject(lamEtaNameShort.str().c_str());
      lamFit->SetParameter( 0, 2.6*tmpHist->GetMaximum() );
      lamFit->SetParameter( 1, lamMassConst );
      lamFit->SetParameter( 2, 0.002 );
      lamFit->SetParameter( 3, 0.75*tmpHist->GetMaximum() );
      lamFit->SetParameter( 4, 1.85*tmpHist->GetMaximum() );

      int fitStatus = tmpHist->Fit("lamFit", "RLE");
      if( !fitStatus ) {
	lamMassEta->SetBinContent( lamEtaCount, lamFit->GetParameter(1) - lamMassConst );
	lamMassEta->SetBinError( lamEtaCount, lamFit->GetParError(1) );
      }
      else {
	lamMassEta->SetBinContent( lamEtaCount, 0.);//lamMassConst );
	lamMassEta->SetBinError( lamEtaCount, 0. );
      }
      lamEtaCount++;
    }

    /*TH1F* lamPhiHists[6];
      lamPhiHists[0] = (TH1F*) gROOT->FindObject("lamMass_phi_-3.14159");
      lamPhiHists[1] = (TH1F*) gROOT->FindObject("lamMass_phi_-2.0944");
      lamPhiHists[2] = (TH1F*) gROOT->FindObject("lamMass_phi_-1.0472");
      lamPhiHists[3] = (TH1F*) gROOT->FindObject("lamMass_phi_-4.44089e-16");
      lamPhiHists[4] = (TH1F*) gROOT->FindObject("lamMass_phi_1.0472");
      lamPhiHists[5] = (TH1F*) gROOT->FindObject("lamMass_phi_2.0944");*/
    double lamMass_phi_binWidth = 2*TMath::Pi()/lamMass_phi_nBins;
    unsigned int lamPhiCount = 1;
    //for( unsigned int ndx = 0; ndx < 6; ndx++ ) {
    for( double lamMass_phi_min = -TMath::Pi();
	 lamMass_phi_min < TMath::Pi() - 0.0001;
	 lamMass_phi_min += lamMass_phi_binWidth ) {
      fin->cd("analyzeLambda");
      ostringstream lamPhiNameShort;
      lamPhiNameShort << "lamMass_phi_" << lamMass_phi_min;
      //TH1F* tmpHist = lamPhiHists[ndx];
      TH1F* tmpHist = (TH1F*) gROOT->FindObject(lamPhiNameShort.str().c_str());
      lamFit->SetParameter( 0, 2.6*tmpHist->GetMaximum() );
      lamFit->SetParameter( 1, lamMassConst );
      lamFit->SetParameter( 2, 0.002 );
      lamFit->SetParameter( 3, 0.75*tmpHist->GetMaximum() );
      lamFit->SetParameter( 4, 1.85*tmpHist->GetMaximum() );

      int fitStatus = tmpHist->Fit("lamFit", "RLE");
      if( !fitStatus ) {
	lamMassPhi->SetBinContent( lamPhiCount, lamFit->GetParameter(1) - lamMassConst);
	lamMassPhi->SetBinError( lamPhiCount, lamFit->GetParError(1) );
      }
      else {
	lamMassPhi->SetBinContent( lamPhiCount, 0.);//lamMassConst );
	lamMassPhi->SetBinError( lamPhiCount, 0. );
      }
      lamPhiCount++;
    }

    TH1F* lamPtHists[7];
    lamPtHists[0] = (TH1F*) gROOT->FindObject("lamMass_pt_0");
    lamPtHists[1] = (TH1F*) gROOT->FindObject("lamMass_pt_0.25");
    lamPtHists[2] = (TH1F*) gROOT->FindObject("lamMass_pt_0.5");
    lamPtHists[3] = (TH1F*) gROOT->FindObject("lamMass_pt_1");
    lamPtHists[4] = (TH1F*) gROOT->FindObject("lamMass_pt_1.5");
    lamPtHists[5] = (TH1F*) gROOT->FindObject("lamMass_pt_2");
    lamPtHists[6] = (TH1F*) gROOT->FindObject("lamMass_pt_3");
    for( unsigned int ndx = 0; ndx < 7; ndx++ ) {
      TH1F* tmpHist = lamPtHists[ndx];
      lamFit->SetParameter( 0, 2.6*tmpHist->GetMaximum() );
      lamFit->SetParameter( 1, lamMassConst );
      lamFit->SetParameter( 2, 0.002 );
      lamFit->SetParameter( 3, 0.75*tmpHist->GetMaximum() );
      lamFit->SetParameter( 4, 1.85*tmpHist->GetMaximum() );

      int fitStatus = tmpHist->Fit("lamFit", "RLE");
      if( !fitStatus ) {
	lamMasspT->SetBinContent( ndx + 1, lamFit->GetParameter(1) - lamMassConst );
	lamMasspT->SetBinError( ndx + 1, lamFit->GetParError(1) );
      }
      else {
	lamMasspT->SetBinContent( ndx + 1, 0.);//lamMassConst );
	lamMasspT->SetBinError( ndx + 1, 0. );
      }
    }
  }

  double ksMass_eta_binWidth = 5. / ksMass_eta_nBins_2D;
  double ksMass_phi_binWidth = 2*TMath::Pi() / ksMass_phi_nBins_2D;

  if(do2D) {
    unsigned int etaCount = 1;
    unsigned int phiCount = 1;
    std::vector<string> theNames;
    std::vector<double> theValues;
    std::vector<int> etaCounts;
    std::vector<int> phiCounts;
    for( double ksMass_eta_min = -2.5;
	 ksMass_eta_min < 2.5;
	 ksMass_eta_min += ksMass_eta_binWidth ) {
      for( double ksMass_phi_min = -TMath::Pi();
	   ksMass_phi_min < TMath::Pi() - 0.0001;
	   ksMass_phi_min += ksMass_phi_binWidth ) {
	fin->cd("analyzeKshort");
	ostringstream shortNames2D;
	//ostringstream longNames2D;
	shortNames2D << "ksMass_2D_" << ksMass_eta_min << "_" << ksMass_phi_min;
	/*longNames2D << "Ks invariant mass, " << ksMass_eta_min << " < #eta < "
	  << ksMass_eta_min + ksMass_eta_binWidth << ", "
	  << ksMass_phi_min << " < #phi < " << ksMass_phi_min + ksMass_phi_binWidth;*/
	cout << shortNames2D.str() << endl;
	TH1F* tmpHist = (TH1F*) gROOT->FindObject(shortNames2D.str().c_str());
	ksFit->SetParameter( 0, 4.4444*tmpHist->GetMaximum() );
	ksFit->SetParameter( 1, ksMassConst );
	ksFit->SetParameter( 2, 0.005 );
	ksFit->SetParameter( 3, 0.1*tmpHist->GetMaximum() );
	ksFit->SetParameter( 4, 1. );
	ksFit->SetParameter( 5, -1. );
	int fitstatus = tmpHist->Fit("ksFit", "RLE");
	//cout << "Fit status: " << fitstatus << endl;
	if( !fitstatus ) {
	  cout << "Mass bias: " << ksFit->GetParameter(1) - ksMassConst
	       << ", Error: " << ksFit->GetParError(1) << endl;
	  ksMassBiasEtaVsPhi->SetBinContent( phiCount, etaCount, ksFit->GetParameter(1) - ksMassConst );
	  ksMassBiasEtaVsPhi->SetBinError( phiCount, etaCount, ksFit->GetParError(1) );
	  theNames.push_back( shortNames2D.str() );
	  theValues.push_back( ksFit->GetParameter(1) - ksMassConst );
	  etaCounts.push_back( etaCount );
	  phiCounts.push_back( phiCount );
	}
	else {
	  cout << "Filling a bad one..." << endl;
	  ksMassBiasEtaVsPhi->SetBinContent( phiCount, etaCount, 0. );
	  ksMassBiasEtaVsPhi->SetBinError( phiCount, etaCount, 0. );
	}
	phiCount++;	  
      }
      phiCount = 1;
      etaCount++;
    }
    /*for( unsigned int kk = 0; kk < theNames.size(); kk++ ) {
      cout << theNames[kk] << ": " << theValues[kk] 
	   << ", eta: " << etaCounts[kk] << ", phi: " << phiCounts[kk] << endl;
	   }*/
    if( drawHistos ) {
      ksMassBiasEtaVsPhi->SetXTitle("#phi of lowest p_{T} daughter track [radians]");
      ksMassBiasEtaVsPhi->SetYTitle("#eta of lowest p_{T} daughter track");
      ksMassBiasEtaVsPhi->Draw("colz");
      //gStyle->SetPadRightMargin(-3.);
      gStyle->SetOptStat(0);
      TPaletteAxis* thePalette = (TPaletteAxis*) c1->FindObject("KsMassBiasEtaVsPhi")->FindObject("palette");
      thePalette->SetX1NDC(0.855);
      thePalette->SetX2NDC(0.9);
      c1->SaveAs("KsMassBiasEtaVsPhi.png");
    }
  }

  fin->Write();
}
