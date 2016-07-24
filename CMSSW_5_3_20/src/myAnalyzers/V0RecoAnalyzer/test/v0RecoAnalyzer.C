void makeVeePlots() {

  M_PI = 3.141592653589793238;
  piMass = 0.13957018;
  protonMass = 0.93827203;
  
  TFile infile("v0analysis_wideMassCuts.root");

  // Get Kshort tree and make histograms
  infile.cd("ksAnalyzer");
  // Trajectory variables
  ksR = new TH1F("ksR", "K^{0}_{s} radial distance from beam line",
		 40, 0., 40.);
  ksEta = new TH1F("ksEta", "K^{0}_{s} momentum #eta",
		   40, -2.5, 2.5);
  ksPhi = new TH1F("ksPhi", "K^{0}_{s} momentum #phi",
		   60, -M_PI, M_PI);
  ksPt = new TH1F("ksPt", "K^{0}_{s} p_{T}",
		  70, 0., 20.);

  ksRFull = new TH1F("ksRFull", "K^{0}_{s} radial distance from beam line",
		 40, 0., 40.);
  ksEtaFull = new TH1F("ksEtaFull", "K^{0}_{s} momentum #eta",
		   40, -2.5, 2.5);
  ksPhiFull = new TH1F("ksPhiFull", "K^{0}_{s} momentum #phi",
		   60, -M_PI, M_PI);
  ksPtFull = new TH1F("ksPtFull", "K^{0}_{s} p_{T}",
		  70, 0., 20.);

  ksRSignal = new TH1F("ksRSignal", "K^{0}_{s} radial distance from beam line",
		 40, 0., 40.);
  ksEtaSignal = new TH1F("ksEtaSignal", "K^{0}_{s} momentum #eta",
		   40, -2.5, 2.5);
  ksPhiSignal = new TH1F("ksPhiSignal", "K^{0}_{s} momentum #phi",
		   60, -M_PI, M_PI);
  ksPtSignal = new TH1F("ksPtSignal", "K^{0}_{s} p_{T}",
		  70, 0., 20.);

  ksRBg = new TH1F("ksRBg", "K^{0}_{s} radial distance from beam line",
		 40, 0., 40.);
  ksEtaBg = new TH1F("ksEtaBg", "K^{0}_{s} momentum #eta",
		   40, -2.5, 2.5);
  ksPhiBg = new TH1F("ksPhiBg", "K^{0}_{s} momentum #phi",
		   60, -M_PI, M_PI);
  ksPtBg = new TH1F("ksPtBg", "K^{0}_{s} p_{T}",
		  70, 0., 20.);

  ksZvsRFull = new TH2F("ksZvsRFull", "K^{0}_{s} z position vs. decay radius",
			120, -120., 120., 40, 0., 40.);
  ksZvsRSignal = new TH2F("ksZvsRSignal", "K^{0}_{s} z position vs. decay radius",
			120, -120., 120., 40, 0., 40.);
  ksZvsRBg = new TH2F("ksZvsRBg", "K^{0}_{s} z position vs. decay radius",
			120, -120., 120., 40, 0., 40.);
  ksZvsRBgSubt = new TH2F("ksZvsRBgSubt", "K^{0}_{s} z position vs. decay radius",
  		      120, -120., 120., 40, 0., 40.);

  ksDaup = new TH1F("ksDaup", "p of K^{0}_{s} daughter tracks",
		    100, 0., 30.);
  ksDaupT = new TH1F("ksDaupT", "p_{T} of K^{0}_{s} daughter tracks",
		     100, 0., 20.);
  ksDauEta = new TH1F("ksDauEta", "#eta of K^{0}_{s} daughter tracks",
		      100, -2.5, 2.5);
  ksDauPhi = new TH1F("ksDauPhi", "#phi of K^{0}_{s} daughter track at vertex position",
		      70, -M_PI, M_PI);
  ksDauNHits = new TH1F("ksDauNHits", "Total number of hits on K^{0}_{s} daughter tracks",
			80, 0., 40.);
  ksDauNValidHits = new TH1F("ksDauNValHits", "Number of valid hits on K^{0}_{s} daughter tracks", 
			     80, 0., 40.);
  ksDauNPixelHits = new TH1F("ksDauNPixelHits", "Number of valid pixel hits on K^{0}_{s} daughter tracks",
			     30 , 0., 15.);
  ksDauNStripHits = new TH1F("ksDauNStripHits", "Number of valid strip hits on K^{0}_{s} daughter tracks",
			     80, 0., 40.);
  ksDauChi2 = new TH1F("ksDauChi2", "#chi^{2} of K^{0}_{s} daughter tracks",
		       100, 0., 100.);
  ksDauNormChi2 = new TH1F("ksDauNormChi2", "#chi^{2}/ndof of K^{0}_{s} daughter tracks",
			   100, 0., 5.);



  //ksm = 0.49767;
  // Mass plots
  double ksXmin = 0.44;
  double ksXmax = 0.56;
  int ksNbins = 60;

  double ksNormXmin = -10.;
  double ksNormXmax = 10.;
  int ksNormNbins = 60;
  ksMass = new TH1F("ksMass", "K^{0}_{s} invariant mass",
		    ksNbins, ksXmin, ksXmax);
  double ksBinWidth = (ksXmax - ksXmin) / (double) ksNbins;
  cout << "ksBinWidth: " << ksBinWidth << endl;
  ksNormMass = new TH1F("ksNormMass", "K^{0}_{s} normalized mass",
			ksNormNbins, ksNormXmin, ksNormXmax);
  double ksNormBinWidth = (ksNormXmax - ksNormXmin) / (double) ksNormNbins;
  ksCandMass = new TH1F("ksCandMass", "K^{0}_{s} mass from Candidate",
			60, 0.44, 0.56);

  ksMass_eta_2__2_5 = new TH1F("ksMass_eta_2__2_5", "K^{0}_{s} invariant mass, 2 < #eta < 2.5",
			       120, 0.1, 0.9);
  ksMass_eta_1_5__2 = new TH1F("ksMass_eta_1_5__2", "K^{0}_{s} invariant mass, 1.5 < #eta < 2",
			       120, 0.1, 0.9);
  ksMass_eta_1__1_5 = new TH1F("ksMass_eta_1__1_5", "K^{0}_{s} invariant mass, 1 < #eta < 1.5",
			       120, 0.1, 0.9);
  ksMass_eta_0_5__1 = new TH1F("ksMass_eta_0_5__1", "K^{0}_{s} invariant mass, 0.5 < #eta < 1",
			       120, 0.1, 0.9);
  ksMass_eta_0__0_5 = new TH1F("ksMass_eta_0__0_5", "K^{0}_{s} invariant mass, 0 < #eta < 0.5",
			       120, 0.1, 0.9);

  ntuple->Project("ksMass_eta_2__2_5", "v0CandMass", 
		  "abs(v0DauEta_lowestpT) > 2.");
  ntuple->Project("ksMass_eta_1_5__2", "v0CandMass", 
		  "abs(v0DauEta_lowestpT) < 2. && abs(v0DauEta_lowestpT) > 1.5");
  ntuple->Project("ksMass_eta_1__1_5", "v0CandMass",
		  "abs(v0DauEta_lowestpT) < 1.5 && abs(v0DauEta_lowestpT) > 1.");
  ntuple->Project("ksMass_eta_0_5__1", "v0CandMass",
		  "abs(v0DauEta_lowestpT) < 1. && abs(v0DauEta_lowestpT) > 0.5");
  ntuple->Project("ksMass_eta_0__0_5", "v0CandMass",
		  "abs(v0DauEta_lowestpT) < 0.5");

  ksMass_eta_2__2_5->Draw();
  c1->SaveAs("ksMass_eta_2__2_5.root");
  ksMass_eta_1_5__2->Draw();
  c1->SaveAs("ksMass_eta_1_5__2.root");
  ksMass_eta_1__1_5->Draw();
  c1->SaveAs("ksMass_eta_1__1_5.root");
  ksMass_eta_0_5__1->Draw();
  c1->SaveAs("ksMass_eta_0_5__1.root");
  ksMass_eta_0__0_5->Draw();
  c1->SaveAs("ksMass_eta_0__0_5.root");

  ntuple->Project("ksRFull", "v0VtxR", "");
  ntuple->Project("ksRSignal", "v0VtxR", 
		  "v0CandMass > (0.49767 - 0.01) && v0CandMass < (0.49767 + 0.01)");
  ntuple->Project("ksRBg", "v0VtxR",
		  "(v0CandMass > 0.445 && v0CandMass < 0.465) || (v0CandMass > 0.55 && v0CandMass < 0.57)");

  ntuple->Project("ksZvsRFull", "v0VtxR:v0VtxZ", "");
  ntuple->Project("ksZvsRSignal", "v0VtxR:v0VtxZ", 
		  "v0CandMass > (0.49767 - 0.01) && v0CandMass < (0.49767 + 0.01)");
  ntuple->Project("ksZvsRBg", "v0VtxR:v0VtxZ",
		  "(v0CandMass > 0.445 && v0CandMass < 0.465) || (v0CandMass > 0.55 && v0CandMass < 0.57)");

  //ntuple->Draw("v0VtxZ:v0VtxR >> ksZvsRFull", "");
  //ntuple->Draw("v0VtxZ:v0VtzR >> ksZvsRSignal",
  //"v0CandMass > (0.49767 - 0.01) && v0CandMass < (0.49767 + 0.01)");
//ntuple->Draw("v0VtxZ:v0VtxR >> ksZvsRBg",
  //"(v0CandMass > 0.445 && v0CandMass < 0.465) || (v0CandMass > 0.55 && v0CandMass < 0.57)");

  ntuple->Project("ksEtaFull", "v0Eta", "");
  ntuple->Project("ksEtaSignal", "v0Eta", 
		  "v0CandMass > (0.49767 - 0.01) && v0CandMass < (0.49767 + 0.01)");
  ntuple->Project("ksEtaBg", "v0Eta",
		  "(v0CandMass > 0.445 && v0CandMass < 0.465) || (v0CandMass > 0.55 && v0CandMass < 0.57)");

  ntuple->Project("ksPhiFull", "v0Phi", "");
  ntuple->Project("ksPhiSignal", "v0Phi", 
		  "v0CandMass > (0.49767 - 0.01) && v0CandMass < (0.49767 + 0.01)");
  ntuple->Project("ksPhiBg", "v0Phi",
		  "(v0CandMass > 0.445 && v0CandMass < 0.465) || (v0CandMass > 0.55 && v0CandMass < 0.57)");

  ntuple->Project("ksPtFull", "v0pT", "");
  ntuple->Project("ksPtSignal", "v0pT", 
		  "v0CandMass > (0.49767 - 0.01) && v0CandMass < (0.49767 + 0.01)");
  ntuple->Project("ksPtBg", "v0pT",
		  "(v0CandMass > 0.445 && v0CandMass < 0.465) || (v0CandMass > 0.55 && v0CandMass < 0.57)");

  ntuple->Project("ksMass", "v0Mass", "");
  ntuple->Project("ksNormMass", "v0Massnorm", "");
  ntuple->Project("ksCandMass", "v0CandMass", "");

  ntuple->Project("ksDaup", "v0Daup", "");
  ntuple->Project("ksDaupT", "v0DaupT", "");
  ntuple->Project("ksDauEta", "v0DauEta", "");
  ntuple->Project("ksDauPhi", "v0DauPhi", "");
  ntuple->Project("ksDauNHits", "v0DauNHits", "");
  ntuple->Project("ksDauNValidHits", "v0DauNValidHits", "");
  ntuple->Project("ksDauNPixelHits", "v0DauNPixelHits", "");
  ntuple->Project("ksDauNStripHits", "v0DauNStripHits", "");
  ntuple->Project("ksDauChi2", "v0DauChi2", "");
  ntuple->Project("ksDauNormChi2", "v0DauNormChi2", "");

  //ksRSignal->Sumw2();
  //ksEtaSignal->Sumw2();
  //ksPtSignal->Sumw2();

  ksR->Add(ksRSignal, ksRBg, 1.0, -0.5);
  ksEta->Add(ksEtaSignal, ksEtaBg, 1.0, -0.5);
  ksPhi->Add(ksPhiSignal, ksPhiBg, 1.0, -0.5);
  ksPt->Add(ksPtSignal, ksPtBg, 1.0, -0.5);
  ksZvsRBgSubt->Add(ksZvsRSignal, ksZvsRBg, 1.0, -0.5);

  ksR->SetXTitle("Transverse K^{0}_{s} decay radius (cm)");
  ksEta->SetXTitle("K^{0}_{s} momentum #eta");
  ksPt->SetXTitle("K^{0}_{s} p_{T} (GeV)");
  ksPhi->SetXTitle("K^{0}_{s} momentum #phi");
  ksZvsRBgSubt->SetXTitle("K^{0}_{s} decay vertex radius (cm)");
  ksZvsRBgSubt->SetYTitle("K^{0}_{s} decay vertex z position (cm)");

  ksR->Draw();
  c1->SaveAs("KsR_backgroundSubtracted.png");
  c1->SaveAs("KsR_backgroundSubtracted.root");

  ksEta->Draw();
  c1->SaveAs("KsEta_backgroundSubtracted.png");
  c1->SaveAs("KsEta_backgroundSubtracted.root");

  ksZvsRBgSubt->Draw();
  c1->SaveAs("KsZvsR_backgroundSubtracted.png");
  c1->SaveAs("KsZvsR_backgroundSubtracted.root");

  ksZvsRFull->Draw();
  c1->SaveAs("KsZvsR.png");
  c1->SaveAs("KsZvsR.root");

  ksPhi->Draw();
  c1->SaveAs("KsPhi_backgroundSubtracted.png");
  c1->SaveAs("KsPhi_backgroundSubtracted.root");

  ksPt->Draw();
  c1->SaveAs("KspT_backgroundSubtracted.png");
  c1->SaveAs("KspT_backgroundSubtracted.root");

  ksRFull->Draw();
  c1->SaveAs("KsR.png");
  c1->SaveAs("KsR.root");

  ksEtaFull->Draw();
  c1->SaveAs("KsEta.png");
  c1->SaveAs("KsEta.root");

  ksPhiFull->Draw();
  c1->SaveAs("KsPhi.png");
  c1->SaveAs("KsPhi.root");

  ksPtFull->Draw();
  c1->SaveAs("KsPt.png");
  c1->SaveAs("KsPt.root");

  //string ksFitStr = (string) ksBinWidth;
  //ksFitStr += string("*gausn(0) + [3] + [4]*(x - 0.49767)");
  string ksFitStr;
  ostringstream ksFitSS;
  //cout << ksBinWidth;
  ksFitSS << ksBinWidth << "*gausn(0) + [3] + [4]*(x-0.49767)";
  cout << ksFitSS.str() << endl;

  //TF1 *ksFit = new TF1("ksFit", "gausn(0) + [3] + [4]*(x - 0.49767)");
  TF1 *ksFit = new TF1("ksFit", ksFitSS.str().c_str());
  ksFit->SetParName( 0, "yield" );
  ksFit->SetParName( 1, "mean" );
  ksFit->SetParName( 2, "sigma" );
  ksFit->SetParName( 3, "const" );
  ksFit->SetParName( 4, "slope" );
  //ksFit->SetParName( 5, "binWidth");

  //double yieldScale = ksMass->GetMaximum();
  //cout << "Maximum is " << yieldScale << endl;

  // For starting sample, yield start is 4 with a maximum of 451.
  // ratio is 2000/450 = 4.4444
  // 4.4444*max = 2000, so make yield = 4.4444 * ksMass->GetMaximum()

  ksFit->SetParameter( 0, 4.4444*ksMass->GetMaximum() );
  //ksFit->SetParameter( 0, 2385 );
  ksFit->SetParameter( 1, 0.49767 );
  ksFit->SetParameter( 2, 0.005 );
  ksFit->SetParameter( 3, 30. );
  ksFit->SetParameter( 4, 50. );
  //ksFit->SetParameter( 5, ksBinWidth );

  //ksFit->Draw();
  //c1->SaveAs("KsFit.png");
  //c1->SaveAs("KsFit.root");

  gStyle->SetOptFit(1);
  ksMass->Draw();
  ksMass->Fit("ksFit", "RL");
  c1->SaveAs("KsMass.png");
  c1->SaveAs("KsMass.root");

  //ksFit->Draw();
  //c1->SaveAs("KsFitPostFit.png");
  //c1->SaveAs("KsFitPostFit.root");

  ostringstream ksNormFitSS;
  ksNormFitSS << ksNormBinWidth << "*gausn(0) + [3] - [4]*(x^2)";
  cout << ksNormFitSS.str() << endl;

  TF1 *ksNormFit = new TF1("ksNormFit", ksNormFitSS.str().c_str(), ksNormXmin, ksNormXmax);
  ksNormFit->SetParName( 0, "yield" );
  ksNormFit->SetParName( 1, "mean" );
  ksNormFit->SetParName( 2, "sigma" );
  ksNormFit->SetParName( 3, "const" );
  ksNormFit->SetParName( 4, "quadInvWidth" );

  ksNormFit->SetParameter( 0, 2500 );
  ksNormFit->SetParameter( 1, 0. );
  ksNormFit->SetParameter( 2, 1. );
  ksNormFit->SetParameter( 3, 80. );
  ksNormFit->SetParameter( 4, 0.17 );

  ostringstream ksNormFitSS2;
  ksNormFitSS2 << ksNormBinWidth << "*gausn(0) + " << ksNormBinWidth << "*gausn(3)";

  TF1 *ksNormFit2 = new TF1("ksNormFit2", ksNormFitSS2.str().c_str(), ksNormXmin, ksNormXmax);
  ksNormFit2->SetParName( 0, "peakYield" );
  ksNormFit2->SetParName( 1, "peakMean" );
  ksNormFit2->SetParName( 2, "peakSigma" );
  ksNormFit2->SetParName( 3, "bgYield" );
  ksNormFit2->SetParName( 4, "bgMean" );
  ksNormFit2->SetParName( 5, "bgSigma" );

  ksNormFit2->SetParameter( 0, 2500 );
  ksNormFit2->SetParameter( 1, 0. );
  ksNormFit2->SetParameter( 2, 1. );
  ksNormFit2->SetParameter( 3, 40*ksNormNbins );
  ksNormFit2->SetParameter( 4, -0.001 );
  ksNormFit2->SetParameter( 5, 15. );

  ksNormMass->Draw();
  ksNormMass->Fit("ksNormFit", "RL");
  //ksNormFit->Draw("same");
  c1->SaveAs("KsNormMass.png");
  c1->SaveAs("KsNormMass.root");

  ksCandMass->Fit("ksFit", "RL");
  ksCandMass->Draw();
  c1->SaveAs("KsCandMass.png");
  c1->SaveAs("KsCandMass.root");

  gStyle->SetOptFit(0);

  ksDaup->Draw();
  c1->SaveAs("KsDaup.png");
  c1->SaveAs("KsDaup.root");

  ksDaupT->Draw();
  c1->SaveAs("KsDaupT.png");
  c1->SaveAs("KsDaupT.root");

  ksDauEta->Draw();
  c1->SaveAs("KsDauEta.png");
  c1->SaveAs("KsDauEta.root");

  ksDauPhi->Draw();
  c1->SaveAs("KsDauPhi.png");
  c1->SaveAs("KsDauPhi.root");

  ksDauNHits->Draw();
  c1->SaveAs("KsDauNHits.png");
  c1->SaveAs("KsDauNHits.root");

  ksDauNValidHits->Draw();
  c1->SaveAs("KsDauNValidHits.png");
  c1->SaveAs("KsDauNValidHits.root");

  ksDauNPixelHits->Draw();
  c1->SaveAs("KsDauNPixelHits.png");
  c1->SaveAs("KsDauNPixelHits.root");

  ksDauNStripHits->Draw();
  c1->SaveAs("KsDauNStripHits.png");
  c1->SaveAs("KsDauNStripHits.root");

  ksDauChi2->Draw();
  c1->SaveAs("KsDauChi2.png");
  c1->SaveAs("KsDauChi2.root");

  ksDauNormChi2->Draw();
  c1->SaveAs("KsDauNormChi2.png");
  c1->SaveAs("KsDauNormChi2.root");

  // Lambda time
  infile.cd("lamAnalyzer");
  // Trajectory variables
  lamR = new TH1F("lamR", "#Lambda^{0} radial distance from beam line",
		 40, 0., 40.);
  lamEta = new TH1F("lamEta", "#Lambda^{0} momentum #eta",
		   40, -2.5, 2.5);
  lamPhi = new TH1F("lamPhi", "#Lambda^{0} momentum #phi",
		    60, -M_PI, M_PI);
  lamPt = new TH1F("lamPt", "#Lambda^{0} p_{T}",
		  70, 0., 20.);

  lamRFull = new TH1F("lamRFull", "#Lambda^{0} radial distance from beam line",
		 40, 0., 40.);
  lamEtaFull = new TH1F("lamEtaFull", "#Lambda^{0} momentum #eta",
		   40, -2.5, 2.5);
  lamPhiFull = new TH1F("lamPhiFull", "#Lambda^{0} momentum #phi",
		    60, -M_PI, M_PI);
  lamPtFull = new TH1F("lamPtFull", "#Lambda^{0} p_{T}",
		  70, 0., 20.);

  lamRSignal = new TH1F("lamRSignal", "#Lambda^{0} radial distance from beam line",
		 40, 0., 40.);
  lamEtaSignal = new TH1F("lamEtaSignal", "#Lambda^{0} momentum #eta",
		   40, -2.5, 2.5);
  lamPhiSignal = new TH1F("lamPhiSignal", "#Lambda^{0} momentum #phi",
		    60, -M_PI, M_PI);
  lamPtSignal = new TH1F("lamPtSignal", "#Lambda^{0} p_{T}",
		  70, 0., 20.);

  lamRBg = new TH1F("lamRBg", "#Lambda^{0} radial distance from beam line",
		 40, 0., 40.);
  lamEtaBg = new TH1F("lamEtaBg", "#Lambda^{0} momentum #eta",
		   40, -2.5, 2.5);
  lamPhiBg = new TH1F("lamPhiBg", "#Lambda^{0} momentum #phi",
		    60, -M_PI, M_PI);
  lamPtBg = new TH1F("lamPtBg", "#Lambda^{0} p_{T}",
		  70, 0., 20.);

  lamZvsRFull = new TH2F("lamZvsRFull", "#Lambda^{0} z position vs. decay radius",
			 120, -120., 120., 40, 0., 40.);
  lamZvsRSignal = new TH2F("lamZvsRSignal", "#Lambda^{0} z position vs. decay radius",
			 120, -120., 120., 40, 0., 40.);
  lamZvsRBg = new TH2F("lamZvsRBg", "#Lambda^{0} z position vs. decay radius",
		       120, -120., 120., 40, 0., 40.);
  lamZvsRBgSubt = new TH2F("lamZvsRBgSubt", "#Lambda^{0} z position vs. decay radius",
			   120, -120., 120., 40, 0., 40.);

  lamDaup = new TH1F("lamDaup", "p of #Lambda^{0} daughter tracks",
		    100, 0., 30.);
  lamDaupT = new TH1F("lamDaupT", "p_{T} of #Lambda^{0} daughter tracks",
		     100, 0., 20.);
  lamDauEta = new TH1F("lamDauEta", "#eta of #Lambda^{0} daughter tracks",
		      100, -2.5, 2.5);
  lamDauPhi = new TH1F("lamDauPhi", "#phi of #Lambda^{0} daughter track at vertex position",
		      70, -M_PI, M_PI);
  lamDauNHits = new TH1F("lamDauNHits", "Total number of hits on #Lambda^{0} daughter tracks",
			80, 0., 40.);
  lamDauNValidHits = new TH1F("lamDauNValHits", "Number of valid hits on #Lambda^{0} daughter tracks", 
			     80, 0., 40.);
  lamDauNPixelHits = new TH1F("lamDauNPixelHits", "Number of valid pixel hits on #Lambda^{0} daughter tracks",
			     30, 0., 15.);
  lamDauNStripHits = new TH1F("lamDauNStripHits", "Number of valid strip hits on #Lambda^{0} daughter tracks",
			     80, 0., 40.);
  lamDauChi2 = new TH1F("lamDauChi2", "#chi^{2} of #Lambda^{0} daughter tracks",
		       100, 0., 100.);
  lamDauNormChi2 = new TH1F("lamDauNormChi2", "#chi^{2}/ndof of #Lambda^{0} daughter tracks",
			   100, 0., 5.);

  //lamm = 0.49767;
  // Mass plots
  double lamXmin = 1.07, lamXmax = 1.16;
  int lamNbins = 60;
  double lamBinWidth = (lamXmax - lamXmin) / (double) lamNbins;

  double lamNormXmin = -10., lamNormXmax = 10.;
  int lamNormNbins = 60;
  double lamNormBinWidth = (lamNormXmax - lamNormXmin) / (double) lamNormNbins;
  lamMass = new TH1F("lamMass", "#Lambda^{0} invariant mass",
		     lamNbins, lamXmin, lamXmax);
  lamNormMass = new TH1F("lamNormMass", "#Lambda^{0} normalized mass",
			 lamNormNbins, lamNormXmin, lamNormXmax);
  lamCandMass = new TH1F("lamCandMass", "#Lambda^{0} mass from Candidate",
			 lamNbins, lamXmin, lamXmax);

  ntuple->Project("lamRFull", "v0VtxR", "");
  ntuple->Project("lamRSignal", "v0VtxR", 
		  "v0CandMass > (1.1156 - 0.01) && v0CandMass < (1.1156 + 0.01)");
  ntuple->Project("lamRBg", "v0VtxR",
		  "(v0CandMass > 1.08 && v0CandMass < 1.1) || (v0CandMass > 1.14 && v0CandMass < 1.16)");
  
  ntuple->Project("lamZvsRFull", "v0VtxR:v0VtxZ", "");
  ntuple->Project("lamZvsRSignal", "v0VtxR:v0VtxZ", 
		  "v0CandMass > (1.1156 - 0.01) && v0CandMass < (1.1156 + 0.01)");
  ntuple->Project("lamZvsRBg", "v0VtxR:v0VtxZ",
		  "(v0CandMass > 1.08 && v0CandMass < 1.1) || (v0CandMass > 1.14 && v0CandMass < 1.16)");

  ntuple->Project("lamEtaFull", "v0Eta", "");
  ntuple->Project("lamEtaSignal", "v0Eta", 
		  "v0CandMass > (1.1156 - 0.01) && v0CandMass < (1.1156 + 0.01)");
  ntuple->Project("lamEtaBg", "v0Eta",
		  "(v0CandMass > 1.08 && v0CandMass < 1.1) || (v0CandMass > 1.14 && v0CandMass < 1.16)");

  ntuple->Project("lamPhiFull", "v0Phi", "");
  ntuple->Project("lamPhiSignal", "v0Phi", 
		  "v0CandMass > (1.1156 - 0.01) && v0CandMass < (1.1156 + 0.01)");
  ntuple->Project("lamPhiBg", "v0Phi",
		  "(v0CandMass > 1.08 && v0CandMass < 1.1) || (v0CandMass > 1.14 && v0CandMass < 1.16)");

  ntuple->Project("lamPtFull", "v0pT", "");
  ntuple->Project("lamPtSignal", "v0pT", 
		  "v0CandMass > (1.1156 - 0.01) && v0CandMass < (1.1156 + 0.01)");
  ntuple->Project("lamPtBg", "v0pT",
		  "(v0CandMass > 1.08 && v0CandMass < 1.1) || (v0CandMass > 1.14 && v0CandMass < 1.16)");
  
  ntuple->Project("lamMass", "v0Mass", "");
  ntuple->Project("lamNormMass", "v0Massnorm", "");
  ntuple->Project("lamCandMass", "v0CandMass", "");

  ntuple->Project("lamDaup", "v0Daup", "");
  ntuple->Project("lamDaupT", "v0DaupT", "");
  ntuple->Project("lamDauEta", "v0DauEta", "");
  ntuple->Project("lamDauPhi", "v0DauPhi", "");
  ntuple->Project("lamDauNHits", "v0DauNHits", "");
  ntuple->Project("lamDauNValidHits", "v0DauNValidHits", "");
  ntuple->Project("lamDauNPixelHits", "v0DauNPixelHits", "");
  ntuple->Project("lamDauNStripHits", "v0DauNStripHits", "");
  ntuple->Project("lamDauChi2", "v0DauChi2", "");
  ntuple->Project("lamDauNormChi2", "v0DauNormChi2", "");

  //lamRSignal->Sumw2();
  //lamEtaSignal->Sumw2();
  //lamPtSignal->Sumw2();

  lamR->Add(lamRSignal, lamRBg, 1.0, -0.5);
  lamEta->Add(lamEtaSignal, lamEtaBg, 1.0, -0.5);
  lamPhi->Add(lamPhiSignal, lamPhiBg, 1.0, -0.5);
  lamPt->Add(lamPtSignal, lamPtBg, 1.0, -0.5);
  lamZvsRBgSubt->Add(lamZvsRSignal, lamZvsRBg, 1.0, -0.5);

  lamR->SetXTitle("Transverse #Lambda^{0} decay radius (cm)");
  lamEta->SetXTitle("#Lambda^{0} momentum #eta");
  lamPt->SetXTitle("#Lambda^{0} p_{T} (GeV)");

  lamR->Draw();
  c1->SaveAs("LamR_backgroundSubtracted.png");
  c1->SaveAs("LamR_backgroundSubtracted.root");

  lamZvsRBgSubt->Draw();
  c1->SaveAs("LamZvsR_backgroundSubtracted.png");
  c1->SaveAs("LamZvsR_backgroundSubtracted.root");

  lamZvsRFull->Draw();
  c1->SaveAs("LamZvsR.png");
  c1->SaveAs("LamZvsR.root");

  lamEta->Draw();
  c1->SaveAs("LamEta_backgroundSubtracted.png");
  c1->SaveAs("LamEta_backgroundSubtracted.root");

  lamPhi->Draw();
  c1->SaveAs("LamPhi_backgroundSubtracted.png");
  c1->SaveAs("LamPhi_backgroundSubtracted.root");

  lamPt->Draw();
  c1->SaveAs("LampT_backgroundSubtracted.png");
  c1->SaveAs("LampT_backgroundSubtracted.root");

  lamRFull->Draw();
  c1->SaveAs("LamR.png");
  c1->SaveAs("LamR.root");

  lamEtaFull->Draw();
  c1->SaveAs("LamEta.png");
  c1->SaveAs("LamEta.root");

  lamPhiFull->Draw();
  c1->SaveAs("LamPhi.png");
  c1->SaveAs("LamPhi.root");

  lamPtFull->Draw();
  c1->SaveAs("LamPt.png");
  c1->SaveAs("LamPt.root");

  //string ksFitStr;

  //y = a*(1.0 + b*x^c)
  ostringstream lamFitSS;
  lamFitSS << lamBinWidth << "*gausn(0) + [3]*(1.0 + [4]*(x-1.1156)^[5])";
  //cout << lamFitSS.str() << endl;
  ostringstream lamFitSS2;
  lamFitSS2 << lamBinWidth << "*gausn(0) + [3]*(x - " << piMass + protonMass 
	    << ")^(1/2) + [4]*(x - " << piMass + protonMass << ")^(3/2)";
  cout << lamFitSS2.str() << endl;

  TF1 *lamFit = new TF1("lamFit", lamFitSS.str().c_str(), 1.07, 1.16);
  lamFit->SetParName( 0, "yield" );
  lamFit->SetParName( 1, "mean" );
  lamFit->SetParName( 2, "sigma" );
  lamFit->SetParName( 3, "const" );
  lamFit->SetParName( 4, "slope" );
  lamFit->SetParName( 5, "exponent");

  lamFit->SetParameter( 0, 4. );
  lamFit->SetParameter( 1, 1.1156 );
  lamFit->SetParameter( 2, 0.005 );
  lamFit->SetParameter( 3, 30. );
  lamFit->SetParameter( 4, 50. );
  lamFit->SetParameter( 5, 0.5 );

  //cout << "LamMassMax: " << lamMass->GetMaximum() << endl;

  // For starting sample, yield start is 700 with a maximum of 172.
  // ratio is 700/172 = 4.0698
  // 4.0698*max = 700, so make yield = 4.0698 * lamMass->GetMaximum()

  //TF1 *lamFit2 = new TF1("lamFit2", lamFitSS2.str().c_str(), lamXmin, lamXmax);
  TF1 *lamFit2 = new TF1("lamFit2", lamFitSS2.str().c_str(), piMass + protonMass, lamXmax);
  lamFit2->SetParName( 0, "yield" );
  lamFit2->SetParName( 1, "mean" );
  lamFit2->SetParName( 2, "sigma" );
  lamFit2->SetParName( 3, "c1" );
  lamFit2->SetParName( 4, "c2" );

  lamFit2->SetParameter( 0, 4.0698*lamMass->GetMaximum() );
  lamFit2->SetParameter( 1, 1.1156 );
  lamFit2->SetParameter( 2, 0.05 );
  lamFit2->SetParameter( 3, 120. );
  lamFit2->SetParameter( 4, 100. );

  gStyle->SetOptFit(1);

  lamMass->Draw();
  lamMass->Fit("lamFit2", "RLE");
  c1->SaveAs("LamMass.png");
  c1->SaveAs("LamMass.root");

  //lamFit2->Draw();
  //c1->SaveAs("LamFit.png");
  //c1->SaveAs("LamFit.root");

  ostringstream lamNormFitSS2;
  lamNormFitSS2 << lamNormBinWidth << "*gausn(0) + " << lamNormBinWidth << "*gausn(3)";

  TF1 *lamNormFit2 = new TF1("lamNormFit2", lamNormFitSS2.str().c_str(), -22., 22.);
  lamNormFit2->SetParName( 0, "peakYield" );
  lamNormFit2->SetParName( 1, "peakMean" );
  lamNormFit2->SetParName( 2, "peakSigma" );
  lamNormFit2->SetParName( 3, "bgYield" );
  lamNormFit2->SetParName( 4, "bgMean" );
  lamNormFit2->SetParName( 5, "bgSigma" );

  lamNormFit2->SetParameter( 0, 650. );
  lamNormFit2->SetParameter( 1, 0. );
  lamNormFit2->SetParameter( 2, 2. );
  lamNormFit2->SetParameter( 3, 40*lamNormNbins) ;
  lamNormFit2->SetParameter( 4, 7. );
  lamNormFit2->SetParameter( 5, 7. );

  lamNormMass->Draw();
  //lamNormMass->Fit("lamNormFit2", "RLE");
  //lamNormFit2->Draw("same");
  c1->SaveAs("LamNormMass.png");
  c1->SaveAs("LamNormMass.root");

  lamCandMass->Draw();
  lamCandMass->Fit("lamFit2", "RLE");
  c1->SaveAs("LamCandMass.png");
  c1->SaveAs("LamCandMass.root");

  lamDaup->Draw();
  c1->SaveAs("LamDaup.png");
  c1->SaveAs("LamDaup.root");

  lamDaupT->Draw();
  c1->SaveAs("LamDaupT.png");
  c1->SaveAs("LamDaupT.root");

  lamDauEta->Draw();
  c1->SaveAs("LamDauEta.png");
  c1->SaveAs("LamDauEta.root");

  lamDauPhi->Draw();
  c1->SaveAs("LamDauPhi.png");
  c1->SaveAs("LamDauPhi.root");

  lamDauNHits->Draw();
  c1->SaveAs("LamDauNHits.png");
  c1->SaveAs("LamDauNHits.root");

  lamDauNValidHits->Draw();
  c1->SaveAs("LamDauNValidHits.png");
  c1->SaveAs("LamDauNValidHits.root");

  lamDauNPixelHits->Draw();
  c1->SaveAs("LamDauNPixelHits.png");
  c1->SaveAs("LamDauNPixelHits.root");

  lamDauNStripHits->Draw();
  c1->SaveAs("LamDauNStripHits.png");
  c1->SaveAs("LamDauNStripHits.root");

  lamDauChi2->Draw();
  c1->SaveAs("LamDauChi2.png");
  c1->SaveAs("LamDauChi2.root");

  lamDauNormChi2->Draw();
  c1->SaveAs("LamDauNormChi2.png");
  c1->SaveAs("LamDauNormChi2.root");

}
