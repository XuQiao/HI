#ifndef lamfit_cxx
#define lamfit_cxx
#include "lamfit.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

int lamfit::Plot(TH1* lamHist) {

    lamHist->Draw();
    if (fittype == 1) {
      lamFit_expbg->SetParameter( 0, lamFit_DG_exp->GetParameter(5) );
      lamFit_expbg->SetParameter( 1, lamFit_DG_exp->GetParameter(6) );
      lamFit_expbg->Draw("same");
    } else if ( fittype == 2 ) {
      lamFit_rootbg->SetParameter( 0, lamFit_DG->GetParameter(5) );
      lamFit_rootbg->SetParameter( 1, lamFit_DG->GetParameter(6) );
      lamFit_rootbg->Draw("same");
    } else if ( fittype == 3 ) {
      lamFit_expppbg->SetParameter( 0, lamFit_DG_highbg->GetParameter(5) );
      lamFit_expppbg->SetParameter( 1, lamFit_DG_highbg->GetParameter(6) );
      lamFit_expppbg->SetParameter( 2, lamFit_DG_highbg->GetParameter(7) );
      lamFit_expppbg->SetParameter( 3, lamFit_DG_highbg->GetParameter(8) );
      lamFit_expppbg->Draw("same");
    }
    labl->Draw("same");
}

int lamfit::Fit(TH1* lamHist, int fitcomm) {

  // These chi2 values are chi2/ndof
  double expChi2 = 0., highbgChi2 = 0., rootsChi2 = 0.;
  bool expIsSG = false, highbgIsSG = false, rootsIsSG = false;
  if ( fitcomm == 1 ) {
    expIsSG = true;
    highbgIsSG = true;
    rootsIsSG = true;
  }
  int expFitStat = FitWithExpBG( lamHist, expChi2 , expIsSG, "QLEM" );
  int highbgFitStat = FitWithHighbgBG( lamHist, highbgChi2, highbgIsSG, "QLEM" );
  int rootsFitStat = FitWithRootsBG( lamHist, rootsChi2, rootsIsSG, "QLEM" );

  // We'll only check the chi2 for those that had a good double gaussian fit,
  //   unless none did.
  bool checkAll = !(expIsSG && highbgIsSG && rootsIsSG) 
    || (expIsSG && highbgIsSG && rootsIsSG);
  bool checkExpAndHighbg = !(expIsSG && highbgIsSG && !rootsIsSG);
  bool checkExpAndRoots = !(expIsSG && !highbgIsSG && rootsIsSG);
  bool checkHighbgAndRoots = !(!expIsSG && highbgIsSG && rootsIsSG);
  bool checkOnlyExp = !(expIsSG && !highbgIsSG && !rootsIsSG);
  bool checkOnlyHighbg = !(!expIsSG && highbgIsSG && !rootsIsSG);
  bool checkOnlyRoots = !(!expIsSG && !highbgIsSG && rootsIsSG);

  // I admit, having to refit again is not ideal, but we need to fill the
  //  signal/noise values, etc. using the good fit.
  if( checkAll ) {
    if( expChi2 < highbgChi2 && expChi2 < rootsChi2 ) {
      norm_chi2 = expChi2;
      return FitWithExpBG( lamHist, expChi2, expIsSG, "LEM" );
    }
    else if( highbgChi2 < expChi2 && highbgChi2 < rootsChi2 ) {
      norm_chi2 = highbgChi2;
      return FitWithHighbgBG( lamHist, highbgChi2, highbgIsSG, "LEM" );
    }
    else if( rootsChi2 < expChi2 && rootsChi2 < highbgChi2 ) {
      norm_chi2 = rootsChi2;
      return FitWithRootsBG( lamHist, rootsChi2, rootsIsSG, "LEM" );
    }
  }
  else if( checkExpAndHighbg ) {
    if( expChi2 < highbgChi2 ) {
      norm_chi2 = expChi2;
      return FitWithExpBG( lamHist, expChi2, expIsSG, "LEM" );
    }
    else {
      norm_chi2 = highbgChi2;
      return FitWithHighbgBG( lamHist, expChi2, expIsSG, "LEM" );
    }
  }
  else if( checkExpAndRoots ) {
    if( expChi2 < rootsChi2 ) {
      norm_chi2 = expChi2;
      return FitWithExpBG( lamHist, expChi2, expIsSG, "LEM" );
    }
    else {
      // Roots was the last to run, so just spit out the FitStat
      norm_chi2 = rootsChi2;
      return rootsFitStat;
      //return FitWithRootsBG( lamHist, rootsChi2, rootsIsSG, "LEM" );
    }
  }
  else if( checkHighbgAndRoots ) {
    if( highbgChi2 < rootsChi2 ) {
      norm_chi2 = highbgChi2;
      return FitWithHighbgBG( lamHist, highbgChi2, highbgIsSG, "LEM" );
    }
    else {
      norm_chi2 = rootsChi2;
      return rootsFitStat;
      //return FitWithRootsBG( lamHist, rootsChi2, rootsIsSG, "LEM" );
    }
  }
  else if (checkOnlyExp ) {
    norm_chi2 = expChi2;
    return FitWithExpBG( lamHist, expChi2, expIsSG, "LEM" );
  }
  else if ( checkOnlyHighbg ) {
    norm_chi2 = highbgChi2;
    return FitWithHighbgBG( lamHist, highbgChi2, highbgIsSG, "LEM" );
  }
  else {
    norm_chi2 = rootsChi2;
    return rootsFitStat;
  }

}



int lamfit::FitWithExpBG(TH1* lamHist, double &chi2, bool &isSG, TString fitopt){

  //cout << "ExpBG: setting yield, mass, bg = " << 4.5*lamHist->GetMaximum() << ", " << lamMassConst << ", " << 0.1*lamHist->GetMaximum() << endl;
  lamFit_SG_exp->SetParameter( 0, 4.0*lamHist->GetMaximum() );
  lamFit_SG_exp->SetParameter( 1, lamMassConst );
  lamFit_SG_exp->SetParameter( 2, 0.002 );
  lamFit_SG_exp->SetParLimits( 2, 0., 0.005 );
  lamFit_SG_exp->SetParameter( 3, 0.01*lamHist->GetMaximum() );
  lamFit_SG_exp->SetParameter( 4, 1.2 );

  // Performing initial SG fit to obtain starting values for DG
  //cout << "Performing initial SG fit..." << endl;
  int fit0 = lamHist->Fit("lamFit_SG_exp", "QLEM");
  fit0 = lamHist->Fit("lamFit_SG_exp", fitopt);
  if (fit0 != 4000) fit0 = lamHist->Fit("lamFit_SG_exp", fitopt);
  if (fit0 != 4000) fit0 = lamHist->Fit("lamFit_SG_exp", fitopt);
  if (fit0 != 4000) fit0 = lamHist->Fit("lamFit_SG_exp", fitopt);

  sigmae = lamFit_SG_exp->GetParError(2);
  int fitstat = -1;
  if ( !isSG ) {
    lamFit_DG_exp->SetParameter( 0, lamFit_SG_exp->GetParameter(0) );
    lamFit_DG_exp->SetParameter( 1, lamFit_SG_exp->GetParameter(1) );
    lamFit_DG_exp->SetParameter( 2, 0.75*lamFit_SG_exp->GetParameter(2) );
    lamFit_DG_exp->SetParLimits( 2, 0., 0.004 );
    lamFit_DG_exp->SetParameter( 3, 0.5 );
    lamFit_DG_exp->SetParameter( 4, min(0.0045,1.3*lamFit_SG_exp->GetParameter(2)) );
    lamFit_DG_exp->SetParLimits( 4, 0., 0.008 );
    lamFit_DG_exp->SetParameter( 5, lamFit_SG_exp->GetParameter(3) );
    lamFit_DG_exp->SetParameter( 6, lamFit_SG_exp->GetParameter(4) );

    fitstat = lamHist->Fit("lamFit_DG_exp", "QLEM");
    fitstat = lamHist->Fit("lamFit_DG_exp", "QLEM");
    fitstat = lamHist->Fit("lamFit_DG_exp", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG_exp", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG_exp", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG_exp", fitopt);

    bool refit = false;
    if( fitstat == 4000 ) {
      double sigma1 = lamFit_DG_exp->GetParameter(2);
      double sigma2 = lamFit_DG_exp->GetParameter(4);
      if( sigma1 > sigma2 ) {
	lamFit_DG_exp->SetParameter(2, sigma2);
	lamFit_DG_exp->SetParameter(4, sigma1);
	double frac = lamFit_DG_exp->GetParameter(3);
	lamFit_DG_exp->SetParameter(3, 1.-frac);
	refit = true;
      }
    } else {
      refit = true;
    }
    if( refit ) fitstat = lamHist->Fit("lamFit_DG_exp", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG_exp", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG_exp", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG_exp", fitopt);
  }

  double height = 0; 
  double heighte = 0;
  sigma = 0;
  yield = 0; yielde = 0;
  if( fitstat >= 0 ) {
    //cout << "Good double gaussian fit, getting parameters." << endl;
    double frac = lamFit_DG_exp->GetParameter(3);
    double sigma1 = lamFit_DG_exp->GetParameter(2);
    double sigma2 = lamFit_DG_exp->GetParameter(4);
    mean = lamFit_DG_exp->GetParameter(1);
    meane = lamFit_DG_exp->GetParError(1);
    sigma = sqrt( frac*sigma1*sigma1 + (1.0-frac)*sigma2*sigma2 );
    if (sigmae <= 0) sigmae = frac*lamFit_DG_exp->GetParError(2) + (1.0-frac)*lamFit_DG_exp->GetParError(2);
    height = lamFit_DG_exp->GetParameter(5);
    heighte = lamFit_DG_exp->GetParError(5);
    yield = lamFit_DG_exp->GetParameter(0);
    yielde = lamFit_DG_exp->GetParError(0);
    fittype = 1;
    chi2 = lamFit_DG_exp->GetChisquare()/lamFit_DG_exp->GetNDF();
  } 
  else {
    //cout << "WARNING: Bad double gaussian fit, refitting with a single gaussian." << endl;
    fitstat = lamHist->Fit("lamFit_SG_exp", fitopt);
    mean = lamFit_SG_exp->GetParameter(1);
    meane = lamFit_SG_exp->GetParError(1);
    sigma = lamFit_SG_exp->GetParameter(2);
    sigmae = lamFit_SG_exp->GetParError(2);
    height = lamFit_SG_exp->GetParameter(3);
    heighte = lamFit_SG_exp->GetParError(3);
    yield = lamFit_SG_exp->GetParameter(0);
    yielde = lamFit_SG_exp->GetParError(0);
    fittype = -1;
    chi2 = lamFit_SG_exp->GetChisquare()/lamFit_SG_exp->GetNDF();
    isSG = true;
  }

  background = 4.0*height*sigma/lamMassBinWidth;
  backgrounde = background*heighte/height;
  //cout << "Returning from lamfit::Fit " << fitstat << endl;

  return fitstat;
}



int lamfit::FitWithRootsBG(TH1* lamHist, double &chi2, bool &isSG, TString fitopt){

  //cout << "RootsBG: setting yield, mass, bg = " << 4.5*lamHist->GetMaximum() << ", " << lamMassConst << ", " << 0.1*lamHist->GetMaximum() << endl;
  lamFit_SG->SetParameter( 0, 4.0*lamHist->GetMaximum() );
  lamFit_SG->SetParameter( 1, lamMassConst );
  lamFit_SG->SetParameter( 2, 0.002 );
  lamFit_SG->SetParLimits( 2, 0., 0.005 );
  lamFit_SG->SetParameter( 3, 0.002*lamHist->GetMaximum() );
  lamFit_SG->SetParameter( 4, 0.002*lamHist->GetMaximum() );

  // Performing initial SG fit to obtain starting values for DG
  //cout << "Performing initial SG fit..." << endl;
  int fit0 = lamHist->Fit("lamFit_SG", "QLEM");
  fit0 = lamHist->Fit("lamFit_SG", fitopt);
  if (fit0 != 4000) fit0 = lamHist->Fit("lamFit_SG", fitopt);
  if (fit0 != 4000) fit0 = lamHist->Fit("lamFit_SG", fitopt);
  if (fit0 != 4000) fit0 = lamHist->Fit("lamFit_SG", fitopt);

  sigmae = lamFit_SG->GetParError(2);
  int fitstat = -1;
  if ( !isSG ) {
    lamFit_DG->SetParameter( 0, lamFit_SG->GetParameter(0) );
    lamFit_DG->SetParameter( 1, lamFit_SG->GetParameter(1) );
    lamFit_DG->SetParameter( 2, 0.75*lamFit_SG->GetParameter(2) );
    lamFit_DG->SetParLimits( 2, 0., 0.004 );
    lamFit_DG->SetParameter( 3, 0.5 );
    lamFit_DG->SetParLimits( 4, 0., 0.008 );
    lamFit_DG->SetParameter( 4, min(0.0045,1.3*lamFit_SG->GetParameter(2)) );
    lamFit_DG->SetParameter( 5, lamFit_SG->GetParameter(3) );
    lamFit_DG->SetParameter( 6, lamFit_SG->GetParameter(4) );

    //cout << "Fitting double Gaussian\n";
    fitstat = lamHist->Fit("lamFit_DG", "QLEM");
    fitstat = lamHist->Fit("lamFit_DG", "QLEM");
    fitstat = lamHist->Fit("lamFit_DG", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG", fitopt);

    bool refit = false;
    if( fitstat == 4000 ) {
      double sigma1 = lamFit_DG->GetParameter(2);
      double sigma2 = lamFit_DG->GetParameter(4);
      if( sigma1 > sigma2 ) {
	lamFit_DG->SetParameter(2, sigma2);
	lamFit_DG->SetParameter(4, sigma1);
	double frac = lamFit_DG->GetParameter(3);
	lamFit_DG->SetParameter(3, 1.-frac);
	refit = true;
      }
    } else {
      refit = true;
    }
    if( refit ) fitstat = lamHist->Fit("lamFit_DG", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG", fitopt);
  }

  double height = 0; 
  double heighte = 0;
  sigma = 0;
  yield = 0; yielde = 0;
  if( fitstat >= 0 ) {
    //cout << "Good double gaussian fit, getting parameters." << endl;
    double frac = lamFit_DG->GetParameter(3);
    double sigma1 = lamFit_DG->GetParameter(2);
    double sigma2 = lamFit_DG->GetParameter(4);
    mean = lamFit_DG->GetParameter(1);
    meane = lamFit_DG->GetParError(1);
    sigma = sqrt( frac*sigma1*sigma1 + (1.0-frac)*sigma2*sigma2 );
    if (sigmae <= 0) sigmae = frac*lamFit_DG->GetParError(2) + (1.0-frac)*lamFit_DG->GetParError(2);
    height = lamFit_DG->GetParameter(5);
    heighte = lamFit_DG->GetParError(5);
    yield = lamFit_DG->GetParameter(0);
    yielde = lamFit_DG->GetParError(0);
    fittype = 2;
    chi2 = lamFit_DG->GetChisquare()/lamFit_DG->GetNDF();
  } 
  else {
    //cout << "WARNING: Bad double gaussian fit, refitting with a single gaussian." << endl;
    fitstat = lamHist->Fit("lamFit_SG", fitopt);
    mean = lamFit_SG->GetParameter(1);
    meane = lamFit_SG->GetParError(1);
    sigma = lamFit_SG->GetParameter(2);
    sigmae = lamFit_SG->GetParError(2);
    height = lamFit_SG->GetParameter(3);
    heighte = lamFit_SG->GetParError(3);
    yield = lamFit_SG->GetParameter(0);
    yielde = lamFit_SG->GetParError(0);
    fittype = -2;
    chi2 = lamFit_SG->GetChisquare()/lamFit_SG->GetNDF();
    isSG = true;
  }

  background = 4.0*height*sigma/lamMassBinWidth;
  backgrounde = background*heighte/height;
  //cout << "Returning from lamfit::Fit " << fitstat << endl;

  return fitstat;
}



int lamfit::FitWithHighbgBG(TH1* lamHist, double &chi2, bool &isSG, TString fitopt){

  //cout << "HighbgBG: setting yield, mass, bg = " << 4.5*lamHist->GetMaximum() << ", " << lamMassConst << ", " << 0.1*lamHist->GetMaximum() << endl;
  lamFit_SG_highbg->SetParameter( 0, 4.0*lamHist->GetMaximum() );
  lamFit_SG_highbg->SetParameter( 1, lamMassConst );
  lamFit_SG_highbg->SetParameter( 2, 0.002 );
  lamFit_SG_highbg->SetParLimits( 2, 0., 0.005 );
  lamFit_SG_highbg->SetParameter( 3, 0.01*lamHist->GetMaximum() );
  lamFit_SG_highbg->SetParameter( 4, 0.4 );
  lamFit_SG_highbg->SetParameter( 5, 1.0 );
  lamFit_SG_highbg->SetParameter( 6, -1.0 );

  // Performing initial SG fit to obtain starting values for DG
  //cout << "Performing initial SG fit..." << endl;
  int fit0 = lamHist->Fit("lamFit_SG_highbg", "QLEM");
  int fit0 = lamHist->Fit("lamFit_SG_highbg", "QLEM");
  int fit0 = lamHist->Fit("lamFit_SG_highbg", "QLEM");
  int fit0 = lamHist->Fit("lamFit_SG_highbg", fitopt);
  if (fit0 != 4000) int fit0 = lamHist->Fit("lamFit_SG_highbg", fitopt);
  if (fit0 != 4000) int fit0 = lamHist->Fit("lamFit_SG_highbg", fitopt);
  if (fit0 != 4000) int fit0 = lamHist->Fit("lamFit_SG_highbg", fitopt);

  sigmae = lamFit_SG_highbg->GetParError(2);
  int fitstat = -1;
  if ( !isSG ) {
    lamFit_DG_highbg->SetParameter( 0, lamFit_SG_highbg->GetParameter(0) );
    lamFit_DG_highbg->SetParameter( 1, lamFit_SG_highbg->GetParameter(1) );
    lamFit_DG_highbg->SetParameter( 2, 0.75*lamFit_SG_highbg->GetParameter(2) );
    lamFit_DG_highbg->SetParLimits( 2, 0., 0.004 );
    lamFit_DG_highbg->SetParameter( 3, 0.5 );
    lamFit_DG_highbg->SetParameter( 4, min(0.0045,1.3*lamFit_SG_highbg->GetParameter(2)) );
    lamFit_DG_highbg->SetParLimits( 4, 0., 0.008 );
    lamFit_DG_highbg->SetParameter( 5, lamFit_SG_highbg->GetParameter(3) );
    lamFit_DG_highbg->SetParameter( 6, lamFit_SG_highbg->GetParameter(4) );
    lamFit_DG_highbg->SetParameter( 7, lamFit_SG_highbg->GetParameter(5) );

    fitstat = lamHist->Fit("lamFit_DG_highbg", "QLEM");
    fitstat = lamHist->Fit("lamFit_DG_highbg", "QLEM");
    fitstat = lamHist->Fit("lamFit_DG_highbg", "QLEM");
    fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
    if (fitstat != 4000) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);

    bool refit = false;
    if( fitstat == 4000 ) {
      double sigma1 = lamFit_DG_highbg->GetParameter(2);
      double sigma2 = lamFit_DG_highbg->GetParameter(4);
      if( sigma1 > sigma2 ) {
	lamFit_DG_highbg->SetParameter(2, sigma2);
	lamFit_DG_highbg->SetParameter(4, sigma1);
	double frac = lamFit_DG_highbg->GetParameter(3);
	lamFit_DG_highbg->SetParameter(3, 1.-frac);
	refit = true;
      }
    } else {
      refit = true;
    }
    if ( refit ) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
    if ( fitstat != 4000 ) fitstat = lamHist->Fit("lamFit_DG_highbg", fitopt);
  }

  double height = 0; 
  double heighte = 0;
  sigma = 0;
  yield = 0; yielde = 0;
  if( fitstat >= 0 ) {
    //cout << "Good double gaussian fit, getting parameters." << endl;
    double frac = lamFit_DG_highbg->GetParameter(3);
    double sigma1 = lamFit_DG_highbg->GetParameter(2);
    double sigma2 = lamFit_DG_highbg->GetParameter(4);
    mean = lamFit_DG_highbg->GetParameter(1);
    meane = lamFit_DG_highbg->GetParError(1);
    sigma = sqrt( frac*sigma1*sigma1 + (1.0-frac)*sigma2*sigma2 );
    if (sigmae <= 0) sigmae = frac*lamFit_DG_highbg->GetParError(2) + (1.0-frac)*lamFit_DG_highbg->GetParError(2);
    height = lamFit_DG_highbg->GetParameter(5);
    heighte = lamFit_DG_highbg->GetParError(5);
    yield = lamFit_DG_highbg->GetParameter(0);
    yielde = lamFit_DG_highbg->GetParError(0);
    fittype = 3;
    chi2 = lamFit_DG_highbg->GetChisquare()/lamFit_DG_highbg->GetNDF();
  } 
  else {
    //cout << "WARNING: Bad double gaussian fit, refitting with a single gaussian." << endl;
    fitstat = lamHist->Fit("lamFit_SG_highbg", fitopt);
    mean = lamFit_SG_highbg->GetParameter(1);
    meane = lamFit_SG_highbg->GetParError(1);
    sigma = lamFit_SG_highbg->GetParameter(2);
    sigmae = lamFit_SG_highbg->GetParError(2);
    height = lamFit_SG_highbg->GetParameter(3);
    heighte = lamFit_SG_highbg->GetParError(3);
    yield = lamFit_SG_highbg->GetParameter(0);
    yielde = lamFit_SG_highbg->GetParError(0);
    fittype = -3;
    chi2 = lamFit_SG_highbg->GetChisquare()/lamFit_SG_highbg->GetNDF();
    isSG = true;
  }

  background = 4.0*height*sigma/lamMassBinWidth;
  backgrounde = background*heighte/height;

  return fitstat;
}

#endif
