#ifndef ksfit_cxx
#define ksfit_cxx
#include "ksfit.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

const double ksMassConst = 0.497614;

int ksfit::Plot(TH1* ksHist) {

    ksHist->Draw();
    if (fittype == 1) {
      ksFit_expbg->SetParameter( 0, ksFit_DG_exp->GetParameter(1) );
      ksFit_expbg->SetParameter( 1, ksFit_DG_exp->GetParameter(5) );
      ksFit_expbg->SetParameter( 2, ksFit_DG_exp->GetParameter(6) );
      ksFit_expbg->Draw("same");
    } else if ( fittype == 2 ) {
      ksFit_linbg->SetParameter( 0, ksFit_DG->GetParameter(1) );
      ksFit_linbg->SetParameter( 1, ksFit_DG->GetParameter(5) );
      ksFit_linbg->SetParameter( 2, ksFit_DG->GetParameter(6) );
      ksFit_linbg->Draw("same");
    } else if ( fittype == 3 ) {
      ksFit_quadbg->SetParameter( 0, ksFit_DG_quad->GetParameter(1) );
      ksFit_quadbg->SetParameter( 1, ksFit_DG_quad->GetParameter(5) );
      ksFit_quadbg->SetParameter( 2, ksFit_DG_quad->GetParameter(6) );
      ksFit_quadbg->SetParameter( 3, ksFit_DG_quad->GetParameter(7) );
      ksFit_quadbg->Draw("same");
    }
    labl->Draw("same");
}

int ksfit::Fit(TH1* ksHist, int fitcomm) {

  // These chi2 values will be chi2/ndof
  double expChi2 = 0., quadChi2 = 0., constChi2 = 0.;
  bool expIsSG = false, quadIsSG = false, constIsSG = false;
  if ( fitcomm == 1 ) {
    expIsSG = true;
    quadIsSG = true;
    constIsSG = true;
  }
  int expFitStat = FitWithExpBG( ksHist, expChi2 , expIsSG, "QLEM" );
  int quadFitStat = FitWithQuadBG( ksHist, quadChi2, quadIsSG, "QLEM" );
  int constFitStat = FitWithConstBG( ksHist, constChi2, constIsSG, "QLEM" );

  // We'll only check the chi2 for those that had a good double gaussian fit,
  //   unless none did.
  bool checkAll = !(expIsSG && quadIsSG && constIsSG) 
    || (expIsSG && quadIsSG && constIsSG);
  bool checkExpAndQuad = !(expIsSG && quadIsSG && !constIsSG);
  bool checkExpAndConst = !(expIsSG && !quadIsSG && constIsSG);
  bool checkQuadAndConst = !(!expIsSG && quadIsSG && constIsSG);
  bool checkOnlyExp = !(expIsSG && !quadIsSG && !constIsSG);
  bool checkOnlyQuad = !(!expIsSG && quadIsSG && !constIsSG);
  bool checkOnlyConst = !(!expIsSG && !quadIsSG && constIsSG);

  // I admit, having to refit again is not ideal, but we need to fill the
  //  signal/noise values, etc. using the good fit.
  if( checkAll ) {
    if( expChi2 < quadChi2 && expChi2 < constChi2 ) {
      norm_chi2 = expChi2;
      return FitWithExpBG( ksHist, expChi2, expIsSG, "LEM" );
    }
    else if( quadChi2 < expChi2 && quadChi2 < constChi2 ) {
      norm_chi2 = quadChi2;
      return FitWithQuadBG( ksHist, quadChi2, quadIsSG, "LEM" );
    }
    else if( constChi2 < expChi2 && constChi2 < quadChi2 ) {
      norm_chi2 = constChi2;
      return FitWithConstBG( ksHist, constChi2, constIsSG, "LEM" );
    }
  }
  else if( checkExpAndQuad ) {
    if( expChi2 < quadChi2 ) {
      norm_chi2 = expChi2;
      return FitWithExpBG( ksHist, expChi2, expIsSG, "LEM" );
    }
    else {
      norm_chi2 = quadChi2;
      return FitWithQuadBG( ksHist, expChi2, expIsSG, "LEM" );
    }
  }
  else if( checkExpAndConst ) {
    if( expChi2 < constChi2 ) {
      norm_chi2 = expChi2;
      return FitWithExpBG( ksHist, expChi2, expIsSG, "LEM" );
    }
    else {
      // Const was the last to run, so just spit out the FitStat
      norm_chi2 = constChi2;
      return constFitStat;
      //return FitWithConstBG( ksHist, constChi2, constIsSG, "LEM" );
    }
  }
  else if( checkQuadAndConst ) {
    if( quadChi2 < constChi2 ) {
      norm_chi2 = quadChi2;
      return FitWithQuadBG( ksHist, quadChi2, quadIsSG, "LEM" );
    }
    else {
      norm_chi2 = constChi2;
      return constFitStat;
      //return FitWithConstBG( ksHist, constChi2, constIsSG, "LEM" );
    }
  }
  else if (checkOnlyExp ) {
    norm_chi2 = expChi2;
    return FitWithExpBG( ksHist, expChi2, expIsSG, "LEM" );
  }
  else if ( checkOnlyQuad ) {
    norm_chi2 = quadChi2;
    return FitWithQuadBG( ksHist, quadChi2, quadIsSG, "LEM" );
  }
  else {
    norm_chi2 = constChi2;
    return constFitStat;
  }

}


int ksfit::FitWithExpBG(TH1* ksHist, double &chi2, bool &isSG, TString fitopt){

  //cout << "ExpBG: setting yield, mass, bg = " << 4.5*ksHist->GetMaximum() << ", " << ksMassConst << ", " << 0.1*ksHist->GetMaximum() << endl;
  ksFit_SG_exp->SetParameter( 0, 4.5*ksHist->GetMaximum() );
  ksFit_SG_exp->SetParameter( 1, ksMassConst );
  ksFit_SG_exp->SetParameter( 2, 0.005 );
  ksFit_SG_exp->SetParLimits( 2, 0., 0.020 );
  ksFit_SG_exp->SetParameter( 3, 0.1*ksHist->GetMaximum() );
  ksFit_SG_exp->SetParameter( 4, -1. );

  // Performing initial SG fit to obtain starting values for DG
  //cout << "Performing initial SG fit..." << endl;
  int fit0 = ksHist->Fit("ksFit_SG_exp", "QLEM");
  fit0 = ksHist->Fit("ksFit_SG_exp", "QLEM");
  fit0 = ksHist->Fit("ksFit_SG_exp", fitopt);

  sigmae = ksFit_SG_exp->GetParError(2);
  int fitstat = -1;
  if ( !isSG ) {
    ksFit_DG_exp->SetParameter( 0, 1.018*ksFit_SG_exp->GetParameter(0) );
    ksFit_DG_exp->SetParameter( 1, ksFit_SG_exp->GetParameter(1) );
    ksFit_DG_exp->SetParameter( 2, 0.85*ksFit_SG_exp->GetParameter(2) );
    ksFit_DG_exp->SetParLimits( 2, 0., 0.010 );
    ksFit_DG_exp->SetParameter( 3, 0.75 );
    ksFit_DG_exp->SetParameter( 4, 1.7*ksFit_SG_exp->GetParameter(2) );
    ksFit_DG_exp->SetParLimits( 4, 0., 0.020 );
    ksFit_DG_exp->SetParameter( 5, ksFit_SG_exp->GetParameter(3) );
    ksFit_DG_exp->SetParameter( 6, ksFit_SG_exp->GetParameter(4) );

    fitstat = ksHist->Fit("ksFit_DG_exp", "QLEM");
    fitstat = ksHist->Fit("ksFit_DG_exp", "QLEM");
    fitstat = ksHist->Fit("ksFit_DG_exp", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG_exp", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG_exp", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG_exp", fitopt);

    bool refit = false;
    if( fitstat == 4000 ) {
      double sigma1 = ksFit_DG_exp->GetParameter(2);
      double sigma2 = ksFit_DG_exp->GetParameter(4);
      if( sigma1 > sigma2 ) {
	ksFit_DG_exp->SetParameter(2, sigma2);
	ksFit_DG_exp->SetParameter(4, sigma1);
	double frac = ksFit_DG_exp->GetParameter(3);
	ksFit_DG_exp->SetParameter(3, 1.-frac);
	refit = true;
      }
    } else {
      refit = true;
    }
    if ( refit ) fitstat = ksHist->Fit("ksFit_DG_exp", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG_exp", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG_exp", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG_exp", fitopt);
  }
  double height = 0; 
  double heighte = 0;
  sigma = 0;
  yield = 0; yielde = 0;
  if( fitstat == 4000 ) {
    //cout << "Good double gaussian fit, getting parameters." << endl;
    double frac = ksFit_DG_exp->GetParameter(3);
    double sigma1 = ksFit_DG_exp->GetParameter(2);
    double sigma2 = ksFit_DG_exp->GetParameter(4);
    mean = ksFit_DG_exp->GetParameter(1);
    meane = ksFit_DG_exp->GetParError(1);
    sigma = sqrt( frac*sigma1*sigma1 + (1.0-frac)*sigma2*sigma2 );
    if (sigmae < 0) sigmae = frac*ksFit_DG_exp->GetParError(2) + (1.0-frac)*ksFit_DG_exp->GetParError(2);
    height = ksFit_DG_exp->GetParameter(5);
    heighte = ksFit_DG_exp->GetParError(5);
    yield = ksFit_DG_exp->GetParameter(0);
    yielde = ksFit_DG_exp->GetParError(0);
    fittype = 1;
    chi2 = ksFit_DG_exp->GetChisquare()/ksFit_DG_exp->GetNDF();
  } 
  else {
    //cout << "WARNING: Bad double gaussian fit, refitting with a single gaussian." << endl;
    fitstat = ksHist->Fit("ksFit_SG_exp", fitopt);
    mean = ksFit_SG_exp->GetParameter(1);
    meane = ksFit_SG_exp->GetParError(1);
    sigma = ksFit_SG_exp->GetParameter(2);
    sigmae = ksFit_SG_exp->GetParError(2);
    height = ksFit_SG_exp->GetParameter(3);
    heighte = ksFit_SG_exp->GetParError(3);
    yield = ksFit_SG_exp->GetParameter(0);
    yielde = ksFit_SG_exp->GetParError(0);
    fittype = -1;
    chi2 = ksFit_SG_exp->GetChisquare()/ksFit_SG_exp->GetNDF();
    isSG = true;
  }

  background = 4.0*height*sigma/ksMassBinWidth;
  backgrounde = background*heighte/height;
  //cout << "Returning from ksfit::Fit " << fitstat << endl;

  return fitstat;
}




int ksfit::FitWithConstBG(TH1* ksHist, double &chi2, bool &isSG, TString fitopt){

  //cout << "ConstBG: setting yield, mass, bg = " << 4.5*ksHist->GetMaximum() << ", " << ksMassConst << ", " << 0.1*ksHist->GetMaximum() << endl;
  ksFit_SG->SetParameter( 0, 4.5*ksHist->GetMaximum() );
  ksFit_SG->SetParameter( 1, ksMassConst );
  ksFit_SG->SetParameter( 2, 0.005 );
  ksFit_SG->SetParLimits( 2, 0., 0.020 );
  ksFit_SG->SetParameter( 3, 0.1*ksHist->GetMaximum() );
  ksFit_SG->SetParameter( 4, 1. );

  // Performing initial SG fit to obtain starting values for DG
  //cout << "Performing initial SG fit..." << endl;
  int fit0 = ksHist->Fit("ksFit_SG", "QLEM");
  fit0 = ksHist->Fit("ksFit_SG", "QLEM");
  fit0 = ksHist->Fit("ksFit_SG", fitopt);

  sigmae = ksFit_SG->GetParError(2);
  int fitstat = -1;
  if ( !isSG ) {
    ksFit_DG->SetParameter( 0, 1.018*ksFit_SG->GetParameter(0) );
    ksFit_DG->SetParameter( 1, ksFit_SG->GetParameter(1) );
    ksFit_DG->SetParameter( 2, 0.85*ksFit_SG->GetParameter(2) );
    ksFit_DG->SetParLimits( 2, 0., 0.010 );
    ksFit_DG->SetParameter( 3, 0.75 );
    ksFit_DG->SetParameter( 4, 1.7*ksFit_SG->GetParameter(2) );
    ksFit_DG->SetParLimits( 4, 0., 0.020 );
    ksFit_DG->SetParameter( 5, ksFit_SG->GetParameter(3) );
    ksFit_DG->SetParameter( 6, ksFit_SG->GetParameter(4) );

    fitstat = ksHist->Fit("ksFit_DG", "QLEM");
    fitstat = ksHist->Fit("ksFit_DG", "QLEM");
    fitstat = ksHist->Fit("ksFit_DG", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG", fitopt);

    bool refit = false;
    if( fitstat == 4000 ) {
      double sigma1 = ksFit_DG->GetParameter(2);
      double sigma2 = ksFit_DG->GetParameter(4);
      if( sigma1 > sigma2 ) {
	ksFit_DG->SetParameter(2, sigma2);
	ksFit_DG->SetParameter(4, sigma1);
	double frac = ksFit_DG->GetParameter(3);
	ksFit_DG->SetParameter(3, 1.-frac);
	refit = true;
      }
    } else {
      refit = true;
    }
    if ( refit ) fitstat = ksHist->Fit("ksFit_DG", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG", fitopt);
  }
  
  double height = 0; 
  double heighte = 0;
  sigma = 0;
  yield = 0; yielde = 0;
  if( fitstat == 4000 ) {
    //cout << "Good double gaussian fit, getting parameters." << endl;
    double frac = ksFit_DG->GetParameter(3);
    double sigma1 = ksFit_DG->GetParameter(2);
    double sigma2 = ksFit_DG->GetParameter(4);
    mean = ksFit_DG->GetParameter(1);
    meane = ksFit_DG->GetParError(1);
    sigma = sqrt( frac*sigma1*sigma1 + (1.0-frac)*sigma2*sigma2 );
    if (sigmae < 0) sigmae = frac*ksFit_DG->GetParError(2) + (1.0-frac)*ksFit_DG->GetParError(2);
    height = ksFit_DG->GetParameter(5);
    heighte = ksFit_DG->GetParError(5);
    yield = ksFit_DG->GetParameter(0);
    yielde = ksFit_DG->GetParError(0);
    fittype = 2;
    chi2 = ksFit_DG->GetChisquare()/ksFit_DG->GetNDF();
  } 
  else {
    //cout << "WARNING: Bad double gaussian fit, refitting with a single gaussian." << endl;
    fitstat = ksHist->Fit("ksFit_SG", fitopt);
    mean = ksFit_SG->GetParameter(1);
    meane = ksFit_SG->GetParError(1);
    sigma = ksFit_SG->GetParameter(2);
    sigmae = ksFit_SG->GetParError(2);
    height = ksFit_SG->GetParameter(3);
    heighte = ksFit_SG->GetParError(3);
    yield = ksFit_SG->GetParameter(0);
    yielde = ksFit_SG->GetParError(0);
    fittype = -2;
    chi2 = ksFit_SG->GetChisquare()/ksFit_SG->GetNDF();
    isSG = true;
  }

  background = 4.0*height*sigma/ksMassBinWidth;
  backgrounde = background*heighte/height;
  //cout << "Returning from ksfit::Fit " << fitstat << endl;

  return fitstat;
}

int ksfit::FitWithQuadBG(TH1* ksHist, double &chi2, bool &isSG, TString fitopt){

  //cout << "QuadBG: setting yield, mass, bg = " << 4.5*ksHist->GetMaximum() << ", " << ksMassConst << ", " << 0.1*ksHist->GetMaximum() << endl;
  ksFit_SG_quad->SetParameter( 0, 4.5*ksHist->GetMaximum() );
  ksFit_SG_quad->SetParameter( 1, ksMassConst );
  ksFit_SG_quad->SetParameter( 2, 0.005 );
  ksFit_SG_quad->SetParLimits( 2, 0., 0.020 );
  ksFit_SG_quad->SetParameter( 3, 0.1*ksHist->GetMaximum() );
  ksFit_SG_quad->SetParameter( 4, 1. );
  ksFit_SG_quad->SetParameter( 5, -1. );

  // Performing initial SG fit to obtain starting values for DG
  //cout << "Performing initial SG fit..." << endl;
  int fit0 = ksHist->Fit("ksFit_SG_quad", "QLEM");
  fit0 = ksHist->Fit("ksFit_SG_quad", "QLEM");
  fit0 = ksHist->Fit("ksFit_SG_quad", "QLEM");
  fit0 = ksHist->Fit("ksFit_SG_quad", fitopt);

  sigmae = ksFit_SG_quad->GetParError(2);
  int fitstat = -1;
  if ( !isSG ) {
    ksFit_DG_quad->SetParameter( 0, 1.018*ksFit_SG_quad->GetParameter(0) );
    ksFit_DG_quad->SetParameter( 1, ksFit_SG_quad->GetParameter(1) );
    ksFit_DG_quad->SetParameter( 2, 0.85*ksFit_SG_quad->GetParameter(2) );
    ksFit_DG_quad->SetParLimits( 2, 0., 0.010 );
    ksFit_DG_quad->SetParameter( 3, 0.75 );
    ksFit_DG_quad->SetParameter( 4, 1.7*ksFit_SG_quad->GetParameter(2) );
    ksFit_DG_quad->SetParLimits( 4, 0., 0.020 );
    ksFit_DG_quad->SetParameter( 5, ksFit_SG_quad->GetParameter(3) );
    ksFit_DG_quad->SetParameter( 6, ksFit_SG_quad->GetParameter(4) );
    ksFit_DG_quad->SetParameter( 7, ksFit_SG_quad->GetParameter(5) );

    fitstat = ksHist->Fit("ksFit_DG_quad", "QLEM");
    fitstat = ksHist->Fit("ksFit_DG_quad", "QLEM");
    fitstat = ksHist->Fit("ksFit_DG_quad", "QLEM");
    fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
    if (fitstat != 4000) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);

    bool refit = false;
    if( fitstat == 4000 ) {
      double sigma1 = ksFit_DG_quad->GetParameter(2);
      double sigma2 = ksFit_DG_quad->GetParameter(4);
      if( sigma1 > sigma2 ) {
	ksFit_DG_quad->SetParameter(2, sigma2);
	ksFit_DG_quad->SetParameter(4, sigma1);
	double frac = ksFit_DG_quad->GetParameter(3);
	ksFit_DG_quad->SetParameter(3, 1.-frac);
	refit = true;
      }
    } else {
      refit = true;
    }
    if ( refit ) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
    if ( fitstat != 4000 ) fitstat = ksHist->Fit("ksFit_DG_quad", fitopt);
  }
  
  double height = 0; 
  double heighte = 0;
  sigma = 0;
  yield = 0; yielde = 0;
  if( fitstat == 4000 ) {
    //cout << "Good double gaussian fit, getting parameters." << endl;
    double frac = ksFit_DG_quad->GetParameter(3);
    double sigma1 = ksFit_DG_quad->GetParameter(2);
    double sigma2 = ksFit_DG_quad->GetParameter(4);
    mean = ksFit_DG_quad->GetParameter(1);
    meane = ksFit_DG_quad->GetParError(1);
    sigma = sqrt( frac*sigma1*sigma1 + (1.0-frac)*sigma2*sigma2 );
    if (sigmae < 0) sigmae = frac*ksFit_DG_quad->GetParError(2) + (1.0-frac)*ksFit_DG_quad->GetParError(2);
    height = ksFit_DG_quad->GetParameter(5);
    heighte = ksFit_DG_quad->GetParError(5);
    yield = ksFit_DG_quad->GetParameter(0);
    yielde = ksFit_DG_quad->GetParError(0);
    fittype = 3;
    chi2 = ksFit_DG_quad->GetChisquare()/ksFit_DG_quad->GetNDF();
  } 
  else {
    //cout << "WARNING: Bad double gaussian fit, refitting with a single gaussian." << endl;
    fitstat = ksHist->Fit("ksFit_SG_quad", fitopt);
    mean = ksFit_SG_quad->GetParameter(1);
    meane = ksFit_SG_quad->GetParError(1);
    sigma = ksFit_SG_quad->GetParameter(2);
    sigmae = ksFit_SG_quad->GetParError(2);
    height = ksFit_SG_quad->GetParameter(3);
    heighte = ksFit_SG_quad->GetParError(3);
    yield = ksFit_SG_quad->GetParameter(0);
    yielde = ksFit_SG_quad->GetParError(0);
    fittype = -3;
    chi2 = ksFit_SG_quad->GetChisquare()/ksFit_SG_quad->GetNDF();
    isSG = true;
  }

  background = 4.0*height*sigma/ksMassBinWidth;
  backgrounde = background*heighte/height;

  return fitstat;
}

#endif
