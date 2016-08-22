#ifndef lamfit_h
#define lamfit_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

  // Constants
  const double piMass = 0.13957018;
  const double protonMass = 0.93827203;
  const double ksMassConst = 0.497614;
  const double lamMassConst = 1.115683;

class lamfit {
public :
  double yield, yielde;
  double background, backgrounde;
  double sigma, sigmae;
  double mean, meane;
  double norm_chi2;
  int fittype;

  TF1* lamFit_SG;
  TF1* lamFit_SG_exp;
  TF1* lamFit_SG_highbg;
  TF1* lamFit_DG;
  TF1* lamFit_DG_exp;
  TF1* lamFit_DG_highbg;
  TF1* lamFit_rootbg;
  TF1* lamFit_expbg;
  TF1* lamFit_expppbg;

  //TH1* lamHist;

  //lamfit(TH1* lamHist);
  lamfit();
  virtual ~lamfit();
  int Fit(TH1* lamHist, int fitcomm = 0);
  int Plot(TH1* lamHist);
  double sig();
  double sig_err();
  double bg();
  double bg_err();
  double res();
  double res_err();
  double mass();
  double mass_err();
  double nChi2();
  int fitbg();

 private:
  int FitWithExpBG(TH1* lamHist, double &expChi2, bool &expIsSG);
  int FitWithConstBG(TH1* lamHist, double &constChi2, bool &constIsSG);
  int FitWithHighbgBG(TH1* lamHist, double &highbgChi2, bool &highbgIsSG);
};

#endif

#ifdef lamfit_cxx
lamfit::lamfit() {

  Double_t lamMassXmin = 1.08;
  Double_t lamMassXmax = 1.18;
  int lamMassNbins = 200;
  Double_t lamMassBinWidth = (lamMassXmax - lamMassXmin) / lamMassNbins;
  Double_t lamMassBinWidth_DG = lamMassBinWidth / sqrt(2*3.1415926536);

  ostringstream lamfitoss;
   lamfitoss << lamMassBinWidth << "*gausn(0) + "
	    << "[3]*sqrt(x - "  << piMass + protonMass << ") + "
	    << "[4]*(x - "  << piMass + protonMass << ")^1.5";

  ostringstream lamfitoss_exp;
  lamfitoss_exp << lamMassBinWidth << "*gausn(0) + "
	    << "[3]*(x - "  << piMass + protonMass << ")^[4]";

  ostringstream lamfitoss_highbg;
  lamfitoss_highbg << lamMassBinWidth << "*gausn(0) + [3]*(x-" << piMass+protonMass<<")^[4]*" 
		   << "(exp([5]*(x-"<<piMass+protonMass<<")+[6]*(x-"<<piMass+protonMass<<")^2))";

  ostringstream lamfitoss_rootbg;
   lamfitoss_rootbg << "[0]*sqrt(x - "  << piMass + protonMass << ") + "
	    << "[1]*(x - "  << piMass + protonMass << ")^1.5";

  ostringstream lamfitoss_expbg;
  lamfitoss_expbg << "[0]*(x - "  << piMass + protonMass << ")^[1]";

  ostringstream lamfitoss_expppbg;
  lamfitoss_expppbg << "[0]*(x-" << piMass+protonMass<<")^[1]*" 
		   << "(exp([2]*(x-"<<piMass+protonMass<<")+[3]*(x-"<<piMass+protonMass<<")^2))";

  ostringstream lamfitoss_DG;
  lamfitoss_DG << lamMassBinWidth_DG
	       << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
	       << lamMassBinWidth_DG 
	       << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
	       << "[5]*sqrt(x - "  << piMass + protonMass << ") + "
	       << "[6]*(x - "  << piMass + protonMass << ")^1.5";
  ostringstream lamfitoss_DG_exp;
  lamfitoss_DG_exp << lamMassBinWidth_DG
		   << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		   << lamMassBinWidth_DG 
		   << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		   << "[5]*(x - "  << piMass + protonMass << ")^[6]";
  ostringstream lamfitoss_DG_highbg;
  lamfitoss_DG_highbg << lamMassBinWidth_DG
		      << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		      << lamMassBinWidth_DG 
		      << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		      << "[5]*(x-" << piMass+protonMass<<")^[6]*" 
		      << "(exp([7]*(x-"<<piMass+protonMass<<")+[8]*(x-"<<piMass+protonMass<<")^2))";

  lamFit_SG = new TF1("lamFit_SG", 
		     lamfitoss.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_SG->SetParName( 0, "yield" );
  lamFit_SG->SetParName( 1, "mean" );
  lamFit_SG->SetParName( 2, "sigma" );
  lamFit_SG->SetParName( 3, "sqrt" );
  lamFit_SG->SetParName( 4, "cubesqrt" );
  lamFit_SG->SetLineColor(2);
  lamFit_SG->SetLineWidth(1);

  lamFit_SG_exp = new TF1("lamFit_SG_exp", 
			 lamfitoss_exp.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_SG_exp->SetParName( 0, "yield" );
  lamFit_SG_exp->SetParName( 1, "mean" );
  lamFit_SG_exp->SetParName( 2, "sigma" );
  lamFit_SG_exp->SetParName( 3, "amp" );
  lamFit_SG_exp->SetParName( 4, "exp" );
  lamFit_SG_exp->SetLineColor(2);
  lamFit_SG_exp->SetLineWidth(1);

  lamFit_SG_highbg = new TF1("lamFit_SG_highbg", 
			  lamfitoss_highbg.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_SG_highbg->SetParName( 0, "yield" );
  lamFit_SG_highbg->SetParName( 1, "mean" );
  lamFit_SG_highbg->SetParName( 2, "sigma" );
  lamFit_SG_highbg->SetParName( 3, "amp" );
  lamFit_SG_highbg->SetParName( 4, "exp" );
  lamFit_SG_highbg->SetParName( 5, "expp1" );
  lamFit_SG_highbg->SetParName( 6, "expp2" );
  lamFit_SG_highbg->SetLineColor(2);
  lamFit_SG_highbg->SetLineWidth(1);

  lamFit_DG = new TF1("lamFit_DG",
		     lamfitoss_DG.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_DG->SetParName( 0, "yield" );
  lamFit_DG->SetParName( 1, "mean" );
  lamFit_DG->SetParName( 2, "sigma1" );
  lamFit_DG->SetParName( 3, "fraction" );
  lamFit_DG->SetParLimits( 3, 0., 1. );
  lamFit_DG->SetParName( 4, "sigma2" );
  lamFit_DG->SetParName( 5, "sqrt" );
  lamFit_DG->SetParName( 6, "cubesqrt" );
  lamFit_DG->SetLineColor(2);
  lamFit_DG->SetLineWidth(1);

  lamFit_DG_exp = new TF1("lamFit_DG_exp",
			 lamfitoss_DG_exp.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_DG_exp->SetParName( 0, "yield" );
  lamFit_DG_exp->SetParName( 1, "mean" );
  lamFit_DG_exp->SetParName( 2, "sigma1" );
  lamFit_DG_exp->SetParName( 3, "fraction" );
  lamFit_DG_exp->SetParLimits( 3, 0., 1. );
  lamFit_DG_exp->SetParName( 4, "sigma2" );
  lamFit_DG_exp->SetParName( 5, "amp" );
  lamFit_DG_exp->SetParName( 6, "exp" );
  lamFit_DG_exp->SetLineColor(2);
  lamFit_DG_exp->SetLineWidth(1);

  lamFit_DG_highbg = new TF1("lamFit_DG_highbg", 
			  lamfitoss_DG_highbg.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_DG_highbg->SetParName( 0, "yield" );
  lamFit_DG_highbg->SetParName( 1, "mean" );
  lamFit_DG_highbg->SetParName( 2, "sigma1" );
  lamFit_DG_highbg->SetParName( 3, "fraction" );
  lamFit_DG_highbg->SetParLimits( 3, 0., 1. );
  lamFit_DG_highbg->SetParName( 4, "sigma2" );
  lamFit_DG_highbg->SetParName( 5, "amp" );
  lamFit_DG_highbg->SetParName( 6, "exp" );
  lamFit_DG_highbg->SetParName( 7, "expp1" );
  lamFit_DG_highbg->SetParName( 8, "expp2" );
  lamFit_DG_highbg->SetLineColor(2);
  lamFit_DG_highbg->SetLineWidth(1);


  lamFit_rootbg = new TF1("lamFit_rootbg", 
		     lamfitoss_rootbg.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_rootbg->SetParName( 0, "sqrt" );
  lamFit_rootbg->SetParName( 1, "cubesqrt" );
  lamFit_rootbg->SetLineColor(2);
  lamFit_rootbg->SetLineWidth(1);

  lamFit_expbg = new TF1("lamFit_expbg", 
			 lamfitoss_expbg.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_expbg->SetParName( 0, "amp" );
  lamFit_expbg->SetParName( 1, "exp" );
  lamFit_expbg->SetLineColor(2);
  lamFit_expbg->SetLineWidth(1);

  lamFit_expppbg = new TF1("lamFit_expppbg", 
			  lamfitoss_expppbg.str().c_str(), lamMassXmin, lamMassXmax);
  lamFit_expppbg->SetParName( 0, "amp" );
  lamFit_expppbg->SetParName( 1, "exp" );
  lamFit_expppbg->SetParName( 2, "expp1" );
  lamFit_expppbg->SetParName( 3, "expp2" );
  lamFit_expppbg->SetLineColor(2);
  lamFit_expppbg->SetLineWidth(1);

}

lamfit::~lamfit()
{
}

double lamfit::sig(){ return yield; }
double lamfit::sig_err(){ return yielde; }
double lamfit::bg(){ return background; }
double lamfit::bg_err(){ return backgrounde; }
double lamfit::res(){ return sigma; }
double lamfit::res_err(){ return sigmae; }
double lamfit::mass(){ return mean; }
double lamfit::mass_err(){ return meane; }
double lamfit::nChi2(){ return norm_chi2; }
int lamfit::fitbg(){ return fittype; }

#endif
