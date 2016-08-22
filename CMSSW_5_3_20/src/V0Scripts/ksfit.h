#ifndef ksfit_h
#define ksfit_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

class ksfit {
public :
  double yield, yielde;
  double background, backgrounde;
  double sigma, sigmae;
  double mean, meane;
  double norm_chi2;
  int fittype;

  TF1* ksFit_SG;
  TF1* ksFit_SG_exp;
  TF1* ksFit_SG_quad;
  TF1* ksFit_DG;
  TF1* ksFit_DG_exp;
  TF1* ksFit_DG_quad;
  TF1* ksFit_linbg;
  TF1* ksFit_expbg;
  TF1* ksFit_quadbg;

  //TH1* ksHist;

  //ksfit(TH1* ksHist);
  ksfit();
  virtual ~ksfit();
  int Fit(TH1* ksHist, int fitcomm = 0);
  int Plot(TH1* ksHist);
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
  int FitWithExpBG(TH1* ksHist, double &expChi2, bool &expIsSG);
  int FitWithConstBG(TH1* ksHist, double &constChi2, bool &constIsSG);
  int FitWithQuadBG(TH1* ksHist, double &quadChi2, bool &quadIsSG);
};

#endif

#ifdef ksfit_cxx
//ksfit::ksfit(TH1* tmpksHist){
ksfit::ksfit() {

  //ksHist = tmpksHist;
  Double_t ksMassXmin = 0.340;
  Double_t ksMassXmax = 0.740;
  int ksMassNbins = 400;
  Double_t ksMassBinWidth = (ksMassXmax - ksMassXmin) / ksMassNbins;
  Double_t ksMassBinWidth_DG = ksMassBinWidth / sqrt(2*3.1415926536);

  ostringstream ksfitoss;
  ksfitoss << ksMassBinWidth << "*gausn(0) + [3] + [4]*(x-[1])";

  ostringstream ksfitoss_exp;
  ksfitoss_exp << ksMassBinWidth << "*gausn(0) + [3]*exp([4]*(x-[1]))";

  ostringstream ksfitoss_quad;
  ksfitoss_quad << ksMassBinWidth << "*gausn(0) + [3] + [4]*(x-[1]) + [5]*(x-[1])^2";

  ostringstream ksfitoss_linbg;
  ksfitoss_linbg << "[1] + [2]*(x-[0])";

  ostringstream ksfitoss_expbg;
  ksfitoss_expbg << "[1]*exp([2]*(x-[0]))";

  ostringstream ksfitoss_quadbg;
  ksfitoss_quadbg << "[1] + [2]*(x-[0]) + [3]*(x-[0])^2";

  ostringstream ksfitoss_DG;
  ksfitoss_DG << ksMassBinWidth_DG
	      << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
	      << ksMassBinWidth_DG 
	      << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
	      << "[5] + [6]*(x-[1])";

  ostringstream ksfitoss_DG_exp;
  ksfitoss_DG_exp << ksMassBinWidth_DG
		  << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		  << ksMassBinWidth_DG 
		  << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		  << "[5]*exp([6]*(x-[1]))";

  ostringstream ksfitoss_DG_quad;
  ksfitoss_DG_quad << ksMassBinWidth_DG
		   << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		   << ksMassBinWidth_DG 
		   << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		   << "[5] + [6]*(x-[1]) + [7]*(x-[1])^2";

  ksFit_SG = new TF1("ksFit_SG", 
		     ksfitoss.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_SG->SetParName( 0, "yield" );
  ksFit_SG->SetParName( 1, "mean" );
  ksFit_SG->SetParName( 2, "sigma" );
  ksFit_SG->SetParName( 3, "const" );
  ksFit_SG->SetParName( 4, "slope" );
  ksFit_SG->SetLineColor(2);
  ksFit_SG->SetLineWidth(1);

  ksFit_SG_exp = new TF1("ksFit_SG_exp", 
			 ksfitoss_exp.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_SG_exp->SetParName( 0, "yield" );
  ksFit_SG_exp->SetParName( 1, "mean" );
  ksFit_SG_exp->SetParName( 2, "sigma" );
  ksFit_SG_exp->SetParName( 3, "econst" );
  ksFit_SG_exp->SetParName( 4, "eslope" );
  ksFit_SG_exp->SetLineColor(2);
  ksFit_SG_exp->SetLineWidth(1);

  ksFit_SG_quad = new TF1("ksFit_SG_quad", 
			  ksfitoss_quad.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_SG_quad->SetParName( 0, "yield" );
  ksFit_SG_quad->SetParName( 1, "mean" );
  ksFit_SG_quad->SetParName( 2, "sigma" );
  ksFit_SG_quad->SetParName( 3, "const" );
  ksFit_SG_quad->SetParName( 4, "slope" );
  ksFit_SG_quad->SetParName( 5, "quadconst" );
  ksFit_SG_quad->SetLineColor(2);
  ksFit_SG_quad->SetLineWidth(1);

  ksFit_DG = new TF1("ksFit_DG",
		     ksfitoss_DG.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_DG->SetParName( 0, "yield" );
  ksFit_DG->SetParName( 1, "mean" );
  ksFit_DG->SetParName( 2, "sigma1" );
  ksFit_DG->SetParName( 3, "fraction" );
  ksFit_DG->SetParLimits( 3, 0., 1. );
  ksFit_DG->SetParName( 4, "sigma2" );
  ksFit_DG->SetParName( 5, "const" );
  ksFit_DG->SetParName( 6, "slope" );
  ksFit_DG->SetLineColor(2);
  ksFit_DG->SetLineWidth(1);

  ksFit_DG_exp = new TF1("ksFit_DG_exp",
			 ksfitoss_DG_exp.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_DG_exp->SetParName( 0, "yield" );
  ksFit_DG_exp->SetParName( 1, "mean" );
  ksFit_DG_exp->SetParName( 2, "sigma1" );
  ksFit_DG_exp->SetParName( 3, "fraction" );
  ksFit_DG_exp->SetParLimits( 3, 0., 1. );
  ksFit_DG_exp->SetParName( 4, "sigma2" );
  ksFit_DG_exp->SetParName( 5, "econst" );
  ksFit_DG_exp->SetParName( 6, "eslope" );
  ksFit_DG_exp->SetLineColor(2);
  ksFit_DG_exp->SetLineWidth(1);

  ksFit_DG_quad = new TF1("ksFit_DG_quad", 
			  ksfitoss_DG_quad.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_DG_quad->SetParName( 0, "yield" );
  ksFit_DG_quad->SetParName( 1, "mean" );
  ksFit_DG_quad->SetParName( 2, "sigma1" );
  ksFit_DG_quad->SetParName( 3, "fraction" );
  ksFit_DG_quad->SetParLimits( 3, 0., 1. );
  ksFit_DG_quad->SetParName( 4, "sigma2" );
  ksFit_DG_quad->SetParName( 5, "const" );
  ksFit_DG_quad->SetParName( 6, "slope" );
  ksFit_DG_quad->SetParName( 7, "quadconst" );
  ksFit_DG_quad->SetLineColor(2);
  ksFit_DG_quad->SetLineWidth(1);

  ksFit_linbg = new TF1("ksFit_linbg", 
		     ksfitoss_linbg.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_linbg->SetParName( 0, "center" );
  ksFit_linbg->SetParName( 1, "const" );
  ksFit_linbg->SetParName( 2, "slope" );
  ksFit_linbg->SetLineColor(2);
  ksFit_linbg->SetLineStyle(2);
  ksFit_linbg->SetLineWidth(1);

  ksFit_expbg = new TF1("ksFit_expbg", 
			 ksfitoss_expbg.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_expbg->SetParName( 0, "center" );
  ksFit_expbg->SetParName( 1, "econst" );
  ksFit_expbg->SetParName( 2, "eslope" );
  ksFit_expbg->SetLineColor(2);
  ksFit_expbg->SetLineStyle(2);
  ksFit_expbg->SetLineWidth(1);

  ksFit_quadbg = new TF1("ksFit_quadbg", 
			  ksfitoss_quadbg.str().c_str(), ksMassXmin, ksMassXmax);
  ksFit_quadbg->SetParName( 0, "center" );
  ksFit_quadbg->SetParName( 1, "const" );
  ksFit_quadbg->SetParName( 2, "slope" );
  ksFit_quadbg->SetParName( 3, "quadconst" );
  ksFit_quadbg->SetLineColor(2);
  ksFit_quadbg->SetLineStyle(2);
  ksFit_quadbg->SetLineWidth(1);

}

ksfit::~ksfit()
{
}

double ksfit::sig(){ return yield; }
double ksfit::sig_err(){ return yielde; }
double ksfit::bg(){ return background; }
double ksfit::bg_err(){ return backgrounde; }
double ksfit::res(){ return sigma; }
double ksfit::res_err(){ return sigmae; }
double ksfit::mass(){ return mean; }
double ksfit::mass_err(){ return meane; }
double ksfit::nChi2(){ return norm_chi2; }
int ksfit::fitbg(){ return fittype; }

#endif
