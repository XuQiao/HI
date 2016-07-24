#include <string>
#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TObject.h"

class Corrector
{
  public:
    Corrector();
    Corrector(TString fileName);
    double getEventWeightNtrk( int M);
    double getEventWeightHFPlus4( double HFEn, bool isRebin);
   
virtual ~Corrector();

  private:
	TFile* f_corr;
	TH1F* h_corrNtrk; 
	TH1F* h_corrHFPlus4; 
	TH1F* reh_corrHFPlus4; 
};

Corrector::Corrector(){
f_corr=new TFile("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/official/output/weight/pPbHist_Hijing_DS.root");

 h_corrNtrk = (TH1F*) f_corr->Get("hEff_Ntrk");
 h_corrHFPlus4 = (TH1F*) f_corr->Get("hEff_HFEnergyPlus4");
 reh_corrHFPlus4 = (TH1F*) f_corr->Get("rehEff_HFEnergyPlus4");
}

Corrector::Corrector(TString fileName )
{
 f_corr = new TFile(fileName);
 h_corrNtrk = (TH1F*) f_corr->Get("hEff_Ntrk");
 h_corrHFPlus4 = (TH1F*) f_corr->Get("hEff_HFEnergyPlus4");
 reh_corrHFPlus4 = (TH1F*) f_corr->Get("rehEff_HFEnergyPlus4");

}

Corrector::~Corrector(){
delete f_corr;
}

double Corrector::getEventWeightNtrk(int M){
if(M<1 || M>29) return 1;
double corrfac = 1./h_corrNtrk->GetBinContent(h_corrNtrk->GetXaxis()->FindBin(M));
return corrfac;
}

double Corrector::getEventWeightHFPlus4(double HFEn, bool isRebin){
if(HFEn <= 0 || HFEn >= 20.) return 1;
double corrfac;
if(isRebin){
if(reh_corrHFPlus4->GetBinContent(reh_corrHFPlus4->GetXaxis()->FindBin(HFEn))==0) return 1;
corrfac = 1./reh_corrHFPlus4->GetBinContent(reh_corrHFPlus4->GetXaxis()->FindBin(HFEn));}
else{
if(h_corrHFPlus4->GetBinContent(h_corrHFPlus4->GetXaxis()->FindBin(HFEn))==0) return 1;
corrfac = 1./h_corrHFPlus4->GetBinContent(h_corrHFPlus4->GetXaxis()->FindBin(HFEn));}
return corrfac;
}
