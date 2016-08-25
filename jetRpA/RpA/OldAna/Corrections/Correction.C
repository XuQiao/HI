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
    double getEventWeightNtrkbak( int M);
    double getEventWeightHFPlus4( double HFEn, bool isRebin);
    double getEventWeightHFPlus4bak( double HFEn, bool isRebin);
   
virtual ~Corrector();

  private:
	TFile* f_corr;
	TH1D* h_corrNtrk; 
	TH1D* h_corrHFPlus4; 
	TH1D* reh_corrHFPlus4; 
};

Corrector::Corrector(){
f_corr=TFile::Open("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/RecoHI/pPbHijing_EffCorr.root");

h_corrNtrk = (TH1D*) f_corr->Get("hEffMult");
h_corrHFPlus4 = (TH1D*) f_corr->Get("hEffHFp4");
reh_corrHFPlus4 = (TH1D*) f_corr->Get("rehEffHFp4");
}

Corrector::Corrector(TString fileName )
{
 f_corr = TFile::Open(fileName);
 h_corrNtrk = (TH1D*) f_corr->Get("hEffMult");
 h_corrHFPlus4 = (TH1D*) f_corr->Get("hEffHFp4");
 reh_corrHFPlus4 = (TH1D*) f_corr->Get("rehEffHFp4");
}

Corrector::~Corrector(){
//delete f_corr;
//delete h_corrNtrk;
//delete h_corrHFPlus4;
//delete reh_corrHFPlus4;
}

double Corrector::getEventWeightNtrk(int M){
if(M<1 || M>29) return 1;
double corrfac = 1./h_corrNtrk->GetBinContent(h_corrNtrk->GetXaxis()->FindBin(M));
return corrfac;
}

double Corrector::getEventWeightNtrkbak(int M){
//if(M<1) return 1;
if(h_corrNtrk->GetBinContent(h_corrNtrk->GetXaxis()->FindBin(M))==0) return 1;
double corrfac = 1./h_corrNtrk->GetBinContent(h_corrNtrk->GetXaxis()->FindBin(M));
return corrfac;
}

double Corrector::getEventWeightHFPlus4(double HFEn, bool isRebin){
if(HFEn >= 20.) return 1;
double corrfac;
if(isRebin){
if(reh_corrHFPlus4->GetBinContent(reh_corrHFPlus4->GetXaxis()->FindBin(HFEn))==0) return 1;
corrfac = 1./reh_corrHFPlus4->GetBinContent(reh_corrHFPlus4->GetXaxis()->FindBin(HFEn));}
else{
if(h_corrHFPlus4->GetBinContent(h_corrHFPlus4->GetXaxis()->FindBin(HFEn))==0) return 1;
corrfac = 1./h_corrHFPlus4->GetBinContent(h_corrHFPlus4->GetXaxis()->FindBin(HFEn));}
return corrfac;
}

double Corrector::getEventWeightHFPlus4bak(double HFEn, bool isRebin){
double corrfac;
if(isRebin){
if(reh_corrHFPlus4->GetBinContent(reh_corrHFPlus4->GetXaxis()->FindBin(HFEn))==0) return 1;
corrfac = 1./reh_corrHFPlus4->GetBinContent(reh_corrHFPlus4->GetXaxis()->FindBin(HFEn));}
else{
if(h_corrHFPlus4->GetBinContent(h_corrHFPlus4->GetXaxis()->FindBin(HFEn))==0) return 1;
corrfac = 1./h_corrHFPlus4->GetBinContent(h_corrHFPlus4->GetXaxis()->FindBin(HFEn));}
return corrfac;
}
