#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TObject.h"
#include <iostream.h>
#include <fstream.h>
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"

bool isMC=kFALSE;
bool isvar;
void ptetaphi_temp(){
if(histonameID.Contains("Eta"))
TString varname="jet #eta";
else
TString varname="jet #phi";
bool var1=kFALSE;
using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
/*
gStyle->SetPadBottomMargin(0.12);
gStyle->SetPadTopMargin     (0.025);
gStyle->SetPadLeftMargin   (0.15);
gStyle->SetPadRightMargin (0.025);
gStyle->SetPadTickX           (1);
gStyle->SetPadTickY           (1);
*/

//double xrange_var[2]={-1.465,0.535};
//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//double binbound_pt[]={100,110,120,140,160,180,200,220,260,300,350,500};
if(varname=="jet #eta"){
if(!var1){
double binbound_var[]={-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0};}
else{
double binbound_var[]={-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};}
int Nbin_var=sizeof(binbound_var)/sizeof(double)-1;
}
else{
const int Nbin_var=20;
double binbound_var[Nbin_var+1];
for(int i=0;i<Nbin_var+1;i++)
binbound_var[i]=-TMath::Pi()+i*TMath::Pi()*2/Nbin_var;
}
//for(int i=17;i<18;i++){
double xrange_var[2]={binbound_var[5]+1e-4,binbound_var[7]-1e-4};

if(!isMC)
TString Titley="Event Fraction";
else TString Titley = "#frac{d#sigma}{dp_{T}}";

TH1D* hdata_var=hdata2F->ProjectionY("hdata_var",hdata2F->GetXaxis()->FindBin(xrange_pt[0]),hdata2F->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hdata_var_all=hdata2F->ProjectionY("hdata_var_all",hdata2F->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hdata2F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hdata_pt=hdata2F->ProjectionX("hdata_pt",hdata2F->GetYaxis()->FindBin(xrange_var[0]),hdata2F->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hdata_pt_all=hdata2F->ProjectionX("hdata_pt_all",hdata2F->GetYaxis()->FindBin(binbound_var[0]+1e-4),hdata2F->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));

TH1D* hPPb_var=hPPb2F->ProjectionY("hPPb_var",hPPb2F->GetXaxis()->FindBin(xrange_pt[0]),hPPb2F->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hPPb_var_all=hPPb2F->ProjectionY("hPPb_var_all",hPPb2F->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hPPb2F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hPPb_pt=hPPb2F->ProjectionX("hPPb_pt",hPPb2F->GetYaxis()->FindBin(xrange_var[0]),hPPb2F->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hPPb_pt_all=hPPb2F->ProjectionX("hPPb_pt_all",hPPb2F->GetYaxis()->FindBin(binbound_var[0]+1e-4),hPPb2F->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));

TH1D* hPP_var=hPP2F->ProjectionY("hPP_var",hPP2F->GetXaxis()->FindBin(xrange_pt[0]),hPP2F->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hPP_var_all=hPP2F->ProjectionY("hPP_var_all",hPP2F->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hPP2F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hPP_pt=hPP2F->ProjectionX("hPP_pt",hPP2F->GetYaxis()->FindBin(xrange_var[0]),hPP2F->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hPP_pt_all=hPP2F->ProjectionX("hPP_pt_all",hPP2F->GetYaxis()->FindBin(binbound_var[0]+1e-4),hPP2F->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));
if(var1){
TH1D* hdata_pt=(TH1D*)hdata->Clone("hdata_pt");
TH1D* hPPb_pt=(TH1D*)hPPb->Clone("hPPb_pt");
TH1D* hPP_pt=(TH1D*)hPP->Clone("hPP_pt");
}
TH1D* re_hdata_pt=(TH1D*)hdata_pt->Rebin(Nbin_pt,"re_hdata_pt",binbound_pt);
TH1D* re_hPP_pt=(TH1D*)hPP_pt->Rebin(Nbin_pt,"re_hPP_pt",binbound_pt);
TH1D* re_hPPb_pt=(TH1D*)hPPb_pt->Rebin(Nbin_pt,"re_hPPb_pt",binbound_pt);
TH1D* re_hdata_var=(TH1D*)hdata_var->Rebin(Nbin_var,"re_hdata_var",binbound_var);
TH1D* re_hPP_var=(TH1D*)hPP_var->Rebin(Nbin_var,"re_hPP_var",binbound_var);
TH1D* re_hPPb_var=(TH1D*)hPPb_var->Rebin(Nbin_var,"re_hPPb_var",binbound_var);

normalizeByBinWidth(re_hdata_pt);
normalizeByBinWidth(re_hPP_pt);
normalizeByBinWidth(re_hPPb_pt);
normalizeByBinWidth(re_hdata_var);
normalizeByBinWidth(re_hPP_var);
normalizeByBinWidth(re_hPPb_var);

if(!isMC){
re_hdata_pt->Scale(1/re_hdata_pt->Integral(re_hdata_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hdata_pt->GetXaxis()->FindBin(xrange_pt[1])));
re_hdata_var->Scale(1/re_hdata_var->Integral(re_hdata_var->GetXaxis()->FindBin(xrange_var[0]),re_hdata_var->GetXaxis()->FindBin(xrange_var[1])));
re_hPPb_pt->Scale(1/re_hPPb_pt->Integral(re_hPPb_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hPPb_pt->GetXaxis()->FindBin(xrange_pt[1])));
re_hPPb_var->Scale(1/re_hPPb_var->Integral(re_hPPb_var->GetXaxis()->FindBin(xrange_var[0]),re_hPPb_var->GetXaxis()->FindBin(xrange_var[1])));
re_hPP_pt->Scale(1/re_hPP_pt->Integral(re_hPP_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hPP_pt->GetXaxis()->FindBin(xrange_pt[1])));
re_hPP_var->Scale(1/re_hPP_var->Integral(re_hPP_var->GetXaxis()->FindBin(xrange_var[0]),re_hPP_var->GetXaxis()->FindBin(xrange_var[1])));
}

	TH1D* ratio_DatavsPPb_pt=(TH1D*)re_hdata_pt->Clone("ratio_DatavsPPb_pt");	ratio_DatavsPPb_pt->Divide(re_hPPb_pt);
	TH1D* ratio_DatavsPP_pt=(TH1D*)re_hdata_pt->Clone("ratio_DatavsPP_pt");	ratio_DatavsPP_pt->Divide(re_hPP_pt);
TH1D* ratio_DatavsPPb_pt_weight=(TH1D*)ratio_DatavsPPb_pt->Clone("ratio_DatavsPPb_pt_weight"); 
TH1D* ratio_DatavsPP_pt_weight=(TH1D*)ratio_DatavsPP_pt->Clone("ratio_DatavsPP_pt_weight");
ratio_DatavsPPb_pt_weight->Multiply(re_hdata_pt);
ratio_DatavsPP_pt_weight->Multiply(re_hdata_pt);
cout<<xrange_pt[0]<<" to "<<xrange_pt[1]<<" pt agreement PPb: "<<re_hdata_pt->Integral(re_hdata_pt->FindBin(xrange_pt[0]),re_hdata_pt->FindBin(xrange_pt[1]))/re_hPPb_pt->Integral(re_hPPb_pt->FindBin(xrange_pt[0]),re_hPPb_pt->FindBin(xrange_pt[1]))*100<<"% "<<endl;
cout<<"before weighting: "<<ratio_DatavsPPb_pt->Integral(ratio_DatavsPPb_pt->FindBin(xrange_pt[0]),ratio_DatavsPPb_pt->FindBin(xrange_pt[1]))<<'\t';
cout<<"after data spetrum weighting: "<<ratio_DatavsPPb_pt_weight->Integral(ratio_DatavsPPb_pt_weight->FindBin(xrange_pt[0]),ratio_DatavsPPb_pt_weight->FindBin(xrange_pt[1]))<<endl;
cout<<xrange_pt[0]<<" to "<<xrange_pt[1]<<" pt agreement PP: "<<re_hdata_pt->Integral(re_hdata_pt->FindBin(xrange_pt[0]),re_hdata_pt->FindBin(xrange_pt[1]))/re_hPP_pt->Integral(re_hPP_pt->FindBin(xrange_pt[0]),re_hPP_pt->FindBin(xrange_pt[1]))*100<<"%"<<endl;
cout<<"before weighting: "<<ratio_DatavsPP_pt->Integral(ratio_DatavsPP_pt->FindBin(xrange_pt[0]),ratio_DatavsPP_pt->FindBin(xrange_pt[1]))<<'\t';
cout<<"after data spectrum weighting: "<<ratio_DatavsPP_pt_weight->Integral(ratio_DatavsPP_pt_weight->FindBin(xrange_pt[0]),ratio_DatavsPP_pt_weight->FindBin(xrange_pt[1]))<<endl;

	TH1D* ratio_DatavsPPb_var=(TH1D*)re_hdata_var->Clone("ratio_DatavsPPb_var");	ratio_DatavsPPb_var->Divide(re_hPPb_var);
	TH1D* ratio_DatavsPP_var=(TH1D*)re_hdata_var->Clone("ratio_DatavsPP_var");	ratio_DatavsPP_var->Divide(re_hPP_var);
TH1D* ratio_DatavsPPb_var_weight=(TH1D*)ratio_DatavsPPb_var->Clone("ratio_DatavsPPb_var_weight");
TH1D* ratio_DatavsPP_var_weight=(TH1D*)ratio_DatavsPP_var->Clone("ratio_DatavsPP_var_weight");
ratio_DatavsPPb_var_weight->Multiply(re_hdata_var);
ratio_DatavsPP_var_weight->Multiply(re_hdata_var);
cout<<xrange_var[0]<<" to "<<xrange_var[1]<<varname<<" agreement PPb: "<<re_hdata_var->Integral(re_hdata_var->FindBin(xrange_var[0]),re_hdata_var->FindBin(xrange_var[1]))/re_hPPb_var->Integral(re_hPPb_var->FindBin(xrange_var[0]),re_hPPb_var->FindBin(xrange_var[1]))*100<<"% "<<endl;
cout<<"before weighting: "<<ratio_DatavsPPb_var->Integral(ratio_DatavsPPb_var->FindBin(xrange_var[0]),ratio_DatavsPPb_var->FindBin(xrange_var[1]))<<'\t';
cout<<"after data spectrum weighting: "<<ratio_DatavsPPb_var_weight->Integral(ratio_DatavsPPb_var_weight->FindBin(xrange_var[0]),ratio_DatavsPPb_var_weight->FindBin(xrange_var[1]))<<endl;
cout<<xrange_var[0]<<" to "<<xrange_var[1]<<varname<<"  agreement PP: "<<re_hdata_var->Integral(re_hdata_var->FindBin(xrange_var[0]),re_hdata_var->FindBin(xrange_var[1]))/re_hPP_var->Integral(re_hPP_var->FindBin(xrange_var[0]),re_hPP_var->FindBin(xrange_var[1]))*100<<"%"<<endl;
cout<<"before weighting: "<<ratio_DatavsPP_var->Integral(ratio_DatavsPP_var->FindBin(xrange_var[0]),ratio_DatavsPP_var->FindBin(xrange_var[1]))<<'\t';
cout<<"after data spetrum weighting: "<<ratio_DatavsPP_var_weight->Integral(ratio_DatavsPP_var_weight->FindBin(xrange_var[0]),ratio_DatavsPP_var_weight->FindBin(xrange_var[1]))<<endl;
}
