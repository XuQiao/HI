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
bool Save=kFALSE;
TString histoname[2]={"jetptEta","jetptphi"};
TString histonamedata[2]={"jetptEtaCombinedSpectra","jetptphiCombinedSpectra"};
void ptetaphi(){
gStyle->SetPadTickY(1);
using namespace std;
c1 = new TCanvas("c1"," ",1000,650);
makeMultiPanelCanvas(c1,3,2,0,0,0.22,0.16,0.12);
for(int ivar=0;ivar<2;ivar++){
TH2F* hPPb=(TH2F*)fPPb->Get(histoname[ivar]);
TH2F* hdata=(TH2F*)fdata->Get(histonamedata[ivar]);
TH2F* hPP=(TH2F*)fPP->Get(histoname[ivar]);
TString varname;
if(histoname[ivar].Contains("Eta"))	varname="jet #eta_{lab}";
else	varname="jet #phi";
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

if(varname=="jet #eta_{lab}"){
double binbound_var[]={-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0};
int Nbin_var=sizeof(binbound_var)/sizeof(double)-1;
}
else{
int Nbin_var=20;
double binbound_var[21];
for(int i=0;i<Nbin_var+1;i++)
binbound_var[i]=-TMath::Pi()+i*TMath::Pi()*2/Nbin_var;
}

double xrange_var[2]={binbound_var[0]+1e-3,binbound_var[Nbin_var]-1e-3};
double xrange_pt[2]={100+1e-3,600-1e-3};
double xrange_pt_pre[2]={100,600};
/*if(xrange_pt[0]<110){
if(varname=="jet #eta_{lab}"){
double ratioptmin=0.79, ratioptmax=1.19;
double ratiovarmin=0.79, ratiovarmax=1.19;
}
else{
double ratioptmin=0.39, ratioptmax=1.39;
double ratiovarmin=0.39, ratiovarmax=1.39;
}
else{
//if(varname=="jet #eta_{lab}"){
*/
double ratioptmin=0.5, ratioptmax=1.35;
double ratiovarmin=0.5, ratiovarmax=1.35;
/*}
else{
double ratioptmin=0.69, ratioptmax=1.29;
doble ratiovarmin=0.69, ratiovarmax=1.29;*/
//}

if(!isMC)
TString Titley="Event Fraction";
else TString Titley = "#frac{d#sigma}{dp_{T}}";

TH1D* hdata_var=hdata->ProjectionY("hdata_var",hdata->GetXaxis()->FindBin(xrange_pt[0]),hdata->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hdata_var_all=hdata->ProjectionY("hdata_var_all",hdata->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hdata->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hdata_pt=hdata->ProjectionX("hdata_pt",hdata->GetYaxis()->FindBin(xrange_var[0]),hdata->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hdata_pt_all=hdata->ProjectionX("hdata_pt_all",hdata->GetYaxis()->FindBin(binbound_var[0]+1e-4),hdata->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));

TH1D* hPPb_var=hPPb->ProjectionY("hPPb_var",hPPb->GetXaxis()->FindBin(xrange_pt[0]),hPPb->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hPPb_var_all=hPPb->ProjectionY("hPPb_var_all",hPPb->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hPPb->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hPPb_pt=hPPb->ProjectionX("hPPb_pt",hPPb->GetYaxis()->FindBin(xrange_var[0]),hPPb->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hPPb_pt_all=hPPb->ProjectionX("hPPb_pt_all",hPPb->GetYaxis()->FindBin(binbound_var[0]+1e-4),hPPb->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));

TH1D* hPP_var=hPP->ProjectionY("hPP_var",hPP->GetXaxis()->FindBin(xrange_pt[0]),hPP->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hPP_var_all=hPP->ProjectionY("hPP_var_all",hPP->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hPP->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hPP_pt=hPP->ProjectionX("hPP_pt",hPP->GetYaxis()->FindBin(xrange_var[0]),hPP->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hPP_pt_all=hPP->ProjectionX("hPP_pt_all",hPP->GetYaxis()->FindBin(binbound_var[0]+1e-4),hPP->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));
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

re_hdata_pt->SetMarkerStyle(20);
re_hdata_pt->SetMarkerSize(1);
re_hdata_pt->SetMarkerColor(1);
re_hdata_pt->SetLineColor(1);

re_hPP_pt->SetLineColor(2);
if(!isMC){
re_hPP_pt->SetFillColor(2);
re_hPP_pt->SetFillStyle(3004);
re_hPP_pt->SetMarkerStyle(0);
}
else{
re_hPP_pt->SetMarkerSize(1.2);
re_hPP_pt->SetMarkerColor(2);
re_hPP_pt->SetMarkerStyle(24);
}

re_hPPb_pt->SetMarkerStyle(28);
re_hPPb_pt->SetMarkerSize(1.2);
re_hPPb_pt->SetMarkerColor(4);
re_hPPb_pt->SetLineColor(4);

re_hPP_var->SetLineColor(2);
if(!isMC){
re_hPP_var->SetFillStyle(3004);
re_hPP_var->SetMarkerStyle(0);
re_hPP_var->SetFillColor(2);
}
else{
re_hPP_var->SetMarkerSize(1);
re_hPP_var->SetMarkerColor(2);
re_hPP_var->SetMarkerStyle(24);
}

re_hdata_var->SetMarkerStyle(20);
re_hdata_var->SetMarkerSize(1);
re_hdata_var->SetMarkerColor(1);
re_hdata_var->SetLineColor(1);

re_hPPb_var->SetMarkerStyle(28);
re_hPPb_var->SetMarkerSize(1);
re_hPPb_var->SetMarkerColor(4);
re_hPPb_var->SetLineColor(4);

TLatex T;
T.SetNDC();
T.SetTextAlign(12);
T.SetTextSize(0.06);
T.SetTextColor(1);
T.SetTextFont(42);
TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1.9,2.6);
TLine *l1 =new TLine(xrange_pt_pre[0],1,xrange_pt_pre[1],1);
l1->SetLineStyle(2);
l1->SetLineColor(1);
TLegend *leg1=new TLegend(0.45,0.64,0.60,0.84);
leg1->SetBorderSize(0);
leg1->SetFillStyle(0);
leg1->SetTextSize(0.05);
leg1->SetFillColor(0);
TLegend *leg2=new TLegend(0.35,0.22,0.65,0.42);
leg2->SetBorderSize(0);
leg2->SetFillStyle(0);
leg2->SetTextSize(0.05);
leg2->SetFillColor(0);

if(ivar==0){
c1->cd(1)->SetLogy();
hFrame->GetYaxis()->SetTitle(Titley);
hFrame->GetXaxis()->SetTitle("");
//hFrame->GetXaxis()->SetRangeUser(27,692);
hFrame->GetXaxis()->SetLimits(xrange_pt_pre[0],xrange_pt_pre[1]);
if(!isMC){
if(xrange_pt[0]>110)
hFrame->GetYaxis()->SetRangeUser(1.01e-5,1);
else
hFrame->GetYaxis()->SetRangeUser(1.01e-3,1);
hFrame->DrawCopy();
re_hPP_pt->Draw("HIST E1same");
re_hdata_pt->Draw("E1same");
}
else{
if(xrange_pt[0]<110)
hFrame->GetYaxis()->SetRangeUser(1.01e-10,2e-2);
else
hFrame->GetYaxis()->SetRangeUser(1.01e-11,1e-5);
//hFrame->GetYaxis()->SetRangeUser(1.01e-10,1e-3);
hFrame->DrawCopy();
re_hPP_pt->Draw("E1same");
}
re_hPPb_pt->Draw("E1same");
T.DrawLatex(0.44,0.72,Form("%s: %.2f--%.2f ",varname.Data(),xrange_var[0],xrange_var[1]));
T.DrawLatex(0.44,0.80,Form("p_{T}^{jet}: %.f--%.f (GeV/c)",xrange_pt[0],xrange_pt[1]));
c1->cd(4);
	hFrame->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
	hFrame->GetYaxis()->SetRangeUser(ratioptmin,ratioptmax);
if(isMC){ 	// for MC PPb vs Pb check
	hFrame->GetYaxis()->SetTitle("#frac{PYTHIA+HIJING}{PYTHIA}");
	hFrame->DrawCopy();
	TH1D* ratio_pt=(TH1D*)re_hPPb_pt->Clone("ratio_pt");
	ratio_pt->Divide(re_hPP_pt);
	ratio_pt->SetMarkerStyle(20);
	ratio_pt->SetMarkerSize(1);
	ratio_pt->Draw("E1same");
//	ratio_pt->Fit("pol0");
	l1->Draw("same");
if(Save){
	c1->Print(Form("pic/DataMC/PPbCompPP_pt_%.fto%.f.pdf",xrange_pt[0],xrange_pt[1]));
	c1->Print(Form("pic/DataMC/PPbCompPP_pt_%.fto%.f.gif",xrange_pt[0],xrange_pt[1]));
}
	}

else{	// for Data vs MC study
	TH1D* ratio_DatavsPPb_pt=(TH1D*)re_hdata_pt->Clone("ratio_DatavsPPb_pt");
	ratio_DatavsPPb_pt->Divide(re_hPPb_pt);
	hFrame->GetYaxis()->SetTitle("#frac{Data}{MC}");
	hFrame->DrawCopy();
	ratio_DatavsPPb_pt->SetMarkerStyle(28);
	ratio_DatavsPPb_pt->SetLineColor(4);
	ratio_DatavsPPb_pt->SetMarkerColor(4);
	ratio_DatavsPPb_pt->SetMarkerSize(1);
	ratio_DatavsPPb_pt->Draw("E1same");

	TH1D* ratio_DatavsPP_pt=(TH1D*)re_hdata_pt->Clone("ratio_DatavsPP_pt");
	ratio_DatavsPP_pt->Divide(re_hPP_pt);
	ratio_DatavsPP_pt->SetMarkerStyle(20);
	ratio_DatavsPP_pt->SetMarkerColor(2);
	ratio_DatavsPP_pt->SetLineColor(2);
	ratio_DatavsPP_pt->SetMarkerSize(1);
	leg2->AddEntry(ratio_DatavsPPb_pt,"Data/PYTHIA+HIJING","lp");
	leg2->AddEntry(ratio_DatavsPP_pt,"Data/PYTHIA","lp");
	ratio_DatavsPP_pt->Draw("E1same");
	leg2->Draw("same");
	l1->Draw("same");
}
}

TH1F* hFrame1=new TH1F("","",800,-4,4);
fixedFontHist(hFrame1,1.9,2.6);
c1->cd(ivar+2);
hFrame1->GetYaxis()->SetLabelSize(0);
hFrame1->GetXaxis()->SetRangeUser(xrange_var[0]-1e-1,xrange_var[1]+1e-1);
hFrame1->GetYaxis()->SetTitle("");
if(!isMC){
hFrame1->GetYaxis()->SetRangeUser(1e-4,1e-1);
if(varname=="jet #eta_{lab}")
hFrame1->DrawCopy();
else{
fixedFontHist(hFrame1,1.9,2.6);
hFrame1->DrawCopy("Y+");
}
re_hPP_var->DrawCopy("HIST E1same");
re_hdata_var->DrawCopy("E1same");
leg1->AddEntry(re_hdata_var,"Data","lp");
leg1->AddEntry(re_hPP_var,"PYTHIA","lpf");
//Gaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),  gPad->GetUxmax(), gPad->GetUymax(),1e-4,1e-1,510,"+LG"); axis->Draw("same");	
}
else{
///if(xrange_pt[0]<110)
hFrame1->GetYaxis()->SetRangeUser(1.01e-3,3e-2);

if(varname!="jet #eta_{lab}"){
fixedFontHist(hFrame1,1.9,1.2);
hFrame1->DrawCopy("Y+");
}
else
//hFrame->GetYaxis()->SetRangeUser(3e-6,1e-4);
hFrame1->DrawCopy();
re_hPP_var->DrawCopy("E1same");
leg1->AddEntry(re_hPP_var,"PYTHIA","lpf");
}
leg1->AddEntry(re_hPPb_var,"PYTHIA+HIJING","lp");
re_hPPb_var->DrawCopy("E1same");
if(varname!="jet #eta_{lab}")
leg1->Draw("same");

c1->cd(ivar+5);
hFrame1->GetXaxis()->SetTitle(Form("%s",varname.Data()));
hFrame1->GetYaxis()->SetRangeUser(ratiovarmin,ratiovarmax);
hFrame1->GetYaxis()->SetLabelSize(0);
TLine *l2 =new TLine(xrange_var[0],1,xrange_var[1],1);
l2->SetLineStyle(2);
l2->SetLineColor(1);
if(isMC){	//for MC PPb vs PP check
	TH1D* ratio_var=(TH1D*)re_hPPb_var->Clone("ratio_var");
	ratio_var->Divide(re_hPP_var);
	if(varname=="jet #eta_{lab}")
        hFrame1->DrawCopy();
        else{
        fixedFontHist(hFrame1,1.9,1.2);
        hFrame1->DrawCopy("Y+");
        }
	ratio_var->SetMarkerStyle(20);
	ratio_var->SetMarkerSize(1);
	ratio_var->DrawCopy("E1same");
	l2->Draw("same");

}
else{	// for Data vs MC study	
	TH1D* ratio_DatavsPPb_var=(TH1D*)re_hdata_var->Clone("ratio_DatavsPPb_var");
	ratio_DatavsPPb_var->Divide(re_hPPb_var);
	hFrame1->GetYaxis()->SetTitle("");
	if(varname=="jet #eta_{lab}") 
	hFrame1->DrawCopy();	
	else{
	fixedFontHist(hFrame1,1.9,1.2);
	hFrame1->DrawCopy("Y+");	
	}
//	if(varname!="jet #eta_{lab}"){   TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),  gPad->GetUxmax(), gPad->GetUymax(),ratiovarmin,ratiovarmax,510,"+L"); axis->Draw("same");	}
	ratio_DatavsPPb_var->SetMarkerStyle(28);
	ratio_DatavsPPb_var->SetMarkerSize(1);
	ratio_DatavsPPb_var->SetMarkerColor(4);
	ratio_DatavsPPb_var->SetLineColor(4);       
	ratio_DatavsPPb_var->DrawCopy("E1same");
	
	TH1D* ratio_DatavsPP_var=(TH1D*)re_hdata_var->Clone("ratio_DatavsPP_var");
	ratio_DatavsPP_var->Divide(re_hPP_var);
	ratio_DatavsPP_var->SetMarkerStyle(20);
	ratio_DatavsPP_var->SetMarkerSize(1);
	ratio_DatavsPP_var->SetMarkerColor(2);
	ratio_DatavsPP_var->SetLineColor(2);
	ratio_DatavsPP_var->DrawCopy("E1same");
	l2->Draw("same");
}
}

if(Save){
if(isMC){
c1->Print(Form("pic/PPbCompPP_%.fto%.f.gif",xrange_pt[0],xrange_pt[1]));
c1->Print(Form("pic/PPbCompPP_%.fto%.f.pdf",xrange_pt[0],xrange_pt[1]));
}
else{
c1->Print(Form("pic/DataPPbCompPP_%.fto%.f.gif",xrange_pt[0],xrange_pt[1]));
c1->Print(Form("pic/DataPPbCompPP_%.fto%.f.pdf",xrange_pt[0],xrange_pt[1]));
}
}

}
