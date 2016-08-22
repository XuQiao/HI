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
#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

bool Save=kTRUE;
TString histonameMC[2]={"jetptEta_0-100%","jetptphi_0-100%"};
//TString histoname[2]={"jetptEtaCombinedSpectra","jetptphiCombinedSpectra"};
TString histoname[2]={"jetptEta","jetptphi"};
TString histoname_woRes[2]={"jetptEta_woRes","jetptphi"};
void ptetaphi_PbP(){
gStyle->SetPadTickY(1);
using namespace std;
c1 = new TCanvas("c1"," ",1000,650);
makeMultiPanelCanvas(c1,3,2,0,0,0.22,0.16,0.12);
//TFile *fDataPPb = new TFile("/home/xuq/Documents/HI/RpA/output/TreeAna/DataPPbcombined.root");
//TFile *fMCPPb = new TFile("/home/xuq/Documents/HI/RpA/output/akPu3PF/MCPPbNoGpluscW_akPu3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root");
for(int ivar=0;ivar<2;ivar++){
TH2F* hPPb=(TH2F*)fMCPbP->Get(histoname[ivar]);
TH2F* hPP=(TH2F*)fMCPPxSec->Get(histoname[ivar]);
TH2F* hdata=(TH2F*)fDataPbPOld->Get(histoname[ivar]);
TH2F* hdata_woRes=(TH2F*)fDataPbPOld->Get(histoname_woRes[ivar]);
TString varname;
if(histoname[ivar].Contains("Eta"))	varname="#eta_{lab}^{jet}";
else	varname="#phi^{jet}";
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

if(varname=="#eta_{lab}^{jet}"){
double binbound_var[]={-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0};
int Nbin_var=sizeof(binbound_var)/sizeof(double)-1;
}
else{
int Nbin_var=20;
double binbound_var[21];
for(int i=0;i<Nbin_var+1;i++)
binbound_var[i]=-TMath::Pi()+i*TMath::Pi()*2/Nbin_var;
}
if(varname=="#eta_{lab}^{jet}"){
//double xrange_var[2]={-1-0.465+1e-3,1-0.465-1e-3};
double xrange_var[2]={-3+1e-3,3-1e-3};
double xrange_var_pre[2]={-3,3};
}
else{
double xrange_var[2]={binbound_var[0]+1e-3,binbound_var[Nbin_var]-1e-3};}
double xrange_pt_pre[2]={47.5+1e-3,480-1e-3};
double xrange_pt_rebin[2]={50+1e-3,600-1e-3};
double ratioptmin=0.7, ratioptmax=1.52;
double ratiovarmin=0.7, ratiovarmax=1.52;

TString Titley="Event Fraction";

TH1D* hdata_var=hdata->ProjectionY("hdata_var",hdata->GetXaxis()->FindBin(xrange_pt_rebin[0]),hdata->GetXaxis()->FindBin(xrange_pt_rebin[1]));
//TH1D* hdata_var_all=hdata->ProjectionY("hdata_var_all",hdata->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hdata->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hdata_pt=hdata_woRes->ProjectionX("hdata_pt",hdata->GetYaxis()->FindBin(xrange_var[0]),hdata->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hdata_pt_all=hdata->ProjectionX("hdata_pt_all",hdata->GetYaxis()->FindBin(binbound_var[0]+1e-4),hdata->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));
//TH1D* hdata_pt =(TH1D*)fDataPPb->Get("jetpt");

TH1D* hPPb_var=hPPb->ProjectionY("hPPb_var",hPPb->GetXaxis()->FindBin(xrange_pt_rebin[0]),hPPb->GetXaxis()->FindBin(xrange_pt_rebin[1]));
//TH1D* hPPb_var_all=hPPb->ProjectionY("hPPb_var_all",hPPb->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hPPb->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hPPb_pt=hPP->ProjectionX("hPPb_pt",hPPb->GetYaxis()->FindBin(xrange_var[0]),hPPb->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hPPb_pt_all=hPPb->ProjectionX("hPPb_pt_all",hPPb->GetYaxis()->FindBin(binbound_var[0]+1e-4),hPPb->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));
//TH1D* hPPb_pt =(TH1D*)fMCPPbxSec->Get("jetpt");

TH1D* re_hdata_pt=(TH1D*)hdata_pt->Rebin(Nbin_pt,"re_hdata_pt",binbound_pt);
TH1D* re_hPPb_pt=(TH1D*)hPPb_pt->Rebin(Nbin_pt,"re_hPPb_pt",binbound_pt);
TH1D* re_hdata_var=(TH1D*)hdata_var->Rebin(Nbin_var,"re_hdata_var",binbound_var);
TH1D* re_hPPb_var=(TH1D*)hPPb_var->Rebin(Nbin_var,"re_hPPb_var",binbound_var);

normalizeByBinWidth(re_hdata_pt);
normalizeByBinWidth(re_hPPb_pt);
normalizeByBinWidth(re_hdata_var);
normalizeByBinWidth(re_hPPb_var);

re_hdata_pt->Scale(1./re_hdata_pt->Integral(re_hdata_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hdata_pt->GetXaxis()->FindBin(xrange_pt[1])));
re_hdata_var->Scale(1./re_hdata_var->Integral(re_hdata_var->GetXaxis()->FindBin(xrange_var_pre[0]),re_hdata_var->GetXaxis()->FindBin(xrange_var_pre[1])));
re_hPPb_pt->Scale(1./re_hPPb_pt->Integral(re_hPPb_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hPPb_pt->GetXaxis()->FindBin(xrange_pt[1])));
re_hPPb_var->Scale(1./re_hPPb_var->Integral(re_hPPb_var->GetXaxis()->FindBin(xrange_var_pre[0]),re_hPPb_var->GetXaxis()->FindBin(xrange_var_pre[1])));
double scale = re_hdata_pt->Integral(re_hdata_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hdata_pt->GetXaxis()->FindBin(xrange_pt[1]))/re_hdata_pt->Integral(re_hdata_pt->GetXaxis()->FindBin(xrange_pt_pre[0]),re_hdata_pt->GetXaxis()->FindBin(xrange_pt_pre[1]));
//re_hdata_pt->Scale((double)scale);
//re_hPPb_pt->Scale((double)scale);

re_hdata_pt->SetMarkerStyle(20);
re_hdata_pt->SetMarkerSize(1);
re_hdata_pt->SetMarkerColor(1);
re_hdata_pt->SetLineColor(1);

re_hPPb_pt->SetFillStyle(3004);
re_hPPb_pt->SetFillColor(2);
re_hPPb_pt->SetLineColor(2);

re_hdata_var->SetMarkerStyle(20);
re_hdata_var->SetMarkerSize(1);
re_hdata_var->SetMarkerColor(1);
re_hdata_var->SetLineColor(1);

re_hPPb_var->SetFillStyle(3004);
re_hPPb_var->SetFillColor(2);
re_hPPb_var->SetLineColor(2);

TLatex T;
T.SetNDC();
T.SetTextAlign(12);
T.SetTextSize(0.06);
T.SetTextColor(1);
T.SetTextFont(42);
TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1.2,1.5);
hFrame->GetXaxis()->SetTitleSize(0.06);
hFrame->GetYaxis()->SetTitleSize(0.06);
hFrame->GetXaxis()->SetLabelSize(0.06);
hFrame->GetYaxis()->SetLabelSize(0.06);
TLine *l1 =new TLine(xrange_pt_pre[0],1,xrange_pt_pre[1],1);
l1->SetLineStyle(2);
l1->SetLineColor(1);
TLegend *leg1=new TLegend(0.30,0.64,0.45,0.84);
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
hFrame->GetXaxis()->SetRangeUser(xrange_pt_pre[0],xrange_pt_pre[1]);
hFrame->GetYaxis()->SetRangeUser(1.01e-8,1);
hFrame->DrawCopy();
re_hdata_pt->Draw("E1same");
re_hPPb_pt->Draw("E1 HISTsame");
if(varname=="#eta_{lab}^{jet}")
T.DrawLatex(0.44,0.80,Form("%.1f < %s < %.1f",varname.Data(),xrange_var[0],xrange_var[1]));
else
T.DrawLatex(0.44,0.80,Form("%.2f < %s < %.2f",varname.Data(),xrange_var[0],xrange_var[1]));
c1->cd(4);
	hFrame->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
	hFrame->GetYaxis()->SetRangeUser(ratioptmin,ratioptmax);
	TH1D* ratio_DatavsPPb_pt=(TH1D*)re_hdata_pt->Clone("ratio_DatavsPPb_pt");
	ratio_DatavsPPb_pt->Divide(re_hPPb_pt);
	hFrame->GetYaxis()->SetTitle("#frac{Data}{MC}");
	hFrame->DrawCopy();
	ratio_DatavsPPb_pt->SetMarkerStyle(20);
	ratio_DatavsPPb_pt->SetLineColor(1);
	ratio_DatavsPPb_pt->SetMarkerColor(1);
	ratio_DatavsPPb_pt->SetMarkerSize(1);
	ratio_DatavsPPb_pt->Draw("E1same");
l1->Draw("same");
}

TH1F* hFrame1=new TH1F("","",800,-4,4);
hFrame1->GetYaxis()->SetLabelSize(0.06);
hFrame1->GetYaxis()->SetTitleSize(0.06);
hFrame1->GetXaxis()->SetTitleSize(0.06);
hFrame1->GetXaxis()->SetLabelSize(0.06);
c1->cd(ivar+2);
hFrame1->GetYaxis()->SetLabelSize(0);
hFrame1->GetXaxis()->SetRangeUser(xrange_var_pre[0]-2e-1,xrange_var_pre[1]+2e-1);
hFrame1->GetYaxis()->SetTitle("");
hFrame1->GetYaxis()->SetRangeUser(1e-4,1e-1);
fixedFontHist(hFrame,1.4,2.0);
if(varname=="#eta_{lab}^{jet}"){
hFrame1->DrawCopy();
T.DrawLatex(0.34,0.80,Form("%.f < p_{T}^{jet} < %.f (GeV/c)",xrange_pt_rebin[0],xrange_pt_rebin[1]));}
else{
hFrame1->DrawCopy("Y+");
}
re_hdata_var->DrawCopy("E1same");
leg1->AddEntry(re_hdata_var,"Proton going positive Data","lp");
/*
if(varname!="#eta_{lab}^{jet}"){
fixedFontHist(hFrame1,1.9,1.2);
hFrame1->DrawCopy("Y+");
}
hFrame1->DrawCopy();
*/
leg1->AddEntry(re_hPPb_var,"HIJING+PYTHIA","lp");
re_hPPb_var->DrawCopy("E1HISTsame");
if(varname!="#eta_{lab}^{jet}")
leg1->Draw("same");

c1->cd(ivar+5);
hFrame1->GetXaxis()->SetTitle(Form("%s",varname.Data()));
hFrame1->GetYaxis()->SetRangeUser(ratiovarmin,ratiovarmax);
hFrame1->GetYaxis()->SetLabelSize(0);
hFrame1->GetXaxis()->CenterTitle();
hFrame1->GetXaxis()->CenterTitle();
TLine *l2 =new TLine(xrange_var_pre[0],1,xrange_var_pre[1],1);
l2->SetLineStyle(2);
l2->SetLineColor(1);
	TH1D* ratio_DatavsPPb_var=(TH1D*)re_hdata_var->Clone("ratio_DatavsPPb_var");
	ratio_DatavsPPb_var->Divide(re_hPPb_var);
	hFrame1->GetYaxis()->SetTitle("");
	if(varname=="#eta_{lab}^{jet}") 
	hFrame1->DrawCopy();	
	else{
	hFrame1->DrawCopy("Y+");
	}
	ratio_DatavsPPb_var->SetMarkerStyle(20);
	ratio_DatavsPPb_var->SetMarkerSize(1);
	ratio_DatavsPPb_var->SetMarkerColor(1);
	ratio_DatavsPPb_var->SetLineColor(1);       
	ratio_DatavsPPb_var->DrawCopy("E1same");	
	l2->Draw("same");
}
if(Save){
c1->Print(Form("DataPbPComp.png"));
c1->Print(Form("DataPbPComp.pdf"));
}

}
