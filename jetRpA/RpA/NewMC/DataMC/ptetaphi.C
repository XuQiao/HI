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
#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
#include "/home/xuq7/HI/jetRpA/RpA/NewMC/produceandcheck/file.h"

bool Save=kTRUE;
TString histoname[2]={"jetptEta","jetptphi"};
TString histoname_woRes[2]={"jetptEta_woRes","jetptphi"};
void ptetaphi(){
TString coll="PbP";
gStyle->SetPadTickY(1);
using namespace std;
c1 = new TCanvas("c1"," ",1000,650);
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
makeMultiPanelCanvas(c1,3,2,0,0,0.22,0.16,0.12);
//TFile *fMCMC = new TFile("/home/maoy/Public/PPb_UnfoPriorGen_akPu3PFOfficialMCNoIDCut_MCJECv8_jtpt20_EtaBin-10_10_Inc_v6.root");
//TFile *fDataMC = new TFile("/home/maoy/Public/DATAPPbakPu3PFJetSpectraKurtCombineJetTriggerEtaWeight7EtabinJetIDCutRecoPt.root");
for(int ivar=0;ivar<2;ivar++){
//TH2F* hPP=(TH2F*)fMCPPxSec->Get(histoname[ivar]);
if(coll=="PPb"){
TH2F* hMC=(TH2F*)fMCPPb->Get(histoname[ivar]);
TH2F* hdata=(TH2F*)fDataPPb->Get(histoname[ivar]);
TH2F* hdata_woRes=(TH2F*)fDataPPb->Get(histoname_woRes[ivar]);
}
if(coll=="PbP"){
TH2F* hMC=(TH2F*)fMCPbP->Get(histoname[ivar]);
TH2F* hdata=(TH2F*)fDataPbP->Get(histoname[ivar]);
TH2F* hdata_woRes=(TH2F*)fDataPbP->Get(histoname_woRes[ivar]);
}
TString varname;
if(histoname[ivar].Contains("Eta"))	varname="#eta_{lab}^{jet}";
else	varname="#phi^{jet}";

if(varname=="#eta_{lab}^{jet}"){
double binbound_var[]={-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0};
int Nbin_var=sizeof(binbound_var)/sizeof(double)-1;
TString varname_ = "eta";
}
else{
int Nbin_var=20;
double binbound_var[21];
for(int i=0;i<Nbin_var+1;i++)
binbound_var[i]=-TMath::Pi()+i*TMath::Pi()*2/Nbin_var;
TString varname_ = "phi";
}
if(varname=="#eta_{lab}^{jet}"){
//double xrange_var[2]={-1-0.465+1e-3,1-0.465-1e-3};
double xrange_var[2]={-2.4+1e-3,2.4-1e-3};
double xrange_var_pre[2]={-2.6,2.7};
}
else{
double xrange_var[2]={binbound_var[0]+1e-3,binbound_var[Nbin_var]-1e-3};
double xrange_var_pre[2]={-3.4,3.4};
}
double ratiomin=0.7, ratiomax=1.52;

TString Titley="Event Fraction";
//TString Titley="Jet cross section d#sigma/dp_{T}";

TH1D* hdata_var=hdata->ProjectionY(Form("hdata_%s",varname_.Data()),hdata->GetXaxis()->FindBin(xrange_pt[0]),hdata->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hdata_var_all=hdata->ProjectionY("hdata_var_all",hdata->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hdata->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hdata_pt=hdata->ProjectionX("hdata_pt",hdata->GetYaxis()->FindBin(xrange_var[0]),hdata->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hdata_pt=hdata->ProjectionX("hdata_pt",0,-1);
//TH1D* hdata_pt_all=hdata->ProjectionX("hdata_pt_all",hdata->GetYaxis()->FindBin(binbound_var[0]+1e-4),hdata->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));
//TH1D* hdata_pt =(TH1D*)fDataPPb->Get("jetpt");
//TH1D* hdata_pt = (TH1D*)fDataMC->Get("jetptEtaBin-10_10");

TH1D* hMC_var=hMC->ProjectionY(Form("hMC_%s",varname_.Data()),hMC->GetXaxis()->FindBin(xrange_pt[0]),hMC->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hMC_var_all=hMC->ProjectionY("hMC_var_all",hMC->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hMC->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hMC_pt=hMC->ProjectionX("hMC_pt",hMC->GetYaxis()->FindBin(xrange_var[0]),hMC->GetYaxis()->FindBin(xrange_var[1]));
//TH1D* hMC_pt=hMC->ProjectionX("hMC_pt",0,-1);
//TH1D* hMC_pt_all=hMC->ProjectionX("hMC_pt_all",hMC->GetYaxis()->FindBin(binbound_var[0]+1e-4),hMC->GetYaxis()->FindBin(binbound_var[Nbin_var]-1e-4));
//TH1D* hMC_pt =(TH1D*)fMCPPb->Get("jetpt");
//TH1D* hMC_pt = (TH1D*)fMCMC->Get("hMeas0");
//        TH1F *hPP_pt1 = (TH1F*)fMCPP->Get(Form("ak3GenJetSpectrum_QCD10001_%s/JetSpectrum_Fine",etabinnamesym[0].Data()));
//        TH1F *hPP_pt2 = (TH1F*)fMCPP->Get(Form("ak3GenJetSpectrum_QCD10001_%s/JetSpectrum_Fine",etabinnamesym[1].Data()));
//        TH1F *hMC_pt = (TH1F*)hPP_pt1->Clone();
//        hMC_pt->Add(hPP_pt2);

TH1D* re_hdata_pt=(TH1D*)hdata_pt->Rebin(Nbin_pt,"re_hdata_pt",binbound_pt);
TH1D* re_hMC_pt=(TH1D*)hMC_pt->Rebin(Nbin_pt,"re_hMC_pt",binbound_pt);
TH1D* re_hdata_var=(TH1D*)hdata_var->Rebin(Nbin_var,"re_hdata_var",binbound_var);
TH1D* re_hMC_var=(TH1D*)hMC_var->Rebin(Nbin_var,"re_hMC_var",binbound_var);

normalizeByBinWidth(re_hdata_pt);
normalizeByBinWidth(re_hMC_pt);
normalizeByBinWidth(re_hdata_var);
normalizeByBinWidth(re_hMC_var);

re_hdata_pt->Scale(1./re_hdata_pt->Integral(re_hdata_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hdata_pt->GetXaxis()->FindBin(xrange_pt[1])));
re_hdata_var->Scale(1./re_hdata_var->Integral(re_hdata_var->GetXaxis()->FindBin(xrange_var[0]),re_hdata_var->GetXaxis()->FindBin(xrange_var[1])));
re_hMC_pt->Scale(1./re_hMC_pt->Integral(re_hMC_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hMC_pt->GetXaxis()->FindBin(xrange_pt[1])));
re_hMC_var->Scale(1./re_hMC_var->Integral(re_hMC_var->GetXaxis()->FindBin(xrange_var[0]),re_hMC_var->GetXaxis()->FindBin(xrange_var[1])));


//re_hdata_var->Scale(1.0/20.4*1e-6/208);
//re_hdata_pt->Scale(1.0/20.4*1e-6/208);
//double scale = re_hdata_pt->Integral(re_hdata_pt->GetXaxis()->FindBin(xrange_pt[0]),re_hdata_pt->GetXaxis()->FindBin(xrange_pt[1]))/re_hdata_pt->Integral(re_hdata_pt->GetXaxis()->FindBin(xrange_pt_pre[0]),re_hdata_pt->GetXaxis()->FindBin(xrange_pt_pre[1]));
//re_hdata_pt->Scale((double)scale);
//re_hMC_pt->Scale((double)scale);

re_hdata_pt->SetMarkerStyle(20);
re_hdata_pt->SetMarkerSize(1);
re_hdata_pt->SetMarkerColor(1);
re_hdata_pt->SetLineColor(1);

re_hMC_pt->SetFillStyle(3004);
re_hMC_pt->SetFillColor(2);
re_hMC_pt->SetLineColor(2);

re_hdata_var->SetMarkerStyle(20);
re_hdata_var->SetMarkerSize(1);
re_hdata_var->SetMarkerColor(1);
re_hdata_var->SetLineColor(1);

re_hMC_var->SetFillStyle(3004);
re_hMC_var->SetFillColor(2);
re_hMC_var->SetLineColor(2);

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
//hFrame->GetXaxis()->SetLabelSize(1);
//hFrame->GetYaxis()->SetLabelSize(1);
TLine *l1 =new TLine(xrange_pt_pre[0],1,xrange_pt_pre[1],1);
l1->SetLineStyle(2);
l1->SetLineColor(1);
TLegend *leg1=new TLegend(0.35,0.64,0.50,0.84);
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
hFrame->GetXaxis()->SetNdivisions(505);
fixedFontHist(hFrame,1.6,2.2);
//hFrame->GetXaxis()->SetRangeUser(27,692);
hFrame->GetXaxis()->SetRangeUser(xrange_pt_pre[0],xrange_pt_pre[1]);
//hFrame->GetYaxis()->SetRangeUser(1.01e-10,1e-2);
hFrame->GetYaxis()->SetRangeUser(1.01e-7,1e1);
hFrame->DrawCopy();
re_hdata_pt->Draw("E1same");
re_hMC_pt->Draw("E1 HISTsame");
if(varname=="#eta_{lab}^{jet}")
T.DrawLatex(0.44,0.80,Form("%.1f < %s < %.1f",varname.Data(),xrange_var[0],xrange_var[1]));
else
T.DrawLatex(0.44,0.80,Form("%.2f < %s < %.2f",varname.Data(),xrange_var[0],xrange_var[1]));
c1->cd(4);
	hFrame->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
	hFrame->GetYaxis()->SetRangeUser(ratiomin,ratiomax);
	TH1D* ratio_DatavsMC_pt=(TH1D*)re_hdata_pt->Clone("ratio_DatavsMC_pt");
	ratio_DatavsMC_pt->Divide(re_hMC_pt);
	hFrame->GetYaxis()->SetTitle("#frac{Data}{MC}");
	hFrame->DrawCopy();
	ratio_DatavsMC_pt->SetMarkerStyle(20);
	ratio_DatavsMC_pt->SetLineColor(1);
	ratio_DatavsMC_pt->SetMarkerColor(1);
	ratio_DatavsMC_pt->SetMarkerSize(1);
	ratio_DatavsMC_pt->Draw("E1same");
l1->Draw("same");
}

TH1F* hFrame1=new TH1F("","",800,-4,4);
//hFrame1->GetYaxis()->SetLabelSize(0.06);
hFrame1->GetYaxis()->SetTitleSize(0.06);
hFrame1->GetXaxis()->SetTitleSize(0.06);
//hFrame1->GetXaxis()->SetLabelSize(0.06);
c1->cd(ivar+2);
hFrame1->GetYaxis()->SetLabelSize(0.06);
hFrame1->GetXaxis()->SetRangeUser(xrange_var_pre[0],xrange_var_pre[1]);
hFrame1->GetYaxis()->SetTitle("");
//hFrame1->GetYaxis()->SetRangeUser(1e-6,1e-2);
//hFrame1->GetYaxis()->SetRangeUser(1.01e-4,1e-3);
hFrame1->GetYaxis()->SetRangeUser(0.008,1e-1);
fixedFontHist(hFrame1,1.9,1.2);
if(varname=="#eta_{lab}^{jet}"){
hFrame1->DrawCopy();
T.DrawLatex(0.34,0.80,Form("%.f < p_{T}^{jet} < %.f (GeV/c)",xrange_pt_pre[0],xrange_pt_pre[1]));
}
else{
hFrame1->DrawCopy("Y+");
}
re_hdata_var->DrawCopy("E1same");
leg1->AddEntry(re_hdata_var,"Pb going positive Data","lp");

//if(varname!="#eta_{lab}^{jet}"){
//fixedFontHist(hFrame1,1.9,1.2);
//hFrame1->DrawCopy("Y+");
//}
//hFrame1->DrawCopy();

leg1->AddEntry(re_hMC_var,"PYTHIA+HIJING","lp");
re_hMC_var->DrawCopy("E1HISTsame");
if(varname!="#eta_{lab}^{jet}")
leg1->Draw("same");

c1->cd(ivar+5);
hFrame1->GetXaxis()->SetTitle(Form("%s",varname.Data()));
hFrame1->GetYaxis()->SetRangeUser(ratiomin,ratiomax);
hFrame1->GetYaxis()->SetLabelSize(0);
hFrame1->GetXaxis()->CenterTitle();
hFrame1->GetXaxis()->CenterTitle();
TLine *l2 =new TLine(xrange_var_pre[0],1,xrange_var_pre[1],1);
l2->SetLineStyle(2);
l2->SetLineColor(1);
	TH1D* ratio_DatavsMC_var=(TH1D*)re_hdata_var->Clone("ratio_DatavsMC_var");
	ratio_DatavsMC_var->Divide(re_hMC_var);
//	hFrame1->GetYaxis()->SetTitle("");
	if(varname=="#eta_{lab}^{jet}") 
	hFrame1->DrawCopy();	
	else{
	hFrame1->DrawCopy("Y+");
	}
	ratio_DatavsMC_var->SetMarkerStyle(20);
	ratio_DatavsMC_var->SetMarkerSize(1);
	ratio_DatavsMC_var->SetMarkerColor(1);
	ratio_DatavsMC_var->SetLineColor(1);       
	ratio_DatavsMC_var->DrawCopy("E1same");	
	l2->Draw("same");
}
if(Save){
c1->Print(Form("Data%sComp.png",coll.Data()));
c1->Print(Form("Data%sComp.pdf",coll.Data()));
}

}
