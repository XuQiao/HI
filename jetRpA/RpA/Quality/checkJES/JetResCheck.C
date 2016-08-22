#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TObject.h"
#include "iostream"
#include "fstream"
#include "root_setting.h"
#include "file.h"

TH1D* myRebin(TH1D* histo,int Nbin, double binbound[]){
//rehisto->Scale(1/rehisto->Integral());

if(binbound[Nbin]==1000){	//rebin pt
//TH1D* rehisto=(TH1D*)histo->Rebin(Nbin,histo->GetName(),binbound);
TH1D* rehisto=(TH1D*)histo->Clone(histo->GetName());
rehisto->GetXaxis()->SetRangeUser(binbound[5],binbound[Nbin]);
}

fixedFontHist(rehisto,1.2,1.6);

return rehisto;
}

void JetResCheck(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);

//------------------------------------------------------------Get Histograms---------------------------------------------
//refptpt, ptrefpt, ptJES, rawptJES, RESrefptpt, rawptwtJES

TString str="rawptJES";


if(!str.Contains("ref")){
TH2F *hdata2D= (TH2F*)fdata->Get(Form("jet%s_0-100%%",str.Data()));
TProfile *hProfdata=(TProfile*)hdata2D->ProfileX("hProfdata");
rehisto_hProfdata=myRebin(hProfdata,Nbin_pt,binbound_pt);
}
TH2F *hPP2D= (TH2F*)fPP->Get(Form("jet%s_0-100%%",str.Data()));
TH2F *hPPb2D= (TH2F*)fPPb->Get(Form("jet%s_0-100%%",str.Data()));
TProfile *hProfPP=(TProfile*)hPP2D->ProfileX("hProfPP");
TProfile *hProfPPb=(TProfile*)hPPb2D->ProfileX("hProfPPb");

//---------------------------------------------------Rebin and Format histograms-----------------------------------------

rehisto_hProfPP=myRebin(hProfPP,Nbin_pt,binbound_pt);
rehisto_hProfPPb=myRebin(hProfPPb,Nbin_pt,binbound_pt);

if(!str.Contains("ref")){
rehisto_hProfdata->SetTitle("");
rehisto_hProfdata->GetYaxis()->SetTitle(hdata2D->GetYaxis()->GetTitle());
rehisto_hProfdata->SetMaximum(TMath::Max(rehisto_hProfPP->GetMaximum(),rehisto_hProfdata->GetMaximum())*1.1);
rehisto_hProfdata->SetMinimum(TMath::Min(rehisto_hProfPP->GetMinimum(),rehisto_hProfdata->GetMinimum())/1.2);
rehisto_hProfdata->SetMarkerSize(1.3);
rehisto_hProfdata->SetMarkerStyle(20);
rehisto_hProfdata->SetMarkerColor(1);
}
rehisto_hProfPP->GetYaxis()->SetTitle(hPP2D->GetYaxis()->GetTitle());
rehisto_hProfPP->SetMaximum(TMath::Max(rehisto_hProfPP->GetMaximum(),rehisto_hProfPPb->GetMaximum())*1.1);
rehisto_hProfPP->SetMinimum(TMath::Min(rehisto_hProfPP->GetMinimum(),rehisto_hProfPPb->GetMinimum())/1.1);
rehisto_hProfPP->SetMarkerSize(0);
rehisto_hProfPP->SetMarkerStyle(0);
rehisto_hProfPP->SetMarkerColor(2);
rehisto_hProfPP->SetFillStyle(3004);
rehisto_hProfPP->SetFillColor(2);
rehisto_hProfPP->SetTitle(Form(""));
rehisto_hProfPPb->SetTitle(Form(""));
rehisto_hProfPPb->SetMarkerSize(1.3);
rehisto_hProfPPb->SetMarkerStyle(30);
rehisto_hProfPPb->SetMarkerColor(4);


//--------------------------------------------------PP&PPb rawptpt,pT Profile histograms-----------------------------------

c1 = new TCanvas("c1"," ",600,600);
makeMultiPanelCanvas(c1,1,1,0.01,0.02,0.14,0.1,0.03);
if(!str.Contains("ref"))
rehisto_hProfdata->Draw();
rehisto_hProfPP->Draw("Hist same");
rehisto_hProfPPb->Draw("same");
TLegend *leg1=new TLegend(0.55,0.7,0.75,0.9);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.05);
if(!str.Contains("ref"))
leg1->AddEntry(rehisto_hProfdata,"Data","lp");
leg1->AddEntry(rehisto_hProfPP,"PYTHIA","lpf");
leg1->AddEntry(rehisto_hProfPPb,"PYTHIA+HIJING","lp");
leg1->Draw("same");
TLine *l =new TLine(binbound_pt[5],1,binbound_pt[Nbin_pt],1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
if(str.Contains("RES") || str.Contains("wt")){
TLatex *T1=new TLatex(0.4,0.92,"Applied Residual");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
}

c1->Print(Form("pic/JetIDcut/Profile_%s.png",str.Data()));
}
