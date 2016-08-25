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
#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

TGraphAsymmErrors* makegraph(TH1* histo1, TH1* histo2){
const int Npoint=1000;
double x[Npoint],y[Npoint],xerrlow[Npoint],xerrhigh[Npoint],yerrlow[Npoint],yerrhigh[Npoint];
TGraphAsymmErrors *g=new TGraphAsymmErrors();
g->BayesDivide(histo1,histo2);
for(int ipoint=0;ipoint<g->GetN();ipoint++){
x[ipoint]=g->GetX()[ipoint];
y[ipoint]=g->GetY()[ipoint];
yerrhigh[ipoint]=g->GetErrorYhigh(ipoint);
yerrlow[ipoint]=g->GetErrorYlow(ipoint);
xerrhigh[ipoint]=0;
xerrlow[ipoint]=0;
}
g0=new TGraphAsymmErrors(g->GetN(),x,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
g0->SetMarkerSize(1.2);
return g0;
}


void DrawFakeFraction(){

using namespace std;
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
//------------------------------------------------------------Get Histograms---------------------------------------------

TH2F* hMC2D[Neta]; TH2F* hMC_real2D[Neta]; TH2F* hMC_fake2D[Neta];
TH1D* hMC[Neta]; TH1D* hMC_real[Neta]; TH1D* hMC_fake[Neta];
TH1F* ratio_hMC_real[Neta]; TH1F* ratio_hMC_fake[Neta];
TGraphAsymmErrors *g_real[Neta]; TGraphAsymmErrors *g_fake[Neta]; 
TH1F* hFrame = new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetRangeUser(47.5,200);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
fixedFontHist(hFrame,1.2,1.6);
hFrame->GetXaxis()->SetLabelSize(0.04);
hFrame->GetYaxis()->SetLabelSize(0.04);
hFrame->GetXaxis()->SetTitleSize(0.04);
hFrame->GetYaxis()->SetTitleSize(0.04);
TCanvas* c1 = new TCanvas("c1"," ",600,600);
TCanvas* c2 = new TCanvas("c2"," ",600,600);
makeMultiPanelCanvas(c1,1,1,-0.1,0,0.14,0.12,0.03);
makeMultiPanelCanvas(c2,1,1,-0.1,0,0.14,0.12,0.03);
TLatex *T1=new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.045);
T1->SetTextColor(1);
T1->SetTextFont(42);
TLegend *leg1=new TLegend(0.68,0.65,0.92,0.88);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetLineWidth(0);
leg1->SetTextSize(0.04);
TLegend *leg2=new TLegend(0.65,0.64,0.88,0.92);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetLineWidth(0);
leg2->SetTextSize(0.04);
TLine *l = new TLine(47.5,1,500,1);
l->SetLineStyle(2);
TString coll="PPb";
if(coll=="PPb")
TFile *fMCJetID = (TFile*)fMCPPbJetID;
else if(coll=="PP")
TFile *fMCJetID = (TFile*)fMCPPxSecJetID;
else if(coll=="PbP")
TFile *fMCJetID = (TFile*)fMCPbPJetIDOld;

int ilist=0;
TString JetIDName=JetIDNameList[ilist];
for(int ieta=0;ieta<Neta;ieta++){
if(coll=="PPb"){
hMC2D[ieta]= (TH2F*)fMCJetID->Get(Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data()));
hMC_real2D[ieta]= (TH2F*)fMCJetID->Get(Form("jetpt%sEtaBin_real%s",JetIDName.Data(),etabinnameswap[ieta].Data()));
hMC_fake2D[ieta]= (TH2F*)fMCJetID->Get(Form("jetpt%sEtaBin_fake%s",JetIDName.Data(),etabinnameswap[ieta].Data()));
TString mc="PYTHIA+HIJING";
}
if(coll=="PbP"){
hMC2D[ieta]= (TH2F*)fMCJetID->Get(Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinname[ieta].Data()));
hMC_real2D[ieta]= (TH2F*)fMCJetID->Get(Form("jetpt%sEtaBin_real%s",JetIDName.Data(),etabinname[ieta].Data()));
hMC_fake2D[ieta]= (TH2F*)fMCJetID->Get(Form("jetpt%sEtaBin_fake%s",JetIDName.Data(),etabinname[ieta].Data()));
TString mc="HIJING+PYTHIA";
}
hMC[ieta]=(TH1D*)hMC2D[ieta]->ProjectionX();
hMC_real[ieta]=(TH1D*)hMC_real2D[ieta]->ProjectionX();
hMC_fake[ieta]=(TH1D*)hMC_fake2D[ieta]->ProjectionX();

hMC[ieta]=(TH1D*)hMC[ieta]->Rebin(Nbin_pt,hMC[ieta]->GetName(),binbound_pt);
hMC_real[ieta]=(TH1D*)hMC_real[ieta]->Rebin(Nbin_pt,hMC_real[ieta]->GetName(),binbound_pt);
hMC_fake[ieta]=(TH1D*)hMC_fake[ieta]->Rebin(Nbin_pt,hMC_fake[ieta]->GetName(),binbound_pt);
}
c1->cd(1);
hFrame->GetYaxis()->SetRangeUser(0.7,1.05);
hFrame->GetYaxis()->SetTitle("Real jet Fraction");
hFrame->DrawCopy();
for(int ieta=0;ieta<Neta;ieta++){
//ratio_hMC_real[ieta]=(TH1F*)hMC_real[ieta]->Clone(Form("ratio_hMC_real_%d",ieta));
//ratio_hMC_real[ieta]->Divide(hMC[ieta]);
//ratio_hMC_real[ieta]->Draw("E1same");
g_real[ieta]=makegraph(hMC_real[ieta],hMC[ieta]);
g_real[ieta]->SetMarkerStyle(markerCode[ieta]);
g_real[ieta]->SetMarkerColor(colorCode[ieta]);
g_real[ieta]->SetLineColor(colorCode[ieta]);
g_real[ieta]->Draw("Psame");
leg1->AddEntry(g_real[ieta],etastring[ieta],"lp");
}
T1->DrawLatex(0.2,0.9,mc);
leg1->Draw("same");
l->Draw("same");

c2->cd(1);
hFrame->GetYaxis()->SetRangeUser(0,0.4);
hFrame->GetYaxis()->SetTitle("Fake jet Fraction");
hFrame->DrawCopy();
for(int ieta=0;ieta<Neta;ieta++){
//ratio_hMC_fake[ieta]=(TH1F*)hMC_fake[ieta]->Clone(Form("ratio_hMC_fake_%d",ieta));
//ratio_hMC_fake[ieta]->Divide(hMC[ieta]);
//ratio_hMC_fake[ieta]->Draw("E1same");
g_fake[ieta]=makegraph(hMC_fake[ieta],hMC[ieta]);
g_fake[ieta]->SetMarkerStyle(markerCode[ieta]);
g_fake[ieta]->SetMarkerColor(colorCode[ieta]);
g_fake[ieta]->SetLineColor(colorCode[ieta]);
g_fake[ieta]->Draw("Psame");
leg2->AddEntry(g_fake[ieta],etastring[ieta],"lp");
}
T1->DrawLatex(0.2,0.9,mc);
leg2->Draw("same");
l->Draw("same");

c1->Print(Form("%srealFraction.png",coll.Data()));
c2->Print(Form("%sfakeFraction.png",coll.Data()));
c2->Print(Form("%sfakeFraction.pdf",coll.Data()));
}

