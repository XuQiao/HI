#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"

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
return g0;
}

void overlayoneJetIDcut(){
gStyle->SetErrorX(0);
gStyle->SetOptStat(kFALSE);
int ilist=12;
int ilistw;
if(ilist== 6)  {    JetID = "Max p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0]=0.05; JetIDcut[1]=1; ilistw=1;}
else if(ilist== 7)  {    JetID = "#Sum p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "Max p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "#Sum p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "Max p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "#Sum p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{jet}";JetIDcut[0]=0; JetIDcut[1]=1.01; ilistw=3;}
else if(ilist== 13) {    JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{raw}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14)   {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})";JetIDcut[0]=0;JetIDcut[1]=0.975;ilistw=2;}

if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};

int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
const int Neta=8;
const TString etabinname[Neta]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12","-10_10"};
const double etabin[Neta]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[Neta]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
int ieta=7;

TString JetIDName=JetIDNameList[ilist];
TString histonameIDData=Form("jetpt%sCombinedSpectraInEtaBin%s",JetIDName.Data(),etabinname[ieta].Data());
TH2F* hdata2F=(TH2F*)fdata->Get(histonameIDData);

const double xrange_JetID[2]={binbound_JetID[0]+1e-4,binbound_JetID[Nbin_JetID]-1e-4};
xrange_JetIDcut[0]=JetIDcut[0]+1e-4;
xrange_JetIDcut[1]=JetIDcut[1]-1e-4;

TGraphAsymmErrors *g = new TGraphAsymmErrors();
TH1D* histo1=(TH1D*)hdata2F->ProjectionX("histo1",hdata2F->GetYaxis()->FindBin(xrange_JetIDcut[0]),hdata2F->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* histo2=(TH1D*)hdata2F->ProjectionX("histo2");
histo1=(TH1D*)histo1->Rebin(Nbin_pt,"histo1",binbound_pt);
histo2=(TH1D*)histo2->Rebin(Nbin_pt,"histo2",binbound_pt);
normalizeByBinWidth(histo1);
normalizeByBinWidth(histo2);
histo1->SetMarkerStyle(20);
histo1->SetMarkerSize(1.5);
histo1->SetMarkerColor(1);
histo1->SetLineColor(1);
histo2->SetMarkerStyle(24);
histo2->SetMarkerSize(1.5);
histo2->SetMarkerColor(2);
histo2->SetLineColor(2);
g=makegraph(histo1,histo2);
g->SetMarkerSize(1.5);

c1 = new TCanvas("c1","",400,700);
makeMultiPanelCanvas(c1,1,2,-0.15,0,0.22,0.16,0.03);
TLatex *T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.03);
T->SetTextColor(1);
T->SetTextFont(42);
T->DrawLatex(0.22,0.98-ilistw*0.05,Form("Cut: %.2f<%s<%.2f",JetID.Data(),xrange_JetIDcut[0],xrange_JetIDcut[1]));
TH1F* hFrame=new TH1F("","",1000,0,1000);
c1->cd(1)->SetLogy();
hFrame->SetLineColor(0);
hFrame->GetXaxis()->SetRangeUser(25,600);
hFrame->GetYaxis()->SetTitle("#frac{dN}{dp_{T}}");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetRangeUser(0.11,1e8);
fixedFontHist(hFrame,1.7,2.8);
hFrame->DrawCopy();
histo1->Draw("same");
histo2->Draw("same");
TLegend *leg1=new TLegend(0.65,0.72,0.85,0.85);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg1->AddEntry(histo2,"Before Cut","lp");
leg1->AddEntry(histo1,"After Cut","lp");
leg1->Draw("same");

c1->cd(2);
g->SetMarkerStyle(20);
g->SetMarkerColor(1);
g->SetLineColor(1);
T->SetTextSize(0.06);
hFrame->GetYaxis()->SetTitle("Cut Efficiency");
hFrame->GetYaxis()->SetRangeUser(0.82,1.025);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->DrawCopy();
g->Draw("Psame");
T->DrawLatex(0.5,0.25,etastring[ieta]);
TLine *l =new TLine(25,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");

c1->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/TreeAna/JetID/pic/overlay_best_JetIDcut.png"));
c1->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/TreeAna/JetID/pic/overlay_best_JetIDcut.pdf"));
}
