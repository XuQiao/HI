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

TString coll="PPb";
int marker[]={20,21,33};
int color[]={1,4,2};
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

TH1* makehisto(int ilist, int ieta, bool isX, bool iscut, int type){
//------------------------------------------------------------Get Histograms---------------------------------------------
double xrange_JetIDcut[2]; double JetIDcut[2];
TString JetIDName = JetIDNameList[ilist];
TString JetID,Unit;
if(ilist== 0)  { JetID = "chargedMax"; JetIDcut[0]=4; JetIDcut[1]=300.;}
else if(ilist== 1)   {   JetID = "chargedSum"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 2)   {   JetID = "neutralMax"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 3)   {   JetID = "neutralSum"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 4)   {   JetID = "photonMax"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 5)   {   JetID = "photonSum"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 6)  {    JetID = "chargedMax/p_{T}"; JetIDcut[0]=0.05; JetIDcut[1]=1;}
else if(ilist== 7)  {    JetID = "chargedSum/p_{T}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "neutralMax/p_{T}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "neutralSum/p_{T}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "photonMax/p_{T}";	JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "photonSum/p_{T}";	JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}",JetIDcut[0]=0; JetIDcut[1]=1.01;}
else if(ilist== 13) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14)   {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})" ;JetIDcut[0]=0;JetIDcut[1]=0.975;}
else if(ilist== 15)     JetID = "charged Multiplicity";
else if(ilist== 16)     JetID = "neutral Multiplicity";
else if(ilist== 17)     JetID = "photon Multiplicity";
else if(ilist== 18)   {  JetID = "#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e"; JetIDcut[0]=0; JetIDcut[1]=1000;}
else{   exit(0);}

if(JetIDName.Contains("pt") || JetIDName=="neuMaxr"){
if(ilist==12 || ilist==13){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else if (ilist==14){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2};}
else if(ilist==8 || ilist==10){
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1.};}

Unit="";
}
else if(JetIDName.Contains("N")){
if(ilist==17){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,30};}
else if(ilist==16){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,10};}
else if(ilist==15){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,24,30,40};}
Unit="";
}
else{
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,300};
Unit="(GeV/c)";
}
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
if(type==0)
TString histoname = Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data());
else if(type==1)
TString histoname = Form("jetpt%sEtaBin_real%s",JetIDName.Data(),etabinnameswap[ieta].Data());
else
TString histoname = Form("jetpt%sEtaBin_fake%s",JetIDName.Data(),etabinnameswap[ieta].Data());

if(coll=="PPb")
TH2F *hMC2D= (TH2F*)fMCPPbJetID->Get(histoname);
if(coll=="PP")
TH2F *hMC2D= (TH2F*)fMCPPxSecJetID->Get(histoname);


const double binbound_pt_coarse[]={0,30,70,100,150,600};
const int Nbin_pt_coarse=sizeof(binbound_pt_coarse)/sizeof(double)-1;
double cuthigh, cutlow;

if(isX){
if(iscut){
cutlow = JetIDcut[0]+1e-4;
cuthigh = JetIDcut[1]-1e-4;
}
else{
cutlow = binbound_JetID[0]+1e-4;
cuthigh = binbound_JetID[Nbin_JetID]-1e-4;
}
TH1D* histo=(TH1D*)hMC2D->ProjectionX("histo",hMC2D->GetYaxis()->FindBin(cutlow),hMC2D->GetYaxis()->FindBin(cuthigh));
histo->SetTitle(Form("%.2f<%s<%.2f",cutlow,JetID.Data(),cuthigh));
//if(type==2)
//histo=(TH1D*)histo->Rebin(Nbin_pt_coarse,"histo",binbound_pt_coarse);
//else
histo=(TH1D*)histo->Rebin(Nbin_pt,"histo",binbound_pt);

}
else{
cutlow = xrange_pt[0];
cuthigh = xrange_pt[1];
TH1D* histo=(TH1D*)hMC2D->ProjectionY("histo",hMC2D->GetXaxis()->FindBin(cutlow),hMC2D->GetXaxis()->FindBin(cuthigh));
histo->SetTitle(Form("%.2f < p_{T} < %.2f [GeV/c]",cutlow,cuthigh));
histo=(TH1D*)histo->Rebin(Nbin_JetID,"histo",binbound_JetID);
}
histo->SetMarkerSize(0.04);
for(int ibin=0;ibin<histo->GetNbinsX();ibin++)
histo->SetBinError(ibin,1e-15);
return histo;
}

void ptJetIDeff(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
//gStyle->SetLabelFont(70);
/*for(int iptbin=0;iptbin<=500;iptbin++)
binbound_pt[iptbin]=iptbin*2;
int Nbin_pt=500;*/

const int Nlist=3;
TH1F* hFrame=new TH1F("","",100000,0,1000);
hFrame->SetTitle("");

//--------------------------------------- MC pt cut efficiency for f/r/i----------------------------------------------------------

TCanvas *c1[Nlist];
for(int ilist=0;ilist<Nlist;ilist++){
c1[ilist] = new TCanvas(Form("c%d",ilist),"",600,600);
makeMultiPanelCanvas(c1[ilist],1,1,-0.15,0,0.17,0.1,0.03);
}
TCanvas* c4 = new TCanvas("c4"," ",600,600);
makeMultiPanelCanvas(c4,1,1,-0.1,0,0.15,0.1,0.03);
TCanvas* c5 = new TCanvas("c5"," ",1200,600);
makeMultiPanelCanvas(c5,4,2,0,0,0.28,0.20,0.03);
hFrame->GetXaxis()->SetTitle("p_{T} [GeV/c]");
hFrame->GetXaxis()->SetNdivisions(505);
hFrame->GetXaxis()->SetLimits(47.5,445);
TH1D* cut_real[Nlist][Neta];TH1D* real[Nlist][Neta];
TH1D* cut_fake[Nlist][Neta];TH1D* fake[Nlist][Neta];
TH1D* cut_inc[Nlist][Neta];TH1D* inc[Nlist][Neta];
TH1D* corr[Nlist][Neta];
TGraphAsymmErrors *geff_real[Nlist][Neta];
TGraphAsymmErrors *geff_fake[Nlist][Neta];
TLegend *leg1=new TLegend(0.02,0.03,0.40,0.35);
TLegend *leg2=new TLegend(0.18,0.65,0.48,0.95);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetLineWidth(0);
leg2->SetLineWidth(0);
leg1->SetTextSize(0.042);
leg2->SetTextSize(0.038);
TLatex *T1=new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.03);
T1->SetTextColor(1);
T1->SetTextFont(42);
TLine *l = new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
int uselist[Nlist]={0,6,14};
for(int ilist=0;ilist<Nlist;ilist++){
for(int ieta=0;ieta<Neta;ieta++){
cut_real[ilist][ieta] = (TH1D*)makehisto(uselist[ilist],ieta,1,1,1);
real[ilist][ieta] = (TH1D*)makehisto(uselist[ilist],ieta,1,0,1);
cut_fake[ilist][ieta] = (TH1D*)makehisto(uselist[ilist],ieta,1,1,2);
fake[ilist][ieta] = (TH1D*)makehisto(uselist[ilist],ieta,1,0,2);
cut_inc[ilist][ieta] = (TH1D*)makehisto(uselist[ilist],ieta,1,1,0);
inc[ilist][ieta] = (TH1D*)makehisto(uselist[ilist],ieta,1,0,0);
}
}

hFrame->GetYaxis()->SetRangeUser(0.78,1.14);
hFrame->GetYaxis()->SetTitle("#frac{1-Fake Fraction}{Real Cut Efficiency}");
for(int ieta=0;ieta<Neta;ieta++){
c5->cd(canvas[ieta]+1);
fixedFontHist(hFrame,1.4,2.0);
hFrame->GetYaxis()->SetTitleSize(0.06);
hFrame->GetYaxis()->SetLabelSize(0.08);
hFrame->GetXaxis()->SetTitleSize(0.06);
hFrame->GetXaxis()->SetLabelSize(0.08);
hFrame->DrawCopy();
for(int ilist=0;ilist<Nlist;ilist++){
corr[ilist][ieta]=(TH1D*)real[ilist][ieta]->Clone(Form("corr_%d_%d",ilist,ieta));
corr[ilist][ieta]->Divide(cut_inc[ilist][ieta]);
corr[ilist][ieta]->SetMarkerSize(1.3);
corr[ilist][ieta]->SetMarkerStyle(marker[ilist]);
corr[ilist][ieta]->SetMarkerColor(color[ilist]);
corr[ilist][ieta]->SetLineColor(color[ilist]);
corr[ilist][ieta]->DrawCopy("same");
if(ieta==3) {
leg1->AddEntry(corr[ilist][ieta],cut_real[ilist][ieta]->GetTitle(),"lp");
leg1->Draw("same");
}
}
T1->SetTextSize(0.065);
T1->DrawLatex(0.58,0.5,etastring[ieta]);
l->Draw("same");
if(ieta==0){
if(coll=="PPb")
T1->DrawLatex(0.35,0.9,"PYTHIA+HIJING");
if(coll=="PP");
T1->DrawLatex(0.35,0.9,"PYTHIA");
}

}
c5->Print(Form("%sCorr_etabin.pdf",coll.Data()));
c5->Print(Form("%sCorr_etabin.png",coll.Data()));
}
/*
fixedFontHist(hFrame,1.1,1.9);
hFrame->GetYaxis()->SetTitleSize(0.04);
hFrame->GetYaxis()->SetLabelSize(0.04);
hFrame->GetXaxis()->SetTitleSize(0.04);
hFrame->GetXaxis()->SetLabelSize(0.04);
hFrame->DrawCopy();
hFrame->GetYaxis()->SetRangeUser(0.85,1.25);
hFrame->GetYaxis()->SetTitle("#frac{1-Fake Fraction}{Real Cut Efficiency}");
for(int ilist=0;ilist<Nlist;ilist++){
c1[ilist]->cd(1);
hFrame->DrawCopy();
for(int ieta=0;ieta<Neta;ieta++){
corr[ilist][ieta]->SetMarkerStyle(markerCode[ieta]);
corr[ilist][ieta]->SetMarkerColor(colorCode[ieta]);
corr[ilist][ieta]->SetLineColor(colorCode[ieta]);
corr[ilist][ieta]->DrawCopy("same");
if(ilist==0)leg2->AddEntry(corr[ilist][ieta],etastring[ieta],"lp");
}
for(int ieta=0;ieta<Neta-1;ieta++) corr[ilist][ieta]->DrawCopy("same");
leg2->Draw("same");
T1->SetTextSize(0.03);
T1->DrawLatex(0.35,0.2,cut_real[ilist][0]->GetTitle());
l->Draw();
}
}*/
/*
c2->cd(1);
hFrame->GetYaxis()->SetTitle("Fake Cut Efficiency");
hFrame->DrawCopy();
for(int ieta=0;ieta<Neta;ieta++){
geff_fake[ieta]=makegraph(cut_fake[ieta],fake[ieta]);
geff_fake[ieta]->SetMarkerStyle(markerCode[ieta]);
geff_fake[ieta]->SetMarkerColor(colorCode[ieta]);
geff_fake[ieta]->SetLineColor(colorCode[ieta]);
geff_fake[ieta]->Draw("Psame");
}
leg1->Draw("same");
T1->DrawLatex(0.2,0.2,cut_real[0]->GetTitle());
l->Draw();

hFrame->GetXaxis()->SetTitle("p_{T} [GeV/c]");
hFrame->GetXaxis()->SetRangeUser(xrange_pt[0],200);
TGraphAsymmErrors *fake_fr[Nlist][Neta];
TGraphAsymmErrors *real_fr[Nlist][Neta];
c3->cd(1);
hFrame->GetYaxis()->SetRangeUser(0.8,1.1);
hFrame->GetYaxis()->SetTitle("Real Fraction After cut");
hFrame->DrawCopy();
for(int ieta=0;ieta<Neta;ieta++){
real_fr[ieta]=makegraph(cut_real[ieta],cut_inc[ieta]);
real_fr[ieta]->SetMarkerStyle(markerCode[ieta]);
real_fr[ieta]->SetMarkerColor(colorCode[ieta]);
real_fr[ieta]->SetLineColor(colorCode[ieta]);
real_fr[ieta]->SetLineColor(colorCode[ieta]);
real_fr[ieta]->Draw("Psame");
}
leg1->Draw("same");
T1->DrawLatex(0.2,0.2,cut_real[0]->GetTitle());
l->Draw();

c4->cd(1);
hFrame->GetYaxis()->SetRangeUser(0,0.2);
hFrame->GetYaxis()->SetTitle("Fake Fraction After cut");
hFrame->DrawCopy();
for(int ieta=0;ieta<Neta;ieta++){
fake_fr[ieta]=makegraph(cut_fake[ieta],cut_inc[ieta]);
fake_fr[ieta]->SetMarkerStyle(markerCode[ieta]);
fake_fr[ieta]->SetMarkerColor(colorCode[ieta]);
fake_fr[ieta]->SetLineColor(colorCode[ieta]);
fake_fr[ieta]->SetLineColor(colorCode[ieta]);
fake_fr[ieta]->Draw("Psame");
}
leg1->Draw("same");
T1->DrawLatex(0.2,0.2,cut_fake[0]->GetTitle());
l->Draw();
*/
