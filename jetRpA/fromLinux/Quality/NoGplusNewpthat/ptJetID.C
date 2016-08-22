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
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"

bool isMC=kTRUE;


TH1* myRebin(TH1* histo,int Nbin, double binbound[]){
//rehisto->Scale(1/rehisto->Integral());

TH1* rehisto=histo->Clone("rehisto");
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
TString histoname=Form("%s",histo->GetName());
if(!histoname.Contains("Prof")){
normalizeByBinWidth(rehisto);
}
if(!binbound[Nbin]==1000)
rehisto->GetXaxis()->SetRangeUser(binbound[0],binbound[Nbin]);
fixedFontHist(rehisto,1.2,1.3);

return rehisto;
}

void ptJetID(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
gStyle->SetLabelFont(70);

//------------------------------------------------------------Get Histograms---------------------------------------------

TFile *fdata=TFile::Open(datafile);
TFile *fPP=TFile::Open(PPfile);
TFile *fPPb=TFile::Open(PPbfile);

double xrange_JetIDcut[2]; double JetIDcut[2];

if(ilist== 0)  { JetID = "chargedMax"; JetIDcut[0]=1.1; JetIDcut[1]=200.;}
else if(ilist== 1)      JetID = "chargedSum";
else if(ilist== 2)      JetID = "neutralMax";
else if(ilist== 3) JetID = "neutralSum";
else if(ilist== 4) JetID = "photonMax";
else if(ilist== 5) JetID = "photonSum";
else if(ilist== 6)  {    JetID = "chargedMax/p_{T}"; if(algo=="ak3PF"){JetIDcut[0]=0; JetIDcut[1]=0.24;} else {JetIDcut[0] = 0;JetIDcut[1]=0.3;}} 
else if(ilist== 7)  {    JetID = "chargedSum/p_{T}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "neutralMax/p_{T}"; JetIDcut[0] = 0; if(algo=="ak3PF")JetIDcut[1]=0.12;  else JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "neutralSum/p_{T}"; JetIDcut[0] = 0; if(algo=="ak3PF")JetIDcut[1]=0.12;  else JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "photonMax/p_{T}";	JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "photonSum/p_{T}";	JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{jet}";JetIDcut[0]=0; JetIDcut[1]=1.01;}
else if(ilist== 13)     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{raw}";
else if(ilist== 14)     JetID = "neutralMax/(neutralMax+chargedMax+photonMax)";
else if(ilist== 15)     JetID = "charged Multiplicity";
else if(ilist== 16)     JetID = "neutral Multiplicity";
else if(ilist== 17)     JetID = "photon Multiplicity";
else{   exit();}

if(JetIDName.Contains("pt") || JetIDName=="neuMaxr"){
if(ilist==12 || ilist==13){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
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
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,200};
Unit="(GeV/c)";
}
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;


TH2F *hdata2D= (TH2F*)fdata->Get(Form("jetpt%s_0-100%%",JetIDName.Data()));
TH2F *hPP2D= (TH2F*)fPP->Get(Form("jetpt%s_0-100%%",JetIDName.Data()));
TH2F *hPPb2D= (TH2F*)fPPb->Get(Form("jetpt%s_0-100%%",JetIDName.Data()));
TH2F *hPP2D_real= (TH2F*)fPP->Get(Form("jetpt%s_real_0-100%%",JetIDName.Data()));
TH2F *hPPb2D_real= (TH2F*)fPPb->Get(Form("jetpt%s_real_0-100%%",JetIDName.Data()));
TH2F *hPP2D_fake= (TH2F*)fPP->Get(Form("jetpt%s_fake_0-100%%",JetIDName.Data()));
TH2F *hPPb2D_fake= (TH2F*)fPPb->Get(Form("jetpt%s_fake_0-100%%",JetIDName.Data()));
TH2F *hPP_Etaphi= (TH2F*)fPP->Get(Form("jetEtaphi_0-100%%"));
TH2F *hPPb_Etaphi= (TH2F*)fPPb->Get(Form("jetEtaphi_0-100%%"));

//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const Double_t TrkBin_Max[]={0,0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 138., 155.2, 172.4, 189.6, 206.8};
//const int nPtBin_Max = sizeof(TrkBin_Max)/sizeof(Double_t)-1 ;

const double binbound_pt_coarse[]={0,30,70,150,200,600};
const int Nbin_pt_coarse=sizeof(binbound_pt_coarse)/sizeof(double)-1;
const double xrange_JetID[2]={binbound_JetID[0]+1e-4,binbound_JetID[Nbin_JetID]-1e-4};
xrange_JetIDcut[0]=JetIDcut[0]+1e-4;
xrange_JetIDcut[1]=JetIDcut[1]-1e-4;
TH1D* hdata_JetID=hdata2D->ProjectionY("hdata_JetID",hdata2D->GetXaxis()->FindBin(xrange_pt[0]),hdata2D->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hdata_pt=hdata2D->ProjectionX("hdata_pt",hdata2D->GetYaxis()->FindBin(xrange_JetID[0]),hdata2D->GetYaxis()->FindBin(xrange_JetID[1]));

TH1D* hPPb_pt=hPPb2D->ProjectionX("hPPb_pt",hPPb2D->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPPb_pt_real=hPPb2D_real->ProjectionX("hPPb_pt_real",hPPb2D_real->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D_real->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPPb_pt_fake=hPPb2D_fake->ProjectionX("hPPb_pt_fake",hPPb2D_fake->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D_fake->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPPb_JetID=hPPb2D->ProjectionY("hPPb_JetID",hPPb2D->GetXaxis()->FindBin(xrange_pt[0]),hPPb2D->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPPb_JetID_real=hPPb2D_real->ProjectionY("hPPb_JetID_real",hPPb2D_real->GetXaxis()->FindBin(xrange_pt[0]),hPPb2D_real->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPPb_JetID_fake=hPPb2D_fake->ProjectionY("hPPb_JetID_fake",hPPb2D_fake->GetXaxis()->FindBin(xrange_pt[0]),hPPb2D_fake->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPPb_JetIDcut_pt=hPPb2D->ProjectionX("hPPb_JetIDcut_pt",hPPb2D->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* hPPb_JetIDcut_pt_real=hPPb2D_real->ProjectionX("hPPb_JetIDcut_pt_real",hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* hPPb_JetIDcut_pt_fake=hPPb2D_fake->ProjectionX("hPPb_JetIDcut_pt_fake",hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[1]));

TH1D* hPP_pt=hPP2D->ProjectionX("hPP_pt",hPP2D->GetYaxis()->FindBin(xrange_JetID[0]),hPP2D->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPP_pt_real=hPP2D_real->ProjectionX("hPP_pt_real",hPP2D_real->GetYaxis()->FindBin(xrange_JetID[0]),hPP2D_real->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPP_pt_fake=hPP2D_fake->ProjectionX("hPP_pt_fake",hPP2D_fake->GetYaxis()->FindBin(xrange_JetID[0]),hPP2D_fake->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPP_JetID=hPP2D->ProjectionY("hPP_JetID",hPP2D->GetXaxis()->FindBin(xrange_pt[0]),hPP2D->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPP_JetID_real=hPP2D_real->ProjectionY("hPP_JetID_real",hPP2D_real->GetXaxis()->FindBin(xrange_pt[0]),hPP2D_real->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPP_JetID_fake=hPP2D_fake->ProjectionY("hPP_JetID_fake",hPP2D_fake->GetXaxis()->FindBin(xrange_pt[0]),hPP2D_fake->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPP_JetIDcut_pt=hPP2D->ProjectionX("hPP_JetIDcut_pt",hPP2D->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPP2D->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* hPP_JetIDcut_pt_real=hPP2D_real->ProjectionX("hPP_JetIDcut_pt_real",hPP2D_real->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* hPP_JetIDcut_pt_fake=hPP2D_fake->ProjectionX("hPP_JetIDcut_pt_fake",hPP2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[1]));

//TProfile *hProfPP=(TProfile*)hPP2D->ProfileX("hProfPP",hPP2D->GetYaxis()->FindBin(xrange_JetID[0]),hPP2D->GetYaxis()->FindBin(xrange_JetID[1]));
TProfile *hProfPP=(TProfile*)hPP2D->ProfileX("hProfPP",1,-1);
TProfile *hProfPP_fake=(TProfile*)hPP2D_fake->ProfileX("hProfPP_fake",1,-1);
TProfile *hProfPP_real=(TProfile*)hPP2D_real->ProfileX("hProfPP_real",1,-1);
TProfile *hProfPPb=(TProfile*)hPPb2D->ProfileX("hProfPPb",1,-1);
TProfile *hProfPPb_fake=(TProfile*)hPPb2D_fake->ProfileX("hProfPPb_fake",1,-1);
TProfile *hProfPPb_real=(TProfile*)hPPb2D_real->ProfileX("hProfPPb_real",1,-1);

//---------------------------------------------------Rebin and Format histograms-----------------------------------------

TH1D* rehisto_hdata_pt=myRebin(hdata_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hdata_JetID=myRebin(hdata_JetID,Nbin_JetID,binbound_JetID);

TH1D* rehisto_hPP_pt=myRebin(hPP_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_pt_real=myRebin(hPP_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_pt_fake=myRebin(hPP_pt_fake,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_JetIDcut_pt=myRebin(hPP_JetIDcut_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_JetIDcut_pt_real=myRebin(hPP_JetIDcut_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_JetIDcut_pt_fake=myRebin(hPP_JetIDcut_pt_fake,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_JetID=myRebin(hPP_JetID,Nbin_JetID,binbound_JetID);
TH1D* rehisto_hPP_JetID_real=myRebin(hPP_JetID_real,Nbin_JetID,binbound_JetID);
TH1D* rehisto_hPP_JetID_fake=myRebin(hPP_JetID_fake,Nbin_JetID,binbound_JetID);

TH1D* rehisto_hPPb_pt=myRebin(hPPb_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_pt_real=myRebin(hPPb_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_pt_fake=myRebin(hPPb_pt_fake,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_JetIDcut_pt=myRebin(hPPb_JetIDcut_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_JetIDcut_pt_real=myRebin(hPPb_JetIDcut_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_JetIDcut_pt_fake=myRebin(hPPb_JetIDcut_pt_fake,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_JetID=myRebin(hPPb_JetID,Nbin_JetID,binbound_JetID);
TH1D* rehisto_hPPb_JetID_real=myRebin(hPPb_JetID_real,Nbin_JetID,binbound_JetID);
TH1D* rehisto_hPPb_JetID_fake=myRebin(hPPb_JetID_fake,Nbin_JetID,binbound_JetID);

TProfile *rehisto_hProfPP=(TProfile*)myRebin(hProfPP,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPP_fake=(TProfile*)myRebin(hProfPP_fake,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPP_real=(TProfile*)myRebin(hProfPP_real,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPPb=(TProfile*)myRebin(hProfPPb,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPPb_fake=(TProfile*)myRebin(hProfPPb_fake,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPPb_real=(TProfile*)myRebin(hProfPPb_real,Nbin_pt,binbound_pt);
cout<<"ptlowbin to pthighbin: "<<'\t'<<"PP fake Entries"<<'\t'<<"PPb fake Entries"<<endl;
for(int i=1;i<rehisto_hProfPP_fake->GetXaxis()->GetNbins();i++){
cout<<rehisto_hProfPP_fake->GetXaxis()->GetBinLowEdge(i)<<" to "<<rehisto_hProfPP_fake->GetXaxis()->GetBinUpEdge(i)<<" : "<<rehisto_hProfPP_fake->GetBinEffectiveEntries(i)<<"\t"<<rehisto_hProfPPb_fake->GetBinEffectiveEntries(i)<<"\t"<<endl;
/*if(rehisto_hProfPP_fake->GetBinEffectiveEntries(i)<5)
rehisto_hProfPP_fake->SetBinEntries(i,0);
if(rehisto_hProfPPb_fake->GetBinEffectiveEntries(i)<5)
rehisto_hProfPPb_fake->SetBinEntries(i,0);*/
}


//TProfile* rehisto_hProfPP=(TProfile*)hProfPP->Clone(hProfPP->GetName());
//TProfile* rehisto_hProfPP_fake=(TProfile*)hProfPP_fake->Clone(hProfPP_fake->GetName());
//TProfile* rehisto_hProfPP_real=(TProfile*)hProfPP_real->Clone(hProfPP_real->GetName());
//TProfile* rehisto_hProfPPb=(TProfile*)hProfPPb->Clone(hProfPPb->GetName());
//TProfile* rehisto_hProfPPb_fake=(TProfile*)hProfPPb_fake->Clone(hProfPPb_fake->GetName());
//TProfile* rehisto_hProfPPb_real=(TProfile*)hProfPPb_real->Clone(hProfPPb_real->GetName());

if(!isMC){
rehisto_hdata_pt->SetMarkerStyle(20);
rehisto_hdata_pt->SetMarkerSize(1.5);
rehisto_hdata_pt->SetMarkerColor(1);
rehisto_hdata_pt->SetLineColor(1);
rehisto_hPP_pt->SetFillStyle(3004);
rehisto_hPP_pt->SetMarkerStyle(0);
rehisto_hPP_pt->SetFillColor(2);
rehisto_hPP_pt->SetMarkerColor(2);
rehisto_hPP_pt->SetLineColor(2);
rehisto_hPP_pt->SetLineColor(2);
}
else{
rehisto_hPP_pt->SetMarkerStyle(24);
rehisto_hPP_pt->SetMarkerSize(1.5);
}
rehisto_hPP_pt->SetTitle(Form("%.2f--%.2f %s",xrange_JetID[0],xrange_JetID[1], Unit.Data()));
rehisto_hPP_pt->GetYaxis()->SetTitle("Probability");
rehisto_hPP_JetIDcut_pt->SetTitle(Form(""));
rehisto_hPP_JetIDcut_pt->GetYaxis()->SetTitle("Probability");

rehisto_hPPb_pt->SetTitle(Form("%.2f--%.2f %s",xrange_JetID[0],xrange_JetID[1], Unit.Data()));
rehisto_hPPb_pt->SetMarkerStyle(24);
rehisto_hPPb_pt->SetMarkerSize(1.5);
rehisto_hPPb_pt->SetMarkerColor(1);
rehisto_hPPb_pt->SetLineColor(1);
rehisto_hPPb_pt->SetLineColor(1);
rehisto_hPPb_pt->GetYaxis()->SetTitle("");
rehisto_hPPb_JetIDcut_pt->SetTitle("");
rehisto_hPPb_JetIDcut_pt->GetYaxis()->SetTitle("");

if(!isMC){
rehisto_hdata_JetID->SetMarkerStyle(20);
rehisto_hdata_JetID->SetMarkerSize(1.5);
rehisto_hdata_JetID->SetMarkerColor(1);
rehisto_hdata_JetID->SetLineColor(1);
rehisto_hPP_JetID->SetFillStyle(3004);
rehisto_hPP_JetID->SetMarkerStyle(0);
rehisto_hPP_JetID->SetFillColor(2);
rehisto_hPP_JetID->SetLineColor(2);
rehisto_hPP_JetID->SetLineColor(2);
}
else{
rehisto_hPP_JetID->SetMarkerStyle(24);
rehisto_hPP_JetID->SetMarkerSize(1.5);
}
rehisto_hPP_JetID->SetTitle("");
rehisto_hPP_JetID->SetMaximum(rehisto_hPP_JetID->GetMaximum()*10);
rehisto_hPP_JetID->SetMinimum(rehisto_hPP_JetID_fake->GetMinimum()/2+1e-11);
rehisto_hPP_JetID->GetYaxis()->SetTitle("Probability");
rehisto_hPP_JetID->GetXaxis()->SetTitle(Form("%s %s",JetID.Data(), Unit.Data()));
rehisto_hPP_JetID_fake->SetMarkerStyle(20);
rehisto_hPP_JetID_fake->GetXaxis()->SetTitle(Form("%s %s",JetID.Data(), Unit.Data()));
rehisto_hPP_JetID_fake->GetYaxis()->SetTitle("Probability");
rehisto_hPP_JetID_fake->SetMarkerColor(2);
rehisto_hPP_JetID_fake->SetLineColor(2);
rehisto_hPP_JetID_fake->SetMarkerSize(1.5);
rehisto_hPP_JetID_real->SetMarkerStyle(29);
rehisto_hPP_JetID_real->SetMarkerColor(4);
rehisto_hPP_JetID_real->SetLineColor(4);
rehisto_hPP_JetID_real->SetMarkerSize(1.5);
rehisto_hPPb_JetID->SetMarkerStyle(24);
rehisto_hPPb_JetID->SetMaximum(rehisto_hPP_JetID->GetMaximum());
rehisto_hPPb_JetID->SetMinimum(rehisto_hPP_JetID->GetMinimum());
rehisto_hPPb_JetID->SetMarkerSize(1.5);
rehisto_hPPb_JetID->SetTitle("");
rehisto_hPPb_JetID->GetYaxis()->SetTitle("");
rehisto_hPPb_JetID->GetXaxis()->SetTitle(Form("%s %s",JetID.Data(), Unit.Data()));
rehisto_hPPb_JetID_fake->SetMarkerStyle(20);
rehisto_hPPb_JetID_fake->GetXaxis()->SetTitle(Form("%s %s",JetID.Data(), Unit.Data()));
rehisto_hPPb_JetID_fake->SetMarkerSize(1.5);
rehisto_hPPb_JetID_fake->SetMarkerColor(2);
rehisto_hPPb_JetID_fake->SetLineColor(2);
rehisto_hPPb_JetID_real->SetMarkerStyle(29);
rehisto_hPPb_JetID_real->SetMarkerSize(1.5);
rehisto_hPPb_JetID_real->SetMarkerColor(4);
rehisto_hPPb_JetID_real->SetLineColor(4);

rehisto_hProfPP->SetMarkerSize(1.3);
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(30.1,599.9);
rehisto_hProfPP_fake->SetMarkerSize(1.3);
rehisto_hProfPP_real->SetMarkerSize(1.3);
rehisto_hProfPP->SetMarkerStyle(24);
rehisto_hProfPP_fake->SetMarkerStyle(20);
rehisto_hProfPP_real->SetMarkerStyle(29);
rehisto_hProfPP->SetMarkerColor(1);
rehisto_hProfPP->SetLineColor(1);
rehisto_hProfPP_fake->SetMarkerColor(2);
rehisto_hProfPP_fake->SetLineColor(2);
rehisto_hProfPP_real->SetMarkerColor(4);
rehisto_hProfPP_real->SetLineColor(4);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");

rehisto_hProfPPb->SetMarkerSize(1.3);
rehisto_hProfPPb_fake->SetMarkerSize(1.3);
rehisto_hProfPPb_real->SetMarkerSize(1.3);
rehisto_hProfPPb->SetMarkerStyle(24);
rehisto_hProfPPb_fake->SetMarkerStyle(20);
rehisto_hProfPPb_real->SetMarkerStyle(29);

rehisto_hProfPPb->SetMarkerColor(1);
rehisto_hProfPPb->SetLineColor(1);
rehisto_hProfPPb_fake->SetMarkerColor(2);
rehisto_hProfPPb_fake->SetLineColor(2);
rehisto_hProfPPb_real->SetMarkerColor(4);
rehisto_hProfPPb_real->SetLineColor(4);

hPP_Etaphi->SetTitle("");
hPP_Etaphi->GetXaxis()->SetRangeUser(-3-1e-4,3+1e-4);
hPPb_Etaphi->SetTitle("");
hPPb_Etaphi->GetYaxis()->SetTitle("");
hPPb_Etaphi->GetXaxis()->SetRangeUser(-3-1e-4,3+1e-4);

fixedFontHist(hFrame,1.2,1.2);
fixedFontHist(rehisto_hProfPP,1.2,1.2);
fixedFontHist(rehisto_hProfPPb,1.2,1.2);
fixedFontHist(hPP_Etaphi,1.2,1.2);
fixedFontHist(hPPb_Etaphi,1.2,1.2);


//--------------------------------------------------------pT distribution------------------------------------------------

c1 = new TCanvas("c1"," ",800,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.08,0.1,0.03);
c1->cd(1)->SetLogy();
rehisto_hPP_pt->SetMarkerStyle(20);
rehisto_hPP_pt->Draw("E1");
if(!isMC)
rehisto_hdata_pt->Draw("E1same");
rehisto_hPPb_pt->Draw("E1same");
TLegend *leg1=new TLegend(0.7,0.7,0.9,0.8);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetLineWidth(0);
leg1->SetTextSize(0.05);
leg1->AddEntry(rehisto_hPP_pt,"PYTHIA","lpf");
if(!isMC)
leg1->AddEntry(rehisto_hdata_pt,"Data","lp");
leg1->AddEntry(rehisto_hPPb_pt,"PYTHIA+HIJING","lp");
leg1->Draw("same");
c1->Print(Form("pic/JetIDcut/%s/%s/jetpt_PPbvsPP.png",algo.Data(),JetIDName.Data()));

//-------------------------------------------------PP&PPb JetID distribution---------------------------------------------

c2 = new TCanvas("c2"," ",1200,600);
makeMultiPanelCanvas(c2,2,1,0,0,0.15,0.1,0.03);
c2->cd(1)->SetLogy();
rehisto_hPP_JetID->Draw("E1");
rehisto_hPP_JetID_fake->Draw("E1same");
rehisto_hPP_JetID_real->Draw("E1same");
TLatex *T1=new TLatex(0.5,0.92,"PYTHIA");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
TLatex T4;
T4.SetNDC();
T4.SetTextSize(0.03);
T4.SetTextFont(42);
T4.DrawLatex(0.65,0.28,Form("p_{T}: %.f-%.f (GeV/c)",xrange_pt[0],xrange_pt[1]));
/*T4.DrawLatex(0.65,0.28,Form("Inc: Mean=%.3f",hPP_JetID->GetMean()));
T4.DrawLatex(0.65,0.24,Form("Fake: Mean=%.3f",rehisto_hPP_JetID_fake->GetMean()));
T4.DrawLatex(0.65,0.2,Form("Fake: RMSError=%.3f",hPP_JetID_fake->GetRMSError()));
T4.DrawLatex(0.65,0.16,Form("Fake: MeanError=%.3f",hPP_JetID_fake->GetMeanError()));
T4.DrawLatex(0.65,0.12,Form("Fake: RMS=%.3f",hPP_JetID_fake->GetRMS()));*/
c2->cd(2)->SetLogy();
rehisto_hPPb_JetID->GetYaxis()->SetLabelSize(0);
rehisto_hPPb_JetID->Draw("E1");
rehisto_hPPb_JetID_fake->Draw("E1same");
rehisto_hPPb_JetID_real->Draw("E1same");
if(ilist==7)
TLegend *leg2=new TLegend(0.7,0.32,0.9,0.45);
else
TLegend *leg2=new TLegend(0.7,0.72,0.9,0.85);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetLineWidth(0);
leg2->SetTextSize(0.05);
leg2->AddEntry(rehisto_hPPb_JetID,"Inclusive","lp");
leg2->AddEntry(rehisto_hPPb_JetID_fake,"fake","lp");
leg2->AddEntry(rehisto_hPPb_JetID_real,"real","lp");
leg2->Draw("same");
TLatex *T2=new TLatex(0.4,0.92,"PYTHIA+HIJING");
T2->SetNDC();
TLatex T7;
T7.SetNDC();
T7.SetTextSize(0.03);
T7.SetTextFont(42);
T7.DrawLatex(0.65,0.28,Form("p_{T}: %.f-%.f (GeV/c)",xrange_pt[0],xrange_pt[1]));
/*T7.DrawLatex(0.65,0.28,Form("Inc: Mean=%.3f",hPPb_JetID->GetMean()));
T7.DrawLatex(0.65,0.24,Form("Fake: Mean=%.3f",hPPb_JetID_fake->GetMean()));
T7.DrawLatex(0.65,0.2,Form("Fake: RMSError=%.3f",hPPb_JetID_fake->GetRMSError()));
T7.DrawLatex(0.65,0.16,Form("Fake: MeanError=%.3f",hPPb_JetID_fake->GetMeanError()));
T7.DrawLatex(0.65,0.12,Form("Fake: RMS=%.3f",hPPb_JetID_fake->GetRMS()));*/
T2->SetTextAlign(12);
T2->SetTextSize(0.05);
T2->SetTextColor(1);
T2->SetTextFont(42);
T2->Draw("same");
c2->Print(Form("pic/JetIDcut/%s/%s/jetid.png",algo.Data(),JetIDName.Data()));

//---------------------------------------PP&PPb JetID f/i&r/i ratio before cut ------------------------------------------ 
c3 = new TCanvas("c3"," ",1200,600);
makeMultiPanelCanvas(c3,2,1,0,0,0.15,0.1,0.03);
c3->cd(1)->SetLogy();
ratio_hPP_JetID_fake=(TH1D*)rehisto_hPP_JetID_fake->Clone("ratio_hPP_JetID_fake");
ratio_hPP_JetID_fake->Divide(rehisto_hPP_JetID);
ratio_hPP_JetID_fake->GetYaxis()->SetTitle("ratio");
ratio_hPP_JetID_fake->SetTitle("");
ratio_hPP_JetID_fake->SetMaximum(5);
ratio_hPP_JetID_fake->SetMinimum(1e-5);
ratio_hPP_JetID_fake->Draw();
ratio_hPP_JetID_real=(TH1D*)rehisto_hPP_JetID_real->Clone("ratio_hPP_JetID_real");
ratio_hPP_JetID_real->Divide(rehisto_hPP_JetID);
ratio_hPP_JetID_real->Draw("same");
T1->Draw("same");
c3->cd(2)->SetLogy();
ratio_hPPb_JetID_fake=(TH1D*)rehisto_hPPb_JetID_fake->Clone("ratio_hPPb_JetID_fake");
ratio_hPPb_JetID_fake->Divide(rehisto_hPPb_JetID);
ratio_hPPb_JetID_fake->SetTitle("");
ratio_hPPb_JetID_fake->SetMaximum(5);
ratio_hPPb_JetID_fake->SetMinimum(1e-5);
ratio_hPPb_JetID_fake->Draw("same");
ratio_hPPb_JetID_real=(TH1D*)rehisto_hPPb_JetID_real->Clone("ratio_hPPb_JetID_real");
ratio_hPPb_JetID_real->Divide(rehisto_hPPb_JetID);
ratio_hPPb_JetID_real->Draw("same");
TLegend *leg3=new TLegend(0.35,0.10,0.65,0.25);
leg3->SetBorderSize(0);
leg3->SetFillColor(0);
leg3->SetLineWidth(0);
leg3->SetTextSize(0.05);
leg3->AddEntry(ratio_hPPb_JetID_fake,"fake Jet/Inc Jet","lp");
leg3->AddEntry(ratio_hPPb_JetID_real,"real Jet/Inc Jet","lp");
leg3->Draw("same");
T2->Draw("same");
c3->Print(Form("pic/JetIDcut/%s/%s/ratio_jetid.png",algo.Data(),JetIDName.Data()));

//---------------------------------------------PPb/Pb JetID distribution------------------------------------------------

TCanvas *c5=new TCanvas("c5");
c5->cd();
TH1D* ratio_JetID=(TH1D*)rehisto_hPPb_JetID->Clone();
ratio_JetID->Divide(rehisto_hPP_JetID);
ratio_JetID->SetTitle(Form("p_{T}^{jet}: %3.f--%3.f (GeV/c)",xrange_pt[0],xrange_pt[1]));
ratio_JetID->GetYaxis()->SetTitle("#frac{PPb}{PP} ratio");
ratio_JetID->SetMarkerStyle(20);
ratio_JetID->SetMarkerSize(1);
ratio_JetID->SetMinimum(0);
ratio_JetID->SetMaximum(2.0);
ratio_JetID->Draw();
/*for(int i=ratio_JetID->GetXaxis()->FindBin(xrange_JetID[0]);i<ratio_JetID->GetXaxis()->FindBin(xrange_JetID[1]);i++){
ymean_JetID+=ratio_JetID->GetBinContent(i)*rehisto_hPP_JetID->GetBinContent(i)*ratio_JetID->GetBinWidth(i)/rehisto_hPP_JetID->Integral("width");
}
ratio_JetID->Fit("pol0");
double ymean_JetID=pol0->GetParameter(0);
TPaveText *pt2=new TPaveText(0.7,0.8,0.9,0.9,"NDC");
pt2->AddText(Form("MeanY=%.4f",ymean_JetID));
pt2->Draw("same");
c5->Update();
*/
TLine *l2 =new TLine(c5->GetUxmin(),1,c5->GetUxmax(),1);
l2->SetLineStyle(2);
l2->SetLineColor(kBlue);
l2->Draw("same");
c5->Print(Form("pic/JetIDcut/%s/%s/ratio_jet%s_PPbvsPP.png",algo.Data(),JetIDName.Data(),JetIDName.Data()));
/*
//---------------------------------------------PPb/Pb pT distribution---------------------------------------------------

TCanvas *c6=new TCanvas("c6");
c6->cd();
TH1D* ratio_pt=(TH1D*)rehisto_hPPb_pt->Clone("ratio_pt");
ratio_pt->Divide(rehisto_hPP_pt);
ratio_pt->SetTitle(Form("%s: %.1f--%.1f %s",JetID.Data(), xrange_JetID[0],xrange_JetID[1], Unit.Data()));
ratio_pt->GetYaxis()->SetTitle("#frac{PPb}{PP} ratio");
ratio_pt->SetMarkerStyle(20);
ratio_pt->SetMarkerSize(1);
ratio_pt->SetMinimum(0);
ratio_pt->SetMaximum(2.0);
ratio_pt->Draw("E1");
ratio_pt->Fit("pol0");
double ymean_pt=pol0->GetParameter(0);
TPaveText *pt1=new TPaveText(0.7,0.8,0.9,0.9,"NDC");
pt1->AddText(Form("MeanY=%.4f",ymean_pt));
pt1->Draw("same");
TLine *l1 =new TLine(xrange_pt[0],1,xrange_pt[1],1);
l1->SetLineStyle(2);
l1->SetLineColor(kBlue);
l1->Draw("same");
c6->Print(Form("pic/JetIDcut/merged/ratio_jetpt_%.3f-%.3f_%.3f-%.3f_PPbvsPP.png",xrange_JetID[0],xrange_JetID[1],xrange_pt[0],xrange_pt[1]));
*/

//----------------------------------------------PP&PPb f&r&i cutvsnocut-------------------------------------------------

c6 = new TCanvas("c6"," ",1200,600);
makeMultiPanelCanvas(c6,2,1,0,0,0.15,0.1,0.03);
c6->cd(1);
TH1D* ratio_hPP_JetIDcutvsnocut_pt=(TH1D*)rehisto_hPP_JetIDcut_pt->Clone("ratio_hPP_JetIDcut_ptcutvsnocut_pt");
ratio_hPP_JetIDcutvsnocut_pt->Divide(rehisto_hPP_pt);
hFrame->GetYaxis()->SetTitle(Form("%s #frac{cut}{no cut}",JetID.Data()));
hFrame->SetMaximum(1.19);
hFrame->SetMinimum(0.01);
hFrame->DrawCopy();
ratio_hPP_JetIDcutvsnocut_pt->SetMarkerSize(1.5);
ratio_hPP_JetIDcutvsnocut_pt->SetMarkerStyle(24);
ratio_hPP_JetIDcutvsnocut_pt->SetMarkerColor(1);
ratio_hPP_JetIDcutvsnocut_pt->SetLineColor(1);
ratio_hPP_JetIDcutvsnocut_pt->Draw("same");
TH1D* ratio_hPP_JetIDcutvsnocut_pt_fake=(TH1D*)rehisto_hPP_JetIDcut_pt_fake->Clone("ratio_hPP_JetIDcutvsnocut_pt_fake");
TH1D* cl_rehisto_hPP_pt_fake=(TH1D*)rehisto_hPP_pt_fake->Clone("cl_rehisto_hPP_pt_fake");
cl_rehisto_hPP_pt_fake=(TH1D*)cl_rehisto_hPP_pt_fake->Rebin(Nbin_pt_coarse,"cl_rehisto_hPP_pt_fake",binbound_pt_coarse);
ratio_hPP_JetIDcutvsnocut_pt_fake=(TH1D*)ratio_hPP_JetIDcutvsnocut_pt_fake->Rebin(Nbin_pt_coarse,"ratio_hPP_JetIDcutvsnocut_pt_fake",binbound_pt_coarse);
normalizeByBinWidth(ratio_hPP_JetIDcutvsnocut_pt_fake);
normalizeByBinWidth(cl_rehisto_hPP_pt_fake);
ratio_hPP_JetIDcutvsnocut_pt_fake->Divide(cl_rehisto_hPP_pt_fake);
ratio_hPP_JetIDcutvsnocut_pt_fake->GetYaxis()->SetTitle("");
ratio_hPP_JetIDcutvsnocut_pt_fake->SetMarkerSize(1.5);
ratio_hPP_JetIDcutvsnocut_pt_fake->SetMarkerStyle(20);
ratio_hPP_JetIDcutvsnocut_pt_fake->SetMarkerColor(2);
ratio_hPP_JetIDcutvsnocut_pt_fake->SetLineColor(2);
ratio_hPP_JetIDcutvsnocut_pt_fake->Draw("same");
TH1D* ratio_hPP_JetIDcutvsnocut_pt_real=(TH1D*)rehisto_hPP_JetIDcut_pt_real->Clone("ratio_hPP_JetIDcutvsnocut_pt_real");
ratio_hPP_JetIDcutvsnocut_pt_real->Divide(rehisto_hPP_pt_real);
ratio_hPP_JetIDcutvsnocut_pt_real->GetYaxis()->SetTitle("");
ratio_hPP_JetIDcutvsnocut_pt_real->SetMarkerSize(1.5);
ratio_hPP_JetIDcutvsnocut_pt_real->SetMarkerStyle(29);
ratio_hPP_JetIDcutvsnocut_pt_real->SetMarkerColor(4);
ratio_hPP_JetIDcutvsnocut_pt_real->SetLineColor(4);
ratio_hPP_JetIDcutvsnocut_pt_real->Draw("same");
TLatex *T5=new TLatex(0.65,0.2,Form("cut: JetID = %.2f - %.2f", JetIDcut[0], JetIDcut[1]));
T5->SetNDC();
T5->SetTextAlign(12);
T5->SetTextSize(0.03);
T5->SetTextColor(1);
T5->SetTextFont(42);
T1->Draw("same");
T5->Draw("same");
TLine *l3 =new TLine(30,1,600,1);
l3->SetLineStyle(2);
l3->SetLineColor(kBlack);
l3->Draw("same");

c6->cd(2);
hFrame->GetYaxis()->SetTitle("");
hFrame->GetYaxis()->SetLabelSize(0);
hFrame->DrawCopy();
TH1D* ratio_hPPb_JetIDcutvsnocut_pt=(TH1D*)rehisto_hPPb_JetIDcut_pt->Clone("ratio_hPPb_JetIDcut_ptcutvsnocut_pt");
ratio_hPPb_JetIDcutvsnocut_pt->Divide(rehisto_hPPb_pt);
ratio_hPPb_JetIDcutvsnocut_pt->SetMarkerSize(1.5);
ratio_hPPb_JetIDcutvsnocut_pt->SetMarkerStyle(24);
ratio_hPPb_JetIDcutvsnocut_pt->SetMarkerColor(1);
ratio_hPPb_JetIDcutvsnocut_pt->SetLineColor(1);
ratio_hPPb_JetIDcutvsnocut_pt->Draw("same");
TH1D* ratio_hPPb_JetIDcutvsnocut_pt_fake=(TH1D*)rehisto_hPPb_JetIDcut_pt_fake->Clone("ratio_hPPb_JetIDcutvsnocut_pt_fake");
TH1D* cl_rehisto_hPPb_pt_fake=(TH1D*)rehisto_hPPb_pt_fake->Clone("cl_rehisto_hPPb_pt_fake");
cl_rehisto_hPPb_pt_fake=(TH1D*)cl_rehisto_hPPb_pt_fake->Rebin(Nbin_pt_coarse,"cl_rehisto_hPPb_pt_fake",binbound_pt_coarse);
ratio_hPPb_JetIDcutvsnocut_pt_fake=(TH1D*)ratio_hPPb_JetIDcutvsnocut_pt_fake->Rebin(Nbin_pt_coarse,"ratio_hPPb_JetIDcutvsnocut_pt_fake",binbound_pt_coarse);
normalizeByBinWidth(ratio_hPPb_JetIDcutvsnocut_pt_fake);
normalizeByBinWidth(cl_rehisto_hPPb_pt_fake);
ratio_hPPb_JetIDcutvsnocut_pt_fake->Divide(cl_rehisto_hPPb_pt_fake);
ratio_hPPb_JetIDcutvsnocut_pt_fake->GetYaxis()->SetTitle("");
ratio_hPPb_JetIDcutvsnocut_pt_fake->SetMarkerSize(1.5);
ratio_hPPb_JetIDcutvsnocut_pt_fake->SetMarkerStyle(20);
ratio_hPPb_JetIDcutvsnocut_pt_fake->SetMarkerColor(2);
ratio_hPPb_JetIDcutvsnocut_pt_fake->SetLineColor(2);
ratio_hPPb_JetIDcutvsnocut_pt_fake->Draw("same");
TH1D* ratio_hPPb_JetIDcutvsnocut_pt_real=(TH1D*)rehisto_hPPb_JetIDcut_pt_real->Clone("ratio_hPPb_JetIDcutvsnocut_pt_real");
ratio_hPPb_JetIDcutvsnocut_pt_real->Divide(rehisto_hPPb_pt_real);
ratio_hPPb_JetIDcutvsnocut_pt_real->GetYaxis()->SetTitle("");
ratio_hPPb_JetIDcutvsnocut_pt_real->SetMarkerSize(1.5);
ratio_hPPb_JetIDcutvsnocut_pt_real->SetMarkerStyle(29);
ratio_hPPb_JetIDcutvsnocut_pt_real->SetMarkerColor(4);
ratio_hPPb_JetIDcutvsnocut_pt_real->SetLineColor(4);
ratio_hPPb_JetIDcutvsnocut_pt_real->Draw("same");
TLegend *leg6=new TLegend(0.7,0.25,0.9,0.4);
leg6->SetBorderSize(0);
leg6->SetFillColor(0);
leg6->SetLineWidth(0);
leg6->SetTextSize(0.05);
leg6->AddEntry(rehisto_hPPb_JetID,"Inclusive","lp");
leg6->AddEntry(rehisto_hPPb_JetID_fake,"fake","lp");
leg6->AddEntry(rehisto_hPPb_JetID_real,"real","lp");
leg6->Draw("same");
T2->Draw("same");
T5->Draw("same");
l3->Draw("same");
c6->Print(Form("pic/JetIDcut/%s/%s/ratio_jetpt_cutvsnocut.png",algo.Data(),JetIDName.Data()));

//---------------------------------------PP&PPb pT f/i&r/i ratio before and after cut-----------------------------------

c7 = new TCanvas("c7"," ",800,800);
makeMultiPanelCanvas(c7,2,2,0,0,0.15,0.1,0.03);
c7->cd(1)->SetLogy();
TH1D* ratio_hPP_JetIDcut_pt_fakevsInc=(TH1D*)rehisto_hPP_JetIDcut_pt_fake->Clone("ratio_hPP_JetIDcut_pt_fakevsInc");
ratio_hPP_JetIDcut_pt_fakevsInc->Divide(rehisto_hPP_JetIDcut_pt);
ratio_hPP_JetIDcut_pt_fakevsInc->SetTitle(Form("PYTHIA cut: %.2f-%.2f", JetIDcut[0],JetIDcut[1]));
ratio_hPP_JetIDcut_pt_fakevsInc->SetMarkerStyle(20);
ratio_hPP_JetIDcut_pt_fakevsInc->SetMarkerSize(1.5);
ratio_hPP_JetIDcut_pt_fakevsInc->SetMarkerColor(2);
ratio_hPP_JetIDcut_pt_fakevsInc->SetLineColor(2);
ratio_hPP_JetIDcut_pt_fakevsInc->SetMaximum(5);
ratio_hPP_JetIDcut_pt_fakevsInc->SetMinimum(1e-6);
ratio_hPP_JetIDcut_pt_fakevsInc->Draw();

TH1D* ratio_hPP_JetIDcut_pt_realvsInc=(TH1D*)rehisto_hPP_JetIDcut_pt_real->Clone("ratio_hPP_JetIDcut_pt_realvsInc");
ratio_hPP_JetIDcut_pt_realvsInc->Divide(rehisto_hPP_JetIDcut_pt);
ratio_hPP_JetIDcut_pt_realvsInc->SetMarkerStyle(29);
ratio_hPP_JetIDcut_pt_realvsInc->SetMarkerSize(1.5);
ratio_hPP_JetIDcut_pt_realvsInc->SetMarkerColor(4);
ratio_hPP_JetIDcut_pt_realvsInc->SetLineColor(4);
ratio_hPP_JetIDcut_pt_realvsInc->Draw("same");

c7->cd(2)->SetLogy();
TH1D* ratio_hPP_pt_fakevsInc=(TH1D*)rehisto_hPP_pt_fake->Clone("ratio_hPP_pt_fakevsInc");
ratio_hPP_pt_fakevsInc->Divide(rehisto_hPP_pt);
ratio_hPP_pt_fakevsInc->SetTitle(Form("PYTHIA nocut %.2f-%.2f", xrange_JetID[0],xrange_JetID[1]));
ratio_hPP_pt_fakevsInc->SetMarkerStyle(20);
ratio_hPP_pt_fakevsInc->SetMarkerSize(1.5);
ratio_hPP_pt_fakevsInc->SetMarkerColor(2);
ratio_hPP_pt_fakevsInc->SetLineColor(2);
ratio_hPP_pt_fakevsInc->SetMaximum(5);
ratio_hPP_pt_fakevsInc->SetMinimum(1e-6);
ratio_hPP_pt_fakevsInc->Draw();

TH1D* ratio_hPP_pt_realvsInc=(TH1D*)rehisto_hPP_pt_real->Clone("ratio_hPP_pt_realvsInc");
ratio_hPP_pt_realvsInc->Divide(rehisto_hPP_pt);
ratio_hPP_pt_realvsInc->SetMarkerStyle(29);
ratio_hPP_pt_realvsInc->SetMarkerSize(1.5);
ratio_hPP_pt_realvsInc->SetMarkerColor(4);
ratio_hPP_pt_realvsInc->SetLineColor(4);
ratio_hPP_pt_realvsInc->Draw("same");

c7->cd(3)->SetLogy();
TH1D* ratio_hPPb_JetIDcut_pt_fakevsInc=(TH1D*)rehisto_hPPb_JetIDcut_pt_fake->Clone("ratio_hPPb_JetIDcut_pt_fakevsInc");
ratio_hPPb_JetIDcut_pt_fakevsInc->Divide(rehisto_hPPb_JetIDcut_pt);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetTitle(Form("PYTHIA+HIJING"));
ratio_hPPb_JetIDcut_pt_fakevsInc->SetMarkerStyle(20);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetMarkerSize(1.5);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetMarkerColor(2);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetLineColor(2);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetMaximum(5);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetMinimum(1e-6);
ratio_hPPb_JetIDcut_pt_fakevsInc->Draw();

TH1D* ratio_hPPb_JetIDcut_pt_realvsInc=(TH1D*)rehisto_hPPb_JetIDcut_pt_real->Clone("ratio_hPPb_JetIDcut_pt_realvsInc");
ratio_hPPb_JetIDcut_pt_realvsInc->Divide(rehisto_hPPb_JetIDcut_pt);
ratio_hPPb_JetIDcut_pt_realvsInc->SetMarkerStyle(29);
ratio_hPPb_JetIDcut_pt_realvsInc->SetMarkerSize(1.5);
ratio_hPPb_JetIDcut_pt_realvsInc->SetMarkerColor(4);
ratio_hPPb_JetIDcut_pt_realvsInc->SetLineColor(4);
ratio_hPPb_JetIDcut_pt_realvsInc->Draw("same");

c7->cd(4)->SetLogy();
TH1D* ratio_hPPb_pt_fakevsInc=(TH1D*)rehisto_hPPb_pt_fake->Clone("ratio_hPPb_pt_fakevsInc");
ratio_hPPb_pt_fakevsInc->Divide(rehisto_hPPb_pt);
ratio_hPPb_pt_fakevsInc->SetTitle(Form("PYTHIA+HIJING nocut"));
ratio_hPPb_pt_fakevsInc->SetMarkerStyle(20);
ratio_hPPb_pt_fakevsInc->SetMarkerSize(1.5);
ratio_hPPb_pt_fakevsInc->SetMarkerColor(2);
ratio_hPPb_pt_fakevsInc->SetLineColor(2);
ratio_hPPb_pt_fakevsInc->SetMaximum(5);
ratio_hPPb_pt_fakevsInc->SetMinimum(1e-6);
ratio_hPPb_pt_fakevsInc->Draw();

TH1D* ratio_hPPb_pt_realvsInc=(TH1D*)rehisto_hPPb_pt_real->Clone("ratio_hPPb_pt_realvsInc");
ratio_hPPb_pt_realvsInc->Divide(rehisto_hPPb_pt);
ratio_hPPb_pt_realvsInc->SetMarkerStyle(29);
ratio_hPPb_pt_realvsInc->SetMarkerSize(1.5);
ratio_hPPb_pt_realvsInc->SetMarkerColor(4);
ratio_hPPb_pt_realvsInc->SetLineColor(4);
ratio_hPPb_pt_realvsInc->Draw("same");
TLegend *leg7=new TLegend(0.65,0.25,0.9,0.35);
leg7->SetBorderSize(0);
leg7->SetFillColor(0);
leg7->SetLineWidth(0);
leg7->SetTextSize(0.05);
leg7->AddEntry(ratio_hPP_JetIDcut_pt_fakevsInc,"fake Jet/Inc Jet","lp");
leg7->AddEntry(ratio_hPP_JetIDcut_pt_realvsInc,"real Jet/Inc Jet","lp");
leg7->Draw("same");
c7->Print(Form("pic/JetIDcut/%s/%s/ratio_jetpt_%s.png",algo.Data(),JetIDName.Data(),JetIDName.Data()));

//-------------------------------PP&PPb pT distribution before and after cut-------------------------------------------
c8 = new TCanvas("c8"," ",1200,700);
makeMultiPanelCanvas(c8,2,2,0,0,0.15,0.1,0.03);
c8->cd(1)->SetLogy();
rehisto_hPP_pt->SetMarkerStyle(24);
rehisto_hPP_pt->SetMarkerSize(1.5);
rehisto_hPP_pt->SetMaximum(1e-2);
rehisto_hPP_pt->SetMinimum(1e-13);
rehisto_hPP_pt_real->SetMarkerStyle(29);
rehisto_hPP_pt_real->SetMarkerColor(4);
rehisto_hPP_pt_real->SetLineColor(4);
rehisto_hPP_pt_real->SetMarkerSize(1.5);
rehisto_hPP_pt_fake->SetMarkerStyle(20);
rehisto_hPP_pt_fake->SetMarkerColor(2);
rehisto_hPP_pt_fake->SetLineColor(2);
rehisto_hPP_pt_fake->SetMarkerSize(1.5);
rehisto_hPP_pt->Draw();
rehisto_hPP_pt_real->Draw("E1same");
rehisto_hPP_pt_fake->Draw("E1same");

c8->cd(2)->SetLogy();
rehisto_hPP_JetIDcut_pt->SetMarkerStyle(24);
rehisto_hPP_JetIDcut_pt->SetMarkerSize(1.5);
rehisto_hPP_JetIDcut_pt->SetMaximum(1e-2);
rehisto_hPP_JetIDcut_pt->SetMinimum(1e-13);
rehisto_hPP_JetIDcut_pt->GetYaxis()->SetTitle("");
rehisto_hPP_JetIDcut_pt_real->SetMarkerSize(1.5);
rehisto_hPP_JetIDcut_pt_real->SetMarkerStyle(29);
rehisto_hPP_JetIDcut_pt_real->SetMarkerColor(4);
rehisto_hPP_JetIDcut_pt_real->SetLineColor(4);
rehisto_hPP_JetIDcut_pt_fake->SetMarkerStyle(20);
rehisto_hPP_JetIDcut_pt_fake->SetMarkerColor(2);
rehisto_hPP_JetIDcut_pt_fake->SetLineColor(2);
rehisto_hPP_JetIDcut_pt_fake->SetMarkerSize(1.5);
rehisto_hPP_JetIDcut_pt->Draw();
rehisto_hPP_JetIDcut_pt_real->Draw("E1same");
rehisto_hPP_JetIDcut_pt_fake->Draw("E1same");

c8->cd(3)->SetLogy();
rehisto_hPPb_pt->SetMarkerStyle(24);
rehisto_hPPb_pt->SetMarkerSize(1.5);
rehisto_hPPb_pt->SetMaximum(1e-2);
rehisto_hPPb_pt->SetMinimum(1e-13);
rehisto_hPPb_pt->GetYaxis()->SetTitle("Probability");
rehisto_hPPb_pt_real->SetMarkerStyle(29);
rehisto_hPPb_pt_real->SetMarkerColor(4);
rehisto_hPPb_pt_real->SetLineColor(4);
rehisto_hPPb_pt_real->SetMarkerSize(1.5);
rehisto_hPPb_pt_fake->SetMarkerStyle(20);
rehisto_hPPb_pt_fake->SetMarkerColor(2);
rehisto_hPPb_pt_fake->SetLineColor(2);
rehisto_hPPb_pt_fake->SetMarkerSize(1.5);
rehisto_hPPb_pt->Draw();
rehisto_hPPb_pt_real->Draw("E1same");
rehisto_hPPb_pt_fake->Draw("E1same");

c8->cd(4)->SetLogy();
rehisto_hPPb_JetIDcut_pt->SetMarkerStyle(24);
rehisto_hPPb_JetIDcut_pt->SetMarkerSize(1.5);
rehisto_hPPb_JetIDcut_pt->SetMaximum(1e-2);
rehisto_hPPb_JetIDcut_pt->SetMinimum(1e-13);
rehisto_hPPb_JetIDcut_pt_real->SetMarkerSize(1.5);
rehisto_hPPb_JetIDcut_pt_real->SetMarkerStyle(29);
rehisto_hPPb_JetIDcut_pt_real->SetMarkerColor(4);
rehisto_hPPb_JetIDcut_pt_real->SetLineColor(4);
rehisto_hPPb_JetIDcut_pt_fake->SetMarkerStyle(20);
rehisto_hPPb_JetIDcut_pt_fake->SetMarkerColor(2);
rehisto_hPPb_JetIDcut_pt_fake->SetLineColor(2);
rehisto_hPPb_JetIDcut_pt_fake->SetMarkerSize(1.5);
rehisto_hPPb_JetIDcut_pt->Draw();
rehisto_hPPb_JetIDcut_pt_real->Draw("E1same");
rehisto_hPPb_JetIDcut_pt_fake->Draw("E1same");
TLegend *leg8= new TLegend(0.7,0.75,0.9,0.9);
leg8->SetBorderSize(0);
leg8->SetFillColor(0);
leg8->SetTextSize(0.05);
leg8->SetLineWidth(0);
leg8->AddEntry(rehisto_hPP_pt,"Inclusive","lp");
leg8->AddEntry(rehisto_hPP_pt_real,"real","lp");
leg8->AddEntry(rehisto_hPP_pt_fake,"fake","lp");
leg8->Draw("same");
c8->Print(Form("pic/JetIDcut/%s/%s/jetpt.png",algo.Data(),JetIDName.Data()));

//--------------------------------------------------PP&PPb JetID,pT Profile histograms-----------------------------------

c9 = new TCanvas("c9"," ",1200,600);
makeMultiPanelCanvas(c9,2,1,0,0,0.15,0.1,0.03);
c9->cd(1);
hFrame->GetYaxis()->SetLabelSize(20);
hFrame->GetYaxis()->SetTitle(Form("%s %s",JetID.Data(), Unit.Data()));
hFrame->SetMaximum(TMath::Max(rehisto_hProfPP_fake->GetMaximum(),rehisto_hProfPP_real->GetMaximum())*1.5);
hFrame->SetMinimum(0);
hFrame->DrawCopy();
rehisto_hProfPP->Draw("same");
rehisto_hProfPP_fake->Draw("same");
rehisto_hProfPP_real->Draw("same");
T1->Draw("same");
c9->cd(2);
hFrame->GetYaxis()->SetLabelSize(0);
hFrame->GetYaxis()->SetTitle("");
hFrame->DrawCopy();
rehisto_hProfPPb->Draw("same");
rehisto_hProfPPb_fake->Draw("same");
rehisto_hProfPPb_real->Draw("same");
TLegend *leg9=new TLegend(0.7,0.7,0.9,0.85);
leg9->SetBorderSize(0);
leg9->SetFillColor(0);
leg9->SetLineWidth(0);
leg9->SetTextSize(0.05);
leg9->AddEntry(rehisto_hProfPP,"Inclusive","lp");
leg9->AddEntry(rehisto_hProfPP_fake,"fake","lp");
leg9->AddEntry(rehisto_hProfPP_real,"real","lp");
leg9->Draw("same");
T2->Draw("same");

c9->Print(Form("pic/JetIDcut/%s/%s/Profile.png",algo.Data(),JetIDName.Data()));

//-------------------------------------------------PP&PPb JetEtaphi 2D histograms----------------------------------------

c10 = new TCanvas("c10","",1200,600);
makeMultiPanelCanvas(c10,2,1,0.12,0.1,0.1,0.1,0.1);
c10->cd(1);
hPP_Etaphi->Draw("colz");
c10->cd(2);
hPPb_Etaphi->Draw("colz");

c10->Print(Form("pic/JetIDcut/%s/jetEtaphi.png",algo.Data()));

}
