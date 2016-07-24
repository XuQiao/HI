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
#include "JetID.h"

bool isMC=kTRUE;

TH1* myRebin(TH1* histo,int Nbin, double binbound[]){
//rehisto->Scale(1/rehisto->Integral());

if(binbound[Nbin]==1000){	//rebin pt
//TH1* rehisto=histo->Clone("rehisto");
cout<<histo->GetEntries()<<'\t';
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
for(int i=1;i<Nbin+1;i++){
TString histoname=Form("%s",histo->GetName());
if(!histoname.Contains("Prof")){
rehisto->SetBinContent(i,rehisto->GetBinContent(i)/rehisto->GetBinWidth(i));
rehisto->SetBinError(i,rehisto->GetBinError(i)/rehisto->GetBinWidth(i));
}
}
rehisto->GetXaxis()->SetRangeUser(binbound[5],binbound[Nbin-1]);
cout<<rehisto->GetEntries()<<endl;
}
else{
TH1D* rehisto=histo->Clone(histo->GetName());
rehisto->GetXaxis()->SetRangeUser(binbound[0],binbound[Nbin]);
}

fixedFontHist(rehisto);

return rehisto;
}

void check1(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);

//------------------------------------------------------------Get Histograms---------------------------------------------

TFile *fdata=TFile::Open(datafile);
TFile *fPP=TFile::Open(PPfile);
TFile *fPPb=TFile::Open(PPbfile);

//TString JetIDNameList[19]={"chSum", "neuSum", "phoSum", "chSumpt", "neuSumpt", "phoSumpt", "chMax", "neuMax", "phoMax", "chMaxpt", "neuMaxpt", "phoMaxpt","neuMaxr","chN","neuN","phoN","SumSumrawpt","SumSumpt","rawptpt"};
int list=11;
//int list=atoi(getenv("LIST"));
TString JetIDName=JetIDNameList[list];

if(JetIDName=="chMax"){
TString JetID = "chargedMax";
TString Unit = "[GeV/c]";
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
double xrange_JetIDcut[2]={2+1e-4,200-1e-4};
}

else if(JetIDName=="chMaxpt"){
TString JetID = "chargedMax/p_{T}";
TString Unit = "";
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
double xrange_JetIDcut[2]={0.02+1e-4,1-1e-4};
}

else if(JetIDName=="phoMax"){
TString JetID = "photonMax";
TString Unit = "[GeV/c]";
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
double xrange_JetIDcut[2]={1+1e-4,200-1e-4};
}

else if(JetIDName=="phoMaxpt"){
TString JetID = "photonMax/p_{T}";
TString Unit = "";
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
double xrange_JetIDcut[2]={0.02+1e-4,1-1e-4};
}

else if(JetIDName=="neuMax"){
TString JetID = "neutralMax";
TString Unit = "[GeV/c]";
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
double xrange_JetIDcut[2]={0.5+1e-4,200-1e-4};
}

else if(JetIDName=="neuMaxpt"){
TString JetID = "neutralMax/p_{T}";
TString Unit = "";
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
double xrange_JetIDcut[2]={0.01+1e-4,1-1e-4};
}

else if(JetIDName=="chSum"){
TString JetID = "chargedSum";
TString Unit = "[GeV/c]";
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
double xrange_JetIDcut[2]={10+1e-4,200-1e-4};
}

else if(JetIDName=="chSumpt"){
TString JetID = "chargedSum/p_{T}";
TString Unit = "";
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
double xrange_JetIDcut[2]={0.02+1e-4,1-1e-4};
}

else if(JetIDName=="neuSum"){
TString JetID = "neutralSum";
TString Unit = "[GeV/c]";
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
double xrange_JetIDcut[2]={9+1e-4,200-1e-4};
}

else if(JetIDName=="neuSumpt"){
TString JetID = "neutralSum/p_{T}";
TString Unit = "";
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
double xrange_JetIDcut[2]={0.01+1e-4,1-1e-4};
}

else if(JetIDName=="phoSum"){
TString JetID = "photonSum";
TString Unit = "[GeV/c]";
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
double xrange_JetIDcut[2]={10+1e-4,200-1e-4};
}

else if(JetIDName=="phoSumpt"){
TString JetID = "photonSum/p_{T}";
TString Unit = "";
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
double xrange_JetIDcut[2]={0.02+1e-4,1-1e-4};
}

else if(JetIDName=="neuMaxr"){
TString JetID = "neutralMax/(neutralMax+chargedMax+photonMax)";
TString Unit = "";
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};
double xrange_JetIDcut[2]={0.01+1e-4,1-1e-4};
}

else if(JetIDName=="SumSumrawpt"){
TString JetID = "(chargedSum+photonSum+neutralSum+muSum+eSum)/raw p_{T}";
TString Unit = "";
double binbound_JetID[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.};
double xrange_JetIDcut[2]={0.00+1e-4,1.15-1e-4};
}

else if(JetIDName=="SumSumpt"){
TString JetID = "(chargedSum+photonSum+neutralSum+muSum+eSum)/p_{T}";
TString Unit = "";
double binbound_JetID[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.};
double xrange_JetIDcut[2]={0.85+1e-4,2-1e-4};
}

else if(JetIDName=="rawptpt"){
TString JetID = "raw p_{T}/p_{T}";
TString Unit = "";
double binbound_JetID[]={0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.08,1.1,1.12,1.14,1.16,1.18,1.2};
double xrange_JetIDcut[2]={0.8+1e-4,1.2-1e-4};
}

else{
exit();
}

TH2F *hdata2D= (TH2F*)fdata->Get(Form("jetpt%s_0-100%%",JetIDName.Data()));
TH2F *hPP2D= (TH2F*)fPP->Get(Form("jetpt%s_0-100%%",JetIDName.Data()));
TH2F *hPPb2D= (TH2F*)fPPb->Get(Form("jetpt%s_0-100%%",JetIDName.Data()));
TH2F *hPP2D_real= (TH2F*)fPP->Get(Form("jetpt%s_real_0-100%%",JetIDName.Data()));
TH2F *hPPb2D_real= (TH2F*)fPPb->Get(Form("jetpt%s_real_0-100%%",JetIDName.Data()));
TH2F *hPP2D_fake= (TH2F*)fPP->Get(Form("jetpt%s_fake_0-100%%",JetIDName.Data()));
TH2F *hPPb2D_fake= (TH2F*)fPPb->Get(Form("jetpt%s_fake_0-100%%",JetIDName.Data()));
TH2F *hPP_Etaphi= (TH2F*)fPP->Get(Form("jetEtaphi_0-100%%"));

//double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,135,150,165,180,205,220,260,300,350,500};
double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
double xrange_pt[2]={binbound_pt[5]+1e-4,binbound_pt[Nbin_pt]-1e-4};
double xrange_JetID[2]={binbound_JetID[0]+1e-4,binbound_JetID[Nbin_JetID]-1e-4};

TProfile *hProfPP=(TProfile*)hPP2D->ProfileX("hProfPP",0,-1,"e");
TProfile *hProfPP_fake=(TProfile*)hPP2D_fake->ProfileX("hProfPP_fake",0,-1,"e");
TProfile *hProfPP_real=(TProfile*)hPP2D_real->ProfileX("hProfPP_real",0,-1,"e");
TProfile *hProfPPb=(TProfile*)hPPb2D->ProfileX("hProfPPb",0,-1,"e");
TProfile *hProfPPb_fake=(TProfile*)hPPb2D_fake->ProfileX("hProfPPb_fake",0,-1,"e");
TProfile *hProfPPb_real=(TProfile*)hPPb2D_real->ProfileX("hProfPPb_real",0,-1,"e");

TProfile *rehisto_hProfPP=(TProfile*)myRebin(hProfPP,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPP_fake=(TProfile*)myRebin(hProfPP_fake,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPP_real=(TProfile*)myRebin(hProfPP_real,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPPb=(TProfile*)myRebin(hProfPPb,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPPb_fake=(TProfile*)myRebin(hProfPPb_fake,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPPb_real=(TProfile*)myRebin(hProfPPb_real,Nbin_pt,binbound_pt);
/*
TProfile* rehisto_hProfPP=(TProfile*)hProfPP->Clone(hProfPP->GetName());
TProfile* rehisto_hProfPP_fake=(TProfile*)hProfPP_fake->Clone(hProfPP_fake->GetName());
TProfile* rehisto_hProfPP_real=(TProfile*)hProfPP_real->Clone(hProfPP_real->GetName());
TProfile* rehisto_hProfPPb=(TProfile*)hProfPPb->Clone(hProfPPb->GetName());
TProfile* rehisto_hProfPPb_fake=(TProfile*)hProfPPb_fake->Clone(hProfPPb_fake->GetName());
TProfile* rehisto_hProfPPb_real=(TProfile*)hProfPPb_real->Clone(hProfPPb_real->GetName());
*/
rehisto_hProfPP->SetMarkerSize(1.3);
rehisto_hProfPP->SetMaximum(TMath::Max(rehisto_hProfPP_real->GetMaximum(),rehisto_hProfPP_fake->GetMaximum())*1.2);
rehisto_hProfPP->SetMinimum(0);
rehisto_hProfPP_fake->SetMarkerSize(1.3);
rehisto_hProfPP_real->SetMarkerSize(1.3);
rehisto_hProfPP->SetMarkerStyle(24);
rehisto_hProfPP_fake->SetMarkerStyle(20);
rehisto_hProfPP_real->SetMarkerStyle(29);
rehisto_hProfPP->SetMarkerColor(1);
rehisto_hProfPP_fake->SetMarkerColor(2);
rehisto_hProfPP_real->SetMarkerColor(4);
rehisto_hProfPP->SetTitle(Form("PYTHIA"));
rehisto_hProfPP->GetXaxis()->SetTitle("p_{T} [GeV/c]");
rehisto_hProfPP->GetYaxis()->SetTitle(Form("%s %s",JetID.Data(), Unit.Data()));

rehisto_hProfPPb->SetTitle(Form("PYTHIA+HIJING"));
rehisto_hProfPPb->GetXaxis()->SetTitle("p_{T} [GeV/c]");
rehisto_hProfPPb->GetYaxis()->SetTitle("");
rehisto_hProfPPb->SetMaximum(rehisto_hProfPP->GetMaximum());
rehisto_hProfPPb->SetMinimum(0);
rehisto_hProfPPb->SetMarkerSize(1.3);
rehisto_hProfPPb_fake->SetMarkerSize(1.3);
rehisto_hProfPPb_real->SetMarkerSize(1.3);
rehisto_hProfPPb->SetMarkerStyle(24);
rehisto_hProfPPb_fake->SetMarkerStyle(20);
rehisto_hProfPPb_real->SetMarkerStyle(29);

rehisto_hProfPPb->SetMarkerColor(1);
rehisto_hProfPPb_fake->SetMarkerColor(2);
rehisto_hProfPPb_real->SetMarkerColor(4);

fixedFontHist(rehisto_hProfPP);
fixedFontHist(rehisto_hProfPPb);


//--------------------------------------------------PP&PPb JetID,pT Profile histograms-----------------------------------

c9 = new TCanvas("c9"," ",1200,600);
makeMultiPanelCanvas(c9,2,1,0.15,0.1,0.15,0.1,0.03);
c9->cd(1);
rehisto_hProfPP->Draw();
rehisto_hProfPP_fake->Draw("same");
rehisto_hProfPP_real->Draw("same");

c9->cd(2);
rehisto_hProfPPb->Draw();
rehisto_hProfPPb_fake->Draw("same");
rehisto_hProfPPb_real->Draw("same");
TLegend *leg9=new TLegend(0.7,0.8,0.85,0.9);
leg9->SetLineWidth(0);
leg9->AddEntry(rehisto_hProfPP,"Inclusive","lp");
leg9->AddEntry(rehisto_hProfPP_fake,"fake","lp");
leg9->AddEntry(rehisto_hProfPP_real,"real","lp");
leg9->Draw("same");
c9->Print(Form("pic/JetIDcut/%s_%s.png",JetIDName.Data(),ff.Data()));
}
