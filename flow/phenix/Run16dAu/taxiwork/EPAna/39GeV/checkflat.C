#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "RpPar.h"

using namespace std;

int GetTotalRun();
void FillGoodRun();
float GetReso(int , int , int, bool);
TString GetRun(int);
float eps = 1e-4;
float GoodRunFit[ncent][nhar][nsub][10000];

TString choosesub(int isub){
    TString str;
        if(isub==0)
         str = "FVTX1LS";
        else if(isub==1)
         str = "FVTX2LS";
        else if(isub==2)
         str = "FVTX3LS";
        else if(isub==3)
         str = "FVTX4LS";
        else if(isub==5)
         str = "FVTX1S";
        else if(isub==4)
         str = "BBCS";
        else if(isub==6)
         str = "FVTX1p2p3LS";
        else if(isub==7)
         str = "FVTX1p2p4LS";
        else if(isub==8)
         str = "FVTX1N";
        else if(isub==9)
         str = "FVTXtrkS";
        else if(isub==10)
         str = "FVTX1LN";
        else if(isub==11)
         str = "FVTX2LN";
        else if(isub==12)
         str = "FVTX3LN";
        else if(isub==13)
         str = "FVTX4LN";
        else
         str = "ABORT";
    return str;
}

void checkflat(int icent = 0, int ihar=1){
    int color[nsub]={1,2,3,4,5,6,7,8};
    float pi = acos(-1.0);
    TString str;
    TFile *fin;
        TH1::SetDefaultSumw2();
        gStyle->SetOptStat(kFALSE);
        int iharE=0;
     if(nhar==1||nhar==2) iharE=1;
        int n = ihar+1.0+iharE;
    int nrun = GetTotalRun();
    float FitGood[nsub][1000];
    float FitGooderr[nsub][1000];
    float runlist[1000];
        TH1F* hpsi = new TH1F("psi","psi",100,-pi,pi);
        TH1F* hunpsi = new TH1F("unpsi","unpsi",100,-pi,pi);
        TH1F* h = new TH1F("h","",100,-4,4);
     for(int irun=0;irun<nrun;irun++){
        std::ifstream corrs("Run16dAu39GeV.lst");
        int index=0; int run=0;
        for(int jrun=0;jrun<irun+1;jrun++){
        corrs>>index>>run;
        }
        runlist[irun] = run;
         fin = TFile::Open(Form("/store/user/qixu/flow/Run16dAu/39GeV/%s",GetRun(irun).Data()));
        for(int isub=0;isub<nsub;isub++){
        str = choosesub(isub);
        if(str=="ABORT") continue;
        hpsi->Reset();
        hunpsi->Reset();
        for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
          TH1F* hunpsitemp = (TH1F*)fin->Get(Form("psi_%d_%d_%d_%d",icent,ibbcz,ihar,isub));
          TH1F* hpsitemp = (TH1F*)fin->Get(Form("psiFla_%d_%d_%d_%d",icent,ibbcz,ihar,isub));
          hpsi->Add(hpsitemp);
          hunpsi->Add(hunpsitemp);
        }
        TF1 *fun = new TF1("fun","pol0",-pi,pi);
        TCanvas *c1 = new TCanvas("c1","c1",600,450);
        h->GetXaxis()->SetTitle(Form("#Psi_{%d}",n));
        h->GetYaxis()->SetTitle(Form("# of events"));
        h->GetYaxis()->SetRangeUser(hunpsi->GetBinContent(hunpsi->GetMinimumBin())*0.9,hunpsi->GetBinContent(hunpsi->GetMaximumBin())*1.1);
        h->Draw("C");
       // h->GetYaxis()->SetRangeUser(0,hunpsi->GetBinContent(hunpsi->GetMaximumBin())*1.2);
      if(hpsi->GetEntries()>1000){
	hpsi->Fit("fun","QR0");
	float par=fun->GetParameter(0);
	float parerr=fun->GetParError(0);
        FitGood[isub][irun]=fun->GetChisquare()/fun->GetNDF();
        FitGooderr[isub][irun]=0;
        fun->SetLineColor(4);
        fun->Draw("same");
      }
        else{
        FitGood[isub][irun]=0;
        FitGooderr[isub][irun]=0;
        }
	hpsi->SetMarkerStyle(20);
	hpsi->SetMarkerSize(0.8);
	hpsi->SetMarkerColor(2);
	hpsi->SetLineColor(2);
	hunpsi->SetMarkerStyle(20);
	hunpsi->SetMarkerSize(0.8);
	hunpsi->SetMarkerColor(1);
	hunpsi->SetLineColor(1);
        TLegend *leg = new TLegend(0.4,0.7,0.8,0.85);
        leg->SetTextSize(0.04);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hunpsi,Form("Run%d %s",run,str.Data()));
        leg->AddEntry(hunpsi,"before flattening","L");
        leg->AddEntry(hpsi,"After flattening","L");
        hpsi->Draw("Csame");
        hunpsi->Draw("Csame");
        leg->Draw("same");
        c1->Print(Form("run-by-run/v%d/%s/checkflat_%d_%d_%d_%d.png",n,str.Data(),icent,ihar,isub,irun));
        delete c1;
        }
        fin->Close();
     }
     TGraphErrors *grRun[nsub];
    TCanvas *c1 = new TCanvas("c1","c1",600,450);
    TCanvas *c2 = new TCanvas("c2","c2",600,450);
    TH1F* h1 = new TH1F("h1","",10000,450000,460000);
    h1->GetXaxis()->SetTitle(Form("Run Number"));
    h1->GetYaxis()->SetTitle("Fit #chi^{2}/Ndf");
    h1->GetXaxis()->SetRangeUser(runlist[0]-5,runlist[nrun-1]+5);
    h1->GetYaxis()->SetRangeUser(0,2);
        for(int isub=0;isub<nsub;isub++){
    TLegend *leg = new TLegend(0.65,0.75,0.88,0.85);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    c1->cd();
    h1->Draw("0");
            grRun[isub] = new TGraphErrors(nrun,runlist,FitGood[isub],0,FitGooderr[isub]);
            grRun[isub] -> SetMarkerColor(color[0]);
            grRun[isub] -> SetMarkerSize(1.);
            grRun[isub] -> SetMarkerStyle(20);
            grRun[isub] -> SetLineColor(color[0]);
            grRun[isub] -> Draw("Psame");
            leg->AddEntry(grRun[isub],choosesub(isub).Data(),"P");
        leg->Draw("same");
      c1->Print(Form("run-by-run/v%d/fitsum%d.png",n,isub));
        }
}

int GetTotalRun(){
    ifstream frun(dataset+".Lst");
    string RunNumber;
    int nrun = -1;
    while(!frun.eof()){
        frun>>RunNumber;
        nrun++;
    }
    return nrun;
}


TString GetRun(int irun){
    ifstream frun(dataset+".Lst");
    string RunNumbertemp;
    int jrun = 0;
    while(!frun.eof() && jrun<=irun){
        frun>>RunNumbertemp;
        jrun++;
    }
    TString RunNumber(RunNumbertemp);
    return RunNumber;
}
