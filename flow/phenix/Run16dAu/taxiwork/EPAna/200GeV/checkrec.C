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

void checkrec(int icent = 0, int ihar=1){
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
    float Meanx[nsub][1000];
    float Meany[nsub][1000];
    float RMSx[nsub][1000];
    float RMSy[nsub][1000];
    float runlist[1000];
        TH1F* hunqx = new TH1F("unqx","hunqx",820,-4.1,4.1);
        TH1F* hqx = new TH1F("qx","hqx",820,-4.1,4.1);
        TH1F* hunqy = new TH1F("unqy","hunqy",820,-4.1,4.1);
        TH1F* hqy = new TH1F("qy","hqy",820,-4.1,4.1);
        TH1F* h = new TH1F("h","",1200,-6,6);
     for(int irun=0;irun<nrun;irun++){
        std::ifstream corrs("Run16dAu200GeV.lst");
        int index=0; int run=0;
        for(int jrun=0;jrun<irun+1;jrun++){
        corrs>>index>>run;
        }
        runlist[irun] = run;
         fin = TFile::Open(Form("/store/user/qixu/flow/Run16dAu/200GeV/%s",GetRun(irun).Data()));
        for(int isub=0;isub<nsub;isub++){
        str = choosesub(isub);
        if(str=="ABORT") continue;
        hunqx->Reset();
        hqx->Reset();
        hunqy->Reset();
        hqy->Reset();
        for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
          TH1F* hunqxtemp = (TH1F*)fin->Get(Form("q_%d_%d_%d_%d_0",icent,ibbcz,ihar,isub));
          TH1F* hqxtemp = (TH1F*)fin->Get(Form("qRec_%d_%d_%d_%d_0",icent,ibbcz,ihar,isub));
          TH1F* hunqytemp = (TH1F*)fin->Get(Form("q_%d_%d_%d_%d_1",icent,ibbcz,ihar,isub));
          TH1F* hqytemp = (TH1F*)fin->Get(Form("qRec_%d_%d_%d_%d_1",icent,ibbcz,ihar,isub));
          hunqx->Add(hunqxtemp);
          hqx->Add(hqxtemp);
          hunqy->Add(hunqytemp);
          hqy->Add(hqytemp);
        }
        TCanvas *c1 = new TCanvas("c1","c1",600,450);
        c1->SetLogy();
        h->GetXaxis()->SetTitle(Form("Qx_{%d}",n));
        h->GetYaxis()->SetTitle(Form("# of events"));
        h->GetXaxis()->SetRangeUser(-4.1,4.1);
        h->GetYaxis()->SetRangeUser(1,hunqx->GetBinContent(hunqx->GetMaximumBin())*11);
        h->Draw();
	hqx->SetMarkerStyle(20);
	hqx->SetMarkerSize(0.8);
	hqx->SetMarkerColor(2);
	hqx->SetLineColor(2);
	hunqx->SetMarkerStyle(20);
	hunqx->SetMarkerSize(0.8);
	hunqx->SetMarkerColor(1);
	hunqx->SetLineColor(1);
        Meanx[isub][irun] = hqx->GetMean();
        RMSx[isub][irun] = hqx->GetRMS();
        TLegend *leg = new TLegend(0.12,0.75,0.5,0.88);
        leg->SetTextSize(0.035);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hunqx,Form("Run%d %s",run,str.Data()));
        leg->AddEntry(hunqx,"before Recentering","L");
        leg->AddEntry(hqx,"After Recentering","L");
        hunqx->Draw("Csame");
        hqx->Draw("Csame");
        leg->Draw("same");
        c1->Print(Form("run-by-run/v%d/%s/checkrec_%d_%d_%d_%d_qx.png",n,str.Data(),icent,ihar,isub,irun));
        delete c1;
        TCanvas *c2 = new TCanvas("c2","c2",600,450);
        c2->SetLogy();
        h->GetXaxis()->SetTitle(Form("Qy_{%d}",n));
        h->GetYaxis()->SetTitle(Form("# of events"));
        h->GetXaxis()->SetRangeUser(-4.1,4.1);
        h->GetYaxis()->SetRangeUser(1,hunqy->GetBinContent(hunqy->GetMaximumBin())*11);
        h->DrawCopy();
	hqy->SetMarkerStyle(20);
	hqy->SetMarkerSize(0.8);
	hqy->SetMarkerColor(2);
	hqy->SetLineColor(2);
	hunqy->SetMarkerStyle(20);
	hunqy->SetMarkerSize(0.8);
	hunqy->SetMarkerColor(1);
	hunqy->SetLineColor(1);
        Meany[isub][irun] = hqy->GetMean();
        RMSy[isub][irun] = hqy->GetRMS();
        TLegend *leg = new TLegend(0.12,0.75,0.5,0.88);
        leg->SetTextSize(0.035);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hunqy,Form("Run%d %s",run,str.Data()));
        leg->AddEntry(hunqy,"before Recentering","L");
        leg->AddEntry(hqy,"After Recentering","L");
        hunqy->Draw("Csame");
        hqy->Draw("Csame");
        leg->Draw("same");
        c2->Print(Form("run-by-run/v%d/%s/checkrec_%d_%d_%d_%d_qy.png",n,str.Data(),icent,ihar,isub,irun));
        delete c2;
        }
        fin->Close();
     }
     TGraphErrors *grRunqx[nsub];
     TGraphErrors *grRunqy[nsub];
     TGraphErrors *grRunqx1[nsub];
     TGraphErrors *grRunqy1[nsub];
    TCanvas *c1 = new TCanvas("c1","c1",600,450);
    TCanvas *c2 = new TCanvas("c2","c2",600,450);
    TH1F* h1 = new TH1F("h1","",10000,450000,460000);
    h1->GetXaxis()->SetTitle(Form("Run Number"));
    h1->GetYaxis()->SetTitle(Form("Mean Qx_{%d}",n));
    h1->GetXaxis()->SetRangeUser(runlist[0]-5,runlist[nrun-1]+5);
    h1->GetYaxis()->SetRangeUser(-1,2);

        for(int isub=0;isub<nsub;isub++){
    TLegend *leg = new TLegend(0.35,0.65,0.55,0.85);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    c1->cd();
    h1->DrawCopy("0");
            grRunqx[isub] = new TGraphErrors(nrun,runlist,Meanx[isub],0,0);
            grRunqx1[isub] = new TGraphErrors(nrun,runlist,RMSx[isub],0,0);
            grRunqx[isub] -> SetMarkerColor(color[0]);
            grRunqx1[isub] -> SetMarkerColor(color[0]);
            grRunqx[isub] -> SetMarkerSize(1.0);
            grRunqx[isub] -> SetMarkerStyle(20);
            grRunqx1[isub] -> SetMarkerSize(1.0);
            grRunqx1[isub] -> SetMarkerStyle(20);
            grRunqx[isub] -> SetLineColor(color[0]);
            grRunqx1[isub] -> SetLineColor(color[0]);
            grRunqx[isub] -> Draw("Psame");
            grRunqx1[isub] -> Draw("Psame");
            leg->AddEntry(grRunqx[isub],choosesub(isub).Data(),"P");
            leg->Draw("same");
            c1->Print(Form("run-by-run/v%d/fitsumQx%d.png",n,isub));
            }

        for(int isub=0;isub<nsub;isub++){
    TLegend *leg = new TLegend(0.35,0.65,0.55,0.85);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    c2->cd();
    h1->GetYaxis()->SetTitle(Form("Mean Qy_{%d}",n));
    h1->DrawCopy("0");
            grRunqy[isub] = new TGraphErrors(nrun,runlist,Meany[isub],0,0);
            grRunqy1[isub] = new TGraphErrors(nrun,runlist,RMSy[isub],0,0);
            grRunqy[isub] -> SetMarkerColor(color[0]);
            grRunqy1[isub] -> SetMarkerColor(color[0]);
            grRunqy[isub] -> SetMarkerSize(1.0);
            grRunqy[isub] -> SetMarkerStyle(20);
            grRunqy1[isub] -> SetMarkerSize(1.0);
            grRunqy1[isub] -> SetMarkerStyle(20);
            grRunqy[isub] -> SetLineColor(color[0]);
            grRunqy1[isub] -> SetLineColor(color[0]);
            grRunqy[isub] -> Draw("Psame");
            grRunqy1[isub] -> Draw("Psame");
            leg->AddEntry(grRunqy[isub],choosesub(isub).Data(),"P");
            leg->Draw("same");
            c2->Print(Form("run-by-run/v%d/fitsumQy%d.png",n,isub));
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
