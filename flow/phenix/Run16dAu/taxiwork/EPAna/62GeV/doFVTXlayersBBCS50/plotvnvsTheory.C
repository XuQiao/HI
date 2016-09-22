#include "SimplifyLife.C"
#include "RpPar.h"

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
        else
         str = "ABORT";
    return str;
}

void plotvnvsTheory(){
    gStyle->SetOptStat(kFALSE);
    int icent = 0;
    int n = 3;
    int color[6] = {1,2,4,7,8,5};
    int style[12] = {20,21,24,25,26,27,29,30,31,32,33,34};
    TGraphErrors *gr[nsub][3][2];
    TGraphErrors *grraw[nsub][3][2];
    TString CNTEP, dire;
    for(int isub=0;isub<nsub;isub++){
    for(int idire=0;idire<3;idire++){
    for(int iCNTEP=0;iCNTEP<1;iCNTEP++){
        if(iCNTEP==0) CNTEP = "NoUseCNTEP";
        if(iCNTEP==1) CNTEP = "UseCNTEP";
        if(idire==0) dire = "";
        if(idire==1) dire = "_east";
        if(idire==2) dire = "_west";
        TString str = choosesub(isub);
        if(str=="ABORT") continue;
        gr[isub][idire][iCNTEP] = new TGraphErrors(Form("Result/%s/v%d_00_%d%s_%s.dat",CNTEP.Data(),n,icent,dire.Data(),str.Data()),"%lg %lg %lg");
        grraw[isub][idire][iCNTEP] = new TGraphErrors(Form("Result/%s/v%draw_00_%d%s_%s.dat",CNTEP.Data(),n,icent,dire.Data(),str.Data()),"%lg %lg %lg");
        SetStyle(*gr[isub][idire][iCNTEP], 1.2, color[idire+3*iCNTEP],style[isub]);
        SetStyle(*grraw[isub][idire][iCNTEP], 1.2, color[idire+3*iCNTEP],style[isub]);
    }
    }
    }
    TFile *fth = TFile::Open("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/taxiwork/TheoryCurves.root");
    TGraphErrors *fAMPT = (TGraphErrors*)fth->Get(Form("v%d_62GeV_AMPT",n));
    TGraphErrors *fSONIC = (TGraphErrors*)fth->Get(Form("v%d_62GeV_SONIC",n));
    TGraphErrors *fSSONIC = (TGraphErrors*)fth->Get(Form("v%d_62GeV_superSONIC",n));
    TFile *fthIP = TFile::Open("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/xuq/phenix/flow/Run16dAu/taxiwork/IPGlasma_dAu.root");
    if(n==2) TGraph *fIP = (TGraph*)fthIP->Get(Form("gv%d",n));
    else TGraph *fIP=NULL;

  
TH1F* h = new TH1F("h","",50,0,5);
h->GetXaxis()->SetRangeUser(0,3.2);
/*
TCanvas *c1 = new TCanvas("c1","c1",800,450);
iCNTEP = 0;
idire = 0;
c1->Divide(2);
c1->cd(1);
h->SetMinimum(0);
h->SetMaximum(0.3);
h->GetXaxis()->SetRangeUser(0,3.2);
//SetTitle(h,"","p_{T}","v_{2}^{raw}");
SetTitle(h,"","p_{T}","v_{2}");
h->DrawCopy();
TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
leg->SetFillColor(0);
leg->SetBorderSize(0);
for(int ilay = 0;ilay<4;ilay++){
leg->AddEntry(gr[ilay+8][idire][iCNTEP],Form("FVTX layer %d",ilay));
gr[ilay+8][idire][iCNTEP]->Draw("Psame");
}
gr[6][0][1]->Draw("Psame");
leg->AddEntry(gr[6][idire][iCNTEP],Form("FVTX -3.0<#eta<-1.0"));
leg->Draw("Psame");

c1->cd(2);
h->SetMinimum(0.8);
h->SetMaximum(1.2);
SetTitle(h,"","p_{T}","v_{2} ratio");
h->DrawCopy();
for(int ilay = 0;ilay<4;ilay++){
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[ilay+8][idire][iCNTEP],gr[6][idire][iCNTEP]);
SetStyle(*grr,1.2,color[idire+3*iCNTEP],style[ilay+8]);
grr->Draw("Psame");
}
*/
TCanvas *c2 = new TCanvas("c2","c2",450,450);
iCNTEP = 0;
idire = 0;
//c2->Divide(2);
c2->cd(1);
if(n==1){
h->SetMinimum(-0.05);
h->SetMaximum(0.);
}
else if(n==2){
h->SetMinimum(0);
h->SetMaximum(0.3);
}
else if(n==3){
h->SetMinimum(-0.05);
h->SetMaximum(0.2);
}
//SetTitle(h,"","p_{T}","v_{2}^{raw}");
SetTitle(*h,"p_{T} (GeV/c)",Form("v_{%d}",n),"");
c2->SetLeftMargin(0.12);
h->GetXaxis()->CenterTitle();
h->GetYaxis()->CenterTitle();
h->GetYaxis()->SetTitleSize(0.06);
h->DrawCopy();
TLatex t;
t.SetNDC();
t.DrawLatex(0.6,0.75,"d+Au 62GeV");
TLegend *leg = new TLegend(0.2,0.65,0.5,0.88);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextSize(0.04);
SetStyle(*gr[4][idire][iCNTEP], 1.2, 2,24);
SetStyle(*gr[6][idire][iCNTEP], 1.2, 4,24);
fAMPT->Draw("Lsame");
fSONIC->SetLineWidth(fAMPT->GetLineWidth());
fSONIC->Draw("Lsame");
fSSONIC->Draw("Lsame");
if(n==2){
fIP->SetLineColor(1);
fIP->Draw("Lsame");
}
gr[4][idire][iCNTEP]->Draw("Psame");
gr[6][idire][iCNTEP]->Draw("Psame");
leg->AddEntry(gr[4][idire][iCNTEP],Form("BBCs"),"P");
if(n==2)
leg->AddEntry(gr[6][idire][iCNTEP],Form("FVTXs -3.0<#eta<-1.0"),"P");
if(n==3)
leg->AddEntry(gr[6][idire][iCNTEP],Form("FVTXs 1+2+3L -3.0<#eta<-1.0"),"P");
leg->AddEntry(fAMPT,Form("AMPT"),"L");
leg->AddEntry(fSONIC,Form("SONIC"),"L");
leg->AddEntry(fSSONIC,Form("superSONIC"),"L");
if(n==2)
leg->AddEntry(fIP,Form("IPGlasma"),"L");
leg->Draw("Psame");
/*
c2->cd(2);
if(n==2){
h->SetMinimum(0.8);
h->SetMaximum(1.2);
}
if(n==3){
h->SetMinimum(0);
h->SetMaximum(2);
}
SetTitle(h,"","p_{T}","v_{2} ratio BBCs/FVTXs");
h->DrawCopy();
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[4][idire][iCNTEP],gr[6][idire][iCNTEP]);
//SetStyle(*grr,1.2,color[idire+3*iCNTEP],style[2]);
grr->Draw("Psame");
*/
c2->Print(Form("v%dTheory.png",n));
}
