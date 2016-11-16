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
        else if(isub==7)
         str = "FVTX1p2p4LS";
        else if(isub==8)
         str = "FVTX1N";
        else if(isub==9)
         str = "FVTXtrkS";
        else
         str = "ABORT";
    return str;
}

void plotvn(){
    gStyle->SetOptStat(kFALSE);
    int icent = 0;
    int n = 2;
    int color[6] = {1,2,4,7,8,5};
    int style[12] = {20,21,24,25,26,27,28,30,31,32,33,34};
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
        gr[isub][idire][iCNTEP] = new TGraphErrors(Form("Result/%s/v%d_%d%s_%s.dat",CNTEP.Data(),n,icent,dire.Data(),str.Data()),"%lg %lg %lg");
        grraw[isub][idire][iCNTEP] = new TGraphErrors(Form("Result/%s/v%draw_%d%s_%s.dat",CNTEP.Data(),n,icent,dire.Data(),str.Data()),"%lg %lg %lg");
        SetStyle(*gr[isub][idire][iCNTEP], 1.2, color[idire+3*iCNTEP],style[isub]);
        SetStyle(*grraw[isub][idire][iCNTEP], 1.2, color[idire+3*iCNTEP],style[isub]);
    }
    }
    }
    TGraphErrors* grold_fvtx1s = new TGraphErrors("Result/v2_pt_dAu_00_05_sys.dat","%lg %lg %lg %lg");
    SetStyle(*grold_fvtx1s,1.2,1,20);
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
TCanvas *c2 = new TCanvas("c2","c2",800,450);
iCNTEP = 0;
idire = 0;
c2->Divide(2);
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
TLegend *leg = new TLegend(0.2,0.7,0.5,0.88);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextSize(0.04);
SetStyle(*gr[4][idire][iCNTEP], 1.2, 2,style[2]);
SetStyle(*gr[5][idire][iCNTEP], 1.2, 1,style[3]);
SetStyle(*gr[6][idire][iCNTEP], 1.2, 4,style[6]);
SetStyle(*gr[7][idire][iCNTEP], 1.2, 7,style[4]);
SetStyle(*gr[8][idire][iCNTEP], 1.2, 1,style[5]);
SetStyle(*gr[9][idire][iCNTEP], 1.2, 4,style[7]);
if(n==2)
grold_fvtx1s->Draw("Psame");
gr[4][idire][iCNTEP]->Draw("Psame");
gr[5][idire][iCNTEP]->Draw("Psame");
gr[6][idire][iCNTEP]->Draw("Psame");
gr[7][idire][iCNTEP]->Draw("Psame");
gr[8][idire][iCNTEP]->Draw("Psame");
gr[9][idire][iCNTEP]->Draw("Psame");
if(n==2)
leg->AddEntry(grold_fvtx1s,"Run8 200GeV 0-5\%");
leg->AddEntry(gr[4][idire][iCNTEP],Form("Run16 200GeV BBCs"),"P");
leg->AddEntry(gr[5][idire][iCNTEP],Form("Run16 200GeV FVTXs cluster "),"P");
leg->AddEntry(gr[6][idire][iCNTEP],Form("Run16 200GeV FVTXs cluster 1+2+3L"),"P");
leg->AddEntry(gr[7][idire][iCNTEP],Form("Run16 200GeV FVTXs cluster 1+2+4L"),"P");
leg->AddEntry(gr[8][idire][iCNTEP],Form("Run16 200GeV FVTXs cluster North"),"P");
leg->AddEntry(gr[9][idire][iCNTEP],Form("Run16 200GeV FVTXs track #eta>-3.5"),"P");
leg->Draw("Psame");

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
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[4][idire][iCNTEP],gr[5][idire][iCNTEP]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[6][idire][iCNTEP],gr[5][idire][iCNTEP]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[7][idire][iCNTEP],gr[5][idire][iCNTEP]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[8][idire][iCNTEP],gr[5][idire][iCNTEP]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[9][idire][iCNTEP],gr[5][idire][iCNTEP]);
grr->Draw("Psame");
c2->Print(Form("v%dBBCSFVTX.png",n));
/*
TCanvas *c3 = new TCanvas("c3","c3",800,450);
idire = 0;
isub=5;
c3->Divide(2);
c3->cd(1);
h->SetMinimum(0);
h->SetMaximum(0.3);
//SetTitle(h,"","p_{T}","v_{2}^{raw}");
SetTitle(h,"","p_{T}","v_{2}");
h->DrawCopy();
TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
leg->SetFillColor(0);
leg->SetBorderSize(0);
gr[isub][idire][1]->Draw("Psame");
gr[isub][idire][0]->Draw("Psame");
leg->AddEntry(gr[isub][idire][1],Form("FVTXs -3.0<#eta<-1.0"));
leg->AddEntry(gr[isub][idire][1],Form("Using Psi EP"),"P");
leg->AddEntry(gr[isub][idire][0],Form("Using phi EP"),"P");
leg->Draw("Psame");

c3->cd(2);
h->SetMinimum(0.8);
h->SetMaximum(1.2);
SetTitle(h,"","p_{T}","v_{2} ratio FVTXs using phi EP/using Psi EP");
h->DrawCopy();
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[isub][idire][0],gr[isub][idire][1]);
SetStyle(*grr,1.2,color[idire+3*0],style[isub]);
grr->Draw("Psame");
*/
TCanvas *c4 = new TCanvas("c4","c4",800,450);
isub=5;
iCNTEP = 0;
c4->Divide(2);
c4->cd(1);
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
h->SetMaximum(0.1);
}
SetTitle(h,"","p_{T}","v_{2}");
h->DrawCopy();
TLegend *leg = new TLegend(0.2,0.7,0.5,0.88);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextSize(0.04);
gr[isub][0][iCNTEP]->Draw("Psame");
gr[isub][1][iCNTEP]->Draw("Psame");
gr[isub][2][iCNTEP]->Draw("Psame");
leg->AddEntry(gr[isub][0][iCNTEP],Form("FVTXs -3.0<#eta<-1.0"),"P");
leg->AddEntry(gr[isub][0][iCNTEP],Form("inclusive"),"P");
leg->AddEntry(gr[isub][1][iCNTEP],Form("East"),"P");
leg->AddEntry(gr[isub][2][iCNTEP],Form("West"),"P");
leg->Draw("Psame");

c4->cd(2);
h->SetMinimum(0.8);
h->SetMaximum(1.2);
SetTitle(h,"","p_{T}","v_{2} ratio FVTXs");
h->DrawCopy();
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[isub][1][iCNTEP],gr[isub][0][iCNTEP]);
//SetStyle(*grr,1.2,color[1+3*iCNTEP],style[isub]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[isub][2][iCNTEP],gr[isub][0][iCNTEP]);
//SetStyle(*grr,1.2,color[2+3*iCNTEP],style[isub]);
grr->Draw("Psame");
c4->Print(Form("v%dEWFVTX.png",n));

TCanvas *c5 = new TCanvas("c5","c5",800,450);
isub=4;
iCNTEP = 0;
c5->Divide(2);
c5->cd(1);
//SetTitle(h,"","p_{T}","v_{2}^{raw}");
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
h->SetMaximum(0.1);
}
SetTitle(h,"","p_{T}","v_{2}");
h->DrawCopy();
TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
leg->SetTextSize(0.05);
leg->SetFillColor(0);
leg->SetBorderSize(0);
gr[isub][0][iCNTEP]->Draw("Psame");
gr[isub][1][iCNTEP]->Draw("Psame");
gr[isub][2][iCNTEP]->Draw("Psame");
leg->AddEntry(gr[isub][0][iCNTEP],Form("Using BBC event plane"),"");
leg->AddEntry(gr[isub][0][iCNTEP],Form("inclusive"),"P");
leg->AddEntry(gr[isub][1][iCNTEP],Form("East"),"P");
leg->AddEntry(gr[isub][2][iCNTEP],Form("West"),"P");
leg->Draw("Psame");

c5->cd(2);
h->SetMinimum(0);
h->SetMaximum(2.);
SetTitle(h,"","p_{T}","v_{2} ratio");
h->DrawCopy();
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[isub][1][iCNTEP],gr[isub][0][iCNTEP]);
//SetStyle(*grr,1.2,color[1+3*iCNTEP],style[isub]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[isub][2][iCNTEP],gr[isub][0][iCNTEP]);
//SetStyle(*grr,1.2,color[2+3*iCNTEP],style[isub]);
grr->Draw("Psame");
c5->Print(Form("v%dEWFVTX4L.png",n));

TCanvas *c6 = new TCanvas("c6","c6",800,450);
iCNTEP = 0;
c6->Divide(2);
c6->cd(1);
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
h->SetMaximum(0.1);
}
//SetTitle(h,"","p_{T}","v_{2}^{raw}");
SetTitle(h,"","p_{T}","v_{2}");
h->DrawCopy();
TLegend *leg = new TLegend(0.5,0.75,0.8,0.85);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextSize(0.04);
grold_fvtx1s->Draw("Psame");
gr[0][0][iCNTEP]->Draw("Psame");
gr[1][0][iCNTEP]->Draw("Psame");
gr[2][0][iCNTEP]->Draw("Psame");
gr[3][0][iCNTEP]->Draw("Psame");
SetStyle(*gr[0][0][iCNTEP],1.2,color[4],style[4]);
SetStyle(*gr[1][0][iCNTEP],1.2,color[1],style[1]);
SetStyle(*gr[2][0][iCNTEP],1.2,color[2],style[2]);
SetStyle(*gr[3][0][iCNTEP],1.2,color[3],style[3]);
leg->AddEntry(grold_fvtx1s,"Run8 200GeV 0-5\%");
leg->AddEntry(gr[0][0][iCNTEP],Form("FVTX 1LS"),"P");
leg->AddEntry(gr[1][0][iCNTEP],Form("FVTX 2LS"),"P");
leg->AddEntry(gr[2][0][iCNTEP],Form("FVTX 3LS"),"P");
leg->AddEntry(gr[3][0][iCNTEP],Form("FVTX 4LS"),"P");
leg->Draw("same");

c6->cd(2);
h->SetMinimum(0);
h->SetMaximum(2);
SetTitle(h,"","p_{T}","v_{2} ratio over bbc");
h->DrawCopy();
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[0][0][iCNTEP],gr[4][0][iCNTEP]);
//SetStyle(*grr,1.2,color[idire+3*iCNTEP],style[0]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[1][0][iCNTEP],gr[4][0][iCNTEP]);
//SetStyle(*grr,1.2,color[idire+3*iCNTEP],style[1]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[2][0][iCNTEP],gr[4][0][iCNTEP]);
//SetStyle(*grr,1.2,color[idire+3*iCNTEP],style[2]);
grr->Draw("Psame");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[3][0][iCNTEP],gr[4][0][iCNTEP]);
//SetStyle(*grr,1.2,color[idire+3*iCNTEP],style[3]);
grr->Draw("Psame");
c6->Print(Form("v%dLayers.png",n));
}
