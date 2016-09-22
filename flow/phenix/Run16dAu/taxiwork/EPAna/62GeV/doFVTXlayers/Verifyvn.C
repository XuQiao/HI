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
        else
         str = "ABORT";
    return str;
}

void Verifyvn(){
    gStyle->SetOptStat(kFALSE);
    int icent = 1;
    int n = 2;
    int color[12] = {1,2,4,5,7,8,1,2,4,5,7,8};
    int style[12] = {20,21,24,25,26,27,29,30,31,32,33,34};
    TGraphErrors *gr[nsub][3][2];
    TGraphErrors *gr1[nsub][3][2];
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
        //gr1[isub][idire][iCNTEP] = new TGraphErrors(Form("../Originaltest/Result/%s/v%d_%d%s_%s.dat",CNTEP.Data(),n,icent,dire.Data(),str.Data()),"%lg %lg %lg");
        SetStyle(*gr[isub][idire][iCNTEP], 1.2, color[idire+3*iCNTEP+isub+3],style[isub]);
      //  SetStyle(*gr1[isub][idire][iCNTEP], 1.2, color[6-(idire+3*iCNTEP)],style[12-isub]);
    }
    }
    }

    TGraphErrors* grold_fvtx1s = new TGraphErrors("Result/v2_pt_dAu_00_05_sys.dat","%lg %lg %lg %lg");
    SetStyle(*grold_fvtx1s,1.2,1,20);
    TGraphErrors* gr_fvtx1s = new TGraphErrors("../../200GeV/doFVTXlayers/Result/NoUseCNTEP/v2_00_0_FVTX1S.dat","%lg %lg %lg");
    SetStyle(*gr_fvtx1s,1.2,1,20);

TCanvas *c1 = new TCanvas("c1","c1",800,450);
iCNTEP = 0;
idire = 0;
c1->Divide(2);
c1->cd(1);
TH1F* h = new TH1F("h","",50,0,5);
h->SetMinimum(0);
h->SetMaximum(0.3);
h->GetXaxis()->SetRangeUser(0,3.2);
//SetTitle(h,"","p_{T}","v_{2}^{raw}");
SetTitle(h,"","p_{T} (GeV/c)","v_{2}");
h->DrawCopy();
TLegend *leg = new TLegend(0.12,0.75,0.5,0.85);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextSize(0.04);
gr_fvtx1s->Draw("Psame");
gr[5][idire][iCNTEP]->Draw("Psame");
gr[4][idire][iCNTEP]->Draw("Psame");
leg->AddEntry(gr_fvtx1s,"Run16 200 GeV 0-5\%","P");
leg->AddEntry(gr[5][idire][iCNTEP],Form("Run16 62GeV 0-5%% FVTX -3.0<#eta<-1.0 as EP"),"P");
leg->AddEntry(gr[4][idire][iCNTEP],Form("Run16 62GeV 0-5\% BBCS as EP"),"P");
//leg->AddEntry(gr1[6][idire][iCNTEP],Form("FVTX -3.0<#eta<-1.0 Ori"));
leg->Draw("Psame");

c1->cd(2);
h->SetMinimum(0);
h->SetMaximum(2);
SetTitle(h,"","p_{T}","v_{2} ratio Run16 62GeV/Run16 200GeV");
h->DrawCopy();
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[5][idire][iCNTEP],gr_fvtx1s);
SetStyle(*grr,1.2,color[idire+3*iCNTEP+5+3],style[5]);
grr->Draw("Psame");

TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[4][idire][iCNTEP],gr_fvtx1s);
SetStyle(*grr,1.2,color[idire+3*iCNTEP+4+3],style[4]);
grr->Draw("Psame");
c1->Print(Form("Verifyv%d_cent%d.png",n,icent));
}

