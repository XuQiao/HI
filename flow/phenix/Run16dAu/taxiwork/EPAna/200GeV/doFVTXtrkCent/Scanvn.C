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

void Scanvn(int n=2){
    gStyle->SetOptStat(kFALSE);
    int icent = 0;
    int color[6] = {1,2,5,4,7,8};
    int style[12] = {20,21,24,25,26,27,29,30,31,32,33,34};
    TGraphErrors *gr[nangle1][nangle2][nsub][3][2];
    TGraphErrors *grraw[nangle1][nangle2][nsub][3][2];
    TGraphErrors *grratio[nangle1][nangle2][nsub][3][2];
    TString CNTEP, dire;
    TGraph *grAngles[nsub];
    int nangle = nangle1 * nangle2;
    float deltax0[10000];
    float diff[10000];
    float ratio[3];
    for(int isub=0;isub<nsub;isub++){
        TString str = choosesub(isub);
        if(str=="ABORT") continue;
TCanvas *c7 = new TCanvas(Form("c7_%d",isub),"",650,400);
c7->Divide(2);
TH1F* h = new TH1F(Form("h%d",isub),"",50,0,5);
TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
for(int iangle1=0;iangle1<nangle1;iangle1++){
for(int iangle2=0;iangle2<nangle2;iangle2++){
    for(int idire=0;idire<3;idire++){
    for(int iCNTEP=0;iCNTEP<1;iCNTEP++){
        if(iCNTEP==0) CNTEP = "NoUseCNTEP";
        if(iCNTEP==1) CNTEP = "UseCNTEP";
        if(idire==0) dire = "";
        if(idire==1) dire = "_east";
        if(idire==2) dire = "_west";
     //   deltax0[iangle1*nangle2+iangle2] = iangle1 * 0.1 - isub * 0.05;
        gr[iangle1][iangle2][isub][idire][iCNTEP] = new TGraphErrors(Form("Result/%s/v%d_%d%d_%d%s_%s.dat",CNTEP.Data(),n,iangle1,iangle2,icent,dire.Data(),str.Data()),"%lg %lg %lg");
        grraw[iangle1][iangle2][isub][idire][iCNTEP] = new TGraphErrors(Form("Result/%s/v%draw_%d%d_%d%s_%s.dat",CNTEP.Data(),n,iangle1,iangle2,icent,dire.Data(),str.Data()),"%lg %lg %lg");
        SetStyle(*gr[iangle1][iangle2][isub][idire][iCNTEP], 1.2, color[idire+3*iCNTEP],style[isub]);
        SetStyle(*grraw[iangle1][iangle2][isub][idire][iCNTEP], 1.2, color[idire+3*iCNTEP],style[isub]);
        grratio[iangle1][iangle2][isub][idire][iCNTEP] = (TGraphErrors*)DivideTwoGraphs(gr[iangle1][iangle2][isub][idire][iCNTEP],gr[iangle1][iangle2][isub][0][iCNTEP]);
        TF1 *fun = new TF1("fun","pol0",0.3,3.2);
        grratio[iangle1][iangle2][isub][idire][iCNTEP] -> Fit("fun","QR0");
        if(iCNTEP == 0){
        ratio[idire] = fun->GetParameter(0);
        }
        }
    }
   if(ratio[1]!=0) diff[iangle1*nangle2+iangle2] = ratio[2]/ratio[1];

int selangle1 = iangle1;
int selangle2 = iangle2;

TCanvas *c5 = new TCanvas(Form("c%d",isub),"c5",650,400);
h->GetXaxis()->SetRangeUser(0,3.2);
iCNTEP = 0;
c5->Divide(2);
c5->cd(1);
if(n==1){
h->SetMinimum(-0.05);
h->SetMaximum(0.);
}
else if(n==2){
h->SetMinimum(0);
h->SetMaximum(0.2);
}
else if(n==3){
h->SetMinimum(0);
h->SetMaximum(0.1);
}
SetTitle(*h,"p_{T}",Form("v_{%d}",n),"");
//SetTitle(h,"","p_{T}","v_{2}");
h->DrawCopy();
TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
leg->SetTextSize(0.05);
leg->SetFillColor(0);
leg->SetBorderSize(0);
gr[selangle1][selangle2][isub][0][iCNTEP]->Draw("Psame");
gr[selangle1][selangle2][isub][1][iCNTEP]->Draw("Psame");
gr[selangle1][selangle2][isub][2][iCNTEP]->Draw("Psame");
if(isub==4)
leg->AddEntry(gr[selangle1][selangle2][isub][0][iCNTEP],Form("Using BBC event plane"),"");
else if(isub==5)
leg->AddEntry(gr[selangle1][selangle2][isub][0][iCNTEP],Form("Using FVTX event plane"),"");
else
leg->AddEntry(gr[selangle1][selangle2][isub][0][iCNTEP],Form("Using FVTX event plane Layer %d",isub+1),"");
leg->AddEntry(gr[selangle1][selangle2][isub][0][iCNTEP],Form("inclusive"),"P");
leg->AddEntry(gr[selangle1][selangle2][isub][1][iCNTEP],Form("East"),"P");
leg->AddEntry(gr[selangle1][selangle2][isub][2][iCNTEP],Form("West"),"P");
leg->Draw("Psame");

c5->cd(2);
h->SetMinimum(0);
h->SetMaximum(2.);
SetTitle(*h,"p_{T}",Form("v_{%d} ratio",n),"");
h->DrawCopy();
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[selangle1][selangle2][isub][1][iCNTEP],gr[selangle1][selangle2][isub][0][iCNTEP]);
SetStyle(*grr,1.2,color[1+3*iCNTEP],style[isub]);
grr->Draw("Psame");
TF1 *fun = new TF1("fun","pol0",0.3,3.2);
fun->SetLineColor(1);
grr->Fit("fun","R");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[selangle1][selangle2][isub][2][iCNTEP],gr[selangle1][selangle2][isub][0][iCNTEP]);
SetStyle(*grr,1.2,color[2+3*iCNTEP],style[isub]);
grr->Draw("Psame");
TF1 *fun = new TF1("fun","pol0",0.3,3.2);
fun->SetLineColor(2);
grr->Fit("fun","R");

TLatex t;
t.SetNDC();
t.SetTextSize(0.04);
t.DrawLatex(0.2,0.85,Form("ratio of W/E = %.2f", diff[selangle1*nangle2+selangle2]));
//t.DrawLatex(0.2,0.8,Form("deltax0 = %.2f",deltax0[selangle1*nangle2+selangle2]));

c5->Print(Form("angleScan%sv%d_%d%d.png",choosesub(isub).Data(),n,selangle1,selangle2));


h->GetXaxis()->SetRangeUser(0,3.2);
iCNTEP = 0;
if(selangle1 % 3 != 0) continue;
c7->cd(1);
if(n==1){
h->SetMinimum(-0.05);
h->SetMaximum(0.);
}
else if(n==2){
h->SetMinimum(0);
h->SetMaximum(0.2);
}
else if(n==3){
h->SetMinimum(0);
h->SetMaximum(0.1);
}
SetTitle(*h,"p_{T}",Form("v_{%d}",n),"");
//SetTitle(h,"","p_{T}","v_{2}");
if(selangle1==0)
h->DrawCopy();
leg->SetTextSize(0.05);
leg->SetFillColor(0);
leg->SetBorderSize(0);
SetStyle(*gr[iangle1][iangle2][isub][0][iCNTEP], 1.2, color[selangle1/3],style[isub]);
gr[selangle1][selangle2][isub][0][iCNTEP]->Draw("Psame");
if(selangle1==0)
if(isub==4)
leg->AddEntry(gr[selangle1][selangle2][isub][0][iCNTEP],Form("Using BBC event plane"),"");
else if(isub==5)
leg->AddEntry(gr[selangle1][selangle2][isub][0][iCNTEP],Form("Using FVTX event plane"),"");
else
leg->AddEntry(gr[selangle1][selangle2][isub][0][iCNTEP],Form("Using FVTX event plane Layer %d",isub+1),"");
//leg->AddEntry(gr[selangle1][selangle2][isub][0][iCNTEP],Form("inclusive deltax0 = %.2f",deltax0[selangle1*nangle2+selangle2]),"P");
leg->Draw("Psame");

c7->cd(2);
h->SetMinimum(0);
h->SetMaximum(2.);
if(selangle1==0)
h->DrawCopy();
SetTitle(*h,"p_{T}",Form("v_{%d} ratio",n),"");
TGraphErrors *grr = (TGraphErrors*)DivideTwoGraphs(gr[selangle1][selangle2][isub][0][iCNTEP],gr[0][selangle2][isub][0][iCNTEP]);
SetStyle(*grr,1.2,color[selangle1/3],style[isub]);
grr->Draw("Psame");
TF1 *fun = new TF1("fun","pol0",0.3,3.2);
fun->SetLineColor(color[selangle1/3]);
grr->Fit("fun","R");

}
}
    c7->Print(Form("angleScan%sv%dratio.png",choosesub(isub).Data(),n));
//    grAngles[isub] = new TGraph(nangle, deltax0, diff);
}
/*
TCanvas *c6 = new TCanvas("c6","c6",650,600);
c6->SetGridx();
c6->SetGridy();
c6->SetRightMargin(0.15);
c6->SetLeftMargin(0.15);
TH1F* h = new TH1F("h","",100,-5,5);
h->GetXaxis()->SetRangeUser(-0.2,1.2);
h->GetYaxis()->SetRangeUser(0.4,1.6);
h->GetXaxis()->SetTitle("deltax0(cm)");
h->GetYaxis()->SetTitle("difference between E and W (W/E)");
h->DrawCopy();
TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
leg->SetTextSize(0.05);
leg->SetFillColor(0);
leg->SetBorderSize(0);
    for(int isub=0;isub<nsub;isub++){
        TString str = choosesub(isub);
        if(str=="ABORT") continue;
grAngles[isub]->SetTitle("");
grAngles[isub]->GetXaxis()->SetTitleOffset(1.0);
grAngles[isub]->GetYaxis()->SetTitleOffset(1.5);
SetStyle(*grAngles[isub],1.2,color[1],style[isub]);
//grAngles -> Draw("surf1");
//c6->Print("angleScan2D.png");
grAngles[isub]->Draw("Psame");
leg->AddEntry(grAngles[isub],Form("FVTX Layer %d",isub+1),"P");
    }
leg->Draw("same");
c6->Print(Form("angleScanv%d.png",n));
*/
}
