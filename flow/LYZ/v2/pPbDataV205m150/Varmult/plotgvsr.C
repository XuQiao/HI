#include "par.h"

int xbin=0;	//xbin<1
double ptmin=0.3, ptmax=6.0;
const int nnu = 5;
//const int marker[nnu] = {20,24,27};
const int marker[nnu] = {20,24,27,30,29};
//const int color[nnu] = {1,2,3};
const int color[nnu] = {1,2,3,4,6};
//const int plotnu[nnu]={8,9,10};
const int plotnu[nnu]={0,1,2,3,4};
void plotgvsr(){
    gStyle->SetOptStat(kFALSE);
TCanvas *c1 = new TCanvas();
c1->SetLogy();
TH1D* hFrame = new TH1D("","",1000,0,1);
hFrame->GetYaxis()->SetRangeUser(5e-7,1);
hFrame->GetXaxis()->SetRangeUser(0.20,0.40);
hFrame->GetXaxis()->SetTitle("r");
hFrame->GetYaxis()->SetTitle("|G^{#theta}(ir)|^{2}");
hFrame->SetTitle("");
hFrame->Draw();
TGraph *gr[nnu];
double r0_theta[nnu];
double G2_theta[nnu];
TLine *l[nnu];
for(int itheta=0;itheta<nnu;itheta++){
gr[itheta]=plotGF(0,plotnu[itheta],r0_theta+itheta,G2_theta+itheta,marker[itheta],color[itheta]);
gr[itheta]->Draw("Psame");
}
TLegend *tg = new TLegend(0.75,0.70-0.10*nnu,0.90,0.70);
tg->SetFillColor(0);
tg->SetBorderSize(0);
tg->SetTextSize(0.04);
for(int itheta=0;itheta<nnu;itheta++)
tg->AddEntry(gr[itheta],Form("#theta = #frac{%d}{%d}*#frac{#pi}{%d}",plotnu[itheta],ntheta,nn),"lp");
tg->Draw("same");
TLatex *t= new TLatex();
t->SetNDC();
t->SetTextSize(0.03);
t->SetTextFont(42);
t->DrawLatex(0.5,0.2,Form("track, %d < mult <%d, %.1f < p_{T} < %.1f", trkbin[xbin+1],trkbin[xbin],ptmin,ptmax));

for(int itheta=0;itheta<nnu;itheta++){
l[itheta] = new TLine(r0_theta[itheta],0,r0_theta[itheta],G2_theta[itheta]);
l[itheta]->SetLineStyle(2);
l[itheta]->SetLineColor(color[itheta]);
l[itheta]->Draw("same");
}
c1->Print("gvsr_thetas.png");


}

TGraph *plotGF(int isSum, int xtheta, double *r0_theta, double *G2_theta, int marker, int color){

if(isSum)
TFile *f = TFile::Open("mergedV_Sum.root");
else
TFile *f = TFile::Open("mergedV_Prod.root");

TVectorD *vecDr = f->Get(Form("D_%d/r",xbin));
TVectorD *vecDg2 = f->Get(Form("D_%d/D_%d/G2",xbin,xtheta));
TVectorD *vecDr0 = f->Get(Form("D_%d/r0",xbin));
TVectorD *vecDr01 = f->Get(Form("D_%d/r01",xbin));

double *r = vecDr->GetMatrixArray();
double *g2 = vecDg2->GetMatrixArray();
double *r0 = vecDr0->GetMatrixArray();
double *r01 = vecDr01->GetMatrixArray();
(*r0_theta) = r0[xtheta];
for(int ir=0;ir<nstepr;ir++)
    if(r[ir] == r01[xtheta]) break;
(*G2_theta) = g2[ir];
TGraph *gr = new TGraph(nstepr,r,g2);
gr->SetMarkerSize(0.5);
gr->SetMarkerColor(color);
gr->SetMarkerStyle(marker);
f->Close();
return gr;
}

