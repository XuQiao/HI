#include "par.h"

int xbin=0;	//xbin<1
int xtheta=0;
double ptmin=0.3, ptmax=6.0;
const int nnu = 2;
const int marker[nnu] = {20,24};
const int color[nnu] = {1,4};
TF1* gcl = new TF1("gcl","TMath::Power(exp(-[0]*x*x/4)*TMath::BesselJ0([1]*x),2)",0,1);
TLine* linv =  new TLine(0,0,0,0);
void plotgvsr_MC(){
gStyle->SetOptStat(kFALSE);
TCanvas *c1 = new TCanvas();
c1->SetLogy();
TH1D* hFrame = new TH1D("","",1000,0,1);
hFrame->GetYaxis()->SetRangeUser(5e-9,1);
hFrame->GetXaxis()->SetRangeUser(0.15,0.50);
hFrame->GetXaxis()->SetTitle("r");
hFrame->GetYaxis()->SetTitle("|G^{#theta}(ir)|^{2}");
hFrame->SetTitle("");
hFrame->Draw();
TGraph *gr[nnu];
double r0[nnu];
double G2[nnu];
TLine *l[nnu];
for(int i=0;i<nnu;i++){
gr[i]=plotGF(1,xtheta,r0+i,G2+i,marker[i],color[i]);
gr[i]->Draw("Psame");
}
gcl->Draw("same");
TLegend *tg = new TLegend(0.75,0.70-0.10*nnu,0.90,0.70);
tg->SetFillColor(0);
tg->SetBorderSize(0);
tg->SetTextSize(0.04);
for(int i=0;i<nnu;i++)
if(i==0)
tg->AddEntry(gr[i],Form("Prod"),"lp");
else
tg->AddEntry(gr[i],Form("Sum"),"lp");
tg->AddEntry(gcl,Form("Theory"),"lp");
tg->Draw("same");
TLatex *t= new TLatex();
t->SetNDC();
t->SetTextSize(0.05);
t->SetTextFont(42);
if(xtheta==0)
    t->DrawLatex(0.6,0.2,Form("mult = %d, theta = %d", (trkbin[xbin+1]+trkbin[xbin])/2,xtheta));
    else
    t->DrawLatex(0.6,0.2,Form("mult = %d, theta = #frac{%d}{%d}#pi", (trkbin[xbin+1]+trkbin[xbin])/2,xtheta,ntheta*nn));

for(int i=0;i<nnu;i++){
l[i] = new TLine(r0[i],0,r0[i],G2[i]);
l[i]->SetLineStyle(2);
l[i]->SetLineColor(color[i]);
l[i]->Draw("same");
}
linv->SetLineStyle(2);
linv->SetLineColor(2);
linv->SetY1(hFrame->GetMinimum());
linv->Draw("same");
c1->Print("gvsr_thetas_MC.png");


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
TVectorD *vecDsigma2 = f->Get(Form("D_%d/sigma2",xbin));
TVectorD *vecDavgmult = f->Get(Form("avgmult"));

double *r = vecDr->GetMatrixArray();
double *g2 = vecDg2->GetMatrixArray();
double *r0 = vecDr0->GetMatrixArray();
double *r01 = vecDr01->GetMatrixArray();
double *sigma2 = vecDsigma2->GetMatrixArray();
double *avgmult = vecDavgmult->GetMatrixArray();
(*r0_theta) = r0[xtheta];
for(int ir=0;ir<nstepr;ir++)
    if(r[ir] == r01[xtheta]) break;
(*G2_theta) = g2[ir];
TGraph *gr = new TGraph(nstepr,r,g2);
if(isSum){
gcl->SetNpx(10000);
gcl->SetParameters(sigma2[xtheta]-inV2*inV2*avgmult[xbin]*avgmult[xbin],inV2*avgmult[xbin]);
linv->SetX1(j01/avgmult[xbin]/inV2);
linv->SetX2(j01/avgmult[xbin]/inV2);
linv->SetY2(5e-9);
}
gr->SetMarkerSize(0.5);
gr->SetMarkerColor(color);
gr->SetMarkerStyle(marker);
f->Close();
return gr;
}

