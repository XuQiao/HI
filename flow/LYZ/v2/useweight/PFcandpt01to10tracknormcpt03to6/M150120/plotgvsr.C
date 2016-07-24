#include "par.h"

void plotgvsr(){

int xbin=0;	//xbin<1
int xtheta=0;	//xtheta<5
int xpt=0;
double ptmin=0.1, ptmax=10.0;

TFile *f = TFile::Open("mergedV_Prod.root");
TVectorD *vecDr = f->Get(Form("D_%d/r",xbin));
TVectorD *vecDg2 = f->Get(Form("D_%d/D_%d/D_%d/G2",xbin,xpt,xtheta));


double *r = vecDr->GetMatrixArray();
double *g2 = vecDg2->GetMatrixArray();

TCanvas *c1 = new TCanvas;
c1->SetLogy();
TGraph *gr=new TGraph(nstepr,r,g2);
gr->GetYaxis()->SetRangeUser(5e-7,1);
//gr->GetXaxis()->SetRangeUser(0.15,0.4);
gr->GetXaxis()->SetTitle("r");
gr->GetYaxis()->SetTitle("|G^{#theta}(ir)|^{2}");
gr->SetTitle("LYZ Sum");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);
gr->Draw("AP");
TLatex *t= new TLatex();
t->SetNDC();
if(xtheta==0)
t->DrawLatex(0.3,0.8,Form("%d < mult <%d, #theta = %d, %.1f < p_{T} < %.1f", trkbin[xbin+1],trkbin[xbin],xtheta,ptmin,ptmax));
else
t->DrawLatex(0.3,0.8,Form("%d < mult <%d, #theta = #frac{%d}{%d}#pi, %.1f < p_{T} < %.1f", trkbin[xbin+1],trkbin[xbin],xtheta,ntheta*nn,ptmin,ptmax));

c1->Print("gvsr.png");

f->Close();

}

