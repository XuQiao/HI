#include "par.h"
void plotgvsr(){

int xbin=0;	//xbin<1
int xptbin=1;	//xptbin<13
int xtheta=4;	//xtheta<5

TFile *f = TFile::Open("mergedV_Sum.root");
TVectorD *vecDr = f->Get(Form("D_%d/r",xbin));
TVectorD *vecDg2 = f->Get(Form("D_%d/D_%d/D_%d/G2",xbin,xptbin,xtheta));
TVectorD *vecDsigma2 = f->Get(Form("D_%d/D_%d/sigma2",xbin,xptbin));
TVectorD *vecDV = f->Get(Form("D_%d/D_%d/V",xbin,xptbin));
TVectorD *vecDavgmult = f->Get(Form("D_%d/avgmult",xbin));

double *r = vecDr->GetMatrixArray();
double *g2 = vecDg2->GetMatrixArray();
double *sigma2 = vecDsigma2->GetMatrixArray();
double *V = vecDV->GetMatrixArray();
double *avgmult = vecDavgmult->GetMatrixArray();
cout<<sigma2[0]<<endl;
cout<<V[0]<<endl;
TCanvas *c1 = new TCanvas;
c1->SetLogy();
TGraph *gr=new TGraph(nstepr,r,g2);
gr->GetXaxis()->SetTitle("r");
gr->GetYaxis()->SetTitle("|G^{#theta}(ir)|");
if(isSimple)
gr->GetXaxis()->SetRangeUser(0,1);
else
gr->GetXaxis()->SetRangeUser(0.05,0.2);
gr->GetYaxis()->SetRangeUser(1e-8,1);
gr->SetTitle("");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);
gr->Draw("AP");
TLatex *t= new TLatex();
t->SetNDC();
t->SetTextSize(0.03);
if(xtheta==0)
t->DrawLatex(0.2,0.8,Form("%d<mult<%d, %.1f<p_{T}<%.1f, theta = %d", trkbin[xbin+1],trkbin[xbin],ptbin[xptbin],ptbin[xptbin+1],xtheta));
else
t->DrawLatex(0.2,0.8,Form("%d<mult<%d, %.1f<p_{T}<%.1f, theta = #frac{%d}{%d}#pi", trkbin[xbin+1],trkbin[xbin],ptbin[xptbin],ptbin[xptbin+1],xtheta,ntheta*nn));

double inV2 = 0.065;
TF1 *gcl = new TF1("gcl","TMath::Power(exp(-[0]*x*x/4)*TMath::BesselJ0([1]*x),2)",0,1);
//gcl->SetParameters(sigma2[xtheta],V[xtheta]*sigma2[xtheta]);
//gcl->SetParameters(300,0.065*300);
gcl->SetParameters(sigma2[xtheta]+V[xtheta]*V[xtheta]*avgmult[xbin]*avgmult[xbin]-inV2*inV2*avgmult[xbin]*avgmult[xbin],inV2*avgmult[xbin]);
gcl->SetLineColor(2);
TLegend *leg = new TLegend(0.6,0.7,0.8,0.9);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.03);
leg->AddEntry(gr,"Analysis result","P");
leg->AddEntry(gcl,"theoretical prediction","L");
leg->Draw("same");
gcl->Draw("same");
c1->Print("gvsr.png");

f->Close();

}

