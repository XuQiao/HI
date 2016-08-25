#include "par.h"
void plotgvsr(){

int xbin=0;	//xbin<1
int xpt=0;
int xtheta=0;	//xtheta<5

TFile *f = TFile::Open("mergedV_Sum.root");
TFile *fProd = TFile::Open("mergedV_Prod.root");
TVectorD *vecDr = f->Get(Form("D_%d/r",xbin));
TVectorD *vecDg2 = f->Get(Form("D_%d/D_%d/D_%d/G2",xbin,xpt,xtheta));
TVectorD *vecDsigma2 = f->Get(Form("D_%d/D_%d/sigma2",xbin,xpt));
TVectorD *vecDV = f->Get(Form("D_%d/D_%d/V",xbin,xpt));
TVectorD *vecDavgmult = f->Get(Form("D_%d/avgmult",xbin));

TVectorD *vecDr_ = fProd->Get(Form("D_%d/r",xbin));
TVectorD *vecDg2_ = fProd->Get(Form("D_%d/D_%d/D_%d/G2",xbin,xpt,xtheta));
TVectorD *vecDsigma2_ = fProd->Get(Form("D_%d/D_%d/sigma2",xbin,xpt));
TVectorD *vecDV_ = fProd->Get(Form("D_%d/D_%d/V",xbin,xpt));
TVectorD *vecDavgmult_ = fProd->Get(Form("D_%d/avgmult",xbin));

double *r = vecDr->GetMatrixArray();
double *g2 = vecDg2->GetMatrixArray();
double *sigma2 = vecDsigma2->GetMatrixArray();
double *V = vecDV->GetMatrixArray();
double *avgmult = vecDavgmult->GetMatrixArray();

double *r_ = vecDr_->GetMatrixArray();
double *g2_ = vecDg2_->GetMatrixArray();
double *sigma2_ = vecDsigma2_->GetMatrixArray();
double *V_ = vecDV_->GetMatrixArray();
double *avgmult_ = vecDavgmult_->GetMatrixArray();

TCanvas *c1 = new TCanvas;
c1->SetLogy();
TGraph *gr=new TGraph(nstepr,r,g2);
TGraph *grProd=new TGraph(nstepr,r_,g2_);
gr->GetXaxis()->SetTitle("r");
gr->GetYaxis()->SetTitle("|G^{#theta}(ir)|");
gr->GetXaxis()->SetRangeUser(0.1,0.4);
gr->GetYaxis()->SetRangeUser(1e-8,2);
gr->SetTitle("");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
grProd->SetMarkerColor(4);
gr->SetMarkerStyle(20);
grProd->SetMarkerStyle(29);
gr->Draw("AP");
grProd->Draw("Psame");
TLatex *t= new TLatex();
t->SetNDC();
if(xtheta==0)
t->DrawLatex(0.2,0.8,Form("%d < mult <%d, theta = %d", trkbin[xbin+1],trkbin[xbin],xtheta));
else
t->DrawLatex(0.2,0.8,Form("%d < mult <%d, theta = #frac{%d}{%d}#pi", trkbin[xbin+1],trkbin[xbin],xtheta,ntheta*nn));

double inV2 = 55098615530239328e-02;
TF1 *gcl = new TF1("gcl","TMath::Power(exp(-[0]*x*x/4)*TMath::BesselJ0([1]*x),2)",0,1);
TF1 *gclProd = new TF1("gclProd","TMath::Power(exp(-[0]*x*x/4)*TMath::BesselJ0([1]*x),2)",0,1);
//gcl->SetParameters(sigma2[xtheta],V[xtheta]*sigma2[xtheta]);
//gcl->SetParameters(300,0.065*300);
cout<<sigma2[xtheta]<<"\t"<<V[xtheta]<<"\t"<<avgmult[xbin]<<endl;
gcl->SetParameters(sigma2[xtheta]+V[xtheta]*V[xtheta]*avgmult[xbin]*avgmult[xbin]-inV2*inV2*avgmult[xbin]*avgmult[xbin],inV2*avgmult[xbin]);
gclProd->SetParameters(sigma2_[xtheta]+V_[xtheta]*V_[xtheta]*avgmult_[xbin]*avgmult_[xbin]-inV2*inV2*avgmult_[xbin]*avgmult_[xbin],inV2*avgmult_[xbin]);
gcl->SetLineColor(2);
gclProd->SetLineColor(3);
TLegend *leg = new TLegend(0.6,0.45,0.8,0.65);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.04);
leg->AddEntry(gr,"LYZ Sum Analysis result","P");
leg->AddEntry(grProd,"LYZ Prod Analysis result","P");
leg->AddEntry(gcl,"theoretical prediction","L");
leg->Draw("same");
gcl->Draw("same");
c1->Print("gvsr.png");

f->Close();

}

