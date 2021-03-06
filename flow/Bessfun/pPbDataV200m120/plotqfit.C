#include "par.h"
#include <iomanip>
#include "../../../jetRpA/RpA/Quality/root_setting.h"

int xbin=0;	//xbin<1
int xtheta=0;
void plotqfit(){
    gStyle->SetOptStat(1011);
    gStyle->SetOptFit(1111);
TFile *f = TFile::Open("mergedV_Sum.root");
TVectorD *vecDavgmult = f->Get(Form("avgmult"));
double *avgmult = vecDavgmult->GetMatrixArray();
TLatex *t= new TLatex();
t->SetNDC();
t->SetTextSize(0.04);
t->SetTextFont(42);
TH1D* hq = (TH1D*)f->Get(Form("D_%d/D_%d/hq",xbin,xtheta));
TH1D* hqx = (TH1D*)f->Get(Form("D_%d/hqx",xbin));
TH1D* hqy = (TH1D*)f->Get(Form("D_%d/hqy",xbin));
TH1D* hq2 = (TH1D*)f->Get(Form("D_%d/hq2",xbin));
TH1D* hq2nonf = (TH1D*)f->Get(Form("D_%d/hq2nonf",xbin));
         multiplyByBinCenter(hq);
         hq->Scale(1./hq->Integral(0,-1,"width"));
        //normalizeByBinWidth(hq);
        //multiplyByBinCenter(hqx);
        hqx->Scale(1./hqx->Integral(0,-1,"width"));
        //normalizeByBinWidth(hqx);
        multiplyByBinCenter(hqy);
        hqy->Scale(1./hqy->Integral(0,-1,"width"));
        //normalizeByBinWidth(hqy);
        hq2->Scale(1./hq2->Integral(0,-1,"width"));
        hq2nonf->Scale(1./hq2nonf->Integral(0,-1,"width"));
/*ffit = new TF1(Form("ffit"),"1./(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
ffit0 = new TF1(Form("ffit0"),"1./(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
f1fit = new TF1(Form("f1fit"),"x/(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
f1fit0 = new TF1(Form("f1fit0"),"x/(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);*/
ffit = new TF1(Form("ffit"),"1./([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
ffit0 = new TF1(Form("ffit0"),"1./([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
f1fit = new TF1(Form("f1fit"),"x/([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
f1fit0 = new TF1(Form("f1fit0"),"x/([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
//ffit->SetParNames("g2","v2","M");
ffit->SetParNames("#sigma2","v2","M");
ffit->SetParameters(0.1,0.05,avgmult[xbin]);
ffit->FixParameter(2,avgmult[xbin]);
f1fit->SetParNames("sigma2","v2","M");
f1fit->SetParameters(0.1,0.05,avgmult[xbin]);
ffit->FixParameter(1,0);
f1fit->FixParameter(1,0);
//f1fit->FixParameter(0,0);
f1fit->FixParameter(2,avgmult[xbin]);
hq->Fit(Form("f1fit"),"R","P",0,10);
TCanvas *c2 = new TCanvas("c2","c2",500,1000);
c2->Divide(1,2);
c2->cd(1)->SetLogy();
fixedFontHist(hqx,1.6,2.0);
hqx->SetTitle("");
hqx->GetXaxis()->SetTitle("q_{x}");
hqx->GetYaxis()->SetTitle("#frac{dN}{dq_{x}}");
hqx->GetYaxis()->SetRangeUser(1e-10,1);
hqx->Fit(Form("f1fit"),"R","",0,10);
f1fit0->SetParameters(f1fit->GetParameter(0),inV2,f1fit->GetParameter(2));
f1fit0->SetLineColor(1);
f1fit0->SetLineStyle(2);
f1fit0->Draw("same");
c2->cd(2)->SetLogy();
TH1D* hqx_cp = (TH1D*)hqx->Clone("hqx_cp");
fixedFontHist(hqx_cp,1.6,2.0);
divideByBinCenter(hqx_cp);
hqx_cp->GetYaxis()->SetTitle("#frac{dN}{q_{x}dq_{x}}");
hqx_cp->GetYaxis()->SetRangeUser(1e-10,10);
hqx_cp->Fit(Form("ffit"),"R","",0,10);
ffit0->SetParameters(ffit->GetParameter(0),inV2,ffit->GetParameter(2));
ffit0->SetLineColor(1);
ffit0->SetLineStyle(2);
ffit0->Draw("same");
t->DrawLatex(0.5,0.2,Form("mult = %.f,input v_{2} = %.3f", avgmult[xbin], inV2));

TCanvas *c3 = new TCanvas("c3","c3",500,1000);
c3->Divide(1,2);
c3->cd(1);
fixedFontHist(hqy,1.6,2.0);
hqy->SetTitle("");
hqy->GetXaxis()->SetTitle("q_{y}");
hqy->GetYaxis()->SetTitle("#frac{dN}{dq_{y}}");
hqy->GetYaxis()->SetRangeUser(0,1);
hqy->SetMarkerStyle(24);
hqy->SetMarkerSize(0.5);
hqy->Fit(Form("f1fit"),"R","",0,10);
f1fit0->SetParameters(f1fit->GetParameter(0),inV2,f1fit->GetParameter(2));
f1fit0->SetLineColor(1);
f1fit0->SetLineStyle(2);
f1fit0->Draw("same");
hqy->Draw("Psame")
c3->cd(2)->SetLogy();
TH1D* hqy_cp = (TH1D*)hqy->Clone("hqy_cp");
fixedFontHist(hqy_cp,1.6,2.0);
divideByBinCenter(hqy_cp);
hqy_cp->GetYaxis()->SetTitle("#frac{dN}{q_{y}dq_{y}}");
hqy_cp->GetYaxis()->SetRangeUser(1e-10,10);
hqy_cp->SetMarkerStyle(24);
hqy_cp->SetMarkerSize(0.5);
hqy_cp->Fit(Form("ffit"),"R","",0,10);
ffit0->SetParameters(ffit->GetParameter(0),inV2,ffit->GetParameter(2));
ffit0->SetLineColor(1);
ffit0->SetLineStyle(2);
ffit0->Draw("same");
hqy_cp->Draw("Psame");
t->DrawLatex(0.5,0.2,Form("mult = %.f,input v_{2} = %.3f", avgmult[xbin], inV2));

TCanvas *c4 = new TCanvas("c4","c4",500,1000);
c4->Divide(1,2);
c4->cd(1);
fixedFontHist(hq2,1.6,2.0);
hq2->SetTitle("");
hq2->GetXaxis()->SetTitle("q_{2}");
hq2->GetYaxis()->SetTitle("#frac{dN}{dq_{2}}");
hq2->GetYaxis()->SetRangeUser(0,1);
hq2->Fit(Form("f1fit"),"R","",0,10);
f1fit0->SetParameters(f1fit->GetParameter(0),inV2,f1fit->GetParameter(2));
f1fit0->SetLineColor(1);
f1fit0->SetLineStyle(2);
f1fit0->Draw("same");
c4->cd(2)->SetLogy();
TH1D* hq2_cp = (TH1D*)hq2->Clone("hq2_cp");
fixedFontHist(hq2_cp,1.6,2.0);
divideByBinCenter(hq2_cp);
hq2_cp->GetYaxis()->SetTitle("#frac{dN}{q_{2}dq_{2}}");
hq2_cp->GetYaxis()->SetRangeUser(1e-10,10);
hq2_cp->Fit(Form("ffit"),"R","",0,10);
ffit0->SetParameters(ffit->GetParameter(0),inV2,ffit->GetParameter(2));
ffit0->SetLineColor(1);
ffit0->SetLineStyle(2);
ffit0->Draw("same");
t->DrawLatex(0.5,0.2,Form("mult = %.f,input v_{2} = %.3f", avgmult[xbin], inV2));

TCanvas *c5 = new TCanvas("c5","c5",500,1000);
c5->Divide(1,2);
c5->cd(1);
fixedFontHist(hq2nonf,1.6,2.0);
hq2nonf->SetTitle("");
hq2nonf->GetXaxis()->SetTitle("q_{2}");
hq2nonf->GetYaxis()->SetTitle("#frac{dN}{dq_{2}}");
hq2nonf->GetYaxis()->SetRangeUser(0,1);
hq2nonf->SetMarkerStyle(24);
hq2nonf->SetMarkerColor(4);
hq2nonf->SetLineColor(4);
hq2nonf->SetMarkerSize(0.5);
ffit->FixParameter(2,avgmult[xbin]*2);
f1fit->FixParameter(2,avgmult[xbin]*2);
ffit->FixParameter(1,0);
f1fit->FixParameter(1,0);
//ffit->FixParameter(0,1);
//f1fit->FixParameter(0,1);
hq2nonf->Fit(Form("f1fit"),"R","",0,10);
f1fit0->SetParameters(f1fit->GetParameter(0),inV2,f1fit->GetParameter(2));
f1fit0->SetLineColor(1);
f1fit0->SetLineStyle(2);
f1fit0->Draw("same");
hq2nonf->Draw("Psame");
c5->cd(2)->SetLogy();
TH1D* hq2nonf_cp = (TH1D*)hq2nonf->Clone("hq2nonf_cp");
fixedFontHist(hq2nonf_cp,1.6,2.0);
divideByBinCenter(hq2nonf_cp);
hq2nonf_cp->GetYaxis()->SetTitle("#frac{dN}{q_{2}dq_{2}}");
hq2nonf_cp->GetYaxis()->SetRangeUser(1e-10,10);
hq2nonf_cp->SetMarkerStyle(24);
hq2nonf_cp->SetMarkerColor(4);
hq2nonf_cp->SetLineColor(4);
hq2nonf_cp->SetMarkerSize(0.5);
hq2nonf_cp->Fit(Form("ffit"),"R","",0,10);
ffit0->SetParameters(ffit->GetParameter(0),inV2,ffit->GetParameter(2));
ffit0->SetLineColor(1);
ffit0->SetLineStyle(2);
ffit0->Draw("same");
hq2nonf_cp->Draw("Psame");
t->DrawLatex(0.5,0.2,Form("mult = %.f,input v_{2} = %.3f", avgmult[xbin], inV2));

c2->Print("hqx_fit.png");
c3->Print("hqy_fit.png");
c4->Print("hq2_fit.png");
c5->Print("hq2nonf_fit.png");
}

