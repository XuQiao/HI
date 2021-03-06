#include "par.h"
#include <iomanip>
#include "../../../jetRpA/RpA/Quality/root_setting.h"

void plotqfit(int fixv2=0,int fixg2=0){
    const int nbin24 = 12;
    const double avgtrkbin[nbin24]={44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,269.2,301.2};
    const double V24[nbin24]={0.02965,0.03913,0.04832,0.04941,0.04822,0.04955,0.049,0.04805,0.04709,0.04665,0.04772,0.04797};
    int xtheta=0;
    gStyle->SetOptStat(1011);
    gStyle->SetOptFit(1111);
TFile *f = TFile::Open("mergedV_Sum.root");
TFile *fout = TFile::Open("qfitV.root","Update");
TVectorD *vecDNevent = (TVectorD*)f->Get(Form("Nevent"));
TVectorD *vecDavgmult = (TVectorD*)f->Get(Form("avgmultall"));
TVectorD *vecDtotmult = (TVectorD*)f->Get(Form("totmultall"));
TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("avgpt"));
TVectorD *vecDavgpt2 = (TVectorD*)f->Get(Form("avgpt2"));
TVectorD *vecDavgtrk = (TVectorD*)f->Get(Form("avgtrk"));
TVectorD *vecDq22 = (TVectorD*)f->Get(Form("q22"));
TVectorD *vecDq24 = (TVectorD*)f->Get(Form("q24"));
double avgmultall = vecDtotmult->Sum()/vecDNevent->Sum();
cout<<avgmultall<<endl;
double *avgmult = vecDavgmult->GetMatrixArray();
double *avgtrk = vecDavgtrk->GetMatrixArray();
double *avgpt = vecDavgpt->GetMatrixArray();
double *avgpt2 = vecDavgpt2->GetMatrixArray();
double *q22 = vecDq22->GetMatrixArray();
double *q24 = vecDq24->GetMatrixArray();
TLatex *t= new TLatex();
t->SetNDC();
t->SetTextSize(0.04);
t->SetTextFont(42);
/*
if(fixg2){
TH1D* hq2all = new TH1D("hq2all","hq2all",1000,0,10);
for(int ibin=0;ibin<nbin;ibin++){	//ibin<1
TH1D* hq2 = (TH1D*)f->Get(Form("D_%d/hq2",ibin));
hq2all->Add(hq2);
}
TF1 *ffit = new TF1(Form("ffit"),"1./(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
TF1* f1fit = new TF1(Form("f1fit"),"x/(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
TH1D* hq2all_cp = (TH1D*)hq2all->Clone("hq2all_cp");
hq2all->Scale(1./hq2all->Integral(0,-1,"width"));
hq2all->Draw();
divideByBinCenter(hq2all_cp);
ffit->SetParameters(1,0.05,avgmultall);
f1fit->SetParameters(1,0.05,avgmultall);
ffit->FixParameter(2,avgmultall);
f1fit->FixParameter(2,avgmultall);
hq2all->Fit(Form("f1fit"),"R","",0,10);
double g2all = f1fit->GetParameter(0);
}
*/
for(int ibin=0;ibin<nbin;ibin++){	//ibin<1
 //   if(ibin!=100) continue;
TH1D* hq = (TH1D*)f->Get(Form("D_%d/D_%d/hq",ibin,xtheta));
TH1D* hqx = (TH1D*)f->Get(Form("D_%d/hqx",ibin));
TH1D* hqy = (TH1D*)f->Get(Form("D_%d/hqy",ibin));
TH1D* hq2 = (TH1D*)f->Get(Form("D_%d/hq2",ibin));
TH1D* hq2nonf = (TH1D*)f->Get(Form("D_%d/hq2nonf",ibin));
int k;
for(k=0;k<nbin24;k++){
    if(avgtrk[ibin]>avgtrkbin[k]&& avgtrk[ibin]<=avgtrkbin[k+1])
        break;
}
        //multiplyByBinCenter(hq);
        hq->Scale(1./hq->Integral(0,-1,"width"));
        //normalizeByBinWidth(hq);
        //multiplyByBinCenter(hqx);
        hqx->Scale(1./hqx->Integral(0,-1,"width"));
        //normalizeByBinWidth(hqx);
        //multiplyByBinCenter(hqy);
        hqy->Scale(1./hqy->Integral(0,-1,"width"));
        //normalizeByBinWidth(hqy);
        hq2->Scale(1./hq2->Integral(0,-1,"width"));
        hq2nonf->Scale(1./hq2nonf->Integral(0,-1,"width"));
TF1 *ffit = new TF1(Form("ffit"),"1./(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
TF1* f1fit = new TF1(Form("f1fit"),"x/(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
//ffit = new TF1(Form("ffit"),"1./([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
//f1fit = new TF1(Form("f1fit"),"x/([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
ffit->SetParNames("g2","v2","M");
f1fit->SetParNames("g2","v2","M");
//ffit->SetParNames("#sigma2","v2","M");
//f1fit->SetParNames("sigma2","v2","M");
ffit->SetParameters(0.1,0.05,avgmult[ibin]);
f1fit->SetParameters(0.1,0.05,avgmult[ibin]);
if(fixv2){
ffit->FixParameter(1,V24[k]);
f1fit->FixParameter(1,V24[k]);
}
if(fixg2){
    double sigma2 = avgpt2[ibin]/2./avgpt[ibin]/avgpt[ibin];
    double g2all = 2*sigma2-1;
    ffit->FixParameter(0,g2all);
    f1fit->FixParameter(0,g2all);
}

//f1fit->FixParameter(0,0);
//f1fit->FixParameter(0,0);
ffit->FixParameter(2,avgmult[ibin]);
f1fit->FixParameter(2,avgmult[ibin]);
/*
hq->Fit(Form("f1fit"),"R","P",0,10);
TCanvas *c2 = new TCanvas("c2","c2",1000,500);
c2->Divide(2,1);
c2->cd(1)->SetLogy();
fixedFontHist(hqx,1.0,2.0);
hqx->SetTitle("");
hqx->GetXaxis()->SetTitle("q_{x}");
hqx->GetYaxis()->SetTitle("#frac{dN}{dq_{x}}");
hqx->GetYaxis()->SetRangeUser(1e-10,1);
hqx->Fit(Form("f1fit"),"R","",0,10);
c2->cd(2)->SetLogy();
TH1D* hqx_cp = (TH1D*)hqx->Clone("hqx_cp");
fixedFontHist(hqx_cp,1.0,2.0);
divideByBinCenter(hqx_cp);
hqx_cp->GetYaxis()->SetTitle("#frac{dN}{q_{x}dq_{x}}");
hqx_cp->GetYaxis()->SetRangeUser(1e-10,10);
hqx_cp->Fit(Form("ffit"),"R","",0,10);
t->DrawLatex(0.5,0.2,Form("N_{trk}^{offline} = %.2f", avgtrk[ibin]));

TCanvas *c3 = new TCanvas("c3","c3",1000,500);
c3->Divide(2,1);
c3->cd(1);
fixedFontHist(hqy,1.0,2.0);
hqy->SetTitle("");
hqy->GetXaxis()->SetTitle("q_{y}");
hqy->GetYaxis()->SetTitle("#frac{dN}{dq_{y}}");
hqy->GetYaxis()->SetRangeUser(0,1);
hqy->SetMarkerStyle(24);
hqy->SetMarkerSize(0.5);
hqy->Fit(Form("f1fit"),"R","",0,10);
hqy->Draw("Psame");
c3->cd(2)->SetLogy();
TH1D* hqy_cp = (TH1D*)hqy->Clone("hqy_cp");
fixedFontHist(hqy_cp,1.0,2.0);
divideByBinCenter(hqy_cp);
hqy_cp->GetYaxis()->SetTitle("#frac{dN}{q_{y}dq_{y}}");
hqy_cp->GetYaxis()->SetRangeUser(1e-10,10);
hqy_cp->SetMarkerStyle(24);
hqy_cp->SetMarkerSize(0.5);
hqy_cp->Fit(Form("ffit"),"R","",0,10);
hqy_cp->Draw("Psame");
t->DrawLatex(0.5,0.2,Form("N_{trk}^{offline} = %.2f",avgtrk[ibin]));
*/
TCanvas *c4 = new TCanvas("c4","c4",1000,500);
c4->Divide(2,1);
c4->cd(1)->SetLeftMargin(0.18);
fixedFontHist(hq2,1.0,2.0);
hq2->SetTitle("");
hq2->GetXaxis()->SetTitle("q_{2}");
hq2->GetYaxis()->SetTitle("#frac{dN}{dq_{2}}");
hq2->GetYaxis()->SetRangeUser(0,1);
hq2->SetMarkerStyle(24);
hq2->SetMarkerColor(4);
hq2->SetLineColor(4);
hq2->SetMarkerSize(0.5);
hq2->Fit(Form("f1fit"),"R","",0,10);
TVectorD vecr;
vecr.ResizeTo(8);
vecr[0]=f1fit->GetParameter(0);
vecr[1]=f1fit->GetParError(0);
vecr[2]=f1fit->GetParameter(1);
vecr[3]=f1fit->GetParError(1);
vecr[4]=avgmult[ibin];
vecr[5]=avgtrk[ibin];
hq2->Draw("Psame");
c4->cd(2)->SetLogy();
c4->cd(2)->SetLeftMargin(0.2);
TH1D* hq2_cp = (TH1D*)hq2->Clone("hq2_cp");
fixedFontHist(hq2_cp,1.0,2.0);
divideByBinCenter(hq2_cp);
hq2_cp->GetYaxis()->SetTitle("#frac{dN}{q_{2}dq_{2}}");
hq2_cp->GetYaxis()->SetRangeUser(1e-10,10);
hq2_cp->SetMarkerStyle(24);
hq2_cp->SetMarkerColor(4);
hq2_cp->SetLineColor(4);
hq2_cp->SetMarkerSize(0.5);
hq2_cp->Fit(Form("ffit"),"R","",0,10);
hq2_cp->Draw("Psame");
t->DrawLatex(0.5,0.2,Form("N_{trk}^{offline} = %.2f", avgtrk[ibin]));

TCanvas *c5 = new TCanvas("c5","c5",1000,500);
c5->Divide(2,1);
c5->cd(1);
fixedFontHist(hq2nonf,1.0,2.0);
hq2nonf->SetTitle("");
hq2nonf->GetXaxis()->SetTitle("q_{2}");
hq2nonf->GetYaxis()->SetTitle("#frac{dN}{dq_{2}}");
hq2nonf->GetYaxis()->SetRangeUser(0,1);
hq2nonf->SetMarkerStyle(24);
hq2nonf->SetMarkerColor(4);
hq2nonf->SetLineColor(4);
hq2nonf->SetMarkerSize(0.5);
ffit->SetParameters(1.5,0.05,avgmult[ibin]*2);
f1fit->SetParameters(1.5,0.05,avgmult[ibin]*2);
ffit->FixParameter(2,avgmult[ibin]*2);
f1fit->FixParameter(2,avgmult[ibin]*2);
if(fixv2){
ffit->FixParameter(1,V24[k]);
f1fit->FixParameter(1,V24[k]);
}
if(fixg2){
    ffit->FixParameter(0,g2all);
    f1fit->FixParameter(0,g2all);
}
//ffit->FixParameter(0,1);
//f1fit->FixParameter(0,1);
hq2nonf->Fit(Form("f1fit"),"R","",0,10);
TVectorD vecrnonf;
vecrnonf.ResizeTo(8);
vecrnonf[0]=f1fit->GetParameter(0);
vecrnonf[1]=f1fit->GetParError(0);
vecrnonf[2]=f1fit->GetParameter(1);
vecrnonf[3]=f1fit->GetParError(1);
vecrnonf[4]=avgmult[ibin];
vecrnonf[5]=avgtrk[ibin];
hq2nonf->Draw("Psame");
c5->cd(2)->SetLogy();
TH1D* hq2nonf_cp = (TH1D*)hq2nonf->Clone("hq2nonf_cp");
fixedFontHist(hq2nonf_cp,1.0,2.0);
divideByBinCenter(hq2nonf_cp);
hq2nonf_cp->GetYaxis()->SetTitle("#frac{dN}{q_{2}dq_{2}}");
hq2nonf_cp->GetYaxis()->SetRangeUser(1e-10,10);
hq2nonf_cp->SetMarkerStyle(24);
hq2nonf_cp->SetMarkerColor(4);
hq2nonf_cp->SetLineColor(4);
hq2nonf_cp->SetMarkerSize(0.5);
hq2nonf_cp->Fit(Form("ffit"),"R","",0,10);
hq2nonf_cp->Draw("Psame");
t->DrawLatex(0.5,0.2,Form("N_{trk}^{offline} = %.2f", avgtrk[ibin]*2));

double v2calc = TMath::Sqrt(TMath::Sqrt(2*q22[ibin]*q22[ibin]-q24[ibin])/avgmult[ibin]);
double g2calc = q22[ibin]-TMath::Sqrt(2*q22[ibin]*q22[ibin]-q24[ibin])-1;
vecr[6]=v2calc;
vecr[7]=g2calc;

fout->cd();
vecr.Write(Form("r_%d_%d_%d",ibin,fixv2,fixg2),TObject::kOverwrite);
vecrnonf.Write(Form("rnonf_%d_%d_%d",ibin,fixv2,fixg2),TObject::kOverwrite);
if(ibin==15)
c4->Print("hq2_fit_15.png");
}
//c5->Print("hq2nonf_fit.png");
}


