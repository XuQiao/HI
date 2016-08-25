#include "par.h"
#include <iomanip>
#include "../../../jetRpA/RpA/Quality/root_setting.h"

void plotvfit(int fixv2=0,int fixg2=0){
    const int nbin24 = 12;
    const double avgtrkbin[nbin24]={44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,269.2,301.2};
    const double V24[nbin24]={0.02965,0.03913,0.04832,0.04941,0.04822,0.04955,0.049,0.04805,0.04709,0.04665,0.04772,0.04797};
    int xtheta=0;
    gStyle->SetOptStat(1011);
    gStyle->SetOptFit(1111);
TFile *f = TFile::Open("mergedV_Sum.root");
TFile *fout = TFile::Open("vfitV.root","Update");
TVectorD *vecDNevent = (TVectorD*)f->Get(Form("Nevent"));
TVectorD *vecDavgmult = (TVectorD*)f->Get(Form("avgmultall"));
TVectorD *vecDtotmult = (TVectorD*)f->Get(Form("totmultall"));
TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("avgpt"));
TVectorD *vecDavgpt2 = (TVectorD*)f->Get(Form("avgpt2"));
TVectorD *vecDavgtrk = (TVectorD*)f->Get(Form("avgtrk"));
double avgmultall = vecDtotmult->Sum()/vecDNevent->Sum();
double *avgmult = vecDavgmult->GetMatrixArray();
double *avgtrk = vecDavgtrk->GetMatrixArray();
double *avgpt = vecDavgpt->GetMatrixArray();
double *avgpt2 = vecDavgpt2->GetMatrixArray();
TLatex *t= new TLatex();
t->SetNDC();
t->SetTextSize(0.04);
t->SetTextFont(42);
for(int ibin=0;ibin<nbin;ibin++){	//ibin<1
 //   if(ibin!=100) continue;
TH1D* hvq2 = (TH1D*)f->Get(Form("D_%d/hvq2",ibin));
int k;
for(k=0;k<nbin24;k++){
    if(avgtrk[ibin]>avgtrkbin[k]&& avgtrk[ibin]<=avgtrkbin[k+1])
        break;
}
        hvq2->Scale(1./hvq2->Integral(0,-1,"width"));
ffit = new TF1(Form("ffit"),"1./([0])*TMath::Exp(-([1]*[1]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]/[0])",0,1);
f1fit = new TF1(Form("f1fit"),"x/([0])*TMath::Exp(-([1]*[1]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]/[0])",0,1);
ffit->SetParNames("#sigma2","v2");
f1fit->SetParNames("#sigma2","v2");
ffit->SetParameters(0.1,0.05);
f1fit->SetParameters(0.1,0.05);
if(fixv2){
ffit->FixParameter(1,V24[k]);
f1fit->FixParameter(1,V24[k]);
}
if(fixg2){
    double sigma2 = 1./avgmult[ibin]*avgpt2[ibin]/2./avgpt[ibin]/avgpt[ibin];
    ffit->FixParameter(0,sigma2);
    f1fit->FixParameter(0,sigma2);
}

TCanvas *c4 = new TCanvas("c4","c4",1000,500);
c4->Divide(2,1);
c4->cd(1)->SetLeftMargin(0.18);
fixedFontHist(hvq2,1.0,2.0);
hvq2->SetTitle("");
hvq2->GetXaxis()->SetTitle("v_{2}^{obs}");
hvq2->GetYaxis()->SetTitle("#frac{dN}{dv_{2}^{obs}}");
hvq2->GetYaxis()->SetRangeUser(0,10);
hvq2->SetMarkerStyle(24);
hvq2->SetMarkerColor(4);
hvq2->SetLineColor(4);
hvq2->SetMarkerSize(0.5);
hvq2->Fit(Form("f1fit"),"R","",0,10);
TVectorD vecr;
vecr.ResizeTo(8);
vecr[0]=f1fit->GetParameter(0);
vecr[1]=f1fit->GetParError(0);
vecr[2]=f1fit->GetParameter(1);
vecr[3]=f1fit->GetParError(1);
vecr[4]=avgmult[ibin];
vecr[5]=avgtrk[ibin];
hvq2->Draw("Psame");
c4->cd(2)->SetLogy();
c4->cd(2)->SetLeftMargin(0.2);
TH1D* hvq2_cp = (TH1D*)hvq2->Clone("hvq2_cp");
fixedFontHist(hvq2_cp,1.0,2.0);
divideByBinCenter(hvq2_cp);
hvq2_cp->GetYaxis()->SetTitle("#frac{dN}{v_{2}^{obs}dv_{2}^{obs}}");
hvq2_cp->GetYaxis()->SetRangeUser(1e-10,100);
hvq2_cp->SetMarkerStyle(24);
hvq2_cp->SetMarkerColor(4);
hvq2_cp->SetLineColor(4);
hvq2_cp->SetMarkerSize(0.5);
hvq2_cp->Fit(Form("ffit"),"R","",0,10);
hvq2_cp->Draw("Psame");
t->DrawLatex(0.5,0.2,Form("N_{trk}^{offline} = %.2f", avgtrk[ibin]));

fout->cd();
vecr.Write(Form("r_%d_%d_%d",ibin,fixv2,fixg2),TObject::kOverwrite);
//if(ibin==15)
//c4->Print("hvq2_fit_15.png");
}
}


