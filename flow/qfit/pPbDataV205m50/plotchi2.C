#include "par.h"
#include <iomanip>
#include "../../../jetRpA/RpA/Quality/root_setting.h"

int xbin=0;	//xbin<1
int xtheta=0;
void plotchi2(){
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
ffit = new TF1(Form("ffit"),"1./(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
ffit0 = new TF1(Form("ffit0"),"1./(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
f1fit = new TF1(Form("f1fit"),"x/(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
f1fit0 = new TF1(Form("f1fit0"),"x/(0.5*(1+[0]))*TMath::Exp(-([1]*[1]*[2]+x*x)/(1+[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/(0.5*(1+[0])))",0,10);
/*ffit = new TF1(Form("ffit"),"1./([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
ffit0 = new TF1(Form("ffit0"),"1./([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
f1fit = new TF1(Form("f1fit"),"x/([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
f1fit0 = new TF1(Form("f1fit0"),"x/([0])*TMath::Exp(-([1]*[1]*[2]+x*x)/(2*[0]))*TMath::BesselI0(x*[1]*TMath::Sqrt([2])/([0]))",0,10);
*/

ffit->SetParNames("g2","v2","M");
//ffit->SetParNames("#sigma2","v2","M");
ffit->FixParameter(2,avgmult[xbin]);
//f1fit->SetParNames("sigma2","v2","M");
f1fit->SetParNames("g2","v2","M");
//ffit->FixParameter(1,0);
//f1fit->FixParameter(1,0);
//f1fit->FixParameter(0,0);
f1fit->FixParameter(2,avgmult[xbin]);

TCanvas *c4 = new TCanvas("c4","c4",800,800);
//c4->Divide(1,2);
//c4->cd(1);
const int Np = 10000;
int Nv2 = 100;
int Ng2 = 100;
double v2[Np],g2[Np],chi2[Np];
TH1D* hq2fit = new TH1D("hq2fit","",1000,0,10);
int cv2=0;
for(double iv2=inV2-inV2/2;iv2<inV2+inV2/2;iv2+=inV2/Nv2){
    int cg2=0;
    for(double ig2=-0.5;ig2<0.5;ig2+=1./Ng2){
        v2[cg2+cv2*Nv2] = iv2;
        g2[cg2+cv2*Nv2] = ig2;
        ffit->SetParameters(iv2,ig2,avgmult[xbin]);
        hq2fit->Reset();
        hq2fit->Add(ffit);
        chi2[cg2+cv2*Nv2]=hq2->Chi2Test(hq2fit,"CHI2/NDF");
        cg2++;
    }
    cv2++;
}
TGraph2D *gr1 = new TGraph2D(Np,v2,g2,chi2); 
gr1->SetName("gr1");
//c4->cd(2);
//TH1D* hq2_cp = (TH1D*)hq2->Clone("hq2_cp");
//divideByBinCenter(hq2_cp);
t->DrawLatex(0.5,0.2,Form("mult = %.f,input v_{2} = %.3f", avgmult[xbin], inV2));

TCanvas *c5 = new TCanvas("c5","c5",800,800);
//c5->Divide(1,2);
//c5->cd(1);
//c5->cd(2)->SetLogy();
TH1D* hq2nonffit = new TH1D("hq2nonffit","",1000,0,10);
int cv2=0;
for(double iv2=inV2-inV2/2;iv2<inV2+inV2/2;iv2+=inV2/Nv2){
    int cg2=0;
    for(double ig2=0;ig2<1;ig2+=1./Ng2){
        v2[cg2+cv2*Nv2] = iv2;
        g2[cg2+cv2*Nv2] = ig2;
        ffit->SetParameters(iv2,ig2,avgmult[xbin]);
        hq2nonffit->Reset();
        hq2nonffit->Add(ffit);
        chi2[cg2+cv2*Nv2]=hq2->Chi2Test(hq2fit,"CHI2/NDF");
        cg2++;
    }
    cv2++;
}
TGraph2D *gr2 = new TGraph2D(Np,v2,g2,chi2); 
gr2->SetName("gr2");
//c4->cd(2);

//TH1D* hq2nonf_cp = (TH1D*)hq2nonf->Clone("hq2nonf_cp");
//divideByBinCenter(hq2nonf_cp);
t->DrawLatex(0.5,0.2,Form("mult = %.f,input v_{2} = %.3f", avgmult[xbin], inV2));
TFile *fout = new TFile("TGraphchi2.root","Recreate");
gr1->Write();
gr2->Write();
fout->Close();
}

