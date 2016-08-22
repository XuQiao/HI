#include "par.h"
void plotV2vstheta(){
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    TLatex T;
    T.SetNDC();
    T.SetTextSize(0.045);
    int xbin = 0;
    TFile *f = new TFile("mergedV_Prod.root");
    TVectorD *vecV = (TVectorD*)f->Get(Form("D_%d/V"));
    TVectorD *vecdeltaV = (TVectorD*)f->Get(Form("D_%d/deltaV"));
    TVectorD *vecVmean = (TVectorD*)f->Get(Form("Vmean"));
    double *V = vecV->GetMatrixArray();
    double *deltaV = vecdeltaV->GetMatrixArray();
    double *Vmean = vecVmean->GetMatrixArray();
    double theta[ntheta];
    for(int itheta=0;itheta<ntheta;itheta++){
        theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
    }
    TH1D *hFrame = new TH1D("","",300,-1,2);
    hFrame->GetXaxis()->SetTitle("#theta");
    hFrame->GetYaxis()->SetTitle("referenceV_{2}");	
    hFrame->GetXaxis()->SetTitleSize(0.04);
    hFrame->GetYaxis()->SetTitleSize(0.04);
    hFrame->GetXaxis()->SetRangeUser(-0.1,1.5);
    hFrame->SetMaximum(0.060);
    hFrame->SetMinimum(0.040);
    hFrame->Draw();
    TGraphErrors *gr = new TGraphErrors(ntheta,theta,V,0,deltaV);
    gr->SetMarkerSize(1.2);
    gr->SetMarkerStyle(20);
    gr->Draw("Psame");
    T.DrawLatex(0.7,0.3,Form("mult = %d",(trkbin[xbin]+trkbin[xbin+1])/2));
    TLine *l = new TLine(0,inV2,1.4,inV2);
    l->SetLineStyle(2);
    l->Draw("same");
    c1->Print("V2vstheta.png");
}
