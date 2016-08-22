#include "TString.h"
#include "TVector.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include <fstream>
#include <iostream>

double MHATg1(double);
double MHATg2(double);
double MHATg4(double);

void MHATfunc(){
    int ifile=0;
    TFile *f = TFile::Open(Form("/cms/store/user/qixu/flow/NewSTEG/pPbDataV205m100/vndata_50k_%d.root",ifile));
    TTree *t = (TTree*)f->Get("tree");
    TFile *fout = new TFile("testV205m100.root","Recreate");
    const int nbinx=50;
    const int nbiny=50;
    //Float_t eta[1000]={1,2.75,3.25,5,6,7};
    Float_t eta[1000];
    Float_t phi[1000];
    int n;
    t->SetBranchAddress("etag",eta);
    t->SetBranchAddress("phig",phi);
    t->SetBranchAddress("n",&n);
    //Int_t N = t->GetEntries();
    Int_t N = 10000;
    double Wg2[nbinx][nbiny] = {};
    double W2g2[nbinx][nbiny] = {};
    double EWg2[nbinx] = {};
    double Wg4[nbinx][nbiny] = {};
    double W2g4[nbinx][nbiny] = {};
    double EWg4[nbinx] = {};
    double a=0.2;
    double amin=0, amax=5, bmin=0, bmax=2*TMath::Pi();
    double astep = (amax-amin)/nbinx;
    double bstep = (bmax-bmin)/nbiny;
 //   TGraph2D* gr2D = new TGraph2D(nbinx*nbiny);
    TH2D* hWg22D = new TH2D("hMHATg2","",nbinx,amin,amax,nbiny,bmin,bmax);
    TH2D* hW2g22D = new TH2D("hMHATg2sq","",nbinx,amin,amax,nbiny,bmin,bmax);
    TH1D* hWg21D = new TH1D("hMHATg21D","",nbiny,bmin,bmax);
//    gr2D->SetName("MHATg2");
//    TGraph2D* grW22D = new TGraph2D(nbinx*nbiny);
    TH2D* hWg42D = new TH2D("hMHATg4","",nbinx,amin,amax,nbiny,bmin,bmax);
    TH2D* hW2g42D = new TH2D("hMHATg4sq","",nbinx,amin,amax,nbiny,bmin,bmax);
    TH1D* hWg41D = new TH1D("hMHATg41D","",nbiny,bmin,bmax);
//    grW22D->SetName("MHATg2sq");
//    TGraph* grEW = new TGraph(nbinx);
    TH1D* hEWg21D = new TH1D("hMHAThg2EW","",nbinx,amin,amax);
    TH1D* hEWg41D = new TH1D("hMHAThg4EW","",nbinx,amin,amax);
//    grEW->SetName("MHATg2E");
    int i=0;
    for(double ia=amin; ia<=amax;ia+=astep){
        if(ia==0)continue;
       // if(ia!=a)continue;
        int j=0;
        for(double ib=bmin; ib<=bmax;ib+=bstep){
        for(int imult = 0;imult<n; imult++){
    //    h->Fill(eta[imult],phi[imult]);
        double x = (phi[imult]-ib)/ia;
    	Wg2[i][j] += 1./n*TMath::Sqrt(1./ia)*MHATg2(x);
    	W2g2[i][j] += TMath::Power(1./n*TMath::Sqrt(1./ia)*MHATg2(x),2);
    	Wg4[i][j] += 1./n*TMath::Sqrt(1./ia)*MHATg4(x);
    	W2g4[i][j] += TMath::Power(1./n*TMath::Sqrt(1./ia)*MHATg4(x),2);
       }
    }
    Wg2[i][j]/=N;
    W2g2[i][j]/=N;
    Wg4[i][j]/=N;
    W2g4[i][j]/=N;
    EWg2[i]+=Wg2[i][j];
    EWg4[i]+=Wg4[i][j];
//    gr2D->SetPoint(j+i*nbiny,ia,ib,W[i][j]);
//    grW22D->SetPoint(j+i*nbiny,ia,ib,W2[i][j]);
    hWg22D->SetBinContent(i+1,j+1,Wg2[i][j]);
    hW2g22D->SetBinContent(i+1,j+1,W2g2[i][j]);
    hWg42D->SetBinContent(i+1,j+1,Wg4[i][j]);
    hW2g42D->SetBinContent(i+1,j+1,W2g4[i][j]);
    if(ia==a){
        hWg21D->SetBinContent(j+1,Wg2[i][j]);
        hWg41D->SetBinContent(j+1,Wg4[i][j]);
    }
    j++;
    }
//    grEW->SetPoint(i,ia,EW[i]);
    hEWg21D->SetBinContent(i+1,EWg2[i]);
    hEWg41D->SetBinContent(i+1,EWg4[i]);
    cout<<i<<endl;
    i++;
    }

    fout->cd();
    hWg22D->Write();
    hW2g22D->Write();
    hEWg21D->Write();
    hWg22D->Write();
    hWg42D->Write();
    hW2g42D->Write();
    hEWg41D->Write();
    hWg21D->Write();
    hWg41D->Write();
    fout->Close();
}

double MHATg1(double x){
    return TMath::Exp(-x*x/2);
}
double MHATg2(double x){
     return (1-x*x)*(TMath::Exp(-x*x/2));
}

double MHATg4(double x){
     return ((6-x*x)*x*x-3)*(TMath::Exp(-x*x/2));
}

