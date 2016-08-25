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

void DrawphiDis(){
    int ifile=0;
    TFile *f = TFile::Open(Form("/cms/store/user/qixu/flow/NewSTEG/pPbDataV205m300/vndata_50k_%d.root",ifile));
    TTree *t = (TTree*)f->Get("tree");
    Float_t eta[1000];
    Float_t phi[1000];
    TH1D* hphi = new TH1D("hphi","hphi",100,-2.4,2.4);
    int n;
    t->SetBranchAddress("etag",eta);
    t->SetBranchAddress("phig",phi);
    t->SetBranchAddress("n",&n);
    Int_t N = t->GetEntries();
    for(int ievt=0;ievt<N;ievt++){
  //      if(ievt!=0) continue;
   // if(ievt%50000==0)        cout<<"Processing "<<ievt<<" events"<<endl;
      t->GetEntry(ievt);
       for(int imult = 0;imult<n; imult++){
           hphi->Fill(eta[imult]);
       }
    }
//    cout<<hphi->Integral("width")/2/TMath::Pi()<<endl;
//    cout<<hphi->GetMean(2)<<endl;
    hphi->Draw();
}

