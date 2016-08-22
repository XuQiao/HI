#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TMath.h"
#include "SetupTree.h"
#include "TDirectory.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "SetSeedOwn.C"
#include "par.h"

using namespace std;

class ridge{
    public:
        ridge(TString);
        ~ridge();
        TVectorD Nevent, tottrk;
        TH2F* s[nbin];
        TH2F* b[nbin];
        treeInt *t;
        Int_t nEvents;
        vector<vector<double> > pVectEta;
        vector<vector<double> > pVectPhi;
        vector<int> pNtrig;
        void beginJob();
        void calcS();
        void calcB();
        void endJob(TString);

    private:
    TString filename;
};

ridge::ridge(TString _filename){
    filename = _filename;
}

ridge::~ridge(){
    delete t;
    pVectEta.clear();
    vector<vector<double> >().swap(pVectEta);
    pVectPhi.clear();
    vector<vector<double> >().swap(pVectPhi);
    pNtrig.clear();
    vector<int>().swap(pNtrig);
    for(int ibin=0;ibin<nbin;ibin++){
    delete s[ibin];
    delete b[ibin];
    }
}

void ridge::beginJob(){
  for(int ibin=0;ibin<nbin;ibin++){
    s[ibin] = new TH2F(Form("signal_%d",ibin),Form("signal_%d",ibin),detastep,detamin,detamax,dphistep,dphimin,dphimax);
    s[ibin]->Sumw2();
    b[ibin] = new TH2F(Form("background_%d",ibin),Form("background_%d",ibin),detastep,detamin,detamax,dphistep,dphimin,dphimax);
    b[ibin]->Sumw2();
  }
  Nevent.ResizeTo(nbin); Nevent.Zero();
  tottrk.ResizeTo(nbin); tottrk.Zero();
  t = new treeInt(Form("%s",filename.Data()));
  t->Setup();
  nEvents = t->GetEntries();
    for(Int_t i=0;i<nEvents;i++){
  //      if(i%10000==0) cout<<"has processed "<<i<<"begin events"<<endl;
        t->GetEntry(i);
        vector<double> pvectEta;
        vector<double> pvectPhi;
        int Ntrig = 0;
        for(int imult=0;imult<t->mult;imult++){
	    if(t->eta[imult]<-2.40||t->eta[imult]>2.40) continue;
            if(t->pt[imult]<ptmin||t->pt[imult]>ptmax) continue; //event selection
            if(t->phi[imult]>0 && t->phi[imult]<TMath::Pi()/2) continue;
            if(t->phi[imult]>TMath::Pi() && t->phi[imult]<3*TMath::Pi()/2) continue;
           // TVector3 pvector;
           // pvector.SetPtEtaPhi(t->pt[imult],t->eta[imult],t->phi[imult]);
            pvectEta.push_back(t->eta[imult]);
            pvectPhi.push_back(t->phi[imult]);
            Ntrig++;
        }
        pVectEta.push_back(pvectEta);
        pVectPhi.push_back(pvectPhi);
        pNtrig.push_back(Ntrig);
    }
}

void ridge::calcS(){
    for(Int_t i=0;i<nEvents;i++){
      //  if(i%10000==0) cout<<"has processed "<<i<<"S events"<<endl;
        if(pNtrig[i]<=2) continue;
        int ntrk = pNtrig[i]; int xbin=-1;
        for(int j=0;j<nbin;j++)
            if(ntrk<trkbin[j]&&ntrk>=trkbin[j+1])
                xbin=j;
            if(xbin<0 || xbin==nbin) continue;
            tottrk[xbin]+=ntrk;
            vector<double> pvectorEta = pVectEta[i];
            vector<double> pvectorPhi = pVectPhi[i];
            int Ntrig = pNtrig[i];
        for(int imult1=0;imult1<Ntrig;imult1++){
              double eta1 = pvectorEta[imult1];
              double phi1 = pvectorPhi[imult1];
            for(int imult2=0;imult2<Ntrig;imult2++){
                if(imult1==imult2) continue;
              double eta2 = pvectorEta[imult2];
              double phi2 = pvectorPhi[imult2];
              s[xbin]->Fill(eta1-eta2,phi1-phi2,1./Ntrig);
            }
        }
        Nevent[xbin]++;
    }
}


void ridge::calcB(){
    for(Int_t i=0;i<nEvents;i++){
   //     if(i%5000==0) cout<<"has processed "<<i<<" B events"<<endl;
        if(pNtrig[i]<=2) continue;
        t->GetEntry(i);
        int ntrk = pNtrig[i]; int xbin=-1;
        for(int k=0;k<nbin;k++)
            if(ntrk<trkbin[k]&&ntrk>=trkbin[k+1])
                xbin=k;
            if(xbin<0 || xbin==nbin) continue;
            for(int ira=0;ira<10;ira++){
            UInt_t iniseed = SetSeedOwn();
            gRandom->SetSeed(iniseed);
            TRandom3* r = new TRandom3(0);
            int j = r->Integer(nEvents);
            if(i==j) continue;
            vector<double> pvectorEta_trig = pVectEta[i];
            vector<double> pvectorPhi_trig = pVectPhi[i];
            vector<double> pvectorEta_ass = pVectEta[j];
            vector<double> pvectorPhi_ass = pVectPhi[j];
            int Ntrig = pNtrig[i];
            int Ntrig_ass = pNtrig[j];
            for(int imult_trig=0;imult_trig<Ntrig;imult_trig++){
                double eta_trig = pvectorEta_trig[imult_trig];
                double phi_trig = pvectorPhi_trig[imult_trig];
                for(int imult_ass=0;imult_ass<Ntrig_ass;imult_ass++){
                double eta_ass = pvectorEta_ass[imult_ass];
                double phi_ass = pvectorPhi_ass[imult_ass];
            b[xbin]->Fill(eta_trig-eta_ass,phi_trig-phi_ass,1./Ntrig);
                }
            }
    }
    }
}

void ridge::endJob(TString outfile){
    TFile *fout = new TFile(Form("%s",outfile.Data()),"Recreate");
    fout->cd();
    Nevent.Write("Nevent");
    tottrk.Write("tottrk");
    for(int ibin=0;ibin<nbin;ibin++){
        s[ibin]->Write();
        b[ibin]->Write();
    }
    fout->Close();
}
    

