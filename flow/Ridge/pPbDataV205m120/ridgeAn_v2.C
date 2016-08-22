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

class ridge{  //use two trees
    public:
        ridge(TString);
        ~ridge();
        TVectorD Nevent, tottrk;
        TH2F* s[nbin];
        TH2F* b[nbin];
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
}

void ridge::calcS(){
  treeInt *t = new treeInt(Form("%s",filename.Data()));
  t->Setup();
  Int_t nEvents = t->GetEntries();
    for(Int_t i=0;i<nEvents;i++){
        if(i%10000==0) cout<<"has processed "<<i<<" S events"<<endl;
        t->GetEntry(i);
        int ntrk = t->mult; int xbin=-1 ;
        for(int j=0;j<nbin;j++)
            if(ntrk<trkbin[j]&&ntrk>=trkbin[j+1])
                xbin=j;
            if(xbin<0 || xbin==nbin) continue;
            tottrk[xbin]+=ntrk;
        int Ntrig = 0;
        for(int imult_trig=0;imult_trig<t->mult;imult_trig++){
	    if(t->eta[imult_trig]<etatrigmin||t->eta[imult_trig]>etatrigmax) continue;
            if(t->pt[imult_trig]<pttrigmin||t->pt[imult_trig]>pttrigmax) continue; //event selection
            Ntrig++;
        }
            if(Ntrig<=2) continue;
        for(int imult_trig=0;imult_trig<t->mult;imult_trig++){
              double eta_trig = t->eta[imult_trig];
              double phi_trig = t->phi[imult_trig];
	    if(eta_trig<etatrigmin||eta_trig>etatrigmax) continue;
            if(t->pt[imult_trig]<pttrigmin||t->pt[imult_trig]>pttrigmax) continue; //event selection
            for(int imult_ass=0;imult_ass<t->mult;imult_ass++){
              double eta_ass = t->eta[imult_ass];
              double phi_ass = t->phi[imult_ass];
	    if(eta_ass<etaassmin||eta_ass>etaassmax) continue;
            if(t->pt[imult_ass]<ptassmin||t->pt[imult_ass]>ptassmax) continue; //event selection
                double deltaeta = eta_trig - eta_ass;
                double deltaphi = phi_trig - phi_ass;
              if(deltaphi>TMath::Pi())
                  deltaphi=deltaphi-2*TMath::Pi();
              if(deltaphi<-TMath::Pi()/2)
                  deltaphi=deltaphi+2*TMath::Pi();
                if(deltaeta == 0 && deltaphi == 0) continue;
              s[xbin]->Fill(deltaeta,deltaphi,1./Ntrig);
            }
        }
        Nevent[xbin]++;
    }
}


void ridge::calcB(){
  treeInt* t_trig = new treeInt(Form("%s",filename.Data()));
  t_trig->Setup();
  treeInt* t_ass = new treeInt(Form("%s",filename.Data()));
  t_ass->Setup();
//  t_ass->SetupExtra();
 // int j[nMix];// index[nMix];
  Int_t nEvents = t_trig->GetEntries();
    for(Int_t i=0;i<nEvents;i++){
        t_trig->GetEntry(i);
        if(i%10000==0) cout<<"has processed "<<i<<" B events"<<endl;
        int ntrk = t_trig->mult; int xbin=-1;
        for(int k=0;k<nbin;k++)
            if(ntrk<trkbin[k]&&ntrk>=trkbin[k+1])
                xbin=k;
            if(xbin<0 || xbin==nbin) continue;
            int Ntrig = 0;
        for(int imult_trig=0;imult_trig<t_trig->mult;imult_trig++){
	    if(t_trig->eta[imult_trig]<etatrigmin||t_trig->eta[imult_trig]>etatrigmax) continue;
            if(t_trig->pt[imult_trig]<pttrigmin||t_trig->pt[imult_trig]>pttrigmax) continue; //event selection
            Ntrig++;
        }
        if(Ntrig<=2) continue;
            /*for(int ira=0;ira<nMix;ira++){
            UInt_t iniseed = SetSeedOwn();
            gRandom->SetSeed(iniseed);
            TRandom3* r = new TRandom3(0);
            j[ira] = r->Integer(nEvents);
            if(i==j[ira]) {ira--; continue;}
            }*/
            for(int ira=0;ira<nMix;ira++){
                int j = i+1+ira;
                if(j>=nEvents) j=j-nEvents;
            t_ass->GetEntry(j);
            for(int imult_trig=0;imult_trig<t_trig->mult;imult_trig++){
                double eta_trig = t_trig->eta[imult_trig];
                double phi_trig = t_trig->phi[imult_trig];
	        if(eta_trig<etatrigmin||eta_trig>etatrigmax) continue;
                if(t_trig->pt[imult_trig]<pttrigmin||t_trig->pt[imult_trig]>pttrigmax) continue; //event selection
                for(int imult_ass=0;imult_ass<t_ass->mult;imult_ass++){
                double eta_ass = t_ass->eta[imult_ass];
                double phi_ass = t_ass->phi[imult_ass];
	        if(eta_ass<etaassmin||eta_ass>etaassmax) continue;
                if(t_ass->pt[imult_ass]<ptassmin||t_ass->pt[imult_ass]>ptassmax) continue; //event selection
                double deltaeta = eta_trig - eta_ass;
                double deltaphi = phi_trig - phi_ass;
              if(deltaphi>TMath::Pi())
                  deltaphi=deltaphi-2*TMath::Pi();
              if(deltaphi<-TMath::Pi()/2)
                  deltaphi=deltaphi+2*TMath::Pi();
                if(deltaeta == 0 && deltaphi == 0) continue;
                b[xbin]->Fill(deltaeta,deltaphi,1.);
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
    

