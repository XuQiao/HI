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
#include <iomanip>
#include "SetSeedOwn.C"
#include "par.h"

using namespace std;

double phicov(double);

class ridge{  //use vector, combine calcS() and calcB() in one event loop
    public:
        ridge(TString);
        ~ridge();
        TVectorD Nevent, tottrk;
        TH1D* hMult[nbin];
        TH2F* s[nbin];
        TH2F* s_sym[nbin];
        TH2F* b[nbin];
        TH2F* b_sym[nbin];
        treeInt *t;
        Int_t nEvents;
        vector<vector<double> > pVectEta_trig;
        vector<vector<double> > pVectPhi_trig;
        vector<vector<double> > pVectEta_ass;
        vector<vector<double> > pVectPhi_ass;
        vector<int> pntrk;
        void beginJob();
        void calc();
        void endJob(TString);

    private:
    TString filename;
};

ridge::ridge(TString _filename){
    filename = _filename;
}

ridge::~ridge(){
    delete t;
    pVectEta_trig.clear();
    vector<vector<double> >().swap(pVectEta_trig);
    pVectPhi_trig.clear();
    vector<vector<double> >().swap(pVectPhi_trig);
    pVectEta_ass.clear();
    vector<vector<double> >().swap(pVectEta_ass);
    pVectPhi_ass.clear();
    vector<vector<double> >().swap(pVectPhi_ass);
    for(int ibin=0;ibin<nbin;ibin++){
    delete s[ibin];
    delete s_sym[ibin];
    delete b[ibin];
    delete b_sym[ibin];
    }
}

void ridge::beginJob(){
    TH1::SetDefaultSumw2();
  for(int ibin=0;ibin<nbin;ibin++){
    hMult[ibin] = new TH1D(Form("mult_%d",ibin),Form("mult_%d",ibin),250,0,250);
    s[ibin] = new TH2F(Form("signal_%d",ibin),Form("signal_%d",ibin),detastep,detamin,detamax,dphistep,dphimin,dphimax);
    s_sym[ibin] = new TH2F(Form("signal_sym_%d",ibin),Form("signal_sym_%d",ibin),detastep,detamin,detamax,dphistep,dphimin,dphimax);
    b[ibin] = new TH2F(Form("background_%d",ibin),Form("background_%d",ibin),detastep,detamin,detamax,dphistep,dphimin,dphimax);
    b_sym[ibin] = new TH2F(Form("background_sym_%d",ibin),Form("background_sym_%d",ibin),detastep,detamin,detamax,dphistep,dphimin,dphimax);
  }
  Nevent.ResizeTo(nbin); Nevent.Zero();
  tottrk.ResizeTo(nbin); tottrk.Zero();
  t = new treeInt(Form("%s",filename.Data()));
  t->Setup();
  nEvents = 10000;//t->GetEntries();
    for(Int_t i=0;i<nEvents;i++){
  //      if(i%10000==0) cout<<"has processed "<<i<<"begin events"<<endl;
        t->GetEntry(i);
        vector<double> pvectEta_trig;
        vector<double> pvectPhi_trig;
        vector<double> pvectEta_ass;
        vector<double> pvectPhi_ass;
        for(int imult=0;imult<t->mult;imult++){
	    if(t->eta[imult]<etatrigmin||t->eta[imult]>etatrigmax) continue;
            if(t->pt[imult]<pttrigmin||t->pt[imult]>pttrigmax) continue; //event selection
            pvectEta_trig.push_back(t->eta[imult]);
            pvectPhi_trig.push_back(t->phi[imult]);
        }
        for(int imult=0;imult<t->mult;imult++){
	    if(t->eta[imult]<etaassmin||t->eta[imult]>etaassmax) continue;
            if(t->pt[imult]<ptassmin||t->pt[imult]>ptassmax) continue; //event selection
            pvectEta_ass.push_back(t->eta[imult]);
            pvectPhi_ass.push_back(t->phi[imult]);
        }
        pntrk.push_back(t->mult);
        pVectEta_trig.push_back(pvectEta_trig);
        pVectPhi_trig.push_back(pvectPhi_trig);
        pVectEta_ass.push_back(pvectEta_ass);
        pVectPhi_ass.push_back(pvectPhi_ass);
    }
}

void ridge::calc(){
   // UInt_t iniseed = SetSeedOwn();
   // gRandom->SetSeed(iniseed);
    TRandom3* r = new TRandom3(0);
   // int nEvents_trig = (int)pVectEta_trig.size();
    int nEvents_ass = (int)pVectEta_ass.size();
   int j[nMix];
    for(Int_t i=0;i<nEvents;i++){
        if(i%5000==0) cout<<"has processed "<<i<<"S events"<<endl;
            vector<double> pvectorEta_trig = pVectEta_trig[i];
            vector<double> pvectorPhi_trig = pVectPhi_trig[i];
            vector<double> pvectorEta_ass = pVectEta_ass[i];
            vector<double> pvectorPhi_ass = pVectPhi_ass[i];
            int Ntrig = pvectorEta_trig.size();
            int Nass = pvectorEta_ass.size();
        int ntrk = pntrk[i]; int xbin=-1;
        for(int k=0;k<nbin;k++)
            if(ntrk<trkbin[k]&&ntrk>=trkbin[k+1])
                xbin=k;
            if(xbin<0 || xbin==nbin) continue;
            tottrk[xbin]+=ntrk;
            hMult[xbin]->Fill(Ntrig);
            if(Ntrig <= 2) continue;
           for(int ira=0;ira<nMix;ira++){
                j[ira] = r->Integer(nEvents_ass);
                if(j[ira]==i) { ira--; continue;}
           }
            for(int imult_trig=0;imult_trig<Ntrig;imult_trig++){
                double eta_trig = pvectorEta_trig[imult_trig];
                double phi_trig = phicov(pvectorPhi_trig[imult_trig]);
            for(int imult_ass=0;imult_ass<Nass;imult_ass++){
              double eta_ass = pvectorEta_ass[imult_ass];
              double phi_ass = phicov(pvectorPhi_ass[imult_ass]);
              double deltaeta = eta_trig - eta_ass;
              double deltaphi = phi_trig - phi_ass;
              if(deltaphi>TMath::Pi())
                  deltaphi=deltaphi-2*TMath::Pi();
              if(deltaphi<-TMath::Pi()/2)
                  deltaphi=deltaphi+2*TMath::Pi();
              if(deltaeta == 0 && deltaphi == 0) continue;
                    double weight = 1.0/Ntrig;
                 s_sym[xbin]->Fill(fabs(deltaeta),fabs(deltaphi),weight);
                 s_sym[xbin]->Fill(-fabs(deltaeta),fabs(deltaphi),weight);
                 s_sym[xbin]->Fill(fabs(deltaeta),-fabs(deltaphi),weight);
                 s_sym[xbin]->Fill(-fabs(deltaeta),-fabs(deltaphi),weight);
                 s_sym[xbin]->Fill(fabs(deltaeta),2*TMath::Pi()-fabs(deltaphi),weight);
                 s_sym[xbin]->Fill(-fabs(deltaeta),2*TMath::Pi()-fabs(deltaphi),weight);
                 s[xbin]->Fill(deltaeta,deltaphi,weight);
            }
                for(int ira=0;ira<nMix;ira++){
                vector<double> pVectorEta_ass = pVectEta_ass[j[ira]];
                vector<double> pVectorPhi_ass = pVectPhi_ass[j[ira]];
                for(int imult_ass=0;imult_ass<(int)pVectorEta_ass.size();imult_ass++){
                double eta_ass = pVectorEta_ass[imult_ass];
                double phi_ass = phicov(pVectorPhi_ass[imult_ass]);
                double deltaeta = eta_trig - eta_ass;
                double deltaphi = phi_trig - phi_ass;
              if(deltaphi>TMath::Pi())
                  deltaphi=deltaphi-2*TMath::Pi();
              if(deltaphi<-TMath::Pi()/2)
                  deltaphi=deltaphi+2*TMath::Pi();
                if(deltaeta == 0 && deltaphi == 0) continue;
                    double weight = 1.0/Ntrig/pVectorEta_ass.size();
                 b_sym[xbin]->Fill(fabs(deltaeta),fabs(deltaphi),weight);
                 b_sym[xbin]->Fill(-fabs(deltaeta),fabs(deltaphi),weight);
                 b_sym[xbin]->Fill(fabs(deltaeta),-fabs(deltaphi),weight);
                 b_sym[xbin]->Fill(-fabs(deltaeta),-fabs(deltaphi),weight);
                 b_sym[xbin]->Fill(fabs(deltaeta),2*TMath::Pi()-fabs(deltaphi),weight);
                 b_sym[xbin]->Fill(-fabs(deltaeta),2*TMath::Pi()-fabs(deltaphi),weight);
                 b[xbin]->Fill(deltaeta,deltaphi,weight);
                }
            }
            }
        Nevent[xbin]++;
    }
}


void ridge::endJob(TString outfile){
    TFile *fout = new TFile(Form("%s",outfile.Data()),"Recreate");
    fout->cd();
    Nevent.Write("Nevent");
    tottrk.Write("tottrk");
    for(int ibin=0;ibin<nbin;ibin++){
        hMult[ibin]->Write();
        s[ibin]->Write();
        s_sym[ibin]->Write();
        b[ibin]->Write();
        b_sym[ibin]->Write();
    }
    fout->Close();
}

double phicov(double phi){
    double phic;
if(phi<-TMath::Pi())
    phic = phi+2*TMath::Pi();
else if(phi>TMath::Pi())
    phic = phi-2*TMath::Pi();
else 
    phic = phi;
return phic;
}

