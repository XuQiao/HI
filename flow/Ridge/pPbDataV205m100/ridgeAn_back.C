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
        vector<vector<double> > pVectEta_trig;
        vector<vector<double> > pVectPhi_trig;
        vector<vector<double> > pVectEta_ass;
        vector<vector<double> > pVectPhi_ass;
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
    pVectPhi.clear();
    vector<vector<double> >().swap(pVectPhi);
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
        vector<double> pvectPt;
        for(int imult=0;imult<t->mult;imult++){
           // TVector3 pvector;
           // pvector.SetPtEtaPhi(t->pt[imult],t->eta[imult],t->phi[imult]);
            pvectEta.push_back(t->eta[imult]);
            pvectPhi.push_back(t->phi[imult]);
            pvectPt.push_back(t->pt[imult]);
        }
        pVectEta.push_back(pvectEta);
        pVectPhi.push_back(pvectPhi);
        pVectPt.push_back(pvectPt);
    }
}

void ridge::calcS(){
    for(Int_t i=0;i<nEvents;i++){
        if(i%10000==0) cout<<"has processed "<<i<<"S events"<<endl;
            vector<double> pvectorEta = pVectEta[i];
            vector<double> pvectorPhi = pVectPhi[i];
            vector<double> pvectorPt = pVectPt[i];
            int trig_mult = pvectorEta.size();
        int ntrk = trig_mult; int xbin=-1;
        for(int j=0;j<nbin;j++)
            if(ntrk<trkbin[j]&&ntrk>=trkbin[j+1])
                xbin=j;
            if(xbin<0 || xbin==nbin) continue;
            tottrk[xbin]+=ntrk;
            int Ntrig = 0;
        for(int imult1=0;imult1<trig_mult;imult1++){
              double eta1 = pvectorEta[imult1];
              double pt1 = pvectorPt[imult1];
	        if(eta1<etatrigmin||eta1>etatrigmax) continue;
                if(pt1<pttrigmin||pt1>pttrigmax) continue; //event selection
                Ntrig++;
        }
        for(int imult1=0;imult1<trig_mult;imult1++){
              double eta1 = pvectorEta[imult1];
              double phi1 = pvectorPhi[imult1];
              double pt1 = pvectorPt[imult1];
	        if(eta1<etatrigmin||eta1>etatrigmax) continue;
                if(pt1<pttrigmin||pt1>pttrigmax) continue; //event selection
            for(int imult2=0;imult2<trig_mult;imult2++){
                if(imult1==imult2) continue;
              double eta2 = pvectorEta[imult2];
              double phi2 = pvectorPhi[imult2];
              double pt2 = pvectorPt[imult1];
	        if(eta2<etatrigmin||eta2>etatrigmax) continue;
                if(pt2<pttrigmin||pt2>pttrigmax) continue; //event selection
              s[xbin]->Fill(eta1-eta2,phi1-phi2,1./Ntrig);
            }
        }
        Nevent[xbin]++;
    }
}


void ridge::calcB(){
    vector<vector<double> > pvectorEta_ass;
    vector<vector<double> > pvectorPhi_ass;
    vector<vector<double> > pvectorPt_ass;
    for(Int_t i=0;i<nEvents;i++){
        if(i%5000==0) cout<<"has processed "<<i<<" B events"<<endl;
            vector<double> pvectorEta_trig = pVectEta[i];
            vector<double> pvectorPhi_trig = pVectPhi[i];
            vector<double> pvectorPt_trig = pVectPt[i];
            int trig_mult = pvectorEta_trig.size();
        int ntrk = trig_mult; int xbin=-1;
        for(int k=0;k<nbin;k++)
            if(ntrk<trkbin[k]&&ntrk>=trkbin[k+1])
                xbin=k;
            if(xbin<0 || xbin==nbin) continue;
            for(int ira=0;ira<nMix;ira++){
            UInt_t iniseed = SetSeedOwn();
            gRandom->SetSeed(iniseed);
            TRandom3* r = new TRandom3(0);
            int j = r->Integer(nEvents);
            if(i==j){ira--; continue;}
            pvectorEta_ass.push_back(pVectEta[j]);
            pvectorPhi_ass.push_back(pVectPhi[j]);
            pvectorPt_ass.push_back(pVectPt[j]);
            }
            int Ntrig = 0;
            for(int imult_trig=0;imult_trig<trig_mult;imult_trig++){
                double eta_trig = pvectorEta_trig[imult_trig];
                double pt_trig = pvectorPt_trig[imult_trig];
	        if(eta_trig<etatrigmin||eta_trig>etatrigmax) continue;
                if(pt_trig<pttrigmin||pt_trig>pttrigmax) continue; //event selection
                Ntrig++;
            }
            for(int imult_trig=0;imult_trig<trig_mult;imult_trig++){
                double eta_trig = pvectorEta_trig[imult_trig];
                double phi_trig = pvectorPhi_trig[imult_trig];
                double pt_trig = pvectorPt_trig[imult_trig];
	        if(eta_trig<etatrigmin||eta_trig>etatrigmax) continue;
                if(pt_trig<pttrigmin||pt_trig>pttrigmax) continue; //event selection
                for(int ira=0;ira<nMix;ira++){
                    int ass_mult = pvectorEta_ass[ira].size();
                for(int imult_ass=0;imult_ass<ass_mult;imult_ass++){
                double eta_ass = pvectorEta_ass[ira][imult_ass];
                double phi_ass = pvectorPhi_ass[ira][imult_ass];
                double pt_ass = pvectorPt_ass[ira][imult_ass];
	        if(eta_ass<etaassmin||eta_ass>etaassmax) continue;
                if(pt_ass<ptassmin||pt_ass>ptassmax) continue; //event selection
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
    

